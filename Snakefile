import pandas
import os

sample_sheet = pandas.read_csv(config['ss'], sep= '\t', comment= '#')

def get_fastq_basenames(fastq_file_list):
    """fastq_file_list is a list of fastq file names. Extract from each
    filename the basename by stripping directory and suffixes. Return a
    dictionary mapping basename to file.
    """
    fastq_names = [os.path.basename(x) for x in fastq_file_list] # Remove directory path
    fastq_names = [x.replace('.fastq.gz', '') for x in fastq_names] # Remove .fastq.gz suffix
    if len(fastq_names) != len(set(fastq_names)):
        raise Exception('Found duplicate fastq basenames')
    fastq_dict = dict(zip(fastq_names, fastq_file_list))
    return fastq_dict

fastq_file_list = list(sample_sheet['R1']) + list(sample_sheet['R2'])
fastq_dict = get_fastq_basenames(fastq_file_list)

rule final_output:
    # The only porpose of this rule is listing the files we want as final
    # output of the workflow. Snakemake will use the rules after this one to
    # generate this files.
    input:
        'multiqc/fastqc_report.html',
        expand('hisat2/{library_id}.bam', library_id= sample_sheet['library_id']),
        'featureCounts/counts.tsv',
        expand('blast/{library_id}.species.tsv', library_id= sample_sheet['library_id']),

# ------
# NB: With the exception of the first rule, which determines the final output,
# the order of the following rules does not matter. Snakemake will chain them in
# the right order using the rules' input/outputs. If there are gaps in the
# chain or ambiguities, like different rules producing the same output,
# snakemake will notified it with a loud error.
# ------

rule check_read_quality:
    input:
        fastq= lambda wildcard: fastq_dict[wildcard.fastq_name], # Get the fastq file corresponding to wildcard {fastq_name}
    output:
        qc= 'fastqc/{fastq_name}_fastqc.zip',
    shell:
        r"""
        fastqc --noextract --outdir fastqc {input.fastq}
        """

rule combine_read_quality_reports:
    input:
        fastqc_reports= expand('fastqc/{fastq_name}_fastqc.zip', fastq_name= fastq_dict.keys()),
    output:
        'multiqc/fastqc_report.html',
    shell:
        r"""
        multiqc --force --outdir multiqc --filename fastqc_report.html {input.fastqc_reports}
        """

rule download_reference_genome:
    output:
        fasta= 'ref/PlasmoDB-49_PbergheiANKA_Genome.fasta',
    shell:
        r"""
        curl -s https://plasmodb.org/common/downloads/release-49/PbergheiANKA/fasta/data/PlasmoDB-49_PbergheiANKA_Genome.fasta > {output.fasta}
        """

rule prepare_reference_annotation:
    output:
        gff= 'ref/PlasmoDB-49_PbergheiANKA.gff',
    shell:
        # addGeneIdToGff.py is a custom is script that adds to each feature the
        # gene_id from which the feature comes from. This is necessary for
        # later steps.
        r"""
        curl -s https://plasmodb.org/common/downloads/release-49/PbergheiANKA/gff/data/PlasmoDB-49_PbergheiANKA.gff \
        | {workflow.basedir}/scripts/addGeneIdToGff.py -tss '' > {output.gff}
        """

rule index_genome:
    input:
        fasta= 'ref/PlasmoDB-49_PbergheiANKA_Genome.fasta',
    output:
        index= 'ref/PlasmoDB-49_PbergheiANKA_Genome.fasta.8.ht2',
    shell:
        r"""
        hisat2-build --seed 1234 {input.fasta} {input.fasta}
        """

rule align_reads:
    input:
        R1= lambda wildcard: sample_sheet[sample_sheet['library_id'] == wildcard.library_id].R1,
        R2= lambda wildcard: sample_sheet[sample_sheet['library_id'] == wildcard.library_id].R2,
        genome= 'ref/PlasmoDB-49_PbergheiANKA_Genome.fasta',
        index= 'ref/PlasmoDB-49_PbergheiANKA_Genome.fasta.8.ht2',
    output:
        bam= 'hisat2/{library_id}.bam',       
        bai= 'hisat2/{library_id}.bam.bai',
        log= 'hisat2/{library_id}.log',
    shell:
        r"""
        hisat2 --summary-file {output.log} --new-summary \
                --max-intronlen 5000 --threads 1 -x {input.genome} -1 {input.R1} -2 {input.R2} \
        | samtools sort -@ 1 > {output.bam}

        samtools index -@ 1 {output.bam}
        """

rule count_reads_in_genes:
    input:
        bam= expand('hisat2/{library_id}.bam', library_id= sample_sheet['library_id']),
        gff= 'ref/PlasmoDB-49_PbergheiANKA.gff',
    output:
        counts= 'featureCounts/counts.tsv',
    shell:
        r"""
        featureCounts -p -T 1 -Q 10 -s 2 -t exon -g gene_id -a {input.gff} -o {output.counts} {input.bam}
        """

rule analyse_differential_gene_expression:
    input:
        sample_sheet= config['ss'],
        counts= 'featureCounts/counts.tsv', 
    output:
        mds= 'edger/mds_plot.pdf',
        dge_table= 'edger/differential_gene_expression.tsv',
        maplot= 'edger/maplot.png',
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R

library(data.table)
library(edgeR)

ss <- fread('{input.sample_sheet}')
counts <- fread(cmd= 'grep -v "^#" {input.counts}')

setnames(counts, names(counts), sub('.bam', '', basename(names(counts))))
counts <- counts[, c('Geneid', ss$library_id), with= FALSE]

mat <- as.matrix(counts, rownames= 'Geneid')
ss[, group := paste(stage, genotype, sep= '_')]

design <- model.matrix(~0 + ss$group)
colnames(design) <- sub('ss$group', '', colnames(design), fixed= TRUE)

y <- DGEList(counts= mat, group= ss$group)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design)

# MDS plot before and after batch correction

pdf('{output.mds}', width= 14/2.54, height= 14/2.54)
mds <- plotMDS(y, col= as.numeric(as.factor(ss$group)))
dev.off()

# Differential expression
# =======================

# These must be the same as those for ATAC
contrasts <- makeContrasts(
    ES_820_vs_ES_G2098= ES_820 - ES_G2098,
    LS_820_vs_LS_G2098= LS_820 - LS_G2098,
    ES_820_vs_ES_G2107= ES_820 - ES_G2107,
    LS_820_vs_LS_G2107= LS_820 - LS_G2107,
    levels= design)
stopifnot(all(abs(colSums(contrasts)) < 1e-6))

fit <- glmFit(y, design, prior.count= 1)
dge <- list()
for(cntr in colnames(contrasts)){{
    print(cntr)
    lfc <- glmTreat(fit, contrast= contrasts[, cntr], lfc= log2(1.5))
    detable <- topTags(lfc, n= nrow(y))$table
    detable$gene_id <- row.names(detable)
    detable <- data.table(detable)
    detable[, contrast := cntr]
    dge[[length(dge)+1]] <- detable
}}
dge <- rbindlist(dge)
dge[, unshrunk.logFC := NULL]

gz <- gzfile('{output.dge_table}', 'w')
write.table(dge, gz, sep= '\t', row.names= FALSE, quote= FALSE)
close(gz)

xord <- c('inf_mosq_percent', 'is_inf_mosq', 'log2_oocysts_per_mosq', 'is_oocysts_per_mosq', 'log2_exfl_XA_per_ml', 'is_exfl', 'log2_exfl_XA_per_ml_only_pos')
dge[, trait := factor(trait, xord)]
nsig <- dge[, list(n_up= sum(FDR < 0.01 & logFC > 0), n_down= sum(FDR < 0.01 & logFC < 0)), trait]
nsig[, n_up:= sprintf('Up = %s', n_up)]
nsig[, n_down:= sprintf('Down = %s', n_down)]

gg <- ggplot(data= dge, aes(x= logCPM, y= logFC)) +
    geom_point(alpha= 0.5, pch= '.') +
    geom_point(data= dge[FDR < 0.01], alpha= 0.5, colour= 'red', pch= '.') +
    geom_smooth(se= FALSE, col= 'grey60', size= 0.1) +
    geom_hline(yintercept= 0, colour= 'blue') +
    geom_text(data= nsig, x= Inf, y= Inf, aes(label= n_up), vjust= 1.3, hjust= 1.1, size= 3) +
    geom_text(data= nsig, x= Inf, y= -Inf, aes(label= n_down), vjust= -1.2, hjust= 1.1, size= 3) +
    facet_wrap(~trait, scales= 'free_y', ncol= 2) +
    theme_light() +
    theme(strip.text= element_text(colour= 'black'))
ggsave('{output.maplot}', width= 16, height= 20, units= 'cm')

xdge <- dcast(dge, gene_id ~ trait, value.var= c('FDR', 'logFC'))
exfl_keep <- unique(c(
    xdge[FDR_is_exfl < 0.01 & abs(logFC_is_exfl) > 2]$gene_id,
    xdge[FDR_is_exfl < 0.01 & FDR_log2_exfl_XA_per_ml_only_pos < 0.01]$gene_id
))
xdge <- xdge[gene_id %in% exfl_keep][order(logFC_log2_exfl_XA_per_ml)]

dat <- logcpm[gene_id %in% exfl_keep & !is.na(log2_exfl_XA_per_ml) ]
dat[, gene_id := factor(gene_id, xdge$gene_id)]

gg <- ggplot(data= dat, aes(x= log2_exfl_XA_per_ml, y= logcpm)) +
    geom_point(pch= '.') +
    facet_wrap(~ gene_id, ncol= 5, scales= 'free_y') +
    geom_smooth(method= 'lm', size= 0.5) +
    theme_minimal() +
    theme(axis.text= element_text(size= 6))

dge_genes <- merge(xdge, genes, by= 'gene_id', sort= FALSE)[, list(gene_id, logFC_is_exfl, logFC_log2_exfl_XA_per_ml, description)]
dge_genes[, logFC_is_exfl := sprintf('%.2f', logFC_is_exfl)]
dge_genes[, logFC_log2_exfl_XA_per_ml := sprintf('%.2f', logFC_log2_exfl_XA_per_ml)]

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """

