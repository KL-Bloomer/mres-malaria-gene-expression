import pandas
import os

sample_sheet = pandas.read_csv(config['ss'], sep= '\t', comment= '#')
sample_sheet['R1'] = [os.path.join(config['datadir'], fq) for fq in sample_sheet['R1']]
sample_sheet['R2'] = [os.path.join(config['datadir'], fq) for fq in sample_sheet['R2']]

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

PbANKA_chroms = ['PbANKA_01_v3',
    'PbANKA_02_v3',
    'PbANKA_03_v3',
    'PbANKA_04_v3',
    'PbANKA_05_v3',
    'PbANKA_06_v3',
    'PbANKA_07_v3',
    'PbANKA_08_v3',
    'PbANKA_09_v3',
    'PbANKA_10_v3',
    'PbANKA_11_v3',
    'PbANKA_12_v3',
    'PbANKA_13_v3',
    'PbANKA_14_v3',
    'PbANKA_00_v3_archived_contig_1',
    'PbANKA_00_v3_archived_contig_2',
    'PbANKA_00_v3_archived_contig_3',
    'PbANKA_00_v3_archived_contig_4',
    'PbANKA_00_v3_archived_contig_5',
    'PbANKA_API_v3',
    'PbANKA_MIT_v3']

rule final_output:
    # The only purpose of this rule is listing the files we want as final
    # output of the workflow. Snakemake will use the rules after this one to
    # generate this files.
    input:
        'multiqc/fastqc_report.html',
        'featureCounts/counts.tsv',
        expand('blast_species/{library_id}.species.tsv', library_id= sample_sheet['library_id']),
        'idxstats/idxstats.tsv',
        'barplot_libsizes_beforenorm.png',
        'edger/differential_gene_expression.tsv',
        'edger/geneid_desc_table.tsv',
        'edger/logrpkm_long.tsv',
        'edger/MDSplots_concatenated.png',
        'edger/MAplot_all_contrasts.png',
        'edger/Volcano_plot_all_contrasts.png',
        'edger/globalexpression_all_contrasts.png',
        'Heatmap_DE_genes.png',
        'clusters_table.tsv',
        'avergene_expr_clusters.png',
        'Heatmap_AP2_genes.png',
        'Heatmap_AP2_genes_FDR.png',
        'Heatmap_DE_genes_logFC.png',
        'Heatmap_genes.png',
        'gene_expression_changes_keygenes.png',
        'topGO_table_clusters.tsv', 

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

rule download_mouse_reference:
    output:
        fa= 'ref/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa',
    shell:
        r"""
        curl -s ftp://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz \
        | zcat > {output.fa}
        """

rule concatenate_genomes:
    input:
        pb= 'ref/PlasmoDB-49_PbergheiANKA_Genome.fasta',
        mm= 'ref/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa',
    output:
        fa= 'ref/PlasmoDB-49_PbergheiANKA-Mus_musculus_GRCm38.fa',
    shell:
        r"""
        cat {input.pb} {input.mm} > {output.fa}
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
        fasta= 'ref/PlasmoDB-49_PbergheiANKA-Mus_musculus_GRCm38.fa',
    output:
        index= 'ref/PlasmoDB-49_PbergheiANKA-Mus_musculus_GRCm38.fa.8.ht2',
    shell:
        r"""
        hisat2-build -p 16 --seed 1234 {input.fasta} {input.fasta}
        """

rule align_reads:
    input:
        R1= lambda wildcard: sample_sheet[sample_sheet['library_id'] == wildcard.library_id].R1,
        R2= lambda wildcard: sample_sheet[sample_sheet['library_id'] == wildcard.library_id].R2,
        genome= 'ref/PlasmoDB-49_PbergheiANKA-Mus_musculus_GRCm38.fa',
        index= 'ref/PlasmoDB-49_PbergheiANKA-Mus_musculus_GRCm38.fa.8.ht2',
    output:
        bam= 'hisat2/{library_id}.bam',
        bai= 'hisat2/{library_id}.bam.bai',
        log= 'hisat2/{library_id}.log',
    shell:
        r"""
        hisat2 --summary-file {output.log} --new-summary \
               --max-intronlen 5000 --threads 16 -x {input.genome} -1 {input.R1} -2 {input.R2} \
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
        featureCounts -p -T 8 -Q 10 -s 2 -t exon -g gene_id -a {input.gff} -o {output.counts} {input.bam}
        """

rule samtools_idxstats:
    input:
        bam= 'hisat2/{library_id}.bam',
    output:
        idxstats= temp('idxstats/{library_id}.tsv'),
    shell:
        r"""
        samtools idxstats {input.bam} \
        | awk -v library_id='{wildcards.library_id}' -v FS='\t' -v OFS='\t' '{{print $0, library_id}}' > {output.idxstats}
        """

rule concatenate_idxstats:
    input:
        stats= expand('idxstats/{library_id}.tsv', library_id= sample_sheet['library_id']),
    output:
        stats= 'idxstats/idxstats.tsv',
    shell:
        r"""
        echo "chrom length mapped unmapped library_id" | tr ' ' '\t' > {output.stats}
        cat {input.stats} >> {output.stats}
        """

rule filter_bam:
    input:
        bam= 'hisat2/{library_id}.bam',
    output:
        bam= temp('bigwig/{library_id}.bam'),
        bai= temp('bigwig/{library_id}.bam.bai'),
    params:
        chroms= ' '.join(PbANKA_chroms),
    shell:
        r"""
        samtools view -u -q 10 -@ 4 {input.bam} {params.chroms} > {output.bam}
        samtools index {output.bam}
        """

rule bam_to_bigwig:
    input:
        bam= 'bigwig/{library_id}.bam',
        bai= 'bigwig/{library_id}.bam.bai',
    output:
        bw= 'bigwig/{library_id}.bw',
    shell:
        r"""
        bamCoverage -b {input.bam} -o {output.bw} \
            --normalizeUsing BPM \
            --numberOfProcessors 8
        """

rule download_blastdb:
    output:
        'blastdb/taxdb.btd',
    shell:
        r"""
	cd `dirname {output}`
        update_blastdb.pl --decompress nt
        """

rule extract_candidate_contaminating_reads:
    input:
        bam= 'hisat2/{library_id}.bam',
    output:
        fasta= 'fasta/{library_id}.fasta',
    params:
        # This is effectively the number of *lines* in output. Since a read
        # takes two lines (one for the fasta header, one for the sequence), the
        # number of reads in output fasta will be half this number
        n_reads= 5000 * 2,
    shell:
        r"""
        samtools fasta -n -F 128 -f 8 -f 4 {input.bam} \
        | awk 'NR <= {params.n_reads}' > {output.fasta}
        """

rule blast_reads:
    input:
        query= 'fasta/{library_id}.fasta',
        db= 'blastdb/taxdb.btd',
    output:
        out= 'blast/{library_id}.tsv',
    shell:
        r"""
        blastdir=`realpath {input.db}`
        blastdir=`dirname $blastdir`
        export BLASTDB=$blastdir

        header="ssciname scomname stitle sseqid sstart send slen qseqid qlen evalue pident qcovhsp"
        echo $header | tr ' ' '\t' > {output.out}
        blastn -task megablast \
               -query {input.query} \
               -db $blastdir/nt \
               -outfmt "6 $header" \
               -max_target_seqs 1 \
               -subject_besthit \
               -word_size 28 \
               -num_threads 8 >> {output.out}
        """
rule summarise_contamination:
      input:
         query= 'fasta/{library_id}.fasta',
         blast= 'blast/{library_id}.tsv',
      output:
         summary= 'blast_species/{library_id}.species.tsv',
      shell:
         r"""
cat <<'EOF' > {rule}.$$.tmp.R

library(data.table)

blast <- fread('{input.blast}')
query <- fread(cmd= 'grep ">" {input.query}', header= FALSE, col.names= 'qseqid', sep= '\t')
query[, qseqid := sub('>', '', qseqid)]

smry <- blast[, list(n_hits= sum(evalue < 0.01)), by= list(species= ssciname)][order(-n_hits)]
smry[, pct_hits := 100 * n_hits/nrow(query)]
smry[, cum_pct_hits := cumsum(pct_hits)]

write.table(smry, '{output.summary}', sep= '\t', row.names= FALSE, quote= FALSE)

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """

rule library_size_barplot:
    input:
        sample_sheet= config['ss'],
        counts= 'featureCounts/counts.tsv',
    output:
        barplot_libsizes_beforenorm= 'barplot_libsizes_beforenorm.png',
    script:
        os.path.join(workflow.basedir, 'scripts/Script_for_barplot.R')

rule differential_gene_expression:
    input:
        sample_sheet= config['ss'],
        counts= 'featureCounts/counts.tsv',
        gff= 'ref/PlasmoDB-49_PbergheiANKA.gff',
    output:
        dge_table= 'edger/differential_gene_expression.tsv',
        geneid_desc_table= 'edger/geneid_desc_table.tsv',
        logrpkm_table= 'edger/logrpkm_long.tsv',
        MDSplots_concatenated= 'edger/MDSplots_concatenated.png',
        MA_plot= 'edger/MAplot_all_contrasts.png',
        Volcano_plot= 'edger/Volcano_plot_all_contrasts.png',
        globalexpression= 'edger/globalexpression_all_contrasts.png',
    script:
        os.path.join(workflow.basedir, 'scripts/Script_for_RNAseq_plots.R')

rule heatmap_and_clustering:
    input:
        sample_sheet= config['ss'],
        logrpkm_table= 'edger/logrpkm_long.tsv',
        dge_table= 'edger/differential_gene_expression.tsv',
        geneid_desc_table= 'edger/geneid_desc_table.tsv',
        GAF= 'ref/PlasmoDB-49_PbergheiANKA_GO.gaf',
    output:
        Heatmap_DE_genes= 'Heatmap_DE_genes.png',
        clusters_table= 'clusters_table.tsv',
        avergene_expr_clusters= 'avergene_expr_clusters.png',
        Heatmap_AP2_genes= 'Heatmap_AP2_genes.png',
        Heatmap_AP2_genes_FDR= 'Heatmap_AP2_genes_FDR.png',
        Heatmap_DE_genes_logFC= 'Heatmap_DE_genes_logFC.png',
        Heatmap_genes= 'Heatmap_genes.png',
    script:
        os.path.join(workflow.basedir, 'scripts/heatmap.R')

rule interesting_gene_plot:
    input:
        sample_sheet= config['ss'],
        logrpkm_table= 'edger/logrpkm_long.tsv',
        interesting_genes= os.path.join(workflow.basedir, 'Interesting_genes.txt'),
        GAF= 'ref/PlasmoDB-49_PbergheiANKA_GO.gaf',
    output:
        gene_expression_changes_keygenes= 'gene_expression_changes_keygenes.png',
    script:
         os.path.join(workflow.basedir, 'scripts/gene_expression_interesting_genes.R')

rule topGO_clusters:
    input:
        clust= 'clusters_table.tsv',
        gene_id_table= 'edger/geneid_desc_table.tsv',
        GAF= 'ref/PlasmoDB-49_PbergheiANKA_GO.gaf',
    output:
        topGO_table_clusters= 'topGO_table_clusters.tsv',
    script:
        os.path.join(workflow.basedir, 'scripts/GO.R')
