### mres-malaria-gene-expression

#Learning about vim: https://vimhelp.org/usr_02.txt.html#02.7

#Created & cloned git repository
https://www.melbournebioinformatics.org.au/tutorials/tutorials/using_git/Using_Git/
#to backtrack to a different version of the files
git log #take log string
git checkout {log string} test.text
git checkout HEAD text.text

git status
git add
git commit -m "message"
git push
git pull
*don't need to use pull and create a new branch when updating files etc.

#Conda & Bioconda installation
https://bioconda.github.io/user/install.html

#Then create conda environment and install the programmes on the file
bash
conda create --yes --name mres-malaria-gene-expression
conda activate mres-malaria-gene-expression
mamba install --freeze-installed -n mres-malaria-gene-expression --yes --file requirements.txt
Note to self: when installing using mamba from the cluster, need to install in the base environment which will translate to the mres-malaria-gene-expression environment

realpath $CONDA_PREFIX
* where the current activated environment lives

#Download, install & connect to VPN Client
https://www.gla.ac.uk/myglasgow/it/vpn/#downloadvpnclientforyourdevice
https://www.cisco.com/c/en/us/support/docs/smb/routers/cisco-rv-series-small-business-routers/Kmgmt-785-AnyConnect-Linux-Ubuntu.html
From the user interface, for connect, type: gucsasa1.cent.gla.ac.uk, followed by username (2592613b) & password

#To make logging into headnode easier
atom ~/.bashrc #opens a cryptic file with default settings for the shell
alias cluster="ssh -oKexAlgorithms=+diffie-hellman-group1-sha1 2592613b@headnode03.cent.gla.ac.uk"
Therefore, just need to open terminal and type "cluster"
ssh bioinf03.hpc.gla.ac.uk

Installed conda and bioconda as on PC.

#Accessing fastqc files and saving into "mres-malaria-gene-expression" folder
 mkdir -p ~/mres-malaria-gene-expression
cp -r /export/home/db291g/Kirstin_fastq ~/mres-malaria-gene-expression/

* To avoid this issue "Unable to negotiate with 130.209.19.28 port 22: no matching key exchange method found...".
atom ~/.ssh/config

, In the atom textfile
Host headnode03.cent.gla.ac.uk
    KexAlgorithms +diffie-hellman-group1-sha1 #save and close

#How to copy requirements file from PC to cluster
scp -r /path/to/requirements.txt 2592613b@headnode03.cent.gla.ac.uk:~/mres-malaria-gene-expression/

#How to copy from cluster to PC - from PC terminal, type the following:
scp -r 2592613b@headnode03.cent.gla.ac.uk:~/mres-malaria-gene-expression/path/to/file /path/to/file/on/local/computer
* need to be on PC to execute the following commands

#To run fastqc
mkdir -p fastqc

#nano script
for fq in path/to/file/data/*.fq
do
	fastqc --noextract --outdir path/to/directory/fastqc $fq
done

#How to open fastqc.html file
firefox file.html

#Creating the sample sheet
*library_id = library name i.e. Ap20-GFP-SFC-16h-R1 - name with only letters, numbers, dots, underscores, hyphens
*Time
*Strain
*R1 Full path to read 1 fastq file
*R2 Full path to read 2 fastq file

*Prepared with Excel - save as plain text with columns separated by TAB (not comma or space).

#To visualise if there are any spaces; where tabs are replaced by ^I; end of each line = $
cat -vet export/home/2592613b/mres-malaria-gene-expression/sample_sheet.tsv
sed -i 's/ //g' ~/mres-malaria-gene-expression/sample_sheet.tsv #the -i option means it will do the replacement directly on the file; "s" is substitution and "/" refers to delimiters...i.e. sed 's /unix/linux' - unix is search pattern & linux is replacement

#Create snakefile_run; remove --dry-run when you want to run the commands (otherwise only shows what would be executed)
snakemake --dry-run --printshellcmds --jobs 1 \
--config ss=$PWD/sample_sheet.tsv \
--directory output \
--snakefile Snakefile
* For future reference, in the Snakefile, run hisat2 command option "--threads 16" instead of "--threads 1" & for featureCounts command "-T 8" instead of "-T 1".

#IGV
Transferred hisat2 files to PC; installed IGV =2.9.2
mamba install --freeze-installed -n MRes_project_2021 --yes --file requirements.txt
open IGV using the command "igv"
 #took pictures of the top 5 genes from the count data - alt/splicing?
 #take pictures of the top 5 or so genes that are differentially expressed.
 #scale to max tracks

#### Edit & re-run the snakefile - download the M.Musculus genome and concatenate with P.berghei
Evidence of contamination due to a high GC content and bimodal graph from the fastqc/multiqc analyses. Therefore, it was decided to blast the unmapped reads and to re-run the mapping with a M.musulus-P.berghei concatenated genome. This would eliminate reads that map preferentially to the M/musculus genome from the downstream analysis (i.e. read counts).
Note to self: No tabs in the snakefile, only whitespace.
Changed the shell command: cat {input.pb} {output.fa} > {output.fa} to cat {input.pb} {input.mm} > {output.fa}

####Samtools
#Number of reads mapped to mouse or plasmodium genome
 samtools idxstats <my.bam>
 # where output columns are - chromosome name, size, number of alignments/mapped reads, reads unmapped but a   mapped mate
 #open a bam file
 samtools view <my.bam> | less -S
 #convert bam to fasta
 samtools fasta <my.bam>

 broadinstitute/picard/explain flags...
 flag 4 = unmapped

####Ran the over-respresented sequences in blast which matched to M.musculus
####Ran the unmapped reads from the files in blast as recordered in Snakefile:

#Corrupt download when trying to download the blast database, therefore:
# Remove leftover from download (this was not necessary as my download seemed to have failed completely):
rm -r /export/home/2592613b/mres-malaria-gene-expression/output/blastdb

# Symlink * could use cp to make hard copies but this will prevent the need to duplicate large files
ln -s /export/home/db291g/Tritume/blastdb /export/home/2592613b/mres-malaria-gene-expression/output/
https://linuxize.com/post/how-to-create-symbolic-links-in-linux-using-the-ln-command/

#And ultimately, read the blast database from Dario's home directory.
#Note to self: can execute a snakemake wrapper script in the background using nohup

nohup ./snakefile_run &
#Therefore, even if I am logged out of the cluster, snakemake will still be going in the background; to see the output from snakemake = less nohup.out; to kill snakemake = pkill snakemake OR immediate kill(emergency) = pkill -9 snakemake

#Issues being logged out:
snakemake --printshellcmds --jobs 1 \
    --config ss=$PWD/sample_sheet.tsv \
    --directory output \
    --snakefile Snakefile \
    --dry-run \
    --unlock

#Snakemake logs under .snakemake in the working directory:
ls -a
/export/home/2592613b/mres-malaria-gene-expression/output/.snakemake/log/

#To get around the "AmbiguousRuleException" error when running rule summarise_contamination, performed the following:
1) output:
    summary= 'blast_species/{library_id}.species.tsv',
		#and made corrections to rule final_output
2) wildcard_constraints:
    library_id = '|'.join([re.escape(x) for x in sample_sheet['library_id']])
		#add this before the rule final_output - gets rid of "greedy behaviour"

#Re-ran the snakemake file on the female gametocyte data
* included female gametocytes (Strain = RFP & Time = 0h)
* *SRR526055* from:
	* https://www.ebi.ac.uk/ena/browser/view/PRJNA374918
* *21725* from:
	* https://www.ebi.ac.uk/ena/browser/view/PRJNA714084
* edited sample sheet to include female gametocyte data
* If you don't want blast to run on female gametocyte data, edit rule final_output:
expand('blast_species/{library_id}.species.tsv', library_id= sample_sheet[sample_sheet['Time'] != 0]['library_id'])
* Prioritise the gene count output:
snakemake --prioritize count_reads_in_genes --jobs 1 ...

#Created MDS plot  - later removed outliers and rRNA contamination and PIR/fam genes (R script)
#removal of the PIR/fam genes did not affect the MDS plot/MA plot, therefore kept those genes in.
#From the MDS plot, there appears to be a batch effect with the RM samples clustering away from the other samples (i.e. GFPcon & Ap20 are confounded by strain). Therefore corrected for batch effect in all subsequent plots.

	rm(list=ls())
	library(data.table)
	library(edgeR)

	setwd('~/MRes_Malaria_2021/git_repos/mres-malaria-gene-expression/')
	getwd()

	ss_min_outliers <- fread('sample_sheet_minusoutliers.tsv')
	ss_min_outliers[, Time := sprintf('%.2d', Time)]
	#Note to self - reassign "time" as a character vector which can sort accordingly (above function)
	ss_min_outliers[, group := paste(Time)]
	ss_min_outliers$Time <- as.factor(ss_min_outliers$Time)
	class(ss_min_outliers$Time)
	# Sanity check that Time is a factor variable
stopifnot(is.factor(ss_min_outliers$Time))

setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/featureCounts/')
counts <- fread(cmd= 'grep -v "^#" counts.tsv')
setnames(counts, names(counts), sub('.bam', '', basename(names(counts))))
counts_min_outliers <- counts[, c('Geneid', ss_min_outliers$library_id), with= FALSE]

library(dplyr)
library(tidyr)
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/ref/')
GFF <- fread(cmd = 'grep -v "^#" PlasmoDB-49_PbergheiANKA.gff')
rRNA_GFF <- GFF %>%
  filter(V3 == "rRNA")
rRNA_ID <- subset(rRNA_GFF, select= V9)
RNA_ID_sep <- separate(data = rRNA_ID, col = V9, into = c("Geneid_Feature", "Parent", "Description", "Geneid"), sep = "([;])")
RNA_Gene_ID <- subset(RNA_ID_sep, select= Geneid)
class(RNA_Gene_ID)
RNA_Gene_ID <- gsub("gene_id=", "", RNA_Gene_ID$Geneid)
class(RNA_Gene_ID)
counts_min_outliers <- counts_min_outliers[!counts_min_outliers$Geneid %in% RNA_Gene_ID]

mat <- as.matrix(counts_min_outliers, rownames= 'Geneid')
#for batch correction, rename mat to raw_counts
raw_counts <- mat

# Sanity check that sample sheet and count matrix are in the same order
stopifnot(identical(ss_min_outliers$library_id, colnames(raw_counts)))
library(sva)
adj_counts <- ComBat_seq(raw_counts, batch= ss_min_outliers$Batch, group= ss_min_outliers$Time)

#raw_counts is the matrix of counts with rRNA and failed libraries removed;
#ss_min_outliers is the sample sheet with column batch as prepared above and Time as factor as before.
#adj_counts is the adjusted matrix of counts that can be used for everything else
#from now on.

design <- model.matrix(~0 + ss_min_outliers$group)
colnames(design) <- sub('ss_min_outliers$group', '', colnames(design), fixed= TRUE)

y <- DGEList(counts= adj_counts,
               group= ss_min_outliers$group)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design)
y$samples

# MDS plot before and after batch correction

mds_min_outliers <- plotMDS(y, plot=FALSE)
mdsout_min_outliers <- as.data.table(mds_min_outliers$cmdscale.out)
mdsout_min_outliers[, library_id := colnames(y)]

# Add library characteristics from sample sheet
mdsout_min_outliers <- merge(mdsout_min_outliers, ss_min_outliers, by= "library_id")

## Include time in the labels - changing library_id to include time
mdsout_min_outliers$library_id <- mdsout_min_outliers [, gsub("Ap20-GFP-SFC-4h-R[0-9]", "Ap20-4h", mdsout_min_outliers$library_id)]
mdsout_min_outliers$library_id <- mdsout_min_outliers [, gsub("Ap20-GFP-SFC-16h-R[0-9]", "Ap20-16h", mdsout_min_outliers$library_id)]
mdsout_min_outliers$library_id <- mdsout_min_outliers [, gsub("Ap20-GFP-SFC-8h-R[0-9]", "Ap20-8h", mdsout_min_outliers$library_id)]
mdsout_min_outliers$library_id <- mdsout_min_outliers [, gsub("GFPcon-4hr-R1", "GFPcon-4h", mdsout_min_outliers$library_id)]
mdsout_min_outliers$library_id <- mdsout_min_outliers [, gsub("GFPcon-16hr-R[0-9]", "GFPcon-16h", mdsout_min_outliers$library_id)]
mdsout_min_outliers$library_id <- mdsout_min_outliers [, gsub("GFPcon-8hr-R1_S9", "GFPcon-8h", mdsout_min_outliers$library_id)]
mdsout_min_outliers$library_id <- mdsout_min_outliers [, gsub("GFPcon-8hr-R2", "GFPcon-8h", mdsout_min_outliers$library_id)]
mdsout_min_outliers$library_id <- mdsout_min_outliers [, gsub("RM-9_S416", "RM-9-24h", mdsout_min_outliers$library_id)]
mdsout_min_outliers$library_id <- mdsout_min_outliers [, gsub("RM-8_S415", "RM-8-24h", mdsout_min_outliers$library_id)]
mdsout_min_outliers$library_id <- mdsout_min_outliers [, gsub("RM-7_S414", "RM-7-24h", mdsout_min_outliers$library_id)]
mdsout_min_outliers$library_id <- mdsout_min_outliers [, gsub("RM-6_S413", "RM-6-12h", mdsout_min_outliers$library_id)]
mdsout_min_outliers$library_id <- mdsout_min_outliers [, gsub("RM-5_S412", "RM-5-24h", mdsout_min_outliers$library_id)]
mdsout_min_outliers$library_id <- mdsout_min_outliers [, gsub("RM-4_S411", "RM-4-6h", mdsout_min_outliers$library_id)]
mdsout_min_outliers$library_id <- mdsout_min_outliers [, gsub("RM-3_S410", "RM-3-4h", mdsout_min_outliers$library_id)]
mdsout_min_outliers$library_id <- mdsout_min_outliers [, gsub("RM-2_S409", "RM-2-2h", mdsout_min_outliers$library_id)]
mdsout_min_outliers$library_id <- mdsout_min_outliers [, gsub("RM-1_S408", "RM-1-12h", mdsout_min_outliers$library_id)]
mdsout_min_outliers$library_id <- mdsout_min_outliers [, gsub("F-.*", "RFP-0h", mdsout_min_outliers$library_id)]
mdsout_min_outliers$library_id <- mdsout_min_outliers [, gsub("F[1-3].*", "RFP-0h", mdsout_min_outliers$library_id)]

# Finally plot:
library(ggrepel)
ggplot(data= mdsout_min_outliers, aes(x= V1, y= V2, label= library_id, colour= Time, shape= Strain)) +
  geom_point() +
  geom_text_repel(size = 3) +
  ggtitle("MDS plot without outliers following rRNA removal - no batch correction") +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  theme_light()
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/')
ggsave('ggplot_min_outliers_min_rRNA_time0_dim1_dim2_nobatchcorrection.png', width= 20, height= 20, units= 'cm')

  #for dim2 vs dim3
	mds3 <- plotMDS(y, ndim = 3, plot=FALSE) #and then the same as above

	#removing the PIR and fam- genes (to see if female gametocyte data clustered better)
	#this could be performed assuming that a table containing the Geneid and description had been created prior - see below
	drop_ids <- gene_id_descr_table[grepl('PIR | fam-', description)]$Geneid
	counts <- counts[!counts$Geneid %in% drop_ids]

	#Creating Geneid vs description table
	gene <- GFF[V3 == "gene", list(attributes = V9)]
	gene_id_only <- gene[, gene_id := sub(".*;gene_id=", "", attributes)]
	gene_id_descr <- gene_id_only[, description := sub(".*;description=", "", attributes)]
	gene_id_descr <- gene_id_descr[, description := sub(";.*", "", description)]
	gene_id_descr_table <- data.table(gene_id_descr$gene_id, gene_id_descr$description)
	colnames(gene_id_descr_table) <- c("Geneid", "description")

	#decoding description column
  class(gene_id_descr_table)
  gene_id_descr_table[, description := sapply(description, URLdecode)]

	#saving table
	setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/')
	gz <- gzfile('geneid_description_table', 'w')
	write.table(gene_id_descr_table, gz, sep= '\t', row.names= FALSE, quote= FALSE)
	close(gz)

#Created MA plot (R script)
#logCPM = average log counts per million across the libraries
#LogCPM are the log counts per million, which can be understood as measuring expression
#level. LogFC is the log fold-change, which is the log difference between your groups.
#log-fold-change is the average (root-mean-square) of the
#largest absolute log-fold-changes between each pair of samples
#transformation and normalisation...logCPM?

#The glm approach to multiple groups is similar to the classic approach,
#but permits more general comparisons to be made
design_glm <- model.matrix(~0+group, data=y$samples)
#colnames(design_glm) <- levels(y$samples$group)
design_glm
#Here, the0+in the model formula is an instruction
#not to include an intercept column and instead to include a column for each group

# These must be the same as those for ATAC
contrasts <- makeContrasts("h24vs16"= group24 - group16,
                           "h24vs12"= group24 - group12,
                           "h16vs12"= group16 - group12,
                           "h16vs8"= group16 - group08,
                           "h16vs4"= group16 - group04,
                           "h12vs8"= group12 - group08,
                           "h8vs4"= group08 - group04,
                           "h4vs0"= group04 - group00,
                           levels= make.names(colnames(design_glm)))

class(design_glm)
stopifnot(all(abs(colSums(contrasts)) < 1e-6))
#all = given a set of logical vectors, are all the values true?
#abs = computes the absolute value of x
#form column sums for numeric arrays or dataframes

fit <- glmFit(y, design_glm, prior.count= 1)
#glmFit = Fit a negative binomial generalized log-linear model to the read counts
#for each gene.

dge <- list()
#Functions to construct, coerce and check for both kinds of R lists.
for(cntr in colnames(contrasts)){{
  print(cntr)
  lfc <- glmTreat(fit, contrast= contrasts[, cntr], lfc= log2(1.5))
  #numeric scalar specifying the absolute value
  #of the log2-fold change threshold
  #above which differential expression is to be considered.
  detable <- topTags(lfc, n= nrow(y))$table
  #n = max number of genes/tags to return
  detable$gene_id <- row.names(detable)
  detable <- data.table(detable)
  detable[, contrast := cntr]
  dge[[length(dge)+1]] <- detable
}}
dge <- rbindlist(dge)
class(dge)
dge[, unshrunk.logFC := NULL]

#https://www.reneshbedre.com/blog/ma.html - ma plots
# Normalisation - https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html

#merge dge table with the gene description table so one can see on the same line
#the gene_id, the difference between time points and the gene description
dge_descr <- merge(dge, gene_id_descr_table, by.x= 'gene_id', by.y= 'Geneid', all.x= TRUE, sort= FALSE)
dge_descr
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/')
gz <- gzfile('dge_descr_table_time0_additionalcontrasts_batchcorrection_9April', 'w')
write.table(dge_descr, gz, sep= '\t', row.names= FALSE, quote= FALSE)
close(gz)

#Finally plot MA
#plot shows DE vs non-DE genes with log2 FC against normalised mean count
 #the higher the mean count, the lower the log fold change required for the gene to be DE
 #for genes with low mean count, need a higher log FC to be DE

dge_descr
is.data.table(dge_descr) == TRUE

xord <- c('h4vs0', 'h8vs4', 'h12vs8', 'h16vs4', 'h16vs8', 'h16vs12', 'h24vs12', 'h24vs16')
dge_descr[, contrast_order := factor(contrast, xord)]

nsig <- dge_descr[, list(n_up= sum(FDR < 0.01 & logFC > 0), n_down= sum(FDR < 0.01 & logFC < 0)), contrast_order]
nsig[, n_up:= sprintf('Up = %s', n_up)]
nsig[, n_down:= sprintf('Down = %s', n_down)]

library(ggplot2)
gg <- ggplot(data= dge_descr, aes(x= logCPM, y= logFC)) +
  geom_point(alpha= 0.5, pch= '.') +
  geom_point(data= dge_descr[FDR < 0.01], alpha= 0.5, colour= 'red', pch= '.') +
  #to avoid overplotting, geom_smooth() and alpha (make points transparent)
  geom_smooth(se= FALSE, col= 'grey60', size= 0.1) +
  geom_hline(yintercept= 0, colour= 'blue') +
  geom_text(data= nsig, x= Inf, y= Inf, aes(label= n_up), vjust= 1.3, hjust= 1.1, size= 3) +
  geom_text(data= nsig, x= Inf, y= -Inf, aes(label= n_down), vjust= -1.2, hjust= 1.1, size= 3) +
  facet_wrap(~contrast_order, ncol= 2) +
  #free_y - individual axes
  ggtitle("MA plot without batch correction") +
  theme_light() +
  theme(strip.text= element_text(colour= 'black'))
gg

ggsave('MAplot_order_extracontrasts_nobatchcorrection_9Apr.png', width= 16, height= 20, units= 'cm')
#http://manual.omicsbox.biobam.com/user-manual/module-transcriptomics/pairwise-differential-expression-analysis/#PairwiseDifferentialExpressionAnalysis-MAPlot


#Created barplot of library sizes - prior to normalisation & sorted by time

	sizes <- data.frame(LibraryID = ss$library_id, LibrarySize = y$samples$lib.size, Strain = ss$Strain,
                    Time = ss$Time, group = ss$group)

	setDT(sizes)
	#Create a label to print and plot

sizes[, label := paste(LibraryID, Time)]
	#must be unique and informative (no repetition)

print(sizes)
xord <- sizes[order(-Time, LibrarySize)]$label #want to order label by time and then size
print(xord) # The order we want

	# New order of levels for the variable label
sizes[, label := factor(label, xord)]

library(ggplot2)

	#order by time
ggplot(sizes, aes(label, LibrarySize/1000000)) +
  geom_bar(stat="identity", aes(fill=Time)) +
  ggtitle("Barplot showing library size before normalisation") +
  xlab("Library name and time") +
  ylab("Library size (x 1 million)") +
  labs(fill = "Time(h)") +
  coord_flip() +
  theme_bw()

setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/')
ggsave('barplot_libsizes_before norm_time.pdf', width = 10, height = 10)

#idxstats workflow - make data analysis easier (see Snakefile) - did not work as expected

#Volcano plot
#displays significance on the y-axis and fold-change on the x-axis. In this case we use the log2 fold change (logFC) on the x-axis, and on the y-axis we’ll use -log10(FDR). This -log10 transformation
#is commonly used for p-values as it means that more significant genes
#have a higher scale.
dge_volcano <- dge_descr
#the significantly differentially expressed genes are the ones on the upper right and left hand side
#add a column to the dataframe to specify if if they are UP- or DOWN- regulated

#add a column of NA's
dge_volcano$diffexpressed <- "NO"

#if logfold change > 0.58 and p-value less than 0.05, set as "UP"
dge_volcano$diffexpressed[dge_volcano$logFC > 0.58 & dge_volcano$PValue < 0.05] <- "UP"
dge_volcano$diffexpressed[dge_volcano$logFC < -0.58 & dge_volcano$PValue < 0.05] <- "DOWN"


##to automate a bit
mycolours <- c("blue", "red", "black")
names(mycolours) <- c("DOWN", "UP", "NO")

### Names of genes besides the points
dge_volcano$delabel <- NA
dge_volcano$delabel[dge_volcano$diffexpressed != "NO"] <- dge_volcano$gene_id[dge_volcano$diffexpressed != "NO"]

library(ggplot2)
library(ggrepel)
volcano_plot <- ggplot(dge_volcano, mapping = aes(x= logFC, y= -log10(PValue), col = diffexpressed, label=delabel)) +
  geom_point() +
  ylab("-log10 (P-value)") +
  xlab("log2 Fold change") +
  ggtitle("Volcano plot with batch correction") +
  geom_vline(xintercept=c(-0.58,0.58), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  facet_wrap(~contrast_order, ncol= 2) +
  scale_colour_manual(values=mycolours)+ ##change point colour
  theme_minimal()
  setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/')
  ggsave('Volcano_plot_order_batchcorrection_9Apr.png', width= 20, height= 25, units= 'cm')

#Create a workflow plot
snakemake [...] / #assuming that [...] are the same commands used before
--forceall --dag > dag.dot # dag.dot refers to the plotting instructions & you can open in  a text editor to #see. dot.pdf is the actual plot.

 *then run on the command line the following
dot -Tpdf dag.dot > dag.pdf
 #OR
dot -Tpng dag.dot > dag.png

#Create a global expression plot

global <- data.table(dge_descr$gene_id, dge_descr$contrast_order, dge_descr$logFC)
global <- rename(global, gene_id = V1, contrast = V2, logFC = V3)
qq <- quantile(global$logFC, p= 0.995)
global[, col := ifelse(abs(global$logFC) > qq, qq, abs(logFC))]

p <- c(0.005, 0.05, 0.95, 0.99)
qbar <- global[, list(qq= quantile(logFC, p= p), p= p), by= contrast]
library(ggbeeswarm)
gg <- ggplot(data= global, aes(x= contrast, y= logFC, colour= col)) +
  geom_line(data= qbar, aes(y= qq, group= p), colour= 'dodgerblue', alpha= 0.5) +
  ggtitle ("Global gene expression plot with batch correction") +
  geom_quasirandom(width= 0.2, size= 0.25, bandwidth= 0.01) +
  scale_colour_gradient(low= 'grey80', high= 'black', guide= FALSE) +
  xlab('') +
  theme_light()
gg
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/')
ggsave('globalexpression_plot_batchcorrection_9Apr.png', width= length(unique(global$contrast)) * 2, height= 10, units= 'cm')
