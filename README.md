### mres-malaria-gene-expression

#Folder structure on the cluster
/export/home/2592613b/mres-malaria-gene-expression
	miniconda3
	data
		{files}.fastq.gz
		file_key.txt
	output
		fastqc
		ref
		multiqc
		featureCounts
		hisat2
		R.output
	requirements.txt
	sample_sheet.tsv
	scripts
		addGeneIdToGff.py

#Folder structure on PC
/home/Kirstin/
	MRes_Malaria_2021
		data
			Sequencing
				A few fastq.gz files
				file_key.txt
		git_repos
			mres-malaria-gene-expression
				README.md
				requirements.txt
				Snakefile
				KLB_diary.txt
				sample_sheet.tsv
				scripts
					addGeneIdToGff.py
		output
			mres-malaria-gene-expression
				hisat
				ref
				fastqc
				R.output

#Learning about vim: https://vimhelp.org/usr_02.txt.html#02.7

#Created & cloned git repository
https://www.melbournebioinformatics.org.au/tutorials/tutorials/using_git/Using_Git/
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

#Accessing fastqc files and saving into "my-project" folder
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

####Ran the overrespresented sequences in blast which matched to M.musculus
####Ran the unmapped reads from the files in blast as per the pdf document, albeit these modifications:
rule final_output, expand('blast_species/{library_id}.species.tsv', library_id= sample_sheet['library_id']),

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

	rm(list=ls())
	library(data.table)
	library(edgeR)

	setwd('~/MRes_Malaria_2021/git_repos/mres-malaria-gene-expression/')
	ss <- fread('sample_sheet.tsv')

	setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/featureCounts/')

	counts <- fread(cmd= 'grep -v "^#" counts.tsv')
	setnames(counts, names(counts), sub('.bam', '', basename(names(counts))))
	counts <- counts[, c('Geneid', ss$library_id), with= FALSE]

	mat <- as.matrix(counts, rownames= 'Geneid')
	ss[, group := paste(Time, Strain, sep= '_')]

	ss$Time <- as.factor(ss$Time)
	class(ss$Time)

	design <- model.matrix(~0 + ss$group)
	colnames(design) <- sub('ss$group', '', colnames(design), fixed= TRUE)

	y <- DGEList(counts= mat, group= ss$group)
	keep <- filterByExpr(y)
	y <- y[keep, , keep.lib.sizes=FALSE]
	y <- calcNormFactors(y)

	y <- estimateDisp(y, design)
	y$samples

mds <- plotMDS(y, plot=FALSE)
mdsout <- as.data.table(mds$cmdscale.out)
mdsout[, library_id := colnames(y)]

	# Add library characteristics from sample sheet
mdsout <- merge(mdsout, ss, by= "library_id")
mdsout$Time <- as.factor(mdsout$Time)
class(mdsout$Time)  

mdsout
	## Include time in the labels - changing library_id to include time

mdsout$library_id <- mdsout [, gsub("Ap20-GFP-SFC-4h-R[0-9]", "Ap20-4h", mdsout$library_id)]
mdsout$library_id <- mdsout [, gsub("Ap20-GFP-SFC-16h-R[0-9]", "Ap20-16h", mdsout$library_id)]
mdsout$library_id <- mdsout [, gsub("Ap20-GFP-SFC-8h-R[0-9]", "Ap20-8h", mdsout$library_id)]
mdsout$library_id <- mdsout [, gsub("GFPcon-4hr-R[1-2]", "GFPcon-4h", mdsout$library_id)]
mdsout$library_id <- mdsout [, gsub("GFPcon-16hr-R[0-9]", "GFPcon-16h", mdsout$library_id)]
mdsout$library_id <- mdsout [, gsub("GFPcon-8hr-R1_S9", "GFPcon-8h", mdsout$library_id)]
mdsout$library_id <- mdsout [, gsub("GFPcon-8hr-R2", "GFPcon-8h", mdsout$library_id)]
mdsout$library_id <- mdsout [, gsub("RM-10_S417", "RM-10-24h", mdsout$library_id)]
mdsout$library_id <- mdsout [, gsub("RM-9_S416", "RM-9-24h", mdsout$library_id)]
mdsout$library_id <- mdsout [, gsub("RM-8_S415", "RM-8-24h", mdsout$library_id)]
mdsout$library_id <- mdsout [, gsub("RM-7_S414", "RM-7-24h", mdsout$library_id)]
mdsout$library_id <- mdsout [, gsub("RM-6_S413", "RM-6-12h", mdsout$library_id)]
mdsout$library_id <- mdsout [, gsub("RM-5_S412", "RM-5-24h", mdsout$library_id)]
mdsout$library_id <- mdsout [, gsub("RM-4_S411", "RM-4-6h", mdsout$library_id)]
mdsout$library_id <- mdsout [, gsub("RM-3_S410", "RM-3-4h", mdsout$library_id)]
mdsout$library_id <- mdsout [, gsub("RM-2_S409", "RM-2-2h", mdsout$library_id)]
mdsout$library_id <- mdsout [, gsub("RM-1_S408", "RM-1-12h", mdsout$library_id)]

library(ggrepel)
ggplot(data= mdsout, aes(x= V1, y= V2, label= library_id, colour= Time, shape= Strain)) +
  geom_point() +
  geom_text_repel(size = 3) +
  ggtitle("MDS plot with outliers prior to rRNA removal") +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  theme_light()
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/')
ggsave('ggplot_outliers_rRNA_dim1_dim2_time.pdf', width= 25, height= 25, units= 'cm')

	#for dim2 vs dim3
	mds3 <- plotMDS(y, ndim = 3, plot=FALSE) #and then the same as above

	ggplot(data= mdsout3, aes(x= V2, y= V3, label= library_id, colour= Time, shape= Strain)) +
  geom_point() +
  geom_text_repel(size = 3) +
  ggtitle("MDS plot with outliers prior to rRNA removal") +
  xlab("Leading logFC dim 2") +
  ylab("Leading logFC dim 3") +
  theme_light()
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/')
ggsave('ggplot_outliers_rRNA_dim2_dim3.pdf', width= 20, height= 20, units= 'cm')

	#removing the outliers
	#edited the sample sheet by removing the libraries that were outliers (RM10 & GFPcon-4h-R2)
	#the rest of the code is the same as above

	#removing the rRNA contamination
	library(dplyr)
	library(tidyr)
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/ref/')
GFF <- fread(cmd = 'grep -v "^#" PlasmoDB-49_PbergheiANKA.gff')
GFF

rRNA_GFF <- GFF %>%
  filter(V3 == "rRNA")

rRNA_ID <- subset(rRNA_GFF, select= V9)
RNA_ID_sep <- separate(data = rRNA_ID, col = V9, into = c("Geneid_Feature", "Parent", "Description", "Geneid"), sep = "([;])")
RNA_Gene_ID <- subset(RNA_ID_sep, select= Geneid)
RNA_Gene_ID <- gsub("gene_id=", "", RNA_Gene_ID$Geneid)
counts_min_outliers <- counts_min_outliers[!counts_min_outliers$Geneid %in% RNA_Gene_ID]

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

#Created MA plot (R script)
* Note to self - reassign "time" as a character vector which can sort accordingly
sprintf('%.2d', 0:11)
ss_min_outliers[, Time := sprintf('%.2d', Time)]
library(ggplot2)
gg <- ggplot(data= dge, aes(x= logCPM, y= logFC)) +
  geom_point(alpha= 0.5, pch= '.') +
  geom_point(data= dge[FDR < 0.01], alpha= 0.5, colour= 'red', pch= '.') +
  #to avoid overplotting, geom_smooth() and alpha (make points transparent)
  geom_smooth(se= FALSE, col= 'grey60', size= 0.1) +
  geom_hline(yintercept= 0, colour= 'blue') +
  geom_text(data= nsig, x= Inf, y= Inf, aes(label= n_up), vjust= 1.3, hjust= 1.1, size= 3) +
  geom_text(data= nsig, x= Inf, y= -Inf, aes(label= n_down), vjust= -1.2, hjust= 1.1, size= 3) +
  facet_wrap(~contrast, scales= 'free_y', ncol= 2) +
  theme_light() +
  theme(strip.text= element_text(colour= 'black'))
ggsave('MAplot_testing_30Mar.pdf', width= 16, height= 20, units= 'cm')

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
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))

setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/')
ggsave('barplot_libsizes_before norm_time.pdf', width = 10, height = 10)

#idxstats workflow - make data analysis easier
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
* edit rule final_output: 'idxstats/idxstats.tsv'

#Create a volcano plot
#Create a workflow plot
#Create a global expression plot
