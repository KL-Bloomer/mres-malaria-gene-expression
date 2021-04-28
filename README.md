### Investigating gene expression changes at the initial stages of _Plasmodium berghei_ mosquito midgut development using an RNA sequencing pipeline - Master of Research in Biomedical Sciences project at the University of Glasgow

#This project involved exploratory analysis of RNA-seq data from the mosquito midgut stages of the _Plasmodium berghei_ parasite. The aim was to determine the changes in gene expression profiles during the development of an immature, female gametocyte into a mature ookinete.  It was hypothesised that there would be gradual changes in gene expression and that the expression of genes involved in ookinete motility and midgut invasion would increase with time. These data were generated by Dr Katarzyna Modrynska's laboratory by culturing three different _P.berghei_ strains (RM, AP2-O, GFPcon) _in vitro_ or _in vivo_ and isolating RNA at 7 time points during mosquito midgut development - 2, 4, 6, 8, 12, 16 and 24hrs.

[Learning about vim](https://vimhelp.org/usr_02.txt.html#02.7)

[Create & clone a git repository](https://www.melbournebioinformatics.org.au/tutorials/tutorials/using_git/Using_Git/)

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

[Conda & Bioconda installation](https://bioconda.github.io/user/install.html)
#Then create conda environment and install the programmes on the file
bash
conda create --yes --name mres-malaria-gene-expression
conda activate mres-malaria-gene-expression
mamba install --freeze-installed -n MRes_project_2021 --yes --file requirements.txt
* Note to self: when installing using mamba from the cluster, need to install in the base environment which will translate to the conda environment

realpath $CONDA_PREFIX
* where the current activated environment lives

[Download, install & connect to VPN Client](https://www.gla.ac.uk/myglasgow/it/vpn/#downloadvpnclientforyourdevice;https://www.cisco.com/c/en/us/support/docs/smb/routers/cisco-rv-series-small-business-routers/Kmgmt-785-AnyConnect-Linux-Ubuntu.html)
* From the user interface, for connect, type: gucsasa1.cent.gla.ac.uk, followed by username (2592613b) & password

#To make logging into headnode easier
atom ~/.bashrc #opens a cryptic file with default settings for the shell
alias cluster="ssh -oKexAlgorithms=+diffie-hellman-group1-sha1 2592613b@headnode03.cent.gla.ac.uk"
* Therefore, just need to open terminal and type "cluster", followed by:
ssh bioinf03.hpc.gla.ac.uk

* Installed conda and bioconda as on PC - environment = mres-malaria-gene-expression.

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
* library_id = library name i.e. Ap20-GFP-SFC-16h-R1 - name with only letters, numbers, dots, underscores, hyphens
* Time
* Strain
* R1 Full path to read 1 fastq file
* R2 Full path to read 2 fastq file

* Prepared with Excel - save as plain text with columns separated by TAB (__not comma or space__).

#To visualise if there are any spaces; where tabs are replaced by ^I; end of each line = $
cat -vet export/home/2592613b/mres-malaria-gene-expression/sample_sheet.tsv
sed -i 's/ //g' ~/mres-malaria-gene-expression/sample_sheet.tsv
* the -i option means it will do the replacement directly on the file;
* "s" is substitution and "/" refers to delimiters...i.e. sed 's /unix/linux' - unix is search pattern & linux is replacement

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

##Edit & re-run the snakefile - download the M.Musculus genome and concatenate with P.berghei
* Evidence of contamination due to a high GC content and bimodal graph from the fastqc/multiqc analyses. Therefore, it was decided to blast the unmapped reads and to re-run the mapping with a M.musulus-P.berghei concatenated genome. This would eliminate reads that map preferentially to the M/musculus genome from the downstream analysis (i.e. read counts).
* Note to self: __No tabs in the snakefile__, only __whitespace__.
* Changed the shell command: cat {input.pb} {output.fa} > {output.fa} __to__ cat {input.pb} {input.mm} > {output.fa}

##Samtools

samtools idxstats <my.bam>
* number of reads mapped to mouse or plasmodium genome; output columns: __chromosome name__, __size__, __number of alignments/mapped reads__, __reads unmapped but a mapped mate__

samtools view <my.bam> | less -S
* open a bam file

samtools fasta <my.bam>
* convert bam to fasta

[explain flags](https://broadinstitute.github.io/picard/explain-flags.html)

* note: decided to run the over-represented sequences in blast which matched to M.musculus

#Ran the unmapped reads from the files in blast (as recorded in Snakefile)
* Corrupt download when trying to download the blast database, therefore:
  * remove leftover from download (this was not necessary as my download seemed to have failed completely):
    rm -r /export/home/2592613b/mres-malaria-gene-expression/output/blastdb
  * [Symlink](https://linuxize.com/post/how-to-create-symbolic-links-in-linux-using-the-ln-command/)
  * Could use cp to make hard copies but this will prevent the need to duplicate large files
    ln -s /export/home/db291g/Tritume/blastdb /export/home/2592613b/mres-malaria-gene-expression/output/
  * And ultimately, read the blast database from Dario's home directory.

#Working with snakemake
* Note to self: can execute a snakemake wrapper script in the background using nohup (very handy!)
nohup ./snakefile_run &
    * Therefore, even if I am logged out of the cluster, snakemake will still be going in the background;
    * To see the output from snakemake = less nohup.out;
    * to kill snakemake = pkill snakemake
    __OR__
    * immediate kill(emergency) = pkill -9 snakemake

* Issues being logged out:
snakemake --printshellcmds --jobs 1 \
    --config ss=$PWD/sample_sheet.tsv \
    --directory output \
    --snakefile Snakefile \
    --dry-run \
    --unlock

* Snakemake logs under .snakemake in the working directory:
ls -a
/export/home/2592613b/mres-malaria-gene-expression/output/.snakemake/log/

* To get around the "AmbiguousRuleException" error when running rule summarise_contamination, performed the following:
1) output:
    summary= 'blast_species/{library_id}.species.tsv',
		#and made corrections to rule final_output
2) wildcard_constraints:
    library_id = '|'.join([re.escape(x) for x in sample_sheet['library_id']])
		#add this before the rule final_output - gets rid of "greedy behaviour"

##Re-ran the snakemake file on the female gametocyte data
* included female gametocytes (Strain = RFP & Time = 0h)
  * [SRR526055](https://www.ebi.ac.uk/ena/browser/view/PRJNA374918)
  * [21725](https://www.ebi.ac.uk/ena/browser/view/PRJNA714084)
* edited sample sheet to include female gametocyte data
* If you don't want blast to run on female gametocyte data, edit rule final_output:
expand('blast_species/{library_id}.species.tsv', library_id= sample_sheet[sample_sheet['Time'] != 0]['library_id'])
* Prioritise the gene count output:
snakemake --prioritize count_reads_in_genes --jobs 1 ...

#MDS plot - code can be found in Script_for_RNAseq_plots.R.
* Removed outliers and rRNA contamination and PIR/fam genes (R script)
* removal of the PIR/fam genes did not affect the MDS plot/MA plot, therefore kept those genes in.
* Appears to be a batch effect with the RM samples clustering away from the other samples (i.e. GFPcon & Ap20 are confounded by strain). Therefore corrected for batch effect in all subsequent plots.

#Created MA plot - code can be found in Script_for_RNAseq_plots.R
* Note to self: logCPM = average log counts per million across the libraries
* LogCPM are the log counts per million, which can be understood as measuring expression level.
* LogFC __(fold change)__ is the log fold-change, which is the log difference between your groups.
* logFC = average (root-mean-square) of the largest absolute log-fold-changes between each pair of samples
* transformation and normalisation = __logCPM__ - normalise for library size/sequencing depth

[ma plots](https://www.reneshbedre.com/blog/ma.html)
[ma plots](http://manual.omicsbox.biobam.com/user-manual/module-transcriptomics/pairwise-differential-expression-analysis/#PairwiseDifferentialExpressionAnalysis-MAPlot)
[normalisation](https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html)

#Barplot of library sizes - prior to normalisation, sorted by time
#Code can be found in Script_for_barplot.R

#idxstats workflow - make data analysis easier (see Snakefile) - did not work as expected

#Volcano plot & global expression plot - see code in Script_for_RNAseq_plots.R

#Create a workflow plot
snakemake [...] /
--forceall --dag > dag.dot
* assuming that [...] are the same commands as used before
* dag.dot refers to the plotting instructions & you can open in a text editor.
* dot.pdf is the actual plot.
* then run on the command line the following
dot -Tpdf dag.dot > dag.pdf
__OR__
dot -Tpng dag.dot > dag.png

#Clustering and "interesting gene" plots - see code in clustering.R

#Heatmaps - code used in heatmap.R

#Plotting genomic data using Gviz (alternative to IGV) - see code in Gviz.R
* Installation: added most recent package version to requirements file and installed with mamba
