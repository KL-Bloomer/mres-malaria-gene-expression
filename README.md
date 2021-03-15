### mres-malaria-gene-expression

#Folder structure on the cluster
/export/home/2592613b/mres-malaria-gene-expression
	miniconda3
	data
		files.fastq.gz
		file_key.txt
	output
		fastqc
		ref
		multiqc
		featureCounts
		hisat2
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
				KLB_diary
				sample_sheet.tsv
				scripts
					addGeneIdToGff.py
		output
			mres-malaria-gene-expression
				hisat
				ref
				fastqc

#Establishing vpn connection

#Created & cloned git repository
https://www.melbournebioinformatics.org.au/tutorials/tutorials/using_git/Using_Git/ 
git status
git add
git commit -m "message"
git push
git pull

#Conda & Bioconda installation
https://bioconda.github.io/user/install.html

#Then create conda environment and install the programmes on the file
bash
conda create --yes --name mres-malaria-gene-expression
conda activate mres-malaria-gene-expression 
mamba install --freeze-installed -n mres-malaria-gene-expression --yes --file requirements.txt

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

#To avoid this issue "Unable to negotiate with 130.209.19.28 port 22: no matching key exchange method found...".
atom ~/.ssh/config

# In the atom textfile
Host headnode03.cent.gla.ac.uk
    KexAlgorithms +diffie-hellman-group1-sha1 #save and close

#How to copy requirements file from PC to cluster
scp -r /path/to/requirements.txt 2592613b@headnode03.cent.gla.ac.uk:~/mres-malaria-gene-expression/

#How to copy from cluster to PC - from PC terminal, type the following:
scp -r 2592613b@headnode03.cent.gla.ac.uk:~/mres-malaria-gene-expression/path/to/file /path/to/file/on/local/computer

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
*R1 Full path to read 1 fastq file
*R2 Full path to read 2 fastq file

Prepared with Excel - save as plain text with columns separated by TAB
by TAB (not comma or space).

#To visualise if there are any spaces; where tabs are replaced by ^I; end of each line = $
cat -vet export/home/2592613b/mres-malaria-gene-expression/sample_sheet.tsv
sed -i 's/ //g' ~/mres-malaria-gene-expression/sample_sheet.tsv #the -i option means it will do the replacement directly on the file; "s" is substitution and "/" refers to delimiters...i.e. sed 's /unix/linux' - unix is search pattern & linux is replacement

#Create snakefile_run; remove --dry-run when you want to run the commands (otherwise only shows what would be executed)
snakemake --dry-run --printshellcmds --jobs 1 \
--config ss=$PWD/sample_sheet.tsv \
--directory output \
--snakefile Snakefile
#For future reference, in the Snakefile, run hisat2 command option "--threads 16" instead of "--threads 1" & for featureCounts command "-T 8" instead of "-T 1".

#IGV
Transferred hisat2 files to PC; installed IGV =2.9.2
mamba install --freeze-installed -n MRes_project_2021 --yes --file requirements.txt
open IGV using the command "igv"
