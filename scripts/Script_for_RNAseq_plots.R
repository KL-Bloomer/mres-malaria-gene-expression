
### Plots and tables for quality control and data visualisation where rRNA was removed 
#subsequently, outliers were removed and a batch correction was performed for downstream analyses
#These include: MA, volcano and global expression plot, the DGE and gene_id-description table

#Load required packages
library(data.table)
library(edgeR)
library(sva)
library(tidyr)
library(ggrepel)
library(ggplot2)
library(ggbeeswarm)


ss <- snakemake@input[['sample_sheet']]
counts_file <- snakemake@input[['counts']]
GFF_file <- snakemake@input[['gff']]
MDSplots_concatenated <- snakemake@output[['MDSplots_concatenated']]
dge_table_file <- snakemake@output[['dge_table']]
geneid_desc_table <- snakemake@output[['geneid_desc_table']]
logrpkm_table_file <- snakemake@output[['logrpkm_table']]
MAplot <- snakemake@output[['MA_plot']]
#Volcano_plot <- snakemake@output[['Volcano_plot']]
#globalexpression <- snakemake@output[['globalexpression']]

#read the sample sheet
ss_min_outliers <- fread(ss)
ss_min_outliers <- ss_min_outliers[Outliers == FALSE,]
ss_min_outliers <- ss_min_outliers[library_type == "RNA-seq",]
ss_min_outliers <- ss_min_outliers[!Stage == "Asex",]
ss_min_outliers[, Time := sprintf('%.2d', Time)]
ss_min_outliers[, group := paste(Time)]
ss_min_outliers$Time <- as.factor(ss_min_outliers$Time)

#read the counts table
counts <- fread(cmd= paste('grep -v "^#"', counts_file))
setnames(counts, names(counts), sub('.bam', '', basename(names(counts))))
counts <- counts[, c('Geneid', ss_min_outliers$library_id), with= FALSE]

#Removing the rRNA contamination
#read the annotation file
GFF <- fread(cmd=paste('grep -v "^#"', GFF_file))
rRNA_GFF <- GFF[V3 == "rRNA"]
rRNA_ID <- subset(rRNA_GFF, select= V9) 
RNA_ID_sep <- separate(data = rRNA_ID, col = V9, into = c("Geneid_Feature", "Parent", "Description", "Geneid"), sep = "([;])")
RNA_Gene_ID <- subset(RNA_ID_sep, select= Geneid)
class(RNA_Gene_ID)
RNA_Gene_ID <- gsub("gene_id=", "", RNA_Gene_ID$Geneid)
counts <- counts[!counts$Geneid %in% RNA_Gene_ID]

#create a matrix of the counts
mat <- as.matrix(counts, rownames= 'Geneid')
mat_log <- log(mat+1)
design <- model.matrix(~0 + ss_min_outliers$group)
colnames(design) <- sub('ss_min_outliers$group', '', colnames(design), fixed= TRUE)
y <- DGEList(counts= mat_log,
             group= ss_min_outliers$group)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design)
y$samples

# MDS plot no batch correction

mds_min_outliers_nobatch <- plotMDS(y, plot=FALSE)
mdsout_min_outliers_nobatch <- as.data.table(mds_min_outliers_nobatch$cmdscale.out)
mdsout_min_outliers_nobatch[, library_id := colnames(y)]
mdsout_min_outliers_nobatch[, sample_set := 'No batch correction']

### MDS no outliers and with batch corretion

# Another sanity check that Time is a factor variable - maybe redundant
stopifnot(is.factor(ss_min_outliers$Time))

#for batch correction, rename mat to raw_counts
raw_counts <- mat_log

# Sanity check that sample sheet and count matrix are in the same order
#issue with this command...
stopifnot(identical(ss_min_outliers$library_id, colnames(raw_counts)))

#load sva for batch correction
adj_counts <- ComBat_seq(raw_counts, batch= ss_min_outliers$Batch, group= ss_min_outliers$Time)

#raw_counts is the matrix of counts with rRNA and failed libraries removed; 
#ss_min_outliers is the sample sheet with column batch as prepared above and Time as factor as before. 
#adj_counts is the adjusted matrix of counts that can be used for everything else 
#from now on.

design <- model.matrix(~0 + ss_min_outliers$group)
colnames(design) <- sub('ss_min_outliers$group', '', colnames(design), fixed= TRUE)
y <- DGEList(counts= adj_counts, #mat
             group= ss_min_outliers$group)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design)
y$samples

# MDS plot before and after batch correction

mds_min_outliers_batch <- plotMDS(y, plot=FALSE)
mdsout_min_outliers_batch <- as.data.table(mds_min_outliers_batch$cmdscale.out)
mdsout_min_outliers_batch[, library_id := colnames(y)]
mdsout_min_outliers_batch[, sample_set := 'Batch correction']

#Now concatenate the two MDS datasets and add the information from the sample sheet:

mdsdata <- rbind(mdsout_min_outliers_batch, mdsout_min_outliers_nobatch)

# Add library characteristics from sample sheet
mdsdata <- merge(mdsdata, ss_min_outliers, by= "library_id")

gg <- ggplot(data= mdsdata, aes(x= V1, y= V2, label= Time, colour= Time, shape= Strain)) +
  geom_point(size = 3) +
  geom_text_repel(size = 4, segment.colour= "grey80") +
  labs( colour = "Time (h)", shape = "Parasite lines") +
  ggtitle("MDS plots") +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  theme_light() +
  facet_wrap(~sample_set, scales = "free")
gg + theme(axis.text=element_text(size=12),
           axis.title=element_text(size=14,face="bold"),
           strip.text = element_text(size = 16, face="bold"))

# Choose appropriate width and height
ggsave(MDSplots_concatenated, width = 40, height = 15, units = "cm")

# Differential expression
# =======================

#create a glm
design_glm <- model.matrix(~0+group, data=y$samples)

# These must be the same as those for ATAC
contrasts <- makeContrasts("h24vs16"= group24 - group16,
                           "h16vs12"= group16 - group12,
                           "h12vs8"= group12 - group08,
                           "h8vs6"= group08 - group06,
                           "h6vs4"= group08 - group04,
                           "h4vs2"= group04 - group02,
                           "h2vs0"= group02 - group00,
                           #"h0vs72"= group00 - group72,
                           levels= make.names(colnames(design_glm)))

stopifnot(all(abs(colSums(contrasts)) < 1e-6))
fit <- glmFit(y, design_glm, prior.count= 1)

#Create the dge table
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
class(dge)
dge[, unshrunk.logFC := NULL]

#merge dge table with the gene description table so one can see on the same line
#the gene_id, the difference between time points and the gene description

### Create the Gene-ID vs Description table

gene <- GFF[V3 == "gene", list(attributes = V9)]
gene_id_only <- gene[, gene_id := sub(".*;gene_id=", "", attributes)]

gene_id_descr <- gene_id_only[, description := sub(".*;description=", "", attributes)]
gene_id_descr <- gene_id_descr[, description := sub(";.*", "", description)]
gene_id_descr_table <- data.table(gene_id_descr$gene_id, gene_id_descr$description)
colnames(gene_id_descr_table) <- c("gene_id", "description")

#decoding description column
gene_id_descr_table[, description := sapply(description, URLdecode)]
write.table(gene_id_descr_table, geneid_desc_table, sep= '\t', row.names= FALSE, quote= FALSE)

#now merge the decription table with the dge table
dge_table <- merge(dge, gene_id_descr_table, by.x= 'gene_id', by.y= 'gene_id', all.x= TRUE, sort= FALSE)
write.table(dge_table, dge_table_file, sep= '\t', row.names= FALSE, quote= FALSE)

##Creating an MA plot
#order the contrasts appropriately
xord <- c('h2vs0', 'h4vs2', 'h6vs4', 'h8vs6', 'h12vs8', 'h16vs12', 'h24vs16')
dge_table[, contrast_order := factor(contrast, xord)]
nsig <- dge_table[, list(n_up= sum(FDR < 0.01 & logFC > 0), n_down= sum(FDR < 0.01 & logFC < 0)), contrast_order]
nsig[, n_up:= sprintf('Up = %s', n_up)]
nsig[, n_down:= sprintf('Down = %s', n_down)]

ggplot(data= dge_table, aes(x= logCPM, y= logFC)) +
  geom_point(alpha= 0.5, pch= '.') +
  geom_point(data= dge_table[FDR < 0.01], alpha= 0.5, colour= 'red', pch= '.') +
  geom_smooth(se= FALSE, col= 'grey60', size= 0.1) +
  geom_hline(yintercept= 0, colour= 'blue') +
  geom_text(data= nsig, x= Inf, y= Inf, aes(label= n_up), vjust= 1.3, hjust= 1.1, size= 3) +
  geom_text(data= nsig, x= Inf, y= -Inf, aes(label= n_down), vjust= -1.2, hjust= 1.1, size= 3) +
  facet_wrap(~contrast_order, ncol= 2) +
  ggtitle("MA plot with batch correction") +
  theme_light() +
  theme(strip.text= element_text(colour= 'black'))
ggsave(MAplot, width= 12, height= 12, units = "cm")

### Generating logrpkm table in long format

#Need to re-read the counts table to include the column "length"
counts <- fread(cmd=paste('grep -v "^#"', counts_file))
setnames(counts, names(counts), sub('.bam', '', basename(names(counts))))
counts <- counts[, c('Length','Geneid', ss_min_outliers$library_id), with= FALSE]
counts <- counts [!counts$Geneid %in% RNA_Gene_ID]

gene_length <- data.table(counts$Geneid, counts$Length)
setnames(gene_length, "V1", "gene_id")
setnames(gene_length, "V2", "length")
gene_length <- gene_length[keep, , ]

logrpkm <- rpkm(y, gene.length = gene_length$length, log=TRUE)
logrpkm <- data.table(logrpkm, keep.rownames = "gene_id")

#Convert to long format
logrpkm_table <- melt(logrpkm, variable.name = "library_id", id.vars = "gene_id",
                           value.name = "logrpkm")
logrpkm_table <- data.table(logrpkm_table)
write.table(logrpkm_table, file=logrpkm_table_file, row.names = FALSE, sep= '\t', quote= FALSE)
