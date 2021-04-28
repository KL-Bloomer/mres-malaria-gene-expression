### Plots and tables for data visualisation where the outliers were removed, 
#along with the rRNA, following batch correction
#These include: MA, volcano and global expression plot, the DGE and gene_id-description table

rm(list=ls())
library(data.table)
library(edgeR)

#read the sample sheet
setwd('~/MRes_Malaria_2021/git_repos/mres-malaria-gene-expression/')

ss_min_outliers <- fread('sample_sheet_minusoutliers.tsv')
ss_min_outliers[, Time := sprintf('%.2d', Time)]
ss_min_outliers[, group := paste(Time)]
ss_min_outliers$Time <- as.factor(ss_min_outliers$Time)

# Another sanity check that Time is a factor variable - maybe redundant
stopifnot(is.factor(ss_min_outliers$Time))

#read the counts table
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/featureCounts/')
counts <- fread(cmd= 'grep -v "^#" counts.tsv')
counts

class(counts)
setnames(counts, names(counts), sub('.bam', '', basename(names(counts))))
counts_min_outliers <- counts[, c('Geneid', ss_min_outliers$library_id), with= FALSE]

#Removing the rRNA contamination
#load required packages
library(dplyr)
library(tidyr)
#read the annotation file
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/ref/')
GFF <- fread(cmd = 'grep -v "^#" PlasmoDB-49_PbergheiANKA.gff')

rRNA_GFF <- GFF %>% 
  filter(V3 == "rRNA")
rRNA_ID <- subset(rRNA_GFF, select= V9) 
RNA_ID_sep <- separate(data = rRNA_ID, col = V9, into = c("Geneid_Feature", "Parent", "Description", "Geneid"), sep = "([;])")
RNA_Gene_ID <- subset(RNA_ID_sep, select= Geneid) 
class(RNA_Gene_ID)
RNA_Gene_ID <- gsub("gene_id=", "", RNA_Gene_ID$Geneid)
counts_min_outliers <- counts_min_outliers[!counts_min_outliers$Geneid %in% RNA_Gene_ID]

#create a matrix of the counts
mat <- as.matrix(counts_min_outliers, rownames= 'Geneid')

#for batch correction, rename mat to raw_counts
raw_counts <- mat

# Sanity check that sample sheet and count matrix are in the same order
#issue with this command...
stopifnot(identical(ss_min_outliers$library_id, colnames(raw_counts)))

#load sva for batch correction
library(sva)
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

# MDS plot - dim1 vs dim2:
library(ggrepel)
ggplot(data= mdsout_min_outliers, aes(x= V1, y= V2, label= library_id, colour= Time, shape= Strain)) +
  geom_point() +
  geom_text_repel(size = 3) +
  ggtitle("MDS plot without outliers following rRNA removal - no batch correction") +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  theme_light()

####Plotting in 3 dimensions now
mds3_min_outliers <- plotMDS(y, ndim = 3, plot=FALSE)
mdsout3_min_outliers <- as.data.table(mds3_min_outliers$cmdscale.out)
mdsout3_min_outliers[, library_id := colnames(y)]

# Add library characteristics from sample sheet
mdsout3_min_outliers <- merge(mdsout3_min_outliers, ss_min_outliers, by= "library_id")
mdsout3_min_outliers

## Include time in the labels - changing library_id to include time
mdsout3_min_outliers$library_id <- mdsout3_min_outliers [, gsub("Ap20-GFP-SFC-4h-R[0-9]", "Ap20-4h", mdsout3_min_outliers$library_id)]
mdsout3_min_outliers$library_id <- mdsout3_min_outliers [, gsub("Ap20-GFP-SFC-16h-R[0-9]", "Ap20-16h", mdsout3_min_outliers$library_id)]
mdsout3_min_outliers$library_id <- mdsout3_min_outliers [, gsub("Ap20-GFP-SFC-8h-R[0-9]", "Ap20-8h", mdsout3_min_outliers$library_id)]
mdsout3_min_outliers$library_id <- mdsout3_min_outliers [, gsub("GFPcon-4hr-R1", "GFPcon-4h", mdsout_min_outliers$library_id)]
mdsout3_min_outliers$library_id <- mdsout3_min_outliers [, gsub("GFPcon-16hr-R[0-9]", "GFPcon-16h", mdsout_min_outliers$library_id)]
mdsout3_min_outliers$library_id <- mdsout3_min_outliers [, gsub("GFPcon-8hr-R1_S9", "GFPcon-8h", mdsout_min_outliers$library_id)]
mdsout3_min_outliers$library_id <- mdsout3_min_outliers [, gsub("GFPcon-8hr-R2", "GFPcon-8h", mdsout_min_outliers$library_id)]
mdsout3_min_outliers$library_id <- mdsout3_min_outliers [, gsub("RM-9_S416", "RM-9-24h", mdsout3_min_outliers$library_id)]
mdsout3_min_outliers$library_id <- mdsout3_min_outliers [, gsub("RM-8_S415", "RM-8-24h", mdsout3_min_outliers$library_id)]
mdsout3_min_outliers$library_id <- mdsout3_min_outliers [, gsub("RM-7_S414", "RM-7-24h", mdsout3_min_outliers$library_id)]
mdsout3_min_outliers$library_id <- mdsout3_min_outliers [, gsub("RM-6_S413", "RM-6-12h", mdsout3_min_outliers$library_id)]
mdsout3_min_outliers$library_id <- mdsout3_min_outliers [, gsub("RM-5_S412", "RM-5-24h", mdsout3_min_outliers$library_id)]
mdsout3_min_outliers$library_id <- mdsout3_min_outliers [, gsub("RM-4_S411", "RM-4-6h", mdsout3_min_outliers$library_id)]
mdsout3_min_outliers$library_id <- mdsout3_min_outliers [, gsub("RM-3_S410", "RM-3-4h", mdsout3_min_outliers$library_id)]
mdsout3_min_outliers$library_id <- mdsout3_min_outliers [, gsub("RM-2_S409", "RM-2-2h", mdsout3_min_outliers$library_id)]
mdsout3_min_outliers$library_id <- mdsout3_min_outliers [, gsub("RM-1_S408", "RM-1-12h", mdsout3_min_outliers$library_id)]
mdsout3_min_outliers$library_id <- mdsout3_min_outliers [, gsub("F-.*", "RFP-0h", mdsout3_min_outliers$library_id)]
mdsout3_min_outliers$library_id <- mdsout3_min_outliers [, gsub("F[1-3].*", "RFP-0h", mdsout3_min_outliers$library_id)]

# MDS plot - dim2 vs dim3:
library(ggrepel)
ggplot(data= mdsout3_min_outliers, aes(x= V2, y= V3, label= library_id, colour= Time, shape= Strain)) +
  geom_point() +
  geom_text_repel(size = 3) +
  ggtitle("MDS plot without outliers following rRNA removal") +
  xlab("Leading logFC dim 2") +
  ylab("Leading logFC dim 3") +
  theme_light()
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/')
ggsave('ggplot_min_outliers_min_rRNA_time0_dim2_dim3.pdf', width= 20, height= 20, units= 'cm')

# Differential expression
# =======================

#create a glm
design_glm <- model.matrix(~0+group, data=y$samples)
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
#glmFit = Fit a negative binomial generalized log-linear model to the read counts 
#for each gene. 
fit <- glmFit(y, design_glm, prior.count= 1)

#Create the dge table
dge <- list()
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

#merge dge table with the gene description table so one can see on the same line 
#the gene_id, the difference between time points and the gene description

### Create the Gene-ID vs Description table

gene <- GFF[V3 == "gene", list(attributes = V9)]
gene_id_only <- gene[, gene_id := sub(".*;gene_id=", "", attributes)]

gene_id_descr <- gene_id_only[, description := sub(".*;description=", "", attributes)]
gene_id_descr <- gene_id_descr[, description := sub(";.*", "", description)]
gene_id_descr_table <- data.table(gene_id_descr$gene_id, gene_id_descr$description)
colnames(gene_id_descr_table) <- c("Geneid", "description")

#decoding description column
class(gene_id_descr_table)
gene_id_descr_table[, description := sapply(description, URLdecode)]
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/')
gz <- gzfile('gene_id_desc_table', 'w')
write.table(gene_id_descr_table, gz, sep= '\t', row.names= FALSE, quote= FALSE)
close(gz)

#now merge the decription table with the dge table
dge_descr <- merge(dge, gene_id_descr_table, by.x= 'gene_id', by.y= 'Geneid', all.x= TRUE, sort= FALSE)
dge_descr
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/')
gz <- gzfile('dge_descr_table_time0_additionalcontrasts_batchcorrection_9April', 'w')
write.table(dge_descr, gz, sep= '\t', row.names= FALSE, quote= FALSE)
close(gz)

##Creating an MA plot
#plot shows DE vs non-DE genes with log2 FC against normalised mean count
#the higher the mean count, the lower the log fold change required for the gene to be DE
#for genes with low mean count, need a higher log FC to be DE

dge_descr
is.data.table(dge_descr) == TRUE

#order the contrasts appropriately
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

###Volcano plot
#Another common visualisation is the volcano plot which displays a 
#measure of significance on the y-axis and fold-change on the x-axis. 
#In this case we use the log2 fold change (logFC) on the x-axis, and on the 
#y-axis we’ll use -log10(FDR). This -log10 transformation 
#is commonly used for p-values as it means that more significant genes 
#have a higher scale. The significantly differentially expressed genes are the ones 
#on the upper right and left hand side - add a column to the dataframe to specify 
#if they are UP- or DOWN- regulated
#order the contrasts appropriately
xord <- c('h4vs0', 'h8vs4', 'h12vs8', 'h16vs4', 'h16vs8', 'h16vs12', 'h24vs12', 'h24vs16')
dge_descr[, contrast_order := factor(contrast, xord)]

dge_volcano <- dge_descr

#add a column of NA's
dge_volcano$diffexpressed <- "NO"

#if logfold change > 0.58 and p-value less than 0.01, set as "UP"
dge_volcano$diffexpressed[dge_volcano$logFC > 0.58 & dge_volcano$PValue < 0.001] <- "UP"
dge_volcano$diffexpressed[dge_volcano$logFC < -0.58 & dge_volcano$PValue < 0.001] <- "DOWN"

# Automate a bit
mycolours <- c("blue", "red", "black")
names(mycolours) <- c("DOWN", "UP", "NO")

# Names of genes besides the points
dge_volcano$delabel <- NA
dge_volcano$delabel[dge_volcano$diffexpressed != "NO"] <- dge_volcano$gene_id[dge_volcano$diffexpressed != "NO"]

library(ggplot2)
library(ggrepel)
volcano_plot <- ggplot(dge_volcano, mapping = aes(x= logFC, y= -log10(PValue), col = diffexpressed, label=delabel)) +
  geom_point() +
  ylab("-log10 (P-value)") +
  xlab("log2 Fold change") +
  ggtitle("Volcano plot with batch correction, p-value < 0.001") +
  geom_vline(xintercept=c(-0.58,0.58), col="red") +
  geom_hline(yintercept=-log10(0.001), col="red") +
  facet_wrap(~contrast_order, ncol= 2) +
  scale_colour_manual(values=mycolours)+ ##change point colour
  theme_minimal()
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/')
ggsave('Volcano_plot_order_batchcorrection_18Apr_p<0.001.png', width= 20, height= 25, units= 'cm')

### Global expression plot

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

### Table of log rpkm counts

#modify counts table to include length & modify data
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/featureCounts/')
counts <- fread(cmd= 'grep -v "^#" counts.tsv')
setnames(counts, names(counts), sub('.bam', '', basename(names(counts))))
counts_min_outliers_length <- counts[, c('Geneid', 'Length', ss_min_outliers$library_id), with= FALSE]
counts_min_outliers_length <- counts_min_outliers_length [!counts_min_outliers_length$Geneid %in% RNA_Gene_ID]

gene_length <- data.table(counts_min_outliers_length$Geneid, counts_min_outliers_length$Length)
setnames(gene_length, "V1", "Geneid")
setnames(gene_length, "V2", "Length")
gene_length <- gene_length[keep, , ]

logrpkm <- rpkm(y, gene.length = gene_length$Length, log=TRUE)
logrpkm <- data.table(logrpkm, keep.rownames = "Geneid")
write.table(logrpkm, file="logrpkm_table", row.names = FALSE, sep= '\t', quote= FALSE)
