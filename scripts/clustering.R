##Clustering and plot of how interesting genes change with time
#https://www.youtube.com/watch?v=9U4h6pZw6f8
#http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning
#http://www.r-tutor.com/gpu-computing/clustering/hierarchical-cluster-analysis
#http://www.r-tutor.com/gpu-computing/clustering/distance-matrix

# =======================

## clustering, rpkm normalises for gene length and sequencing depth(library size)
#https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/
rm(list=ls())
#Packages required
library(tidyr)
library(data.table)
library(ggplot2)
library(ggdendro)
library(dplyr)
library("gplots")
library("RColorBrewer")
library(reshape2)

#Use rpkm table to cluster all the genes 
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/featureCounts/')
logrpkm_table <- fread('logrpkm_table')

#converting to long format
logrpkm_table_long <- melt(logrpkm_table, variable.name = "Time", id.vars = "Geneid",
                           value.name = "Logrpkm_counts")
logrpkm_table_long <- data.table(logrpkm_table_long)


class(logrpkm_table)
logrpkm_table.mat <- as.matrix(logrpkm_table, rownames = "Geneid")

#run clustering 
dist_matrix_euc <- dist(logrpkm_table) #Euclidean distance; error - introduced NA
which(is.na(dist_matrix_euc))

corr_mat <- cor(t(logrpkm_table.mat)) #creating a correlation
distmat_cor <- as.dist(1-corr_mat) # determining distance from corr matrix

#hclust algorithm
h.cl_complete_cor <- hclust(distmat_cor, method = "complete")
h.cl_median_cor <- hclust(distmat_cor, method = "median")
h.cl_average_cor <- hclust(distmat_cor, method = "average")
h.cl_ward_cor <- hclust(distmat_cor, method = "ward.D")
summary(h.cl_complete_cor)


# Visualization using the default theme named theme_dendro(); ggdendrogram = wrapper
# generate graphs for all the above hclust
q <- ggdendrogram(h.cl_median_cor, leaf_labels = FALSE, rotate = FALSE, size = 2)+
  ggtitle("Median hierarchical clustering using a correlation distance matrix")+
  theme(axis.text.x = element_blank())+
  theme(plot.title = element_text(size = 50, face = "bold"))
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/')
ggsave('dendrogram_corr_medianclust.png', width = 25, height = 25)

## Using cutree to create groups
groups <- cutree(h.cl_complete_cor, k=50) # if 5000 genes, should have about 100 genes per group
groups <- as.data.table(groups, keep.rownames = "Geneid")

setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/')
gene_id_descr_table <- fread('gene_id_desc_table (1)')
gene_id_desc_table_cut <- gene_id_descr_table[gene_id_descr_table$Geneid %in% groups$Geneid]

groups <- merge(groups, gene_id_desc_table_cut, by = "Geneid")
groups <- groups[order(groups)]
write.table(groups, file="genes_group_descr_allgenes.txt", row.names = FALSE, sep= '\t', quote= FALSE)

#### Clustering the DE genes with FDR < 0.01
## using logrpkm table that has been filtered for DG with FDR<0.01
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/')
dge_descr <- fread('dge_descr_table_time0_additionalcontrasts_batchcorrection_9April (1)', header = TRUE)
DE_filter <- dge_descr %>% #requires dplyr library
  filter(FDR < 0.01)

logrpkm_table_DE <- logrpkm_table[logrpkm_table$Geneid %in% DE_filter$gene_id] #filter 

logrpkm_table_DE.mat <- as.matrix(logrpkm_table_DE, rownames = "Geneid")

#run clustering 
dist_matrix_euc <- dist(logrpkm_table_DE.mat) #Euclidean distance; error - introduced NA
which(is.na(dist_matrix_euc))

corr_mat <- cor(t(logrpkm_table_DE.mat)) #creating a correlation
distmat_cor <- as.dist(1-corr_mat) # determining distance from corr matrix

#hclust algorithm - correlation distance matrix
h.cl_complete_cor <- hclust(distmat_cor, method = "complete")
h.cl_median_cor <- hclust(distmat_cor, method = "median")
h.cl_average_cor <- hclust(distmat_cor, method = "average")
h.cl_ward_cor <- hclust(distmat_cor, method = "ward.D")
summary(h.cl_complete_cor)

# Visualization using the default theme named theme_dendro(); ggdendrogram = wrapper
# generate graphs for all the above hclust
q <- ggdendrogram(h.cl_ward_cor, leaf_labels = FALSE, rotate = FALSE, size = 2)+
  ggtitle("Ward.D hierarchical clustering of DE genes with FDR < 0.01 using a correlation distance matrix")+
  theme(axis.text.x = element_blank())+
  theme(plot.title = element_text(size = 30, face = "bold"))
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/')
ggsave('dendrogram_corr_ward_DEclust.png', width = 25, height = 25)

## Using cutree to create groups
groups <- cutree(h.cl_complete_cor, k=25) # if 5000 genes, should have about 100 genes per group
groups <- as.data.table(groups, keep.rownames = "Geneid")

setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/')
gene_id_descr_table <- fread('gene_id_desc_table (1)')
gene_id_desc_table_cut <- gene_id_descr_table[gene_id_descr_table$Geneid %in% groups$Geneid]

groups <- merge(groups, gene_id_desc_table_cut, by = "Geneid")
groups <- groups[order(groups)]
write.table(groups, file="genes_group_descr_DE_FDR<0.01.txt", row.names = FALSE, sep= '\t', quote= FALSE)

### Clustering the logFC
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/')
dge_descr <- fread('dge_descr_table_time0_additionalcontrasts_batchcorrection_9April (1)')
mat <- dcast(data = dge_descr, gene_id ~ contrast, value.var = 'logFC')
mat <- as.matrix(mat, rownames = "gene_id")
class(mat)
#run clustering 
dist_matrix_euc <- dist(mat) #Euclidean distance; error - introduced NA
which(is.na(dist_matrix_euc))

corr_mat <- cor(t(mat)) #creating a correlation
distmat_cor <- as.dist(1-corr_mat) # determining distance from corr matrix

#hclust algorithm - correlation distance matrix
h.cl_complete_cor <- hclust(distmat_cor, method = "complete")
h.cl_median_cor <- hclust(distmat_cor, method = "median")
h.cl_average_cor <- hclust(distmat_cor, method = "average")
h.cl_ward_cor <- hclust(distmat_cor, method = "ward.D")

# Visualization using the default theme named theme_dendro(); ggdendrogram = wrapper
# generate graphs for all the above hclust
q <- ggdendrogram(h.cl_median_cor, leaf_labels = FALSE, rotate = FALSE, size = 2)+
  ggtitle("Median hierarchical clustering of DE genes logFC using a correlation distance matrix")+
  theme(axis.text.x = element_blank())+
  theme(plot.title = element_text(size = 30, face = "bold"))
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/')
ggsave('dendrogram_corr_median_DElogFC.png', width = 25, height = 25)

#group clusters - about 5000 genes so use 50 clusters & 100 clusters
groups <- cutree(h.cl_complete_cor, k=50) 
groups <- as.data.table(groups, keep.rownames = "Geneid")
gene_id_desc_table_cut <- gene_id_descr_table[gene_id_descr_table$Geneid %in% groups$Geneid]

groups <- merge(groups, gene_id_desc_table_cut, by = "Geneid")
groups <- groups[order(groups)]
write.table(groups, file="genes_group_descr_logFC_50groups.txt", row.names = FALSE, sep= '\t', quote= FALSE)


### Plotting the gene expression changes over time of key genes
#Instead of plotting the PBANKA ids it would be more useful to 
#plot the gene names
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/')
genes <- fread('Interesting_genes.txt')
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/ref/')
gaf <- fread(cmd= "grep -v '^!' PlasmoDB-49_PbergheiANKA_GO.gaf", select= c(2, 3, 10), col.names= c('Geneid', 'gene_name', 'description'))
gaf <- unique(gaf)

#Take care that not all genes have a name so you should fill in any 
#missing values semi-manually. Then add these gene names to your dataframe 
#and from then on work with "gene_name" instead of "Geneid"

genes <- merge(genes, gaf, by= 'Geneid', all.x= TRUE)
set(genes, i = 11L, 3L, "Cap93")
set(genes, i=11L, 4L, "oocyst capsule protein Cap93")
set(genes, i=18L, 3L, "WARP")
set(genes, i=18L, 4L, "von Willebrand factor A domain-related protein")

### Interesting genes
selected <- c('PBANKA_0800500', 'PBANKA_1319500', 'PBANKA_1300700', 'PBANKA_1035200', 
              'PBANKA_0417600', 'PBANKA_1315300', 'PBANKA_0204500', 'PBANKA_1037800',
              'PBANKA_1228900', 'PBANKA_0515000', 'PBANKA_0514900', 'PBANKA_1436600',
              'PBANKA_0402600', 'PBANKA_0905900', 'PBANKA_1363700','PBANKA_0620600',
              'PBANKA_1414900', 'PBANKA_1217700', 'PBANKA_1312700', 'PBANKA_1415700',
              'PBANKA_1001800', 'PBANKA_1302800', 'PBANKA_1436100', 'PBANKA_0905200',
              'PBANKA_0615200', 'PBANKA_0412900', 'PBANKA_1227400', 'PBANKA_1432200',
              'PBANKA_0314200', 'PBANKA_1301300')


slct <- logrpkm_table_long[logrpkm_table_long$Geneid %in% selected]
slct <- merge(slct, genes, by= 'Geneid')
slct$Time <- slct [, gsub("Ap20-GFP-SFC-4h-R[0-9]", "Ap20 4h", slct$Time)]
slct$Time <- slct [, gsub("Ap20-GFP-SFC-16h-R[0-9]", "Ap20 16h", slct$Time)]
slct$Time <- slct [, gsub("Ap20-GFP-SFC-8h-R[0-9]", "Ap20 8h", slct$Time)]
slct$Time <- slct[, gsub("GFPcon-4hr-R1", "GFPcon 4h", slct$Time)]
slct$Time <- slct [, gsub("GFPcon-16hr-R[0-9]", "GFPcon 16h", slct$Time)]
slct$Time <- slct [, gsub("GFPcon-8hr-R1_S9", "GFPcon 8h", slct$Time)]
slct$Time <- slct [, gsub("GFPcon-8hr-R2", "GFPcon 8h", slct$Time)]
slct$Time <- slct [, gsub("RM-9_S416", "RM9 24h", slct$Time)]
slct$Time <- slct [, gsub("RM-8_S415", "RM8 24h", slct$Time)]
slct$Time <- slct [, gsub("RM-7_S414", "RM7 24h", slct$Time)]
slct$Time <- slct [, gsub("RM-6_S413", "RM6 12h", slct$Time)]
slct$Time <- slct [, gsub("RM-5_S412", "RM5 24h", slct$Time)]
slct$Time <- slct [, gsub("RM-4_S411", "RM4 6h", slct$Time)]
slct$Time <- slct [, gsub("RM-3_S410", "RM3 4h", slct$Time)]
slct$Time <- slct [, gsub("RM-2_S409", "RM2 2h", slct$Time)]
slct$Time <- slct [, gsub("RM-1_S408", "RM1 12h", slct$Time)]
slct$Time <- slct [, gsub("F-.*", "SRR526055_RFP 0h", slct$Time)]
slct$Time <- slct [, gsub("F[1-3].*", "21725_RFP 0h", slct$Time)]

slct <- separate(data = slct, col = Time, into = c("Strain", "Time"), sep = "([ ])")
slct$Time <- slct [, gsub("h", "", slct$Time)]
#slct$Time <- as.numeric(slct$Time)
#slct$Time <- sprintf('%02d', slct$Time)
slct[, label := paste(Time, Strain)]

average <- slct %>%                            # name of the dataset
  group_by(Time, gene_name) %>%               # grouping the data 
  summarize(m = mean(Logrpkm_counts)) %>%      # calculating the mean
  ungroup()                      # ungroup the data

#order genes by shape of expression profile
average <- as.data.table(average)
average$Time <- as.numeric(average$Time)
ave_mat <- dcast(data= average, gene_name ~ Time, value.var= 'm')
ave_mat <- as.data.table(ave_mat)
class(ave_mat)
ave_mat <- as.matrix(ave_mat, rownames = "gene_name")

#dd <- dist(ave_mat), euclidean distance
str(ave_mat)
corr_mat <- cor(t(ave_mat)) #creating a correlation 
distmat_cor <- as.dist(1-corr_mat) # determining distance from corr matrix

#dd <- dist(ave_mat)
dendro <- as.dendrogram(hclust(distmat_cor))
# Ladderize could give slightly nicer ordering. 
# You need library(dendextend)
#reorganizes the internal structure of the tree to 
#get the ladderized effect when plotted.
library(dendextend)
dendro <- ladderize(dendro)

# This is the new order and we refactor the datasets accordingly
gene_order <- labels(dendro)
slct$gene_name <- factor(slct$gene_name, gene_order)
average$gene_name <- factor(average$gene_name, gene_order)

slct$Time <- as.numeric(slct$Time)

test <- ggplot(data = slct, aes(x=Time, y = Logrpkm_counts, group =1)) +
  geom_point(size = 0.5) + 
  facet_wrap(~gene_name) +
  theme(axis.text.x = element_text(angle = 90))+
  ggtitle("Changes in gene expression of key genes over time") +
  xlab("Time (hr)") +
  ylab("Counts (log(rpkm)") +
  theme(axis.text.x = element_text(size = 7, angle=90), plot.margin=unit(c(1.5,7,1.5,1.5),"cm"))

test + 
  geom_line(data = average, aes(x = Time, y=m, group = 1), colour = "blue")+
  guides(y.sec = guide_axis()) 

setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/')
ggsave('gene_expression_changes_keygenes_extragenes_27Apr.png', width= 30, height= 20, units= 'cm')

#Plotting the average clusters

#Center individual gene expressions by subtracting the mean. 
#In this way we make the average expression of each gene centered on zero. 
#Say you have the data.table (not data.frame) of logrpkm in long format with 
#RPKM in column "rpkm":

#Use rpkm table to cluster all the genes 
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/featureCounts/')
logrpkm_table <- fread('logrpkm_table')

#converting to long format
logrpkm_table_long <- melt(logrpkm_table, variable.name = "Time", id.vars = "Geneid",
                           value.name = "Logrpkm_counts")
logrpkm_table_long <- data.table(logrpkm_table_long)
logrpkm_table_long[, center_rpkm := scale(Logrpkm_counts, center= TRUE, scale= FALSE), by= Geneid]

#Now assign to each gene the cluster it belongs to genes_group_descr_allgenes 
#is a data.table with column Geneid and cluster ID
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/')
genes_group_descr_allgenes <- fread('genes_group_descr_allgenes.txt')
logrpkm_table_long <- merge(logrpkm_table_long, genes_group_descr_allgenes, by= 'Geneid')

#Add time point to each library_id in logrpkm, unless you already have such column:

logrpkm_table_long$Time <- logrpkm_table_long [, gsub("Ap20-GFP-SFC-4h-R[0-9]", "Ap20 4h", logrpkm_table_long$Time)]
logrpkm_table_long$Time <- logrpkm_table_long [, gsub("Ap20-GFP-SFC-16h-R[0-9]", "Ap20 16h", logrpkm_table_long$Time)]
logrpkm_table_long$Time <- logrpkm_table_long [, gsub("Ap20-GFP-SFC-8h-R[0-9]", "Ap20 8h", logrpkm_table_long$Time)]
logrpkm_table_long$Time <- logrpkm_table_long[, gsub("GFPcon-4hr-R1", "GFPcon 4h", logrpkm_table_long$Time)]
logrpkm_table_long$Time <- logrpkm_table_long [, gsub("GFPcon-16hr-R[0-9]", "GFPcon 16h", logrpkm_table_long$Time)]
logrpkm_table_long$Time <- logrpkm_table_long [, gsub("GFPcon-8hr-R1_S9", "GFPcon 8h", logrpkm_table_long$Time)]
logrpkm_table_long$Time <- logrpkm_table_long [, gsub("GFPcon-8hr-R2", "GFPcon 8h", logrpkm_table_long$Time)]
logrpkm_table_long$Time <- logrpkm_table_long[, gsub("RM-9_S416", "RM9 24h", logrpkm_table_long$Time)]
logrpkm_table_long$Time <- logrpkm_table_long[, gsub("RM-8_S415", "RM8 24h", logrpkm_table_long$Time)]
logrpkm_table_long$Time <- logrpkm_table_long [, gsub("RM-7_S414", "RM7 24h", logrpkm_table_long$Time)]
logrpkm_table_long$Time <- logrpkm_table_long[, gsub("RM-6_S413", "RM6 12h", logrpkm_table_long$Time)]
logrpkm_table_long$Time <- logrpkm_table_long[, gsub("RM-5_S412", "RM5 24h", logrpkm_table_long$Time)]
logrpkm_table_long$Time <- logrpkm_table_long [, gsub("RM-4_S411", "RM4 6h", logrpkm_table_long$Time)]
logrpkm_table_long$Time <- logrpkm_table_long[, gsub("RM-3_S410", "RM3 4h", logrpkm_table_long$Time)]
logrpkm_table_long$Time <- logrpkm_table_long [, gsub("RM-2_S409", "RM2 2h", logrpkm_table_long$Time)]
logrpkm_table_long$Time <- logrpkm_table_long [, gsub("RM-1_S408", "RM1 12h", logrpkm_table_long$Time)]
logrpkm_table_long$Time <- logrpkm_table_long [, gsub("F-.*", "SRR526055_RFP 0h", logrpkm_table_long$Time)]
logrpkm_table_long$Time <- logrpkm_table_long [, gsub("F[1-3].*", "21725_RFP 0h", logrpkm_table_long$Time)]

logrpkm_table_long<- separate(data = logrpkm_table_long, col = Time, into = c("Strain", "Time"), sep = "([ ])")
logrpkm_table_long$Time <- logrpkm_table_long[, gsub("h", "", logrpkm_table_long$Time)]

#Average expression within clusters and time points:
avg_clst <- logrpkm_table_long[, list(Logrpkm_counts= mean(Logrpkm_counts), sd= sd(Logrpkm_counts)), by= list(Time, groups)]

#Now you should be able to plot
avg_clst$Time <- as.numeric(avg_clst$Time) #Time = continuous
avg_clst$Time <- sprintf('%.2d', avg_clst$Time) #Time = character

ggplot(data= avg_clst, aes(x= Time, y= Logrpkm_counts, group =1)) +
  geom_line() +
  geom_point(size=0.5) +
  facet_wrap(~groups) +
  geom_errorbar(aes(ymin=Logrpkm_counts - sd, ymax=Logrpkm_counts+sd), width=.2, position=position_dodge(.9)) +
  ggtitle("Changes in average gene expression of clusters over time") +
  xlab("Time (hr)") +
  ylab("Counts (log(rpkm)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 7, angle=90), plot.margin=unit(c(1.5,7,1.5,1.5),"cm")) +
  guides(y.sec = guide_axis()) 

setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/')
ggsave('gene_expression_changes_50clusters_timediscrete.png', width= 30, height= 20, units= 'cm')

setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/Data_for_report/')
interesting_genes <- fread('Interesting_genes_clustnmbr.txt')
interesting_genes <- interesting_genes[order(cluster)]
write.table(interesting_genes, file="Interesting_genes_clustnmbr.txt", row.names = FALSE, sep= '\t', quote= FALSE)
