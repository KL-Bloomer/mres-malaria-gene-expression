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
write.table(logrpkm_table_long, file="logrpkm_table_long.txt", row.names = FALSE, sep= '\t', quote= FALSE)



logrpkm_table.mat <- as.matrix(logrpkm_table, rownames = "Geneid")

#run clustering 
corr_mat <- cor(t(logrpkm_table.mat)) #creating a correlation
distmat_cor <- as.dist(1-corr_mat) # determining distance from corr matrix

#hclust algorithm
h.cl_complete_cor <- hclust(distmat_cor, method = "complete")
h.cl_median_cor <- hclust(distmat_cor, method = "median")
h.cl_average_cor <- hclust(distmat_cor, method = "average")
h.cl_ward_cor <- hclust(distmat_cor, method = "ward.D")

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

#generate correlation distance matrix
corr_mat <- cor(t(mat)) #creating a correlation
distmat_cor <- as.dist(1-corr_mat) # determining distance from corr matrix

#hclust algorithm 
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


#Plotting the average clusters

#Center individual gene expressions by subtracting the mean. 
#In this way we make the average expression of each gene centered on zero. 
#Say you have the data.table (not data.frame) of logrpkm in long format with 
#RPKM in column "rpkm":

#Use rpkm table to cluster DE genes
#setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/featureCounts/')
#logrpkm_table <- fread('logrpkm_table')

#converting to long format
logrpkm_table_long <- melt(logrpkm_table_DE, variable.name = "Time", id.vars = "Geneid",
                           value.name = "Logrpkm_counts")
logrpkm_table_long <- data.table(logrpkm_table_long)
logrpkm_table_long[, center_rpkm := scale(Logrpkm_counts, center= TRUE, scale= TRUE), by= Geneid]

#Now assign to each gene the cluster it belongs to cluster_scaled_cluster8
#is a data.table with column Geneid and cluster ID
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/')
genes_group_descr_DE <- fread('clusters_scaled_8clusters.txt')
logrpkm_table_long <- merge(logrpkm_table_long, genes_group_descr_DE, by= 'Geneid')

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

ggplot(data= avg_clst, aes(x= Time, y= Logrpkm_counts, group =1)) +
  geom_line() +
  geom_point(size=0.5) +
  facet_wrap(~groups) +
  geom_errorbar(aes(ymin=Logrpkm_counts - sd, ymax=Logrpkm_counts+sd), width=.2, position=position_dodge(.9)) +
  ggtitle("Changes in average gene expression of clusters over time - centered and scaled") +
  xlab("Time (hr)") +
  ylab("Normalised expression (log2 rpkm)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 7, angle=90), plot.margin=unit(c(1.5,7,1.5,1.5),"cm")) +
  guides(y.sec = guide_axis()) +
  scale_x_continuous(breaks = c(0, 4, 8, 12, 16, 24))

setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/')
ggsave('gene_expression_changes_8clusters_timecontinuous_DEgenes.png', width= 30, height= 20, units= 'cm')



#Use rpkm table to cluster all the genes 
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/featureCounts/')
logrpkm_table <- fread('logrpkm_table')

#converting to long format
logrpkm_table_long <- melt(logrpkm_table, variable.name = "Time", id.vars = "Geneid",
                           value.name = "Logrpkm_counts")
logrpkm_table_long <- data.table(logrpkm_table_long)
logrpkm_table_long[, center_rpkm := scale(Logrpkm_counts, center= TRUE, scale= TRUE), by= Geneid]

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
  ggtitle("Changes in average gene expression of clusters over time - centered and scaled") +
  xlab("Time (hr)") +
  ylab("Normalised expression (log2 rpkm)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 7, angle=90), plot.margin=unit(c(1.5,7,1.5,1.5),"cm")) +
  guides(y.sec = guide_axis()) +
  scale_x_continuous(breaks = c(0, 4, 8, 12, 16, 24))

setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/')
ggsave('gene_expression_changes_50clusters_timecontinuous_scaled2May.png', width= 30, height= 20, units= 'cm')

setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/Data_for_report/')
interesting_genes <- fread('Interesting_genes_clustnmbr.txt')
interesting_genes <- interesting_genes[order(cluster)]
write.table(interesting_genes, file="Interesting_genes_clustnmbr.txt", row.names = FALSE, sep= '\t', quote= FALSE)



###Testing gene cluster from heatmap object 

#Now assign to each gene the cluster it belongs to cluster_scaled_cluster8
#is a data.table with column Geneid and cluster ID

testing <- as.data.table(hm$carpet, hm$rowMeans, hm$rowSDs)
Means <- hm$rowMeans
names <- hm$carpet
logrpkm_table_long <- melt(names, variable.name = "Time", id.vars = "Geneid",
                           value.name = "Logrpkm_counts")
logrpkm_table_long <- rename(logrpkm_table_long, Time = Var1, Geneid = Var2)
logrpkm_table_long <- data.table(logrpkm_table_long)
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/')
genes_group_descr_DE <- fread('clusters_scaled_8clusters.txt')
logrpkm_table_long <- merge(logrpkm_table_long, genes_group_descr_DE, by= 'Geneid')
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
avg_clst <- logrpkm_table_long[, list(Logrpkm_counts = mean(Logrpkm_counts), sd= sd(Logrpkm_counts)), by= list(Time, groups)]

#Now you should be able to plot
avg_clst$Time <- as.numeric(avg_clst$Time) #Time = continuous

ggplot(data= avg_clst, aes(x= Time, y= Logrpkm_counts, group =1)) +
  geom_line() +
  geom_point(size=0.5) +
  facet_wrap(~groups) +
  geom_errorbar(aes(ymin=Logrpkm_counts - sd, ymax=Logrpkm_counts+sd), width=.2, position=position_dodge(.9)) +
  ggtitle("Changes in average gene expression of clusters over time - scaled") +
  xlab("Time (hr)") +
  ylab("Average gene Z-score") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), plot.margin=unit(c(1.5,7,1.5,1.5),"cm")) +
  guides(y.sec = guide_axis()) +
  scale_x_continuous(breaks = c(0, 4, 8, 12, 16, 24))

setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/')
ggsave('testing_gene_expression_changes_8clusters_timecontinuous_DEgenes.png', width= 30, height= 20, units= 'cm')

