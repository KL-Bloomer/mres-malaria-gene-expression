#Heatmap, identification of gene clusters & average gene expression of clusters plot

#Packages required
library(data.table)
library(gplots)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(ggplot2)

logrpkm_table_long <- snakemake@input[['logrpkm_table']]
ss <- snakemake@input[['sample_sheet']]
dge_table <- snakemake@input[['dge_table']]
geneid_desc_table <- snakemake@input[['geneid_desc_table']]
Heatmap_DE_genes <- snakemake@output[['Heatmap_DE_genes']]
clusters_table <- snakemake@output[['clusters_table']]
avergene_expr_clusters <- snakemake@output[['avergene_expr_clusters']]
Heatmap_DE_genes_logFC <- snakemake@output[['Heatmap_DE_genes_logFC']]
Heatmap_genes <- snakemake@output[['Heatmap_genes']]

#Use rpkm table to cluster all the genes 
logrpkm_table <- fread(logrpkm_table_long)
#convert to wide format
logrpkm_table <- dcast(data = logrpkm_table, gene_id ~ Time, value.var = 'logrpkm')

# Set up colour vector for time variable
ss <- fread(ss)
ss <- ss[Outliers == FALSE, ]
ss[, Time := sprintf('%.2d', Time)]
ss$Time <- as.factor(ss$Time)

# Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
col.cell <- c("purple","orange","green","black", "red", "pink", "blue", "grey")[ss$Time]

# a function to assign colors based on treatment time 
treatment_times <- c(0,2,4,6,8,12,16,24)
treatment_colours_options <- c("purple","orange","green","black", "red", "pink", "blue", "grey")

#### Clustering & heatmap for DE genes with FDR < 0.01
## using logrpkm table that has been filtered for DG with FDR<0.01
dge_descr <- fread(dge_table, header = TRUE)
DE_filter <- dge_descr %>% 
  filter(FDR < 0.01)

logrpkm_table_DE <- logrpkm_table[logrpkm_table$gene_id %in% DE_filter$gene_id] #filter 
logrpkm_table_DE.mat <- as.matrix(logrpkm_table_DE, rownames = "gene_id")

#Heatmap: DEG with FDR<0.01; correlation distance matrix and "complete" hierarchical clustering
png(file=Heatmap_DE_genes, width= 800, height= 750)
par(cex.main = 1.5)
hm <- heatmap.2(logrpkm_table_DE.mat, col=brewer.pal(11,"RdBu"),
                distfun = function(logrpkm_table_DE.mat) as.dist(1-cor(t(logrpkm_table_DE.mat))),
                main="Row Z-score for normalised expression of DEG 
                with FDR < 0.01 at different time points",
                ColSideColors=col.cell,scale="row", cexCol=1.5, 
                key = TRUE, keysize = 1.2, key.title = NULL, 
                density.info = "none", trace="none", labRow = TRUE, 
                margins = c(16,10), xlab = "", ylab = "", lwid = c(5,15), lhei = c(3,15))
legend("topright",legend=paste(treatment_times,"h"),fill=treatment_colours_options,cex=1.0)
dev.off()

#Clustering the genes into 8 clusters
hc <- as.hclust(hm$rowDendrogram)
groups <- cutree(hc, k = 8)
groups <- as.data.table(groups, keep.rownames = "Geneid")

#loading a table which contains the gene ids and their descriptions
gene_id_descr_table <- fread(geneid_desc_table)
gene_id_desc_table_cut <- gene_id_descr_table[gene_id_descr_table$Geneid %in% groups$Geneid]

groups <- merge(groups, gene_id_desc_table_cut, by = "Geneid") #table with geneid, cluster no. & description
groups <- groups[order(groups)]
write.table(groups, file=clusters_table, row.names = FALSE, sep= '\t', quote= FALSE)

### Creating a plot showing the average gene Z-score in each cluster

logrpkm_table_long <- hm$carpet
logrpkm_table_long <- as.data.table(logrpkm_table_long, keep.rownames = "Time")
logrpkm_table_long <- melt(logrpkm_table_long, variable.name = "Geneid", id.vars = "Time",
                           value.name = "logrpkm")
logrpkm_table_long <- data.table(logrpkm_table_long)

logrpkm_table_long <- merge(logrpkm_table_long, groups, by= "Geneid")
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
avg_clst <- logrpkm_table_long[, list(logrpkm = mean(logrpkm), sd= sd(logrpkm)), by= list(Time, groups)]

#Now you should be able to plot
avg_clst$Time <- as.numeric(avg_clst$Time) #Time = continuous

ggplot(data= avg_clst, aes(x= Time, y= logrpkm, group =1)) +
  geom_line() +
  geom_point(size=0.5) +
  facet_wrap(~groups) +
  geom_errorbar(aes(ymin=logrpkm - sd, ymax=logrpkm+sd), width=.2, position=position_dodge(.9)) +
  ggtitle("Changes in average gene expression of clusters over time - scaled") +
  xlab("Time (hr)") +
  ylab("Average gene Z-score") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), plot.margin=unit(c(1.5,7,1.5,1.5),"cm")) +
  guides(y.sec = guide_axis()) +
  scale_x_continuous(breaks = c(0, 4, 8, 12, 16, 24))

ggsave(avergene_expr_clusters, width= 30, height= 20, units= 'cm')


#Heatmap for all the genes - correlation distance matrix and default for hclust = "complete" method
logrpkm_table.mat <- as.matrix(logrpkm_table, rownames = "gene_id")
png(file=Heatmap_genes, width= 800, height= 750)
par(cex.main=1.5)
heatmap.2(logrpkm_table.mat, col=brewer.pal(11,"RdBu"),
          distfun = function(logrpkm_table.mat) as.dist(1-cor(t(logrpkm_table.mat))),
         main="Z-score for normalised expression 
         of all the genes at different time points",
         ColSideColors=col.cell,scale="row", cexCol=1.5,
         key = TRUE, keysize = 1.2,  key.title = "Colour Key",
         density.info = "none", trace="none", labRow = TRUE, 
         margins = c(15,10), lwid = c(5,15), lhei = c(3,15))
legend("topright",legend=paste(treatment_times,"h"),fill=treatment_colours_options,cex=1.0)
dev.off()

#Heatmap showing logFC of DEG
mat <- dcast(data = dge_descr, gene_id ~ contrast, value.var = 'logFC')
mat <- as.matrix(mat, rownames = "gene_id")

png(file=Heatmap_DE_genes_logFC, width= 800, height= 750)
par(cex.axis = 1.5)
heatmap.2(mat, col=brewer.pal(11,"RdBu"),
          distfun = function(mat) as.dist(1-cor(t(mat))),
          main="Z-score for logFC of DEG at different time points",
          scale="row", cexCol=1.5, 
          key = TRUE, keysize = 1.2,  key.title = "Colour Key",
          density.info = "none", trace="none", labRow = TRUE, 
          margins = c(10,5), xlab = "Contrasts", ylab = "Genes", lwid = c(5,12), lhei = c(3,15))
dev.off()



