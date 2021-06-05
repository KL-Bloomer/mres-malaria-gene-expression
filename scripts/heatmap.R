#Heatmap, identification of gene clusters & average gene expression of clusters plot

#Packages required
library(data.table)
library(gplots)
library(RColorBrewer)
library(ggplot2)


logrpkm_table_long <- snakemake@input[['logrpkm_table']]
ss <- snakemake@input[['sample_sheet']]
dge_table <- snakemake@input[['dge_table']]
geneid_desc_table <- snakemake@input[['geneid_desc_table']]
GAF <- snakemake@input[['GAF']]
Heatmap_DE_genes <- snakemake@output[['Heatmap_DE_genes']]
Heatmap_genes <- snakemake@output[['Heatmap_genes']]
clusters_table <- snakemake@output[['clusters_table']]
avergene_expr_clusters <- snakemake@output[['avergene_expr_clusters']]
Heatmap_AP2_genes <- snakemake@output[['Heatmap_AP2_genes']]
Heatmap_AP2_genes_FDR <- snakemake@output[['Heatmap_AP2_genes_FDR']]
Heatmap_DE_genes_logFC <- snakemake@output[['Heatmap_DE_genes_logFC']]


#Use rpkm table to cluster all the genes 
logrpkm_table <- fread(logrpkm_table_long)
#convert to wide format
logrpkm_table <- dcast(data = logrpkm_table, gene_id ~ library_id, value.var = 'logrpkm')
logrpkm_table <- logrpkm_table[, c('gene_id', ss$library_id), with= FALSE]

# Set up colour vector for time variable
ss <- fread(ss)
ss <- ss[Outliers == FALSE, ]
ss[, Time := sprintf('%.2d', Time)]
ss[, group := paste(Time)]
ss$Time <- as.factor(ss$Time)
ss$Time <- factor(c("00", "00", "00","00","00","00","02","04","04","04","04","06","08","08","08","08","12","12","16","16","16","16","24","24","24","24"))

# Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
col.cell <- c("purple","orange","green","black", "red", "pink", "blue", "grey")[ss$Time]

# a function to assign colors based on treatment time 
treatment_times <- c(0,2,4,6,8,12,16,24)
treatment_colours_options <- c("purple","orange","green","black", "red", "pink", "blue", "grey")

#Heatmap key colours
Colors= c("midnightblue", "thistle", "deeppink3")
Colors=colorRampPalette(Colors)(100)

#### Clustering & heatmap for DE genes with FDR < 0.01
## using logrpkm table that has been filtered for DG with FDR<0.01
dge_descr <- fread(dge_table, header = TRUE)
DE_filter <- dge_descr[FDR<0.01] 

logrpkm_table_DE <- logrpkm_table[logrpkm_table$gene_id %in% DE_filter$gene_id] #filter 
logrpkm_table_DE.mat <- as.matrix(logrpkm_table_DE, rownames = "gene_id")
neworder <- c("F-1_SRR5260552", "F-2_SRR5260551", "F-3_SRR5260550","F1_21725_8_4","F2_21725_8_5","F3_21725_8_6",
              "RM-2_S409","RM-3_S410","GFPcon-4hr-R1","Ap20-GFP-SFC-4h-R1","Ap20-GFP-SFC-4h-R2","RM-4_S411",
              "Ap20-GFP-SFC-8h-R1","Ap20-GFP-SFC-8h-R2","GFPcon-8hr-R1_S9","GFPcon-8hr-R2","RM-1_S408","RM-6_S413",
              "Ap20-GFP-SFC-16h-R1","Ap20-GFP-SFC-16h-R2","GFPcon-16hr-R1","GFPcon-16hr-R2","RM-5_S412","RM-7_S414",
              "RM-8_S415","RM-9_S416")
logrpkm_table_DE.mat <- logrpkm_table_DE.mat[,neworder]

#Row side column for clusters
cluster_numbers <- c(1,2,3,4,5,6,7,8)
#run clustering 
col1 <- brewer.pal(8,"Dark2")
corr_mat <- cor(t(logrpkm_table_DE.mat)) #creating a correlation
distmat_cor <- as.dist(1-corr_mat) # determining distance from corr matrix
hclust <- hclust(distmat_cor, method = "complete")
gr.row <- cutree(hclust, k=8)

#Heatmap: DEG with FDR<0.01; correlation distance matrix and "complete" hierarchical clustering
png(file=Heatmap_DE_genes, width= 800, height= 750)
par(cex.main = 1.5)
hm <- heatmap.2(logrpkm_table_DE.mat, col=Colors,
                distfun = function(logrpkm_table_DE.mat) as.dist(1-cor(t(logrpkm_table_DE.mat))),
                main="Gene Z-score for normalised expression of DEG 
                with FDR<0.01 at different time points",
                ColSideColors=col.cell, RowSideColors=col1[gr.row], scale="row", cexCol=1.5, 
                key = TRUE, keysize = 1.2, key.title = NULL, 
                density.info = "none", trace="none", labRow = FALSE, labCol = FALSE,
                margins = c(24,13), xlab = "", ylab = "", lwid = c(5,15), lhei = c(3,15))
legend("topright",legend=paste(treatment_times), title = "Time (h)", fill=treatment_colours_options,cex=1.0)
legend("bottomleft", legend=paste(cluster_numbers), title = "Clusters",fill=col1, cex=1.0)
dev.off()

#Clustering the genes into 8 clusters
hc <- as.hclust(hm$rowDendrogram)
groups <- cutree(hc, k = 8)
groups <- as.data.table(groups, keep.rownames = "gene_id")

#loading a table which contains the gene ids and their descriptions
gene_id_descr_table <- fread(geneid_desc_table)
gene_id_desc_table_cut <- gene_id_descr_table[gene_id_descr_table$gene_id %in% groups$gene_id]

groups <- merge(groups, gene_id_desc_table_cut, by = "gene_id") #table with geneid, cluster no. & description
groups <- groups[order(groups)]
write.table(groups, file=clusters_table, row.names = FALSE, sep= '\t', quote= FALSE)

### Creating a plot showing the average gene Z-score in each cluster - SCALED values

logrpkm_table <- hm$carpet
logrpkm_table <- as.data.table(logrpkm_table, keep.rownames = "Time")
logrpkm_table_long <- melt(logrpkm_table, variable.name = "gene_id", id.vars = "library_id",
                           value.name = "logrpkm")
logrpkm_table_long <- data.table(logrpkm_table_long)
logrpkm_table_long <- merge(logrpkm_table_long, groups, by= "gene_id")

#Add time point to each library_id in logrpkm, unless you already have such column:

# Add Time to each library_id
logrpkm_table_long <- merge(logrpkm_table_long, ss[, list(library_id, Time)], by= 'library_id')

avg_clst <- logrpkm_table_long[, list(
  zscore = mean(zscore), 
  ngenes= length(unique(.SD$gene_id)),  #.SD = internal variable of data.table which contains each combination of "groups" and "time"
  sd= sd(zscore)), by= list(groups, Time)]

#Add a printable panel title
avg_clst[, panel_title := paste('Cluster ', groups, ' | N = ', ngenes, sep= "")]
avg_clst[, Time := as.numeric(as.character(Time))]

# Reorder clusters from bigger to smaller
avg_clst <- avg_clst[order(-ngenes)]
xord <- unique(avg_clst$panel_title)
avg_clst[, panel_title := factor(panel_title, xord)]
#Time is continuous

ggplot(data= avg_clst, aes(x= Time, y= zscore, by = Time)) +
  geom_line() +
  geom_point(size=0.5) +
  facet_wrap(~panel_title, nrow=2) +
  geom_errorbar(aes(ymin=zscore - sd, ymax=zscore+sd), width=.2, position=position_dodge(.9)) +
  ggtitle("Changes in average gene expression of clusters over time - scaled") +
  xlab("Time (hr)") +
  ylab("Average gene Z-score") +
  theme_bw() +
  theme(axis.text = element_text(size = 12), 
        plot.margin=unit(c(1.5,7,1.5,1.5),"cm"),
        strip.text = element_text(size = 14, face = "bold"),
        axis.title=element_text(size=14,face="bold")) +
  guides(y.sec = guide_axis()) +
  scale_x_continuous(breaks = c(0, 4, 8, 12, 16, 24))

ggsave(avergene_expr_clusters, width= 30, height= 20, units= 'cm')

##Heatmap for all the genes correlation distance matrix and default for hclust = "complete" method

logrpkm_table.mat <- as.matrix(logrpkm_table, rownames = "gene_id")
neworder <- c("F-1_SRR5260552", "F-2_SRR5260551", "F-3_SRR5260550","F1_21725_8_4","F2_21725_8_5","F3_21725_8_6",
              "RM-2_S409","RM-3_S410","GFPcon-4hr-R1","Ap20-GFP-SFC-4h-R1","Ap20-GFP-SFC-4h-R2","RM-4_S411",
              "Ap20-GFP-SFC-8h-R1","Ap20-GFP-SFC-8h-R2","GFPcon-8hr-R1_S9","GFPcon-8hr-R2","RM-1_S408","RM-6_S413",
              "Ap20-GFP-SFC-16h-R1","Ap20-GFP-SFC-16h-R2","GFPcon-16hr-R1","GFPcon-16hr-R2","RM-5_S412","RM-7_S414",
              "RM-8_S415","RM-9_S416")
logrpkm_table.mat <- logrpkm_table.mat[,neworder]

png(file=Heatmap_genes, width= 800, height= 750)
par(cex.main=1.5)
heatmap.2(logrpkm_table.mat, col=Colors,
          distfun = function(logrpkm_table.mat) as.dist(1-cor(t(logrpkm_table.mat))),
          dendrogram = "row",
          Colv = FALSE,
          main="Gene Z-score for normalised expression for 
          all genes at different time points",
          ColSideColors=col.cell,scale="row", cexCol=1.5,
          key = TRUE, keysize = 1.2,  key.title = "Colour Key",
          density.info = "none", trace="none", labRow = FALSE, labCol = FALSE, 
          margins = c(10,10), lwid = c(5,15), lhei = c(3,15))
legend("topright",legend=paste(treatment_times),title = "Time(h)", fill=treatment_colours_options,cex=1.0)
dev.off()

## Heatmap for all AP2 transcription factors

#Use rpkm table to cluster all the genes 
logrpkm_table <- fread(logrpkm_table_long)
#convert to wide format
logrpkm_table <- dcast(data = logrpkm_table, gene_id ~ library_id, value.var = 'logrpkm')
logrpkm_table <- logrpkm_table[, c('gene_id', ss$library_id), with= FALSE]

#Filter for AP2 genes 
AP2_filter <- dge_descr[grepl("AP2 domain", description, fixed = TRUE)]

#Gene names
gaf <- fread(cmd= paste('grep -v '^!'', GAF), select= c(2, 3, 10), col.names= c('Geneid', 'gene_name', 'description'))
gaf <- unique(gaf)

logrpkm_table_AP2 <- logrpkm_table[logrpkm_table$gene_id %in% AP2_filter$gene_id] #filter 
logrpkm_table_AP2 <- merge(logrpkm_table_AP2, gaf, by.x= 'gene_id', by.y= 'Geneid', all.x= TRUE, sort= FALSE)
logrpkm_table_AP2 <- logrpkm_table_AP2[, c("gene_name", ss$library_id), with= FALSE]
logrpkm_table_AP2.mat <- as.matrix(logrpkm_table_AP2, rownames = "gene_name")
neworder <- c("F-1_SRR5260552", "F-2_SRR5260551", "F-3_SRR5260550","F1_21725_8_4","F2_21725_8_5","F3_21725_8_6",
              "RM-2_S409","RM-3_S410","GFPcon-4hr-R1","Ap20-GFP-SFC-4h-R1","Ap20-GFP-SFC-4h-R2","RM-4_S411",
              "Ap20-GFP-SFC-8h-R1","Ap20-GFP-SFC-8h-R2","GFPcon-8hr-R1_S9","GFPcon-8hr-R2","RM-1_S408","RM-6_S413",
              "Ap20-GFP-SFC-16h-R1","Ap20-GFP-SFC-16h-R2","GFPcon-16hr-R1","GFPcon-16hr-R2","RM-5_S412","RM-7_S414",
              "RM-8_S415","RM-9_S416")
logrpkm_table_AP2.mat <- logrpkm_table_AP2.mat[,neworder]

#Row side column for clusters
cluster_numbers <- c(1,2,3,4,5,6,7,8,9)
#run clustering 
col1 <- brewer.pal(9,"Paired")
corr_mat <- cor(t(logrpkm_table_AP2.mat)) #creating a correlation
distmat_cor <- as.dist(1-corr_mat) # determining distance from corr matrix
hclust <- hclust(distmat_cor, method = "complete")
gr.row <- cutree(hclust, k=9)

png(file=Heatmap_AP2_genes, width= 800, height= 750)
par(cex.main = 1.5)
hm <- heatmap.2(logrpkm_table_AP2.mat, col=Colors,
                distfun = function(logrpkm_table_DE.mat) as.dist(1-cor(t(logrpkm_table_DE.mat))),
                Colv = FALSE,
                dendrogram = "row",
                main=" Gene Z-score for normalised expression of AP2 TFs
                at different time points",
                ColSideColors=col.cell,scale="row", cexRow = 1.5,
                key = TRUE, keysize = 1.2, key.title = NULL, 
                density.info = "none", trace="none",labCol = FALSE,
                margins = c(25,16), xlab = "", ylab = "", lwid = c(5,15), lhei = c(3,15))
legend("topright",legend=paste(treatment_times), title = "Time (h)", fill=treatment_colours_options,cex=1.0)
legend("bottomleft", legend=paste(cluster_numbers), title = "Clusters",fill=col1, cex=1.0)
dev.off()

### Heatmap for AP2 transcription factors with FDR<0.01 

#Use rpkm table to cluster all the genes 
logrpkm_table <- fread(logrpkm_table_long)
#convert to wide format
logrpkm_table <- dcast(data = logrpkm_table, gene_id ~ library_id, value.var = 'logrpkm')
logrpkm_table <- logrpkm_table[, c('gene_id', ss$library_id), with= FALSE]

#Filter for AP2 genes
AP2_filter <- DE_filter[grepl("AP2 domain", description, fixed = TRUE)]

#Gene names
gaf <- fread(cmd= paste('grep -v '^!'', GAF_file), select= c(2, 3, 10), col.names= c('Geneid', 'gene_name', 'description'))
gaf <- unique(gaf)

logrpkm_table_AP2 <- logrpkm_table[logrpkm_table$gene_id %in% AP2_filter$gene_id] #filter 
logrpkm_table_AP2 <- merge(logrpkm_table_AP2, gaf, by.x= 'gene_id', by.y= 'Geneid', all.x= TRUE, sort= FALSE)
logrpkm_table_AP2 <- logrpkm_table_AP2[, c("gene_name", ss$library_id), with= FALSE]
logrpkm_table_AP2.mat <- as.matrix(logrpkm_table_AP2, rownames = "gene_name")
neworder <- c("F-1_SRR5260552", "F-2_SRR5260551", "F-3_SRR5260550","F1_21725_8_4","F2_21725_8_5","F3_21725_8_6",
              "RM-2_S409","RM-3_S410","GFPcon-4hr-R1","Ap20-GFP-SFC-4h-R1","Ap20-GFP-SFC-4h-R2","RM-4_S411",
              "Ap20-GFP-SFC-8h-R1","Ap20-GFP-SFC-8h-R2","GFPcon-8hr-R1_S9","GFPcon-8hr-R2","RM-1_S408","RM-6_S413",
              "Ap20-GFP-SFC-16h-R1","Ap20-GFP-SFC-16h-R2","GFPcon-16hr-R1","GFPcon-16hr-R2","RM-5_S412","RM-7_S414",
              "RM-8_S415","RM-9_S416")
logrpkm_table_AP2.mat <- logrpkm_table_AP2.mat[,neworder]

#Row side column for clusters
cluster_numbers <- c(1,2,3,4,5,6,7,8,9)
#run clustering 
col1 <- brewer.pal(9,"Paired")
corr_mat <- cor(t(logrpkm_table_AP2.mat)) #creating a correlation
distmat_cor <- as.dist(1-corr_mat) # determining distance from corr matrix
hclust <- hclust(distmat_cor, method = "complete")
gr.row <- cutree(hclust, k=9)

#Heatmap: AP2; correlation distance matrix and "complete" hierarchical clustering
png(file=Heatmap_AP2_genes_FDR, width= 800, height= 750)
par(cex.main = 1.5)
hm <- heatmap.2(logrpkm_table_AP2.mat, col=Colors,
                distfun = function(logrpkm_table_AP2.mat) as.dist(1-cor(t(logrpkm_table_AP2.mat))),
                Colv = FALSE,
                dendrogram = "row",
                main="Gene Z-score for normalised expression of AP2 TFs 
                with FDR<0.01 at different time points",
                ColSideColors=col.cell, RowSideColors =col1[gr.row], scale="row", cexRow = 1.5,
                key = TRUE, keysize = 1.2, key.title = NULL, 
                density.info = "none", trace="none",labCol = FALSE, labRow = TRUE, 
                margins = c(21,19), xlab = "", ylab = "", lwid = c(5,15), lhei = c(3,15))
legend("topright",legend=paste(treatment_times), title = "Time (h)", fill=treatment_colours_options,cex=1.0)
legend("bottomleft", legend=paste(cluster_numbers), title = "Clusters",fill=col1, cex=1.0)
dev.off()


#Heatmap showing logFC of DEG
mat <- dcast(data = dge_descr, gene_id ~ contrast, value.var = 'logFC')
mat <- as.matrix(mat, rownames = "gene_id")

png(file=Heatmap_DE_genes_logFC, width= 800, height= 750)
par(cex.axis = 1.5)
heatmap.2(mat, col=Colors,
          distfun = function(mat) as.dist(1-cor(t(mat))),
          main="Gene Z-score for logFC of DEG 
          at different time points",
          scale="row", cexCol=1.5, 
          key = TRUE, keysize = 1.2,  key.title = "Colour Key",
          density.info = "none", trace="none", labRow = FALSE,
          margins = c(10,5), xlab = "Contrasts", ylab = "Genes", lwid = c(5,12), lhei = c(3,15))
dev.off()





