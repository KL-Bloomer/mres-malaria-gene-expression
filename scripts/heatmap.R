#Heatmap and identification of gene clusters
rm(list=ls())
#Packages required
library(data.table)
library(gplots)
library(RColorBrewer)
library(dplyr)

#Use rpkm table to cluster all the genes 
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/featureCounts/')
logrpkm_table <- fread('logrpkm_table')

# Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)

# Set up colour vector for time variable
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/')
libraryid_and_time <- fread("libraryid_and_time.txt") #created from y$samples
libraryid_and_time$group <- as.factor(libraryid_and_time$group)
col.cell <- c("purple","orange","green","black", "red", "pink", "blue", "grey")[libraryid_and_time$group] #need to run script_for_RNAseq_plots.R for y

# a function to assign colors based on treatment time 
# http://www.rapidtables.com/web/color/RGB_Color.htm
#https://github.com/LeahBriscoe/AdvancedHeatmapTutorial/blob/master/AdvancedHeatmapTutorial.R
treatment_times <- c(0,2,4,6,8,12,16,24)
treatment_colours_options <- c("purple","orange","green","black", "red", "pink", "blue", "grey")

#### Clustering & heatmap for DE genes with FDR < 0.01
## using logrpkm table that has been filtered for DG with FDR<0.01
dge_descr <- fread('dge_descr_table_time0_additionalcontrasts_batchcorrection_9April (1)', header = TRUE)
DE_filter <- dge_descr %>% #requires dplyr library
  filter(FDR < 0.01)

logrpkm_table_DE <- logrpkm_table[logrpkm_table$Geneid %in% DE_filter$gene_id] #filter 
logrpkm_table_DE.mat <- as.matrix(logrpkm_table_DE, rownames = "Geneid")

#Heatmap: DEG with FDR<0.01; correlation distance matrix and "complete" hierarchical clustering
png(file="Heatmap_DE_genes_FDR<0.01_logrpkm_forproject.png", width= 800, height= 750)
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

#Identifying which clusters the genes fall into
hc <- as.hclust(hm$rowDendrogram)
groups <- cutree(hc, k = 8)
groups <- as.data.table(groups, keep.rownames = "Geneid")

#loading a table which contains the gene ids and their descriptions
gene_id_descr_table <- fread('gene_id_desc_table (1)')
gene_id_desc_table_cut <- gene_id_descr_table[gene_id_descr_table$Geneid %in% groups$Geneid]

groups <- merge(groups, gene_id_desc_table_cut, by = "Geneid") #table with geneid, cluster no. & description
groups <- groups[order(groups)]
write.table(groups, file="clusters_scaled_8clusters.txt", row.names = FALSE, sep= '\t', quote= FALSE)

#Heatmap for all the genes - correlation distance matrix and default for hclust = "complete" method
logrpkm_table.mat <- as.matrix(logrpkm_table, rownames = "Geneid")
png(file="Heatmap_genes_logrpkm.png", width= 800, height= 750)
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
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/')
dge_descr <- fread('dge_descr_table_time0_additionalcontrasts_batchcorrection_9April (1)', header = TRUE)
mat <- dcast(data = dge_descr, gene_id ~ contrast, value.var = 'logFC')
mat <- as.matrix(mat, rownames = "gene_id")

png(file="Heatmap_DE_genes_logFC_test.png", width= 800, height= 750)
par(cex.axis = 1.5)
heatmap.2(mat, col=brewer.pal(11,"RdBu"),
          distfun = function(mat) as.dist(1-cor(t(mat))),
          main="Z-score for logFC of DEG at different time points",
          scale="row", cexCol=1.5, 
          key = TRUE, keysize = 1.2,  key.title = "Colour Key",
          density.info = "none", trace="none", labRow = TRUE, 
          margins = c(10,5), xlab = "Contrasts", ylab = "Genes", lwid = c(5,12), lhei = c(3,15))
dev.off()


#For future reference
#gained some insight from this tutorial: https://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html
#https://jcoliver.github.io/learn-r/006-heatmaps.html 

