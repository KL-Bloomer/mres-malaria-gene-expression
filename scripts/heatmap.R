#Heatmap and identification of gene clusters
rm(list=ls())
#Packages required
library(data.table)
library(gplots)
library(RColorBrewer)
library(dplyr)

#input and output files
logrpkm_table <- fread('logrpkm_table')


#Use rpkm table to cluster all the genes 
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/featureCounts/')
logrpkm_table <- fread('logrpkm_table')

# Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)

# Set up colour vector for time variable
setwd('~/MRes_Malaria_2021/git_repos/mres-malaria-gene-expression/')
ss_min_outliers <- fread("sample_sheet.tsv") #created from y$samples
ss_min_outliers <- ss_min_outliers[Outliers == TRUE, ]
ss_min_outliers[, Time := sprintf('%.2d', Time)]
ss_min_outliers$Time <- as.factor(ss_min_outliers$Time)


col.cell <- c("purple","orange","green","black", "red", "pink", "blue", "grey")[ss_min_outliers$Time] #need to run script_for_RNAseq_plots.R for y

# a function to assign colors based on treatment time 
# http://www.rapidtables.com/web/color/RGB_Color.htm
#https://github.com/LeahBriscoe/AdvancedHeatmapTutorial/blob/master/AdvancedHeatmapTutorial.R
treatment_times <- c(0,2,4,6,8,12,16,24)
treatment_colours_options <- c("purple","orange","green","black", "red", "pink", "blue", "grey")

#### Clustering & heatmap for DE genes with FDR < 0.01
## using logrpkm table that has been filtered for DG with FDR<0.01
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/')
dge_descr <- fread('dge_descr_table_time0_additionalcontrasts_batchcorrection_9April (1)', header = TRUE)
DE_filter <- dge_descr %>% #requires dplyr library
  filter(FDR < 0.01)

logrpkm_table_DE <- logrpkm_table[logrpkm_table$Geneid %in% DE_filter$gene_id] #filter 
logrpkm_table_DE.mat <- as.matrix(logrpkm_table_DE, rownames = "Geneid")

#Heatmap: DEG with FDR<0.01; correlation distance matrix and "complete" hierarchical clustering
png(file="AAtestHeatmap_DE_genes_FDR<0.01_logrpkm_forproject.png", width= 800, height= 750)
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
groups <- cutree(hc, k = 6)
groups <- as.data.table(groups, keep.rownames = "Geneid")

#loading a table which contains the gene ids and their descriptions
gene_id_descr_table <- fread('gene_id_desc_table (1)')
gene_id_desc_table_cut <- gene_id_descr_table[gene_id_descr_table$Geneid %in% groups$Geneid]

groups <- merge(groups, gene_id_desc_table_cut, by = "Geneid") #table with geneid, cluster no. & description
groups <- groups[order(groups)]
write.table(groups, file="clusters_scaled_6clusters.txt", row.names = FALSE, sep= '\t', quote= FALSE)

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

###Testing gene cluster from heatmap object 

#Now assign to each gene the cluster it belongs to cluster_scaled_cluster8
#is a data.table with column Geneid and cluster ID

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



