####HEATMAP

# Get some nicer colours

mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for time variable
col.cell <- c("purple","orange","green","black", "red", "pink", "blue", "grey")[y$samples$group]


# a function to assign colors based on treatment time 
# http://www.rapidtables.com/web/color/RGB_Color.htm
#https://github.com/LeahBriscoe/AdvancedHeatmapTutorial/blob/master/AdvancedHeatmapTutorial.R
treatment_times <- c(0,2,4,6,8,12,16,24)
treatment_colours_options <- c("purple","orange","green","black", "red", "pink", "blue", "grey")
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/')

 # determining distance from corr matrix
png(file="Heatmap_genes_logrpkm.png", width= 800, height= 750)

#"complete" hclust = default and distance matrix = correlation matrix - all the genes
par(cex.main=1.5)
heatmap.2(logrpkm_table.mat, col=brewer.pal(11,"RdBu"),
          distfun = function(logrpkm_table.mat) as.dist(1-cor(t(logrpkm_table.mat))),
         main="Heatmap showing the logrpkm of all the genes",
         ColSideColors=col.cell,scale="row", cexCol=1.5,
         key = TRUE, keysize = 1.2,  key.title = "Colour Key",
         density.info = "none", trace="none", labRow = TRUE, 
         margins = c(15,7), xlab = "Library_id", ylab = "Genes", lwid = c(5,15), lhei = c(3,15))
legend("topright",legend=paste(treatment_times,"h"),fill=treatment_colours_options,cex=0.6)
dev.off()

png(file="Heatmap_DE_genes_FDR<0.01_logrpkm.png", width= 800, height= 750)
# Heatmap for the DE genes with FDR < 0.01
heatmap.2(logrpkm_table_DE.mat, col=brewer.pal(11,"RdBu"),
          distfun = function(logrpkm_table_DE.mat) as.dist(1-cor(t(logrpkm_table_DE.mat))),
          main="Heatmap showing the logrpkm of DE genes with FDR < 0.01",
          ColSideColors=col.cell,scale="row", cexCol=1.5, 
          key = TRUE, keysize = 1.2, key.title = NULL, 
          density.info = "none", trace="none", labRow = TRUE, 
          margins = c(17,7), xlab = "Library_id", ylab = "Genes", lwid = c(5,15), lhei = c(3,15))
legend("topright",legend=paste(treatment_times,"h"),fill=treatment_colours_options,cex=0.6)
dev.off()

png(file="Heatmap_DE_genes_logFC.png", width= 800, height= 750)
## Heatmap for all the DE genes showing the logFC 
## note to self, be careful of mat
heatmap.2(mat, col=brewer.pal(11,"RdBu"),
          distfun = function(mat) as.dist(1-cor(t(mat))),
          main="Heatmap showing the logFC of DE genes",
          #ColSideColors=col.cell,
          scale="row", cexCol=1.5, 
          key = TRUE, keysize = 1.2,  key.title = "Colour Key",
          density.info = "none", trace="none", labRow = TRUE, 
          margins = c(10,5), xlab = "Contrasts", ylab = "Genes", lwid = c(5,15), lhei = c(3,15))
dev.off()








 
##This seems to work very nicely
 
library(ggplot2)
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/')
png(file="High_var_genes.heatmap_testing.png")
heatmap.plot <- ggplot(data = logcounts, aes(x = variable, y = accession)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2() +
  theme(axis.text.y = element_text(size = 6))


#used this tutorial: https://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html
#https://jcoliver.github.io/learn-r/006-heatmaps.html 

