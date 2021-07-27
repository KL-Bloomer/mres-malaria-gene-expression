### Plotting temporal changes in the average gene Z-score of AP2 target genes 
#Packages required
library(data.table)
library(ggplot2)
library(ggrepel)

ap2_FG <- snakemake@input[['ap2_FG']]
ap2_O <- snakemake@input[['ap2_O']]
ap2_O3 <- snakemake@input[['ap2_O3']]
ap2_O4 <- snakemake@input[['ap2_O4']]
clust <- snakemake@input[['clust']]
ss_file <- snakemake@input[['sample_sheet']]
zscore_logrpkm <- snakemake@input[['zscore_logrpkm']]
target_plot <- snakemake@output[['target_plot']]

#load AP2 gene target files and concatenate
ap2fg <- fread(ap2_FG, select= 'tss_id')
ap2fg[, target := 'AP2-FG']

ap2o <- fread(ap2_O, select= 'tss_id')
ap2o[, target := 'AP2-O']

ap2o3 <- fread(ap2_O3, select= 'tss_id')
ap2o3[, target := 'AP2-O3']

ap2o4 <- fread(ap2_O4, select= 'tss_id')
ap2o4[, target := 'AP2-O4']

ap2 <- rbindlist(list(ap2fg, ap2o, ap2o3, ap2o4))

#convert tss_id to gene_id by stripping the part after the dot
ap2[, gene_id := sub('\\..*', '', tss_id)]
ap2[, tss_id := NULL] # We don't need this column anymore

#remove any duplicates (may be multiple peaks assigned to the same gene)
ap2 <- unique(ap2)

#Merge table of clusters to the table of AP2 targets
clusters <- fread(clust)
clusters <- clusters[, c("gene_id","groups")]
setnames(clusters, "groups", "cluster_id")
ap2_clust <- merge(clusters, ap2, by = "gene_id", all=TRUE)

### remove NA values 
genes <- na.omit(ap2_clust)

ss <- fread(ss_file)
ss <- ss[Outliers == FALSE,]

#Use zscore logrpkm table to cluster all the genes
logrpkm_table_long <- fread(zscore_logrpkm)

key_genes <- merge(logrpkm_table_long, genes, by= 'gene_id')
key_genes <- merge(key_genes, ss[, list(library_id, Time)], by= 'library_id')

avg <- key_genes[, list(zscore = mean(zscore), sd= sd(zscore), ngenes= length(unique(.SD$gene_id))), 
                 by= list(Time, target)]

#Add a printable panel title
avg[, panel_title := paste('Transcription factor', target, ' | N = ', ngenes, sep = " ")]
avg[, Time := as.numeric(as.character(Time))]

##Creating plot showing the avr zscore of genes 
gg <- ggplot(data= avg, aes(x= Time, y= zscore, by = Time)) +
  geom_line() +
  geom_point(size = 1.0)+
  geom_errorbar(aes(ymin= zscore - sd, ymax= zscore + sd), width=.2, position=position_dodge(.9)) +
  facet_wrap(~panel_title, nrow=2) +
  ggtitle("Temporal changes in average gene Z-score of AP2 transcription factor target genes") +
  xlab("Time (hr)") +
  ylab("Average gene Z-score")+
  theme_linedraw()
gg <- gg +
  scale_x_continuous(breaks= c(0, 4, 8, 12, 16, 24))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12, face = "bold"),
        strip.text = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust=0.5, size = 14, face = "bold"))+
  guides(y.sec = guide_axis())
ggsave(target_plot, width= 25, height= 15, units= 'cm')
