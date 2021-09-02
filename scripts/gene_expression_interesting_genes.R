### Plotting the gene expression changes over time of key genes

#Packages required
library(data.table)
library(ggplot2)
library(ggrepel)

ss_file <- snakemake@input[['sample_sheet']]
logrpkm_table <- snakemake@input[['logrpkm_table']]
gene_expression_changes_keygenes <- snakemake@output[['gene_expression_changes_keygenes']]

#Gene for Isabelle: PBANKA_1447900

genes <- data.table(
gene_id = "PBANKA_1447900")

ss <- fread(ss_file)
ss <- ss[Outliers == FALSE,]

#Use rpkm table to cluster all the genes
logrpkm_table_long <- fread(logrpkm_table)

key_genes <- merge(logrpkm_table_long, genes, by= 'gene_id')
key_genes <- merge(key_genes, ss[, list(library_id, Time)], by= 'library_id')

avg <- key_genes[, list(logrpkm= mean(logrpkm), sd= sd(logrpkm)), by= list(Time, gene_id)]


gg <- ggplot(data= avg, aes(x= Time, y= logrpkm, group = gene_id)) +
  geom_line(size = 0.7) +
  geom_errorbar(aes(ymin= logrpkm - sd, ymax= logrpkm + sd), width=.2, size=0.7) +
  ggtitle("Temporal gene expression changes of PBANKA_1447900") +
  xlab("Time (hr)") +
  ylab("Normalised expression (log2 RPKM)")+
  theme_linedraw()

gg <- gg +
  scale_x_continuous(breaks= unique(avg$Time))+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14, face = "bold"),
        strip.text = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust=0.5, size = 15, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14))
ggsave(gene_expression_changes_keygenes, width= 25, height= 15, units= 'cm')

