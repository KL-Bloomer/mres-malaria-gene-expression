### Plotting the gene expression changes over time of key genes

#Packages required
library(data.table)
library(ggplot2)
library(ggrepel)

interesting_genes <- snakemake@input[['interesting_genes']]
GAF <- snakemake@input[['GAF']]
logrpkm_table <- snakemake@input[['logrpkm_table']]
gene_expression_changes_keygenes <- snakemake@output[['gene_expression_changes_keygenes']]

#Instead of plotting the PBANKA ids it would be more useful to 
#plot the gene names
genes <- fread(interesting_genes)
gaf <- fread(cmd= paste("grep -v '^!'", GAF), select= c(2, 3, 10), col.names= c('Geneid', 'gene_name', 'description'))
gaf <- unique(gaf)

#Take care that not all genes have a name so you should fill in any 
#missing values semi-manually. Then add these gene names to your dataframe 
#and from then on work with "gene_name" instead of "Geneid"

genes <- merge(genes, gaf, by.x= 'gene_id', by.y = "Geneid", all.x = TRUE)
genes <- set(genes, i=16L, 4L, "WARP")
genes <- set(genes, i=16L, 5L, "von Willebrand factor A domain-related protein")

ss <- fread(ss)
ss <- ss[Outliers == FALSE,]

#Use rpkm table to cluster all the genes 
logrpkm_table_long <- fread(logrpkm_table)

key_genes <- merge(logrpkm_table_long, genes, by= 'gene_id')
key_genes <- merge(key_genes, ss[, list(library_id, Time)], by= 'library_id')

avg <- key_genes[, list(logrpkm= mean(logrpkm), sd= sd(logrpkm)), by= list(Time, gene_name, Functional_group)] 

# Position of gene names:
labels <- avg[, list(Time= max(Time)), by= list(gene_name)]
labels <- merge(labels, avg, by= c('gene_name', 'Time'))

gg <- ggplot(data= avg, aes(x= Time, y= logrpkm, group = gene_name)) +
  geom_line(size = 0.5) +
  geom_errorbar(aes(ymin= logrpkm - sd, ymax= logrpkm + sd), width=.2) +
  facet_wrap(~Functional_group) +
  ggtitle("Changes in gene expression of key genes over time") +
  xlab("Time (hr)") +
  ylab("Normalised expression (log2 RPKM)")

gg <- gg + 
  geom_text_repel(data= labels, xlim= c(max(labels$Time), NA), aes(x = Time, y=logrpkm, label= gene_name), inherit.aes = FALSE, hjust= 1, segment.colour= 'grey57') +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.5), add = 0), breaks= unique(avg$Time))+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        strip.text = element_text(size = 14, face="bold"))+
  guides(y.sec = guide_axis()) 
ggsave(gene_expression_changes_keygenes, width= 30, height= 20, units= 'cm')

