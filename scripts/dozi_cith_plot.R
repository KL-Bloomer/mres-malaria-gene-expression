#DOZI/CITH plot

dozi_file <- snakemake@input[['dozi_file']]
clust <- snakemake@input[['clust']]
ss_file <- snakemake@input[['sample_sheet']]
logrpkm_table <- snakemake@input[['logrpkm_table']]
rbp_target_plot <- snakemake@output[['rbp_target_plot']]

library(data.table)
library(ggplot2)
library(ggrepel)

dozi <- fread(dozi_file)
dozi[, gene_id := sprintf('%s0', PBANKA)]
dozi <- dozi[, c("CITH", "DOZI", 'DOZI/CITH', 'gene_id')]
head(dozi)

cido <- dozi[gene_id %in% dozi[(DOZI == 1) & (CITH == 1)]$gene_id]
cido[, target := 'CITH/DOZI']
cido <- cido[, c("gene_id","target")]

DOZI <- dozi[gene_id %in% dozi[(DOZI == 1) & (CITH != 1)]$gene_id]
DOZI[, target := 'DOZI']
DOZI <- DOZI[, c("gene_id","target")]

cith <- dozi[gene_id %in% dozi[(CITH == 1) & (DOZI != 1)]$gene_id]
cith[, target := 'CITH']
cith <- cith[, c("gene_id","target")]

rbp <- rbindlist(list(DOZI, cith, cido))


#Merge table of clusters to the table of AP2 targets
clusters <- fread(clust)
clusters <- clusters[, c("gene_id","groups")]
setnames(clusters, "groups", "cluster_id")
rbp_clust <- merge(clusters, rbp, by = "gene_id", all=TRUE)

### remove NA values
genes <- na.omit(rbp_clust)

logrpkm_table_long <- fread(logrpkm_table)
ss <- fread(ss_file)
ss <- ss[Outliers == FALSE,]

genes <- merge(logrpkm_table_long, genes, by= 'gene_id')
genes <- merge(genes, ss[, list(library_id, Time)], by= 'library_id')

avg <- genes[, list(logrpkm = mean(logrpkm), sd= sd(logrpkm), ngenes= length(unique(.SD$gene_id))),
                           by= list(Time, target)]

# Labels for RNA binding proteins names & number of genes:
avg[, rbp_ngenes := paste(' ', target, ' | N = ', ngenes, sep = " ")]
labels <- avg[, list(Time= max(Time)), by= list(rbp_ngenes)]
labels <- merge(labels, avg, by= c('rbp_ngenes', 'Time'))

##Creating plot showing the avr expression of genes
gg <- ggplot(data= avg, aes(x= Time, y= logrpkm, by = target, colour = target)) +
  geom_errorbar(aes(ymin= logrpkm - sd, ymax= logrpkm + sd), width=.1, size = 0.5, position=position_dodge(.6)) +
  geom_line(size =0.5, position = position_dodge(0.6)) +
  geom_point(size = 1.5, position = position_dodge(0.6))+
  ggtitle("Temporal changes in the average normalised gene expression of
          genes bound by translational repressor proteins") +
  xlab("Time (hr)") +
  ylab("Average normalised gene expression (log2 rpkm)")+
  ylim(0, 13)+
  theme_linedraw()

gg <- gg +
  geom_text_repel(data= labels, xlim= c(max(labels$Time), NA), aes(x = Time, y=logrpkm, label= rbp_ngenes), size = 4, inherit.aes = FALSE, hjust= 1, segment.colour= 'grey57') +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.4), add = 0), breaks= c(0, 4, 8, 12, 16, 24))+
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=13, face = "bold"),
        strip.text = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust=0.5, size = 16, face = "bold"),
        legend.position = "none")+
  guides(y.sec = guide_axis())
ggsave(rbp_target_plot, width= 20, height= 12, units= 'cm')
