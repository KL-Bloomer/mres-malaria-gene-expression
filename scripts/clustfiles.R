#Separating genes into separate files for each cluster and creating a file for genes not within a cluster

library(tidyr)
library(stringr)
library(dplyr)
library(data.table)

clust <- snakemake@input[['clust']]
gene_id_table <- snakemake@input[['gene_id_table']]
GFF_file <- snakemake@input[['gff']]
clst1_file <- snakemake@output[['clst1']]
clst2_file <- snakemake@output[['clst2']]
clst3_file <- snakemake@output[['clst3']]
clst4_file <- snakemake@output[['clst4']]
clst5_file <- snakemake@output[['clst5']]
clst6_file <- snakemake@output[['clst6']]
clst7_file <- snakemake@output[['clst7']]
clst8_file <- snakemake@output[['clst8']]
clstnot_file <- snakemake@output[['clstnot']]


clusters <- fread(clust)
allgenes_table <- fread(gene_id_table)

clst1 <- clusters[groups == 1]
clst2 <- clusters[groups == 2]
clst3 <- clusters[groups == 3]
clst4 <- clusters[groups == 4]
clst5 <- clusters[groups == 5]
clst6 <- clusters[groups == 6]
clst7 <- clusters[groups == 7]
clst8 <- clusters[groups == 8]
clstnot <- allgenes_table[!allgenes_table$gene_id %in% clusters$gene_id]

###Isolating GFF files

GFF<- fread(cmd=paste('grep -v "^#"', GFF_file))
gff <- separate(data = GFF, col = V9, into = c("Geneid_Feature"), sep = "([;])")
gff <- gff %>%
  mutate_at("Geneid_Feature", str_replace, "ID=", "")
gff[, V4 := ifelse(V7 == '-', V5, V4)]
gff[, V5 := V4]

GFF_clst1 <- gff[gff$Geneid_Feature %in% clst1$gene_id]
setnames(GFF_clst1, "Geneid_Feature", "gene_id")
setnames(GFF_clst1, "V1", "#V1")
write.table(GFF_clst1, file=clst1_file, row.names = FALSE, sep= '\t', quote= FALSE)

#Clust2
GFF_clst2 <- gff[gff$Geneid_Feature %in% clst2$gene_id]
setnames(GFF_clst2, "Geneid_Feature", "gene_id")
setnames(GFF_clst2, "V1", "#V1")
write.table(GFF_clst2, file=clst2_file, row.names = FALSE, sep= '\t', quote= FALSE)

#Clust3
GFF_clst3 <- gff[gff$Geneid_Feature %in% clst3$gene_id]
setnames(GFF_clst3, "Geneid_Feature", "gene_id")
setnames(GFF_clst3, "V1", "#V1")
write.table(GFF_clst3, file=clst3_file, row.names = FALSE, sep= '\t', quote= FALSE)

#Clust4
GFF_clst4 <- gff[gff$Geneid_Feature %in% clst4$gene_id]
setnames(GFF_clst4, "Geneid_Feature", "gene_id")
setnames(GFF_clst4, "V1", "#V1")
write.table(GFF_clst4, file=clst4_file, row.names = FALSE, sep= '\t', quote= FALSE)

#Clust5
GFF_clst5 <- gff[gff$Geneid_Feature %in% clst5$gene_id]
setnames(GFF_clst5, "Geneid_Feature", "gene_id")
setnames(GFF_clst5, "V1", "#V1")
write.table(GFF_clst5, file=clst5_file, row.names = FALSE, sep= '\t', quote= FALSE)

#Clust6
GFF_clst6 <- gff[gff$Geneid_Feature %in% clst6$gene_id]
setnames(GFF_clst6, "Geneid_Feature", "gene_id")
setnames(GFF_clst6, "V1", "#V1")
write.table(GFF_clst6, file=clst6_file, row.names = FALSE, sep= '\t', quote= FALSE)

#Clust7
GFF_clst7 <- gff[gff$Geneid_Feature %in% clst7$gene_id]
setnames(GFF_clst7, "Geneid_Feature", "gene_id")
setnames(GFF_clst7, "V1", "#V1")
write.table(GFF_clst7, file=clst7_file, row.names = FALSE, sep= '\t', quote= FALSE)

#Clust8
GFF_clst8 <- gff[gff$Geneid_Feature %in% clst8$gene_id]
setnames(GFF_clst8, "Geneid_Feature", "gene_id")
setnames(GFF_clst8, "V1", "#V1")
write.table(GFF_clst8, file=clst8_file, row.names = FALSE, sep= '\t', quote= FALSE)

#Not in cluster
GFF_not <- gff[gff$Geneid_Feature %in% clstnot$gene_id]
setnames(GFF_not, "Geneid_Feature", "gene_id")
setnames(GFF_not, "V1", "#V1")
write.table(GFF_not, file=clstnot_file, row.names = FALSE, sep= '\t', quote= FALSE)

#### Creating neg/background file for each clusters
all_genes_gff <- gff[gff$Geneid_Feature %in% allgenes_table$gene_id]
setnames(all_genes_gff, "Geneid_Feature", "gene_id")
setnames(all_genes_gff, "V1", "#V1")

neg_clst1 <- all_genes_gff[!all_genes_gff$gene_id %in% clst1$gene_id]
setwd('~/MRes_Malaria_2021/git_repos/mres-malaria-gene-expression/Tables/')
write.table(neg_clst1, file="neg_clst1.gff", row.names = FALSE, sep= '\t', quote= FALSE)

#Clust2
neg_clst2 <- all_genes_gff[!all_genes_gff$gene_id %in% clst2$gene_id]
setwd('~/MRes_Malaria_2021/git_repos/mres-malaria-gene-expression/Tables/')
write.table(neg_clst2, file="neg_clst2.gff", row.names = FALSE, sep= '\t', quote= FALSE)

#Clust3
neg_clst3 <- all_genes_gff[!all_genes_gff$gene_id %in% clst3$gene_id]
setwd('~/MRes_Malaria_2021/git_repos/mres-malaria-gene-expression/Tables/')
write.table(neg_clst3, file="neg_clst3.gff", row.names = FALSE, sep= '\t', quote= FALSE)

#Clust4
neg_clst4 <- all_genes_gff[!all_genes_gff$gene_id %in% clst4$gene_id]
setwd('~/MRes_Malaria_2021/git_repos/mres-malaria-gene-expression/Tables/')
write.table(neg_clst4, file="neg_clst4.gff", row.names = FALSE, sep= '\t', quote= FALSE)

#Clust5
neg_clst5 <- all_genes_gff[!all_genes_gff$gene_id %in% clst5$gene_id]
setwd('~/MRes_Malaria_2021/git_repos/mres-malaria-gene-expression/Tables/')
write.table(neg_clst5, file="neg_clst5.gff", row.names = FALSE, sep= '\t', quote= FALSE)

#Clust6
neg_clst6 <- all_genes_gff[!all_genes_gff$gene_id %in% clst6$gene_id]
setwd('~/MRes_Malaria_2021/git_repos/mres-malaria-gene-expression/Tables/')
write.table(neg_clst6, file="neg_clst6.gff", row.names = FALSE, sep= '\t', quote= FALSE)

#Clust7
neg_clst7 <- all_genes_gff[!all_genes_gff$gene_id %in% clst7$gene_id]
setwd('~/MRes_Malaria_2021/git_repos/mres-malaria-gene-expression/Tables/')
write.table(neg_clst7, file="neg_clst7.gff", row.names = FALSE, sep= '\t', quote= FALSE)

#Clust8
neg_clst8 <- all_genes_gff[!all_genes_gff$gene_id %in% clst8$gene_id]
setwd('~/MRes_Malaria_2021/git_repos/mres-malaria-gene-expression/Tables/')
write.table(neg_clst8, file="neg_clst8.gff", row.names = FALSE, sep= '\t', quote= FALSE)
