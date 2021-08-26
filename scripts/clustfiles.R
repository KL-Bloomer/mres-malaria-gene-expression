##Creating gff files for each cluster and the negative controls

library(tidyr)
library(stringr)
library(dplyr)
library(data.table)

clust <- snakemake@input[['clust']]
dge_table_file <- snakemake@input[['dge_table']]
GFF_file <- snakemake@input[['gff']]
clst_pos_files <- snakemake@output[['clst_pos']]
clst_neg_files <- snakemake@output[['clst_neg']]
clst_out_file <- snakemake@output[['clst_out']]

clusters <- fread(clust)
dge_table <- fread(dge_table_file)
GFF<- fread(cmd=paste('grep -v "^#"', GFF_file))

gff <- separate(data = GFF, col = V9, into = c("Geneid_Feature"), sep = "([;])")
gff <- gff %>%
  mutate_at("Geneid_Feature", str_replace, "ID=", "")
gff[, V4 := ifelse(V7 == '-', V5, V4)]
gff[, V5 := V4]
setnames(gff, "Geneid_Feature", "gene_id")

dge_table <- dge_table[, 'gene_id', with= FALSE]
dge_table <- unique(dge_table)
dge_table_gff <- gff[gff$gene_id %in% dge_table$gene_id]

cluster_ids <- unique(clusters$groups)
for (clst in cluster_ids) {
  clst_genes <- clusters[groups == clst]$gene_id
  sample <- gff[gene_id %in% clst_genes]
  neg_clst_genes <- dge_table_gff[!dge_table_gff$gene_id %in% clst_genes]
#Make a filename for the positive set
  clst_pos_file <- paste('meme/clst_pos', clst, '.gff', sep= '')
#Check the filename we created is in the list of outputs
  stopifnot(clst_pos_file %in% clst_pos_files)
  write.table(sample, file= clst_pos_file, row.names = FALSE, col.names = FALSE, sep= '\t', quote= FALSE)
#Make a filename for the negative set
  clst_neg_file <- paste('meme/clst_neg', clst, '.gff', sep= '')
  stopifnot(clst_neg_file %in% clst_neg_files)
  write.table(neg_clst_genes, file= clst_neg_file, row.names = FALSE, col.names = FALSE, sep= '\t', quote= FALSE)
}
