#Enrichment  AP2 & metabolic
library(data.table)

##AP2 target enrichment
ap2_FG <- snakemake@input[['ap2_FG']]
ap2_O <- snakemake@input[['ap2_O']]
ap2_O3 <- snakemake@input[['ap2_O3']]
ap2_O4 <- snakemake@input[['ap2_O4']]
clust <- snakemake@input[['clust']]
gene_id_table <- snakemake@input[['gene_id_table']]
path <- snakemake@input[['pathways']]
conversion_table <- snakemake@input[['conversion']]
AP2_table <- snakemake@output[['AP2_table']]
path_table <- snakemake@output[['path_table']]

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

#Merge the table of expressed genes to the table of clusters (similar to the GO analysis)
#and then merge to the table of AP2 targets
clusters <- fread(clust)
allgenes_table <- fread(gene_id_table)

clusters <- merge(clusters, allgenes_table, by = "gene_id", all = TRUE)
clusters <- clusters[, c("gene_id","groups")]
setnames(clusters, "groups", "cluster_id")
ap2_clust <- merge(clusters, ap2, by = "gene_id", all=TRUE)
ap2_clust[is.na(ap2_clust)] <- FALSE

n_genes <- length(unique(ap2_clust$gene_id))
cluster_ids <- unique(ap2_clust$cluster_id)
target <- unique(ap2_clust$target)

counts <- list()
for(clst in cluster_ids) {
  for(trg in target) {
    count_both <- length(unique(ap2_clust[cluster_id == clst & target == trg]$gene_id))
    cnt <- data.table(
      cluster_id= clst,
      target= trg,
      count_both= count_both,
      count_only_cluster= length(unique(ap2_clust[cluster_id == clst]$gene_id)) - count_both,
      count_only_target= length(unique(ap2_clust[target == trg]$gene_id)) - count_both
    )
    cnt[, count_none := n_genes - (count_both + count_only_cluster + count_only_target)]
    counts[[length(counts) + 1]] <- cnt
  }
}
counts <- rbindlist(counts)

counts[, p.value := NA]
for(i in 1:nrow(counts)) {
  x <- counts[i]
  p <- fisher.test(matrix(c(x$count_both, x$count_only_cluster, x$count_only_target, x$count_none), nrow=2), conf.int= FALSE)$p.value
  counts$p.value[i] <- p
}
# not sure if this is relevant but let's adjust for multiple testing:
counts[, fdr := p.adjust(p.value, method= 'fdr')]
unts[ , c(3,4)])

write.table(counts, file= AP2_table, row.names = FALSE, sep= '\t', quote= FALSE)

#Metabolic pathway enrichment
#The "universe" of genes should be all the genes in the metabolic pathway table rather
#all the genes in the genome.

#load pathways files and P.berghei and P. falciparum conversion table
pathways <- fread(path, select= c('PFID New Name', 'Map_Name' ))
pathways <- unique(pathways)
conversion <- fread(conversion_table)
conversion[, Orthogroup:= NULL]
conversion <- conversion[!duplicated(conversion$PbergheiANKA)]
conversion <- conversion[!duplicated(conversion$Pfalciparum3D7)]

pathcon <- merge(pathways, conversion, by.x = "PFID New Name", by.y = 'Pfalciparum3D7', all = TRUE)

#remove NA values  -select for Pbergehi genes which are associated with a pathway
pathcon <- na.omit(pathcon)
pathcon[, 'PFID New Name':= NULL]
setnames(pathcon, "Map_Name", "target")

#Merge the table of expressed genes to the table of clusters (similar to the GO analysis)
#and then merge to the table of metabolic pathways
clusters <- fread(clust)
clusters <- clusters[, c("gene_id","groups")]
setnames(clusters, "groups", "cluster_id")
comb <- merge(clusters, pathcon, by.x = "gene_id", by.y = 'PbergheiANKA',all=TRUE)
comb[is.na(comb)] <- FALSE

n_genes <- length(unique(comb$gene_id))
cluster_ids <- unique(comb$cluster_id)
target <- unique(comb$target)

counts <- list()
for(clst in cluster_ids) {
  for(trg in target) {
    count_both <- length(unique(comb[cluster_id == clst & target == trg]$gene_id))
    cnt <- data.table(
      cluster_id= clst,
      target= trg,
      count_both= count_both,
      count_only_cluster= length(unique(comb[cluster_id == clst]$gene_id)) - count_both,
      count_only_target= length(unique(comb[target == trg]$gene_id)) - count_both
    )
    cnt[, count_none := n_genes - (count_both + count_only_cluster + count_only_target)]
    counts[[length(counts) + 1]] <- cnt
  }
}
counts <- rbindlist(counts)

counts[, p.value := NA]
for(i in 1:nrow(counts)) {
  x <- counts[i]
  p <- fisher.test(matrix(c(x$count_both, x$count_only_cluster, x$count_only_target, x$count_none), nrow= 2), conf.int= FALSE)$p.value
  counts$p.value[i] <- p
}

# not sure if this is relevant but let's adjust for multiple testing:
counts[, fdr := p.adjust(p.value, method= 'fdr')]
write.table(counts, file=path_table, row.names = FALSE, sep= '\t', quote= FALSE)

