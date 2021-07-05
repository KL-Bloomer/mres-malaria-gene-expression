#Enrichment  AP2 & metabolic 
rm(list = ls())
library(data.table)

##AP2 target enrichment

#load AP2 gene target files and concatenate
setwd('~/MRes_Malaria_2021/data/Enrichment/')
ap2fg <- fread('AP2-FG.targets.csv', select= 'tss_id')
ap2fg[, target := 'AP2-FG']

ap2o <- fread('AP2-O.targets.csv', select= 'tss_id')
ap2o[, target := 'AP2-O']

ap2o3 <- fread('AP2-O3.targets.csv', select= 'tss_id')
ap2o3[, target := 'AP2-O3']

ap2o4 <- fread('AP2-O4.targets.csv', select= 'tss_id')
ap2o4[, target := 'AP2-O4']

ap2 <- rbindlist(list(ap2fg, ap2o, ap2o3, ap2o4))

#convert tss_id to gene_id by stripping the part after the dot
ap2[, gene_id := sub('\\..*', '', tss_id)]
ap2[, tss_id := NULL] # We don't need this column anymore

#remove any duplicates (may be multiple peaks assigned to the same gene)
ap2 <- unique(ap2) 

#Merge the table of expressed genes to the table of clusters (similar to the GO analysis) 
#and then merge to the table of AP2 targets
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/RNA-seq_tables/')
clusters <- fread('clusters_scaled_8clusters_20May.txt')
allgenes_table <- fread('gene_id_desc_table (1)')

clusters <- merge(clusters, allgenes_table, by = "gene_id", all = TRUE)
clusters <- clusters[, c("gene_id","groups")]
setnames(clusters, "groups", "cluster_id")
ap2 <- merge(clusters, ap2, by = "gene_id", all=TRUE)
n_genes <- length(unique(ap2$gene_id))
cluster_ids <- unique(ap2[!is.na(cluster_id)]$cluster_id)
target <- unique(ap2[!is.na(target)]$target)

counts <- list()
for(clst in cluster_ids) {
  for(trg in target) {
    cnt <- data.table(
      cluster_id= clst,
      target= trg,
      count_both= nrow(ap2[cluster_id == clst & target == trg]),
      count_only_cluster= nrow(ap2[cluster_id == clst & target != trg]),
      count_only_target= nrow(ap2[cluster_id != clst & target == trg])
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
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/RNA-seq_tables/')
write.table(counts, file="counts_ap2.txt", row.names = FALSE, sep= '\t', quote= FALSE)

#Metabolic pathway enrichment

#load pathways files and P.berghei and P. falciparum conversion table
setwd('~/MRes_Malaria_2021/data/Enrichment/')
pathways <- fread('MetabolicPathways.csv', select= c('PFID New Name', 'Map_Name' ))
conversion <- fread('PbergheiANKA__v__Pfalciparum3D7.genes.tsv')

pathcon <- merge(pathways, conversion, by.x = "PFID New Name", by.y = 'Pfalciparum3D7', all = TRUE)

#remove NA values  -select for Pbergehi genes which are associated with a pathway
pathcon <- na.omit(pathcon)
pathcon[, Orthogroup:= NULL]
pathcon[, 'PFID New Name':= NULL]

#Merge the table of expressed genes to the table of clusters (similar to the GO analysis) 
#and then merge to the table of AP2 targets
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/RNA-seq_tables/')
clusters <- fread('clusters_scaled_8clusters_20May.txt')
allgenes_table <- fread('gene_id_desc_table (1)')

clusters <- merge(clusters, allgenes_table, by = "gene_id", all = TRUE)
clusters <- clusters[, c("gene_id","groups")]
setnames(clusters, "groups", "cluster_id")
comb <- merge(clusters, pathcon, by.x = "gene_id", by.y = 'PbergheiANKA',all=TRUE)

n_genes <- length(unique(comb$gene_id))
cluster_ids <- unique(comb[!is.na(cluster_id)]$cluster_id)
target <- unique(comb[!is.na(Map_Name)]$Map_Name)

counts <- list()
for(clst in cluster_ids) {
  for(trg in target) {
    cnt <- data.table(
      cluster_id= clst,
      target= trg,
      count_both= nrow(comb[cluster_id == clst & target == trg]),
      count_only_cluster= nrow(comb[cluster_id == clst & target != trg]),
      count_only_target= nrow(comb[cluster_id != clst & target == trg])
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
counts_pvalue <- counts[p.value < 0.05]
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/RNA-seq_tables/')
write.table(counts, file="counts_path.txt", row.names = FALSE, sep= '\t', quote= FALSE)
