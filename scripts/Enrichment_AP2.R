#Enrichment  AP2 & metabolic
rm(list = ls())
library(data.table)

##AP2 target enrichment

#load AP2 gene target files and concatenate
setwd('~/MRes_Malaria_2021/data/Enrichment_files/')
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
#check that rows Sum to 5245
rowSums(counts[ , c(3,4,5, 6)])
rowSums(counts[ , c(3,4)])
colSums(counts[ , c(3,4)])
colSums(counts[ , 5])
test <- ap2_clust[cluster_id == 3 & target != "AP2-O3"]
test <- ap2_clust[cluster_id != 7 & target == "AP2-O3"]
ap2_clust[cluster_id != 8 & target == "AP2-O3"]
ap2_clust[cluster_id == 7 & target == "FALSE"]

counts_p5 <- counts[p.value < 0.05]
counts_p1 <- counts[p.value < 0.01]
setwd('~/MRes_Malaria_2021/git_repos/mres-malaria-gene-expression/Tables/')
write.table(counts_p1, file="counts_ap2_p01.txt", row.names = FALSE, sep= '\t', quote= FALSE)


#Metabolic pathway enrichment
#The "universe" of genes should be all the genes in the metabolic pathway table rather
#all the genes in the genome.

#load pathways files and P.berghei and P. falciparum conversion table
setwd('~/MRes_Malaria_2021/data/Enrichment_files/')
pathways <- fread('MetabolicPathways.csv', select= c('PFID New Name', 'Map_Name' ))
pathways <- unique(pathways)
conversion <- fread('PbergheiANKA__v__Pfalciparum3D7.genes.tsv')
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
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/RNA-seq_tables/')
clusters <- fread('clusters_scaled_8clusters_20May.txt')

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
counts_p01 <- counts[p.value < 0.01]
counts_p05 <- counts[p.value < 0.05]
setwd('~/MRes_Malaria_2021/git_repos/mres-malaria-gene-expression/Tables/')
write.table(counts_p01, file="counts_path_p01.txt", row.names = FALSE, sep= '\t', quote= FALSE)

## enrichment for conoid/apical proteins

#load AP2 gene target files and concatenate
setwd('~/MRes_Malaria_2021/data/Enrichment_files/')
conoid <- fread('conoid_analysis.csv', header = TRUE, select= c('gene_id'))
#
conoid[, target := 'Conoid/apical']

#remove any duplicates (may be multiple peaks assigned to the same gene)
conoid <- unique(conoid)

#Merge the table of expressed genes to the table of clusters (similar to the GO analysis)
#and then merge to the table of AP2 targets
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/RNA-seq_tables/')
clusters <- fread('clusters_scaled_8clusters_20May.txt')
allgenes_table <- fread('gene_id_desc_table (1)')

clusters <- merge(clusters, allgenes_table, by = "gene_id", all = TRUE)
clusters <- clusters[, c("gene_id","groups")]
setnames(clusters, "groups", "cluster_id")
conoid_clust <- merge(clusters, conoid, by = "gene_id", all=TRUE)
conoid_clust[is.na(conoid_clust)] <- FALSE

n_genes <- length(unique(conoid_clust$gene_id))
cluster_ids <- unique(conoid_clust$cluster_id)
target <- unique(conoid_clust$target)

counts <- list()
for(clst in cluster_ids) {
  for(trg in target) {
    count_both <- length(unique(conoid_clust[cluster_id == clst & target == trg]$gene_id))
    cnt <- data.table(
      cluster_id= clst,
      target= trg,
      count_both= count_both,
      count_only_cluster= length(unique(conoid_clust[cluster_id == clst]$gene_id)) - count_both,
      count_only_target= length(unique(conoid_clust[target == trg]$gene_id)) - count_both
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
#check that rows Sum to 5245
rowSums(counts[ , c(3,4,5, 6)])
rowSums(counts[ , c(3,4)])

setwd('~/MRes_Malaria_2021/git_repos/mres-malaria-gene-expression/Tables/')
write.table(counts, file="counts_conoid.txt", row.names = FALSE, sep= '\t', quote= FALSE)

