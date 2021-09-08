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
dozi_file <- snakemake@input[['dozi_file']]
conoid_file <- snakemake@input[['conoid_file']]
AP2_table <- snakemake@output[['AP2_table']]
path_table <- snakemake@output[['path_table']]
conoid_table <- snakemake@output[['conoid_table']]
cith_table <- snakemake@output[['cith_table']]
dozi_table <- snakemake@output[['dozi_table']]
cith_dozi_table <- snakemake@output[['cith_dozi_table']]

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
counts[, odds.ratio := NA]
for(i in 1:nrow(counts)) {
x <- counts[i]
ft <- fisher.test(matrix(c(x$count_both, x$count_only_cluster, x$count_only_target, x$count_none), nrow=2), conf.int= FALSE)
counts$p.value[i] <- ft$p.value
counts$odds.ratio[i] <- ft$estimate
}
# not sure if this is relevant but let's adjust for multiple testing:
counts[, fdr := p.adjust(p.value, method= 'fdr')]
#counts_p5 <- counts[p.value < 0.05]
#counts_p1 <- counts[p.value < 0.01]


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
counts[, odds.ratio := NA]
for(i in 1:nrow(counts)) {
x <- counts[i]
ft <- fisher.test(matrix(c(x$count_both, x$count_only_cluster, x$count_only_target, x$count_none), nrow=2), conf.int= FALSE)
counts$p.value[i] <- ft$p.value
counts$odds.ratio[i] <- ft$estimate
}

# not sure if this is relevant but let's adjust for multiple testing:
counts[, fdr := p.adjust(p.value, method= 'fdr')]
#counts_p01 <- counts[p.value < 0.01]
#counts_p05 <- counts[p.value < 0.05]
write.table(counts, file=path_table, row.names = FALSE, sep= '\t', quote= FALSE)

############################### enrichment for conoid/apical proteins ############################
conoid <- fread(conoid_file, header = TRUE, select= c('gene_id'))
conoid[, target := 'Conoid/apical']

#remove any duplicates (may be multiple peaks assigned to the same gene)
conoid <- unique(conoid)

#Merge the table of expressed genes to the table of clusters (similar to the GO analysis) and then to the conoid table
clusters <- fread(clust)
allgenes_table <- fread(gene_id_table)
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
counts[, odds.ratio := NA]
for(i in 1:nrow(counts)) {
x <- counts[i]
ft <- fisher.test(matrix(c(x$count_both, x$count_only_cluster, x$count_only_target, x$count_none), nrow=2), conf.int= FALSE)
counts$p.value[i] <- ft$p.value
counts$odds.ratio[i] <- ft$estimate
}

# not sure if this is relevant but let's adjust for multiple testing:
counts[, fdr := p.adjust(p.value, method= 'fdr')]
#check that rows Sum to 5245
rowSums(counts[ , c(3,4,5, 6)])
rowSums(counts[ , c(3,4)])
write.table(counts, file=conoid_table, row.names = FALSE, sep= '\t', quote= FALSE)

#######################################################################################

#DOZI/CITH enrichment

##CITH
dozi <- fread(dozi_file)
dozi[, gene_id := sprintf('%s0', PBANKA)]
dozi <- dozi[, c("CITH", "DOZI", 'DOZI/CITH', 'gene_id')]

head(dozi)

cith <- dozi[gene_id %in% dozi[(CITH == 1) & (DOZI != 1)]$gene_id]
#cith <- dozi[gene_id %in% dozi[(CITH == 1)]$gene_id] #if you want genes bound by CITH or CITH and DOZI
cith[, target := 'cith']
cith <- cith[, c("gene_id","target")]

#Merge the table of expressed genes to the table of clusters (similar to the GO analysis)
#and then merge to the table of AP2 targets
clusters <- fread(clust)
allgenes_table <- fread(gene_id_table)

clusters <- merge(clusters, allgenes_table, by = "gene_id", all = TRUE)
clusters <- clusters[, c("gene_id","groups")]
setnames(clusters, "groups", "cluster_id")
cith_clust <- merge(clusters, cith, by = "gene_id", all =TRUE)

cith_clust[is.na(cith_clust)] <- FALSE
cith_clust <- unique(cith_clust)
#test <- cith_clust[!cith_clust$gene_id %in% allgenes_table$gene_id] #genes that are not in the gene_id_table
cith_clust <- cith_clust[cith_clust$gene_id %in% allgenes_table$gene_id] # remove genes not in gene_id_table

n_genes <- length(unique(cith_clust$gene_id))
cluster_ids <- unique(cith_clust$cluster_id)
target <- unique(cith_clust$target)

counts <- list()
for(clst in cluster_ids) {
  for(trg in target) {
    count_both <- length(unique(cith_clust[cluster_id == clst & target == trg]$gene_id))
    cnt <- data.table(
      cluster_id= clst,
      target= trg,
      count_both= count_both,
      count_only_cluster= length(unique(cith_clust[cluster_id == clst]$gene_id)) - count_both,
      count_only_target= length(unique(cith_clust[target == trg]$gene_id)) - count_both
    )
    cnt[, count_none := n_genes - (count_both + count_only_cluster + count_only_target)]
    counts[[length(counts) + 1]] <- cnt
  }
}
counts <- rbindlist(counts)

counts[, p.value := NA]
counts[, odds.ratio := NA]
for(i in 1:nrow(counts)) {
x <- counts[i]
ft <- fisher.test(matrix(c(x$count_both, x$count_only_cluster, x$count_only_target, x$count_none), nrow=2), conf.int= FALSE)
counts$p.value[i] <- ft$p.value
counts$odds.ratio[i] <- ft$estimate
}

# not sure if this is relevant but let's adjust for multiple testing:
counts[, fdr := p.adjust(p.value, method= 'fdr')]

#check that rows Sum to 5245
rowSums(counts[ , c(3,4,5, 6)])
rowSums(counts[ , c(3,4)])
write.table(counts, file=cith_table, row.names = FALSE, sep= '\t', quote= FALSE)

##########################################################################
#DOZI enrichment

DOZI <- dozi[gene_id %in% dozi[(DOZI == 1) & (CITH != 1)]$gene_id]
DOZI[, target := 'DOZI']
DOZI <- DOZI[, c("gene_id","target")]

#merge DOZI to clusters (merged clust to gene_id_table)
dozi_clust <- merge(clusters, DOZI, by = "gene_id", all =TRUE)

dozi_clust[is.na(dozi_clust)] <- FALSE
dozi_clust <- unique(dozi_clust)
dozi_clust <- dozi_clust[dozi_clust$gene_id %in% allgenes_table$gene_id]

n_genes <- length(unique(dozi_clust$gene_id))
cluster_ids <- unique(dozi_clust$cluster_id)
target <- unique(dozi_clust$target)

counts <- list()
for(clst in cluster_ids) {
  for(trg in target) {
    count_both <- length(unique(dozi_clust[cluster_id == clst & target == trg]$gene_id))
    cnt <- data.table(
      cluster_id= clst,
      target= trg,
      count_both= count_both,
      count_only_cluster= length(unique(dozi_clust[cluster_id == clst]$gene_id)) - count_both,
      count_only_target= length(unique(dozi_clust[target == trg]$gene_id)) - count_both
    )
    cnt[, count_none := n_genes - (count_both + count_only_cluster + count_only_target)]
    counts[[length(counts) + 1]] <- cnt
  }
}
counts <- rbindlist(counts)

counts[, p.value := NA]
counts[, odds.ratio := NA]
for(i in 1:nrow(counts)) {
x <- counts[i]
ft <- fisher.test(matrix(c(x$count_both, x$count_only_cluster, x$count_only_target, x$count_none), nrow=2), conf.int= FALSE)
counts$p.value[i] <- ft$p.value
counts$odds.ratio[i] <- ft$estimate
}

# not sure if this is relevant but let's adjust for multiple testing:
counts[, fdr := p.adjust(p.value, method= 'fdr')]

#check that rows Sum to 5245
rowSums(counts[ , c(3,4,5, 6)])
rowSums(counts[ , c(3,4)])
write.table(counts, file=dozi_table, row.names = FALSE, sep= '\t', quote= FALSE)

########################################################################

#DOZI & CITH enrichment

cido <- dozi[gene_id %in% dozi[(DOZI == 1) & (CITH == 1)]$gene_id]
cido[, target := 'CITH/DOZI']
cido <- cido[, c("gene_id","target")]

cido_clust <- merge(clusters, cido, by = "gene_id", all =TRUE)

cido_clust[is.na(cido_clust)] <- FALSE
cido_clust <- unique(cido_clust)
cido_clust <- cido_clust[cido_clust$gene_id %in% allgenes_table$gene_id]

n_genes <- length(unique(cido_clust$gene_id))
cluster_ids <- unique(cido_clust$cluster_id)
target <- unique(cido_clust$target)

counts <- list()
for(clst in cluster_ids) {
  for(trg in target) {
    count_both <- length(unique(cido_clust[cluster_id == clst & target == trg]$gene_id))
    cnt <- data.table(
      cluster_id= clst,
      target= trg,
      count_both= count_both,
      count_only_cluster= length(unique(cido_clust[cluster_id == clst]$gene_id)) - count_both,
      count_only_target= length(unique(cido_clust[target == trg]$gene_id)) - count_both
    )
    cnt[, count_none := n_genes - (count_both + count_only_cluster + count_only_target)]
    counts[[length(counts) + 1]] <- cnt
  }
}
counts <- rbindlist(counts)

counts[, p.value := NA]
counts[, odds.ratio := NA]
for(i in 1:nrow(counts)) {
x <- counts[i]
ft <- fisher.test(matrix(c(x$count_both, x$count_only_cluster, x$count_only_target, x$count_none), nrow=2), conf.int= FALSE)
counts$p.value[i] <- ft$p.value
counts$odds.ratio[i] <- ft$estimate
}

# not sure if this is relevant but let's adjust for multiple testing:
counts[, fdr := p.adjust(p.value, method= 'fdr')]

#check that rows Sum to 5245
rowSums(counts[ , c(3,4,5, 6)])
rowSums(counts[ , c(3,4)])
write.table(counts, file=cith_dozi_table, row.names = FALSE, sep= '\t', quote= FALSE)

