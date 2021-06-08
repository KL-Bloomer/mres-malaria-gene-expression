#topGO analysis for 8 clusters
library(data.table)
library(topGO)

clust <- snakemake@input[['clust']]
gene_id_table <- snakemake@input[['gene_id_table']]
GAF <- snakemake@input[['GAF']]
topGO_table_clusters <- snakemake@output[['topGO_table_clusters']]

clusters <- fread(clust)
allgenes_table <- fread(gene_id_table)
clusters <- merge(clusters, allgenes_table, by = "gene_id", all = TRUE)
clusters <- clusters[, c("gene_id","groups")]
setnames(clusters, "groups", "cluster_id")

gaf <- fread(GAF, skip=1)

# Here we just reformat the gaf table into a more convenient shape
setnames(gaf, c('V2', 'V5'), c('gene_id', 'go'))
gafu <- as.matrix(unique(gaf[, list(gene_id, go)]))
gene2go <- list()
for(g in unique(gafu[, 'gene_id'])){{
  gg <- which(gafu[, 'gene_id'] == g)
  gene2go[[g]] <- gafu[gg, 'go']
}}

# Here we loop through each contrast and we test each of the 3 main branches of the GO tree (BP, MF, CC)
gout <- list()
for(id in unique(clusters$cluster_id)) {
  allgenes <- clusters$cluster_id == id
  allgenes[is.na(allgenes)] <- FALSE
  allgenes <- as.numeric(allgenes)
  names(allgenes) <- clusters$gene_id

  for(go_set in c('BP', 'MF', 'CC')){
    # Here we prepare the object for topGO
    sampleGOdata <- new("topGOdata",
                        ontology = go_set,
                        allGenes = allgenes,
                        geneSelectionFun= function(x) {return(x == 1)},
                        nodeSize = 5,
                        annot = annFUN.gene2GO, gene2GO= gene2go)

    # This is the actual testing
    test_results <- runTest(sampleGOdata, algorithm = "weight01", statistic = "fisher")
    got <- data.table(GenTable(sampleGOdata, p.value= test_results, topNodes= 100, numChar= 1000))
    got[, groups := id]
    gout[[length(gout)+1]] <- got
  }
}
# Finally, we reshape the output and write it out:
gout <- rbindlist(gout)
gout[, p.value := as.numeric(p.value)]
gout <- gout[order(p.value)]
setnames(gout, c('GO.ID', 'Term'), c('go_id', 'go_term'))

write.table(x= gout[p.value < 0.05], file= topGO_table_clusters, sep= '\t', row.names= FALSE, quote= FALSE)
