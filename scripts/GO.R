#topGO analysis for 8 clusters 

rm(list=ls())
library(data.table)
library(topGO)

# dge_table is the table of differential expression from edgeR. It should have columns logFC, FDR and contrast
# gaf is the annotation file downloaded from PlasmoDB as shown before

setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/RNA-seq_tables/')
clusters <- fread('clusters_scaled_8clusters_20May.txt')


allgenes_table <- fread('gene_id_desc_table (1)')

clusters <- merge(clusters, allgenes_table, by = "gene_id", all = TRUE)
clusters <- clusters[, c("gene_id","groups")]
setnames(clusters, "groups", "cluster_id")

setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/ref/')
gaf <- fread('PlasmoDB-49_PbergheiANKA_GO.gaf', skip=1)

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

setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/RNA-seq_tables/')
write.table(x= gout[p.value < 0.05], file= 'GO_table_clusters.topgo', sep= '\t', row.names= FALSE, quote= FALSE)

setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/RNA-seq_tables/')
GO <- fread('GO_table_clusters.topgo')
GO <- GO[order(GO$groups)]
write.table(x= GO, file= 'GO_table_clusters_order.topgo', sep= '\t', row.names= FALSE, quote= FALSE)

# Here we divide genes into UP and DOWN regulated, we do this for each of the contrasts. 
# I.e. we test for enrichment the up-regulated set and the down-regulated set in each contrast
up <- copy(dge)
up[, contrast := sprintf('%s::up', contrast)]
up[, is_de := ifelse(logFC > 0, FDR, 1)]
down <- copy(dge)
down[, contrast := sprintf('%s::down', contrast)]
down[, is_de := ifelse(logFC < 0, FDR, 1)]
de <- rbind(up, down)

# Here we loop through each contrast and we test each of the 3 main branches of the GO tree (BP, MF, CC)
gout <- list()
for(cntr in unique(de$contrast)) {{
  allgenes <- de[contrast == cntr]$is_de
  names(allgenes) <- de[contrast == cntr]$gene_id
  
for(go_set in c('BP', 'MF', 'CC')){{
    # Here we prepare the object for topGO
    sampleGOdata <- new("topGOdata",
                        ontology = go_set,
                        allGenes = allgenes,
                        geneSelectionFun= function(x) {{return(x < 0.01)}},
                        nodeSize = 5,
                        annot = annFUN.gene2GO, gene2GO= gene2go)
    
    # This is the actual testing
    test_results <- runTest(sampleGOdata, algorithm = "weight01", statistic = "fisher")
    got <- data.table(GenTable(sampleGOdata, p.value= test_results, topNodes= 100, numChar= 1000))
    got[, contrast := cntr]
    gout[[length(gout)+1]] <- got
  }}
}}
