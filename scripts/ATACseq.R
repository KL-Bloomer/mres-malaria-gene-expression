#ATACseq analysis

gam_file <- snakemake@input[['gam_file']]
ook_file <- snakemake@input[['ook_file']]
gam_filtered <- snakemake@output[['gam_filter']]
ook_filtered <- snakemake@output[['ook_filter']]

#Packages required
library(data.table)

#load ATACseq files and filter for peaks with FC/input > 2.0 and tss_distance_rank = 1 (only 1 TSS/peak)
gam <- fread(gam_file)
gam[, summit_start := start + summit]
gam[, summit_end := summit_start + 1] 
setcolorder(gam, c('#chrom', 'summit_start', 'summit_end'))

setnames(gam, "logfc:avg", "logfc")
gam_filter <- gam[logfc>2.0]
gam_filter <- gam_filter[tss_distance_rank == 1]
#gam_filter[ , c(5:14, 16)] <- list(NULL)
write.table(gam_filter, file = gam_filtered, sep= '\t', row.names= FALSE, quote= FALSE)

ook <- fread(ook_file)
ook[, summit_start := start + summit]
ook[, summit_end := summit_start + 1] 
setcolorder(ook, c('#chrom', 'summit_start', 'summit_end'))
setnames(ook, "logfc:avg", "logfc")
ook_filter <- ook[logfc>2.0]
ook_filter <- ook_filter[tss_distance_rank == 1]
write.table(ook_filter, file = ook_filtered, sep= '\t', row.names= FALSE, quote= FALSE)

