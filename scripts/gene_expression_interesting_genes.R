### Plotting the gene expression changes over time of key genes

#Packages required
library(dendextend)
library(data.table)
library(ggplot2)


interesting_genes <- snakemake@input[['interesting_genes']]
gaf <- snakemake@input[['gaf']]
logrpkm_table <- snakemake@input[['logrpkm_table']]
gene_expression_changes_keygenes <- snakemake@output[['gene_expression_changes_keygenes']]

#Instead of plotting the PBANKA ids it would be more useful to 
#plot the gene names

genes <- fread(Interesting_genes)
gaf <- fread(cmd= "grep -v '^!' PlasmoDB-49_PbergheiANKA_GO.gaf", select= c(2, 3, 10), col.names= c('Geneid', 'gene_name', 'description'))
gaf <- unique(gaf)

#Take care that not all genes have a name so you should fill in any 
#missing values semi-manually. Then add these gene names to your dataframe 
#and from then on work with "gene_name" instead of "Geneid"

genes <- merge(genes, gaf, by= 'Geneid', all.x= TRUE)
set(genes, i = 11L, 3L, "Cap93")
set(genes, i=11L, 4L, "oocyst capsule protein Cap93")
set(genes, i=18L, 3L, "WARP")
set(genes, i=18L, 4L, "von Willebrand factor A domain-related protein")

### Interesting genes
selected <- c('PBANKA_0800500', 'PBANKA_1319500', 'PBANKA_1300700', 'PBANKA_1035200', 
              'PBANKA_0417600', 'PBANKA_1315300', 'PBANKA_0204500', 'PBANKA_1037800',
              'PBANKA_1228900', 'PBANKA_0515000', 'PBANKA_0514900', 'PBANKA_1436600',
              'PBANKA_0402600', 'PBANKA_0905900', 'PBANKA_1363700','PBANKA_0620600',
              'PBANKA_1414900', 'PBANKA_1217700', 'PBANKA_1312700', 'PBANKA_1415700',
              'PBANKA_1001800', 'PBANKA_1302800', 'PBANKA_1436100', 'PBANKA_0905200',
              'PBANKA_0615200', 'PBANKA_0412900', 'PBANKA_1227400', 'PBANKA_1432200',
              'PBANKA_0314200', 'PBANKA_1301300')


slct <- logrpkm_table[logrpkm_table$Geneid %in% selected]
slct <- merge(slct, genes, by= 'Geneid')
slct$Time <- slct [, gsub("Ap20-GFP-SFC-4h-R[0-9]", "Ap20 4h", slct$Time)]
slct$Time <- slct [, gsub("Ap20-GFP-SFC-16h-R[0-9]", "Ap20 16h", slct$Time)]
slct$Time <- slct [, gsub("Ap20-GFP-SFC-8h-R[0-9]", "Ap20 8h", slct$Time)]
slct$Time <- slct[, gsub("GFPcon-4hr-R1", "GFPcon 4h", slct$Time)]
slct$Time <- slct [, gsub("GFPcon-16hr-R[0-9]", "GFPcon 16h", slct$Time)]
slct$Time <- slct [, gsub("GFPcon-8hr-R1_S9", "GFPcon 8h", slct$Time)]
slct$Time <- slct [, gsub("GFPcon-8hr-R2", "GFPcon 8h", slct$Time)]
slct$Time <- slct [, gsub("RM-9_S416", "RM9 24h", slct$Time)]
slct$Time <- slct [, gsub("RM-8_S415", "RM8 24h", slct$Time)]
slct$Time <- slct [, gsub("RM-7_S414", "RM7 24h", slct$Time)]
slct$Time <- slct [, gsub("RM-6_S413", "RM6 12h", slct$Time)]
slct$Time <- slct [, gsub("RM-5_S412", "RM5 24h", slct$Time)]
slct$Time <- slct [, gsub("RM-4_S411", "RM4 6h", slct$Time)]
slct$Time <- slct [, gsub("RM-3_S410", "RM3 4h", slct$Time)]
slct$Time <- slct [, gsub("RM-2_S409", "RM2 2h", slct$Time)]
slct$Time <- slct [, gsub("RM-1_S408", "RM1 12h", slct$Time)]
slct$Time <- slct [, gsub("F-.*", "SRR526055_RFP 0h", slct$Time)]
slct$Time <- slct [, gsub("F[1-3].*", "21725_RFP 0h", slct$Time)]

slct <- separate(data = slct, col = Time, into = c("Strain", "Time"), sep = "([ ])")
slct$Time <- slct [, gsub("h", "", slct$Time)]
#slct$Time <- as.numeric(slct$Time)
#slct$Time <- sprintf('%02d', slct$Time)
slct[, label := paste(Time, Strain)]

average <- slct %>%                            # name of the dataset
  group_by(Time, gene_name) %>%               # grouping the data 
  summarize(m = mean(logrpkm)) %>%      # calculating the mean
  ungroup()                      # ungroup the data

#order genes by shape of expression profile
average <- as.data.table(average)
average$Time <- as.numeric(average$Time)
ave_mat <- dcast(data= average, gene_name ~ Time, value.var= 'm')
ave_mat <- as.data.table(ave_mat)
class(ave_mat)
ave_mat <- as.matrix(ave_mat, rownames = "gene_name")

str(ave_mat)
corr_mat <- cor(t(ave_mat)) #creating a correlation 
distmat_cor <- as.dist(1-corr_mat) # determining distance from corr matrix

dendro <- as.dendrogram(hclust(distmat_cor))
# Ladderize could give slightly nicer ordering. 
# You need library(dendextend)
#reorganizes the internal structure of the tree to 
#get the ladderized effect when plotted.
dendro <- ladderize(dendro)

# This is the new order and we refactor the datasets accordingly
gene_order <- labels(dendro)
slct$gene_name <- factor(slct$gene_name, gene_order)
average$gene_name <- factor(average$gene_name, gene_order)

slct$Time <- as.numeric(slct$Time)

test <- ggplot(data = slct, aes(x=Time, y = Logrpkm_counts, group =1)) +
  geom_point(size = 0.5) + 
  facet_wrap(~gene_name) +
  theme(axis.text.x = element_text(angle = 90))+
  ggtitle("Changes in gene expression of key genes over time") +
  xlab("Time (hr)") +
  ylab("Normalised expression (log2 rpkm)") +
  theme(axis.text.x = element_text(size = 7), plot.margin=unit(c(1.5,7,1.5,1.5),"cm"))

test + 
  geom_line(data = average, aes(x = Time, y=m, group = 1), colour = "blue")+
  guides(y.sec = guide_axis()) 

ggsave(gene_expression_changes_keygenes, width= 30, height= 20, units= 'cm')
