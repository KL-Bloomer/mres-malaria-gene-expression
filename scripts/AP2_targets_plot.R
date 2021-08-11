### Plotting temporal changes in the average gene Z-score of AP2 target genes 
#Packages required
library(data.table)
library(ggplot2)
library(ggrepel)
library(emmeans)                  
library(lme4)                  
library(pbkrtest)
library(lmerTest)

ap2_FG <- snakemake@input[['ap2_FG']]
ap2_O <- snakemake@input[['ap2_O']]
ap2_O3 <- snakemake@input[['ap2_O3']]
ap2_O4 <- snakemake@input[['ap2_O4']]
clust <- snakemake@input[['clust']]
ss_file <- snakemake@input[['sample_sheet']]
zscore_logrpkm <- snakemake@input[['zscore_logrpkm']]
logrpkm_table <- snakemake@input[['logrpkm_table']]
AP2_FG_O3_plot <- snakemake@output[['AP2_FG_O3_plot']]
AP2_O_O4_plot <- snakemake@output[['AP2_O_O4_plot']]
AP2_O_O3_plot <- snakemake@output[['AP2_O_O3_plot']]
AP2_target_plot <- snakemake@output[['AP2_target_plot']]
lmer_out_file <- snakemake@output[['lmer_out_file']]

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

#Merge table of clusters to the table of AP2 targets
clusters <- fread(clust)
clusters <- clusters[, c("gene_id","groups")]
setnames(clusters, "groups", "cluster_id")
ap2_clust <- merge(clusters, ap2, by = "gene_id", all=TRUE)

### remove NA values 
genes <- na.omit(ap2_clust)

################################################# Plot zscore_logrpkm table ######################################
logrpkm_table_long <- fread(zscore_logrpkm)

ss <- fread(ss_file)
ss <- ss[Outliers == FALSE,]

key_genes <- merge(logrpkm_table_long, genes, by= 'gene_id')
key_genes <- merge(key_genes, ss[, list(library_id, Time)], by= 'library_id')

avg <- key_genes[, list(zscore = mean(zscore), sd= sd(zscore), ngenes= length(unique(.SD$gene_id))),
                 by= list(Time, target)]

#Add a printable panel title
avg[, panel_title := paste(' ', target, ' | N = ', ngenes, sep = " ")]
avg[, Time := as.numeric(as.character(Time))]

##Creating plot showing the avr zscore of genes
gg <- ggplot(data= avg, aes(x= Time, y= zscore, by = Time)) +
  geom_line() +
  geom_point(size = 1.0)+
  geom_errorbar(aes(ymin= zscore - sd, ymax= zscore + sd), width=.2, position=position_dodge(.9)) +
  facet_wrap(~panel_title, nrow=2) +
  ggtitle("Temporal changes in average gene Z-score of key 
          developmental transcription factor target genes") +
  xlab("Time (hr)") +
  ylab("Average gene Z-score")+
  theme_linedraw()
gg <- gg +
  scale_x_continuous(breaks= c(0, 4, 8, 12, 16, 24))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12, face = "bold"),
        strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust=0.5, size = 16, face = "bold"))+
  guides(y.sec = guide_axis())
ggsave(AP2_target_plot, width= 22, height= 15, units= 'cm')

############################################################

#Generating plots for genes regulated by one or two transcription factors 

#table with number of trxn factors and comma separated list of trxn factors
geneTarget <- genes[, list(N= length(unique(target)), 
                         target= paste(sort(target), collapse= ' & ')), by= gene_id] 


###for AP2-FG and AP2-O3 plot
key_genes <- geneTarget[gene_id %in% geneTarget[(target == "AP2-FG & AP2-O3") | (target == "AP2-FG") |
                                              (target == "AP2-O3")]$gene_id]
key_genes <- merge(logrpkm_table_long, key_genes, by= 'gene_id')
key_genes <- merge(key_genes, ss[, list(library_id, Time)], by= 'library_id')

avg <- key_genes[, list(zscore = mean(zscore), sd= sd(zscore), ngenes= length(unique(.SD$gene_id))), 
                 by= list(Time, target)]

#Add a printable panel title
avg[, panel_title := paste(' ', target, ' | N = ', ngenes, sep = " ")]
avg[, Time := as.numeric(as.character(Time))]

##Creating plot showing the avr zscore of genes 
gg <- ggplot(data= avg, aes(x= Time, y= zscore, by = Time)) +
  geom_line() +
  geom_point(size = 1.0)+
  geom_errorbar(aes(ymin= zscore - sd, ymax= zscore + sd), width=.2, position=position_dodge(.9)) +
  facet_wrap(~panel_title, nrow=2) +
  ggtitle("Temporal changes in average gene Z-score of 
          AP2-FG & AP2-O3 transcription factor target genes") +
  xlab("Time (hr)") +
  ylab("Average gene Z-score")+
  theme_linedraw()
gg <- gg +
  scale_x_continuous(breaks= c(0, 4, 8, 12, 16, 24))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12, face = "bold"),
        strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust=0.5, size = 16, face = "bold"))+
  guides(y.sec = guide_axis())
ggsave(AP2_FG_O3_plot, width= 22, height= 15, units= 'cm')

################################################################

#for AP2-O and AP2-O4 plots

key_genes <- geneTarget[gene_id %in% geneTarget[(target == "AP2-O & AP2-O4") | (target == "AP2-O") |
                                              (target == "AP2-O4")]$gene_id]
key_genes <- merge(logrpkm_table_long, key_genes, by= 'gene_id')
key_genes <- merge(key_genes, ss[, list(library_id, Time)], by= 'library_id')

avg <- key_genes[, list(zscore = mean(zscore), sd= sd(zscore), ngenes= length(unique(.SD$gene_id))),
                 by= list(Time, target)]

#Add a printable panel title
avg[, panel_title := paste(' ', target, ' | N = ', ngenes, sep = " ")]
avg[, Time := as.numeric(as.character(Time))]

##Creating plot showing the avr zscore of genes
gg <- ggplot(data= avg, aes(x= Time, y= zscore, by = Time)) +
  geom_line() +
  geom_point(size = 1.0)+
  geom_errorbar(aes(ymin= zscore - sd, ymax= zscore + sd), width=.2, position=position_dodge(.9)) +
  facet_wrap(~panel_title, nrow=2) +
  ggtitle("Temporal changes in average gene Z-score of
          AP2-O & AP2-O4 transcription factor target genes") +
  xlab("Time (hr)") +
  ylab("Average gene Z-score")+
  theme_linedraw()
gg <- gg +
  scale_x_continuous(breaks= c(0, 4, 8, 12, 16, 24))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12, face = "bold"),
        strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust=0.5, size = 16, face = "bold"))+
  guides(y.sec = guide_axis())
ggsave(AP2_O_O4_plot, width= 22, height= 15, units= 'cm')

##############################################################

### For AP2-O and AP2-O3 plots

key_genes <- geneTarget[gene_id %in% geneTarget[(target == "AP2-O & AP2-O3") | (target == "AP2-O") |
                                              (target == "AP2-O3")]$gene_id]
key_genes <- merge(logrpkm_table_long, key_genes, by= 'gene_id')
key_genes <- merge(key_genes, ss[, list(library_id, Time)], by= 'library_id')

avg <- key_genes[, list(zscore = mean(zscore), sd= sd(zscore), ngenes= length(unique(.SD$gene_id))),
                 by= list(Time, target)]

#Add a printable panel title
avg[, panel_title := paste(' ', target, ' | N = ', ngenes, sep = " ")]
avg[, Time := as.numeric(as.character(Time))]

##Creating plot showing the avr zscore of genes
gg <- ggplot(data= avg, aes(x= Time, y= zscore, by = Time)) +
  geom_line() +
  geom_point(size = 1.0)+
  geom_errorbar(aes(ymin= zscore - sd, ymax= zscore + sd), width=.2, position=position_dodge(.9)) +
  facet_wrap(~panel_title, nrow=2) +
  ggtitle("Temporal changes in average gene Z-score of
          AP2-O & AP2-O3 transcription factor target genes") +
  xlab("Time (hr)") +
  ylab("Average gene Z-score")+
  theme_linedraw()
gg <- gg +
  scale_x_continuous(breaks= c(0, 4, 8, 12, 16, 24))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12, face = "bold"),
        strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust=0.5, size = 16, face = "bold"))+
  guides(y.sec = guide_axis())
ggsave(AP2_O_O3_plot, width= 22, height= 15, units= 'cm')

#############################################################################################

######################## Plot normalised expression logrpkm #################################
### need to work on taking into account Dario's suggestions

logrpkm_table_long <- fread(logrpkm_table)



# Calculating the slope and the mean normalised expression values (log2 rpkm) for the targets of
#the transcription factors, alone and in combination using a linear mixed effects model which
# accounts for multiple genes and library replicates in each transcription factor group:

sink(lmer_out_file) # Send output to this file instead of to the Terminal

cat ('# AP2-FG & AP2-O3 \n\n')
key_genes <- geneTarget[gene_id %in% geneTarget[(target == "AP2-FG & AP2-O3") | (target == "AP2-FG") |
                                              (target == "AP2-O3")]$gene_id]
key_genes <- merge(logrpkm_table_long, key_genes, by= 'gene_id')
key_genes <- merge(key_genes, ss[, list(library_id, Time)], by= 'library_id')

cat('# Model fit \n\n')

fitr <- lmer(logrpkm ~ Time * target + (1 + Time|gene_id), data= key_genes, REML= FALSE, control = lmerControl(optimizer ="Nelder_Mead"))
summary(fitr)

cat('\n\n# Compare slopes \n\n')

trendr <- emtrends(fitr, 'target', var="Time", pbkrtest.limit = 10000, lmerTest.limit = 10000)
trendr
pairs(trendr)

cat('\n\n# Compare means \n\n')

means <- emmeans(fitr, 'target', by= 'Time', pbkrtest.limit = 10000, lmerTest.limit = 10000)
means
pairs(means)

cat('\n# ================ \n')

cat ('# AP2-O & AP2-O4 \n\n')

key_genes <- geneTarget[gene_id %in% geneTarget[(target == "AP2-O & AP2-O4") | (target == "AP2-O") |
                                              (target == "AP2-O4")]$gene_id]
key_genes <- merge(logrpkm_table_long, key_genes, by= 'gene_id')
key_genes <- merge(key_genes, ss[, list(library_id, Time)], by= 'library_id')

cat('# Model fit \n\n')

fitr <- lmer(logrpkm ~ Time * target + (1 + Time|gene_id), data= key_genes, REML= FALSE, control = lmerControl(optimizer ="Nelder_Mead"))
summary(fitr)

cat('\n\n# Compare slopes \n\n')

trendr <- emtrends(fitr, 'target', var="Time", pbkrtest.limit = 10000, lmerTest.limit = 10000)
trendr
pairs(trendr)

cat('\n\n# Compare means \n\n')

means <- emmeans(fitr, 'target', by= 'Time', pbkrtest.limit = 10000, lmerTest.limit = 10000)
means
pairs(means)

cat('\n# ================ \n')

cat ('# AP2-O & AP2-O3 \n\n')

key_genes <- geneTarget[gene_id %in% geneTarget[(target == "AP2-O & AP2-O3") | (target == "AP2-O") |
                                              (target == "AP2-O3")]$gene_id]
key_genes <- merge(logrpkm_table_long, key_genes, by= 'gene_id')
key_genes <- merge(key_genes, ss[, list(library_id, Time)], by= 'library_id')

cat('# Model fit \n\n')

fitr <- lmer(logrpkm ~ Time * target + (1 + Time|gene_id), data= key_genes, REML= FALSE, control = lmerControl(optimizer ="Nelder_Mead"))
summary(fitr)

cat('\n\n# Compare slopes \n\n')

trendr <- emtrends(fitr, 'target', var="Time", pbkrtest.limit = 10166, lmerTest.limit = 10166)
trendr
pairs(trendr)

cat('\n\n# Compare means \n\n')

means <- emmeans(fitr, 'target', by= 'Time',pbkrtest.limit = 10166, lmerTest.limit = 10166)
means
pairs(means)

cat('\n# ================ \n')

sink() # Stop sending output to file

