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
logrpkm_table <- snakemake@input[['logrpkm_table']]
tf_plot_trendline <- snakemake@output[['tf_plot_trendline']]
#lmer_out_file <- snakemake@output[['lmer_out_file']]
trendr_file <- snakemake@output[['trendr_file']]
pairs_trendr_file <- snakemake@output[['pairs_trendr_file']]
means_file <- snakemake@output[['means_file']]
pairs_means_file <- snakemake@output[['pairs_means_file']]


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

######################## Plot normalised expression logrpkm #################################
#Use logrpkm raw scores to cluster all genes as scaling makes the standard deviation 1 regardless of input.
#If you use raw logrpkm's you can also test whether some genes bound by one
#and/or another transcription factor tend to have different mean of expression.
logrpkm_table_long <- fread(logrpkm_table)

ss <- fread(ss_file)
ss <- ss[Outliers == FALSE,]

#table with number of trxn factors and comma separated list of trxn factors
geneTarget <- genes[, list(N= length(unique(target)),
                         target= paste(sort(target), collapse= ' & ')), by= gene_id]

#for AP2-FG and AP2-O3 plot
ap2fg_ap2o3 <- geneTarget[gene_id %in% geneTarget[(target == "AP2-FG & AP2-O3") | (target == "AP2-FG") |
                                              (target == "AP2-O3")]$gene_id]

ap2fg_ap2o3 <- merge(logrpkm_table_long, ap2fg_ap2o3, by= 'gene_id')
ap2fg_ap2o3 <- merge(ap2fg_ap2o3, ss[, list(library_id, Time)], by= 'library_id')

fitr <- lmer(logrpkm ~ Time * target + (1 + Time|gene_id), data= ap2fg_ap2o3, REML= FALSE, control = lmerControl(optimizer ="Nelder_Mead"))
ap2fg_ap2o3$predict <- predict(fitr)

avr_ap2fg_ap2o3 <- ap2fg_ap2o3[, list(logrpkm = mean(logrpkm), sd= sd(logrpkm), predict = mean(predict), ngenes= length(unique(.SD$gene_id))),
                 by= list(Time, target)]
avr_ap2fg_ap2o3[,tf_group := "AP2-FG and AP2-O3"]


#for AP2-O and AP2-O4
ap2o_ap2o4 <- geneTarget[gene_id %in% geneTarget[(target == "AP2-O & AP2-O4") | (target == "AP2-O") |
                                              (target == "AP2-O4")]$gene_id]

ap2o_ap2o4 <- merge(logrpkm_table_long, ap2o_ap2o4, by= 'gene_id')
ap2o_ap2o4 <- merge(ap2o_ap2o4, ss[, list(library_id, Time)], by= 'library_id')

fitr <- lmer(logrpkm ~ Time * target + (1 + Time|gene_id), data= ap2o_ap2o4, REML= FALSE, control = lmerControl(optimizer ="Nelder_Mead"))
ap2o_ap2o4$predict <- predict(fitr)

avr_ap2o_ap2o4 <- ap2o_ap2o4[, list(logrpkm = mean(logrpkm), sd= sd(logrpkm), predict = mean(predict), ngenes= length(unique(.SD$gene_id))),
                           by= list(Time, target)]

avr_ap2o_ap2o4[,tf_group := "AP2-O and AP2-O4"]

#for AP2-O and AP2-O3 plot

ap2o_ap2o3 <- geneTarget[gene_id %in% geneTarget[(target == "AP2-O & AP2-O3") | (target == "AP2-O") |
                                              (target == "AP2-O3")]$gene_id]

ap2o_ap2o3 <- merge(logrpkm_table_long, ap2o_ap2o3, by= 'gene_id')
ap2o_ap2o3 <- merge(ap2o_ap2o3, ss[, list(library_id, Time)], by= 'library_id')

fitr <- lmer(logrpkm ~ Time * target + (1 + Time|gene_id), data= ap2o_ap2o3, REML= FALSE, control = lmerControl(optimizer ="Nelder_Mead"))
ap2o_ap2o3$predict <- predict(fitr)

avr_ap2o_ap2o3 <- ap2o_ap2o3[, list(logrpkm = mean(logrpkm), sd= sd(logrpkm), predict = mean(predict), ngenes= length(unique(.SD$gene_id))),
                         by= list(Time, target)]

avr_ap2o_ap2o3[, tf_group := "AP2-O and AP2-O3"]

#merge the average data
avg <- rbindlist(list(avr_ap2fg_ap2o3, avr_ap2o_ap2o4, avr_ap2o_ap2o3))

# Labels with transcription factor combinations & the number of genes:
avg[, tf_ngenes := paste(' ', target, ' | N = ', ngenes, sep = " ")]
labels <- avg[, list(Time= max(Time)), by= list(tf_ngenes)]
labels <- merge(labels, avg, by= c('tf_ngenes', 'Time'))

##Creating plot showing the avr expression of transcription factor target genes with fitted lines

gg <- ggplot(data= avg, aes(x= Time, y= logrpkm, by = target, colour = target)) +
  geom_errorbar(aes(ymin= logrpkm - sd, ymax= logrpkm + sd), width=.1, size = 0.5, position=position_dodge(.6)) +
  geom_line(size =0.5, position = position_dodge(0.6)) +
  geom_point(size = 1.5, position = position_dodge(0.6))+
  geom_line(aes(y=predict), size =0.2) +
  facet_wrap(~tf_group, nrow=2) +
  ggtitle("Temporal changes in the average normalised gene expression of
          transcription factor target genes") +
  xlab("Time (hr)") +
  ylab("Average normalised gene expression (log2 rpkm)")+
  theme_linedraw()

gg <- gg +
  geom_text_repel(data= labels, xlim= c(max(labels$Time), NA), aes(x = Time, y=logrpkm, label= tf_ngenes), size = 3, inherit.aes = FALSE, hjust= 1, segment.colour= 'grey57') +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.7), add = 0), breaks= c(0, 4, 8, 12, 16, 24))+
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=13, face = "bold"),
        strip.text = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust=0.5, size = 16, face = "bold"),
        legend.position = "none")+
        guides(y.sec = guide_axis())
ggsave(tf_plot_trendline, width= 25, height= 15, units= 'cm')

# Calculating the slope and the mean normalised expression values (log2 rpkm) for the targets of
#the transcription factors, alone and in combination using a linear mixed effects model which
# accounts for multiple genes and library replicates in each transcription factor group:

#sink(lmer_out_file) # Send output to this file instead of to the Terminal

cat ('# AP2-FG & AP2-O3 \n\n')

cat('# Model fit \n\n')

fitr <- lmer(logrpkm ~ Time * target + (1 + Time|gene_id), data= ap2fg_ap2o3, REML= FALSE, control = lmerControl(optimizer ="Nelder_Mead"))
summary(fitr)

cat('\n\n# Compare slopes \n\n')

trendr <- emtrends(fitr, 'target', var="Time", pbkrtest.limit = 10000, lmerTest.limit = 10000)
fg_o3_trendr <- as.data.table(trendr)
fg_o3_pairs_trendr <- as.data.table(pairs(trendr))

cat('\n\n# Compare means \n\n')

means <- emmeans(fitr, 'target', by= 'Time', pbkrtest.limit = 10000, lmerTest.limit = 10000)
fg_o3_means <- as.data.table(means)
fg_o3_pairs_means <- as.data.table(pairs(means))

cat('\n# ================ \n')

cat ('# AP2-O & AP2-O4 \n\n')

cat('# Model fit \n\n')

fitr <- lmer(logrpkm ~ Time * target + (1 + Time|gene_id), data= ap2o_ap2o4, REML= FALSE, control = lmerControl(optimizer ="Nelder_Mead"))
summary(fitr)

cat('\n\n# Compare slopes \n\n')

trendr <- emtrends(fitr, 'target', var="Time", pbkrtest.limit = 10000, lmerTest.limit = 10000)
o_o4_trendr <- as.data.table(trendr)
o_o4_pairs_trendr <- as.data.table(pairs(trendr))

cat('\n\n# Compare means \n\n')

means <- emmeans(fitr, 'target', by= 'Time', pbkrtest.limit = 10000, lmerTest.limit = 10000)
o_o4_means <- as.data.table(means)
o_o4_pairs_means <- as.data.table(pairs(means))

cat('\n# ================ \n')

cat ('# AP2-O & AP2-O3 \n\n')

cat('# Model fit \n\n')

fitr <- lmer(logrpkm ~ Time * target + (1 + Time|gene_id), data= ap2o_ap2o3, REML= FALSE, control = lmerControl(optimizer ="Nelder_Mead"))
summary(fitr)

cat('\n\n# Compare slopes \n\n')

trendr <- emtrends(fitr, 'target', var="Time", pbkrtest.limit = 10166, lmerTest.limit = 10166)
o_o3_trendr <- as.data.table(trendr)
o_o3_pairs_trendr <- as.data.table(pairs(trendr))

cat('\n\n# Compare means \n\n')

means <- emmeans(fitr, 'target', by= 'Time',pbkrtest.limit = 10166, lmerTest.limit = 10166)
o_o3_means <- as.data.table(means)
o_o3_pairs_means <- as.data.table(pairs(means))

cat('\n# ================ \n')

cat('# concatenate trendr \n\n')

trendr <- rbindlist(list(fg_o3_trendr, o_o4_trendr, o_o3_trendr))
write.table(trendr, file=trendr_file, row.names = FALSE, sep= '\t', quote= FALSE)

cat('# concatenate pairs(trendr) \n\n')

pairs_trendr <- rbindlist(list(fg_o3_pairs_trendr, o_o4_pairs_trendr, o_o3_pairs_trendr))
write.table(pairs_trendr, file=pairs_trendr_file, row.names = FALSE, sep= '\t', quote= FALSE)

cat('# concatenate means \n\n')

means <- rbindlist(list(fg_o3_means, o_o4_means, o_o3_means))
write.table(means, file=means_file, row.names = FALSE, sep= '\t', quote= FALSE)

cat('# concatenate pairs(means) \n\n')

pairs_means <- rbindlist(list(fg_o3_pairs_means, o_o4_pairs_means, o_o3_pairs_means))
write.table(pairs_means, file=pairs_means_file, row.names = FALSE, sep= '\t', quote= FALSE)

#sink() # Stop sending output to file

