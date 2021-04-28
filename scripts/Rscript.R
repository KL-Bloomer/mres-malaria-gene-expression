R script - mds plot

rm(list=ls())
library(data.table)
library(edgeR)

setwd('~/mres-malaria-gene-expression')
getwd()

ss <- fread('sample_sheet.tsv')

setwd('~/mres-malaria-gene-expression/output/featureCounts/')

counts <- fread(cmd= 'grep -v "^#" counts.tsv')

setnames(counts, names(counts), sub('.bam', '', basename(names(counts))))
counts <- counts[, c('Geneid', ss$library_id), with= FALSE]

mat <- as.matrix(counts, rownames= 'Geneid')

ss[, group := paste(Time, Strain, sep= '_')]

design <- model.matrix(~0 + ss$group)
colnames(design) <- sub('ss$group', '', colnames(design), fixed= TRUE)

y <- DGEList(counts= mat, group= ss$group)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design)

# MDS plot before and after batch correction

pdf('mds.pdf', width= 14/2.54, height= 14/2.54)
mds <- plotMDS(y, col= as.numeric(as.factor(ss$group)))
pdf('mds.pdf')

pdf('mds_mon.pdf', width= 14/2.54, height= 14/2.54)
> mds <- plotMDS(y, xlab="PC1", ylab="PC2", col=as.numeric(as.factor(ss$group)), pch=as.numeric(as.factor(ss$group)))
> dev.off()

pdf('mds_mon_2.pdf')
>  mds <- plotMDS(y, xlab="Dimension 1", ylab="Dimension 2", main="MDS plot", col=as.numeric(as.factor(ss$group)), pch=as.numeric(as.factor(ss$group)))
> mds <- plotMDS(y, xlab="Dimension 1", ylab="Dimension 2", main="MDS plot", col=as.numeric(as.factor(ss$group)), pch=19)

pdf('mds_mon_gg.pdf')
mds <- plotMDS(y)


testing <- data.frame(mds$cmdscale.out, ss$group)
colnames(testing)[2] <- "MDS2"
> colnames(testing)[1] <- "MDS1"
> colnames(testing)[3] <- "group"
pdf("mds_ggplot.pdf")
> mds_plot <- ggplot(data=testing, aes(MDS1, MDS2, colour = group)) + geom_point() + theme_bw()
> dev.off()
