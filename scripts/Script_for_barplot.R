### Barplot for visualisation of library size - before library size normalisation, 

#Packages required
library(data.table)
library(edgeR)
library(RColorBrewer)
library(ggplot2)

# Input and output files:
ss <- snakemake@input[['sample_sheet']]
counts_file <- snakemake@input[['counts']]
barplot_libsizes_beforenorm <- snakemake@output[['barplot_libsizes_beforenorm']]

#read the sample sheet
ss <- fread(ss)
ss <- ss[library_type == "RNA-seq",]
ss[, Time := sprintf('%.2d', Time)]
ss[, group := paste(Time)]
ss$Time <- as.factor(ss$Time)

# Another sanity check that Time is a factor variable - maybe redundant
stopifnot(is.factor(ss$Time))

#read the counts table
counts <- fread(cmd= paste('grep -v "^#"', counts_file))
setnames(counts, names(counts), sub('.bam', '', basename(names(counts))))
counts <- counts[, c('Geneid', ss$library_id), with= FALSE]

#create a matrix of the counts
mat <- as.matrix(counts, rownames= 'Geneid')

design <- model.matrix(~0 + ss$group)
colnames(design) <- sub('ss$group', '', colnames(design), fixed= TRUE)
y <- DGEList(counts= mat,
             group= ss$group)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]

## Creating bar plot of lib sizes; prior to normalisation etc. 
# Create a new data frame
sizes <- data.frame(LibraryID = ss$library_id, LibrarySize = y$samples$lib.size, Strain = ss$Strain,
                    Time = ss$Time, group = ss$group)

#Check if sizes is a data table
setDT(sizes)

#Make a nice label to print and plot; unique and informative
sizes[, label := paste(LibraryID, Time)]
print(sizes)
xord <- sizes[order(-Time, LibrarySize)]$label
print(xord) # The order we want

# New order of levels for the variable label
sizes[, label := factor(label, xord)]

#plot barplot; ordered by time
gg <- ggplot(sizes, aes(label, LibrarySize/1000000)) +
  geom_bar(stat="identity", aes(fill=Time)) +
  ggtitle("Library sizes before normalisation") +
  xlab("Library name and time") +
  ylab("Library size (x 1 million)") +
  labs(fill = "Time(h)") + 
  coord_flip() 
  #theme_bw()
gg + theme(axis.title=element_text(size=10,face="bold"), legend.title = element_text(size = 10, face="bold"),
           plot.title=element_text(size=15, face="bold"),axis.text = element_text(face="bold"))
ggsave(barplot_libsizes_beforenorm, width = 10, height = 8)



