### Barplot for visualisation of library size - before library size normalisation, 
#removal of outliers, rRNA and batch correction

rm(list=ls())
library(data.table)

#read sample sheet with outliers
setwd('~/MRes_Malaria_2021/git_repos/mres-malaria-gene-expression/')
ss <- fread('sample_sheet.tsv')
ss$Time <- as.factor(ss$Time)

#read counts with outliers
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/featureCounts/')
counts <- fread(cmd= 'grep -v "^#" counts.tsv')
setnames(counts, names(counts), sub('.bam', '', basename(names(counts))))
counts <- counts[, c('Geneid', ss$library_id), with= FALSE]

## Creating bar plot of lib sizes; prior to normalisation etc. 
# Create a new dataframe
sizes <- data.frame(LibraryID = ss$library_id, LibrarySize = y$samples$lib.size, Strain = ss$Strain,
                    Time = ss$Time, group = ss$group)

setDT(sizes)
class(sizes$LibrarySize)
#OR
is.data.table(sizes) == TRUE

#Make a nice label to print and plot; unique and informative
sizes[, label := paste(LibraryID, Time)]
print(sizes)
xord <- sizes[order(-Time, LibrarySize)]$label
print(xord) # The order we want

# New order of levels for the variable label
sizes[, label := factor(label, xord)]

library(RColorBrewer)
library(ggplot2)
#plot barplot; ordered by time
ggplot(sizes, aes(label, LibrarySize/1000000)) +
  geom_bar(stat="identity", aes(fill=Time)) +
  ggtitle("Barplot showing library size before normalisation") +
  xlab("Library name and time") +
  ylab("Library size (x 1 million)") +
  labs(fill = "Time(h)") + 
  coord_flip() +
  theme_bw()
setwd('~/MRes_Malaria_2021/output/mres-malaria-gene-expression/R_output/')
ggsave('barplot_libsizes_before norm_time.png', width = 10, height = 10)

