suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
#library(ggtree)
args <- commandArgs(T)

feature_table_file <- args[1]
outprefix <- args[2]

result_file  <- paste(outprefix, ".result.txt", sep="")
barstat_file <- paste(outprefix, ".barstat.txt", sep="")
outpng <- paste(outprefix, ".barplot.png", sep="")
outpdf <- paste(outprefix, ".barplot.pdf", sep="")

#读取Species/featur-table.tsv
data <- read.table(
    feature_table_file,
    header = T,
    sep = "\t",
    row.names = 1,
    comment.char = "",
    skip = 1,
    check.names = F
)

# 统计每个物种分类组成level计数
result <- data.frame(
    Kigdom = colSums(data[grepl("d__", rownames(data)), ]),
    Phylum = colSums(data[grepl("p__", rownames(data)), ]),
    Class  = colSums(data[grepl("c__", rownames(data)), ]),
    Order  = colSums(data[grepl("o__", rownames(data)), ]),
    Family = colSums(data[grepl("f__", rownames(data)), ]),
    Genus  = colSums(data[grepl("g__", rownames(data)), ]),
    Species = colSums(data[grepl("s__", rownames(data)), ]),
    Unassigned = colSums(data[grepl("Unassigned", rownames(data)), ])
)
head(result)
write.table(
    data.frame(SampleID=rownames(result), result),
    result_file,
    quote = F, sep = "\t", row.names = F, col.names = T
    )

barstat <- data.frame(
    Kingdom = result$Kigdom - result$Phylum,
    Phylum  = result$Phylum - result$Class,
    Class   = result$Class  - result$Order,
    Order   = result$Order  - result$Family,
    Family  = result$Family - result$Genus,
    Genus   = result$Genus  - result$Species,
    Species = result$Species- result$Unassigned,
    Unassigned = result$Unassigned
)
barstat <- data.frame(SampleID =rownames(result), barstat)
head(barstat)
write.table(
    barstat,
    barstat_file,
    quote = F, sep = "\t", row.names = F, col.names = T    
)

barplotdata <- melt(barstat, id.vars=c("SampleID"), variable.name = "Levels")
head(barplotdata)

p <- ggplot(barplotdata, aes(x=SampleID, y=value, fill=Levels)) + 
    geom_bar(stat = "identity", width=0.7) +
    theme_classic() +
    theme(
        axis.text = element_text(family="serif"),
        axis.text.x = element_text(angle = 90, size = 6),
    )

ggsave(outpdf, p, height=7, width =14)
ggsave(outpng, p, height=7, width =14)

