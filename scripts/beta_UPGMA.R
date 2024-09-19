#UPGMA分析
suppressMessages(library("phangorn"))
suppressMessages(library(ggtree))
#suppressMessages(library(patchwork))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))

args <- commandArgs(T)

matrixfile <- args[1]
metafile <- args[2]
groupby <- args[3]
outprefix <- args[4]
circle_outpng <- paste0(outprefix, ".circular.png")
circle_outpdf <- paste0(outprefix, ".circular.pdf")
rect_outpng <- paste0(outprefix, ".rectangular.png")
rect_outpdf <- paste0(outprefix, ".rectangular.pdf")

#读取距离矩阵文件
MatrixData <- read.table(matrixfile,
                         header = T,
                         row.names = 1,
                         sep = "\t",
                         check.names = F)

#读取分组信息表
MetaData <- read.table(metafile,
                       row.names = 1,
                       header = T,
                       sep = "\t",
                       comment = "",
                       check.names = F)

MetaData$Groups <- MetaData[[groupby]]

bool_idx <- rownames(MatrixData) %in% rownames(MetaData)
MatrixData <- MatrixData[bool_idx, bool_idx]

#按照分组列名把样本分组
#GroupInfo <- split(rownames(MetaData), MetaData[[groupby]])
GroupInfo <- split(rownames(MetaData), MetaData["Groups"])
print(GroupInfo)

#使用upgma聚类算法
tree_upgma <- upgma(MatrixData)
tree_upgma_root <- root(tree_upgma, out = 1) #有根树

# 给树添加分组信息
tree_upgma <- groupOTU(tree_upgma, GroupInfo)
print(tree_upgma)

#标签字体大小
#leaf_text_size <- min(c(5, 360/length(tree_upgma$tip)))
leaf_text_size <- min(c(5, 200 / length(tree_upgma$tip)))

print(leaf_text_size)
p1 <- ggtree(tree_upgma,
             branch.length = "none", #"branch.length",
             layout = "circular",
             linetype = 1,
             size = 0.65,
             aes(color = group),
             ladderize = F) +
  geom_tiplab(size = leaf_text_size) +
  labs(title = "UPGMA Tree(circular)") +
  theme(legend.position = "top",
        axis.text = element_text(size = 8),
        plot.title = element_text(hjust = 0.5))

p2 <- ggtree(tree_upgma,
             branch.length = "none", #"branch.length",
             layout = "rectangular",
             linetype = 1,
             size = 0.65,
             aes(color = group),
             ladderize = F) +
  geom_tiplab(size = 0.5 * leaf_text_size) +
  labs(title = "UPGMA Tree(rectangular)") +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5))

#p <- p1 + p2 +
#    plot_layout(widths = c(1,1), ncol=2, nrow=1)

ggsave(circle_outpng, p1)
ggsave(circle_outpdf, p1)
ggsave(rect_outpng, p2)
ggsave(rect_outpdf, p2)
