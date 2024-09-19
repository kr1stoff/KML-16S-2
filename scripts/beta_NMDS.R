#NMDS
suppressMessages(library(ggplot2))
suppressMessages(library(vegan))
suppressMessages(library(dplyr))

args <- commandArgs(T)

metafile <- args[1] #样本分组信息文件
matrixfile <- args[2] #qiime2的距离矩阵文件
groupby <- args[3] #分组表头
outprefix <- args[4] #输出图片文件前缀
outpng <- paste0(outprefix, ".png")
outpdf <- paste0(outprefix, ".pdf")

#读取样本信息分组文件
MetaData <- read.table(metafile,
                       header = T,
                       sep = "\t",
                       check.names = F,
                       row.names = 1,
                       comment.char = "")

MetaData$Groups <- MetaData[[groupby]]

#读取距离矩阵文件
MatrixData <- read.table(matrixfile,
                         header = T,
                         sep = "\t",
                         check.names = F,
                         row.names = 1,
                         comment.char = "")

# 样本筛选
bool_idx <- rownames(MatrixData) %in% rownames(MetaData)
MatrixData <- MatrixData[bool_idx, bool_idx]

#NMSS计算
dfNmds <- metaMDS(MatrixData, trace = F)

#绘图数据
plotdata <- data.frame(dfNmds$points)

#data数据顺序获取group分组信息
grouplist <- NULL
for (spname in rownames(plotdata)) {
  #groupname <- MetaData[spname, groupby]
  groupname <- MetaData[spname, "Groups"]
  grouplist <- c(grouplist, groupname)
}
plotdata$group <- grouplist
print(head(plotdata))

#绘图
p <- ggplot(plotdata,
            aes(x = MDS1, y = MDS2, color = group, group = group, fill = group)) +
  geom_point(size = 2) +
  theme_classic() +
  stat_ellipse( #添加置信区间
    geom = "polygon",
    level = 0.95,
    alpha = 0.3
  ) +
  geom_text(  #添加文本标签
    aes(label = rownames(plotdata)),
    vjust = 1.5,
    size = 2,
    color = "black"
  ) +
  labs( #副标题处添加stress
    subtitle = paste0("stress=", round(dfNmds$stress, 2)))

#ggsave(outpng, p, dpi=300)
#ggsave(outpdf, p, height = 10, width = 10, dpi=300)
ggsave(outpng, p)
ggsave(outpdf, p)
