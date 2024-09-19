### 绘图异常需要修改
#PLS_DA分析
suppressMessages(library(mixOmics))
suppressMessages(library(ggplot2))
suppressMessages(library(ggalt))
suppressMessages(library(dplyr))

args <- commandArgs(T)

otufile <- args[1]
metafile <- args[2]
groupby <- args[3]
outprefix <- args[4]
ptitle <- args[5]

outfile <- paste0(outprefix, ".txt")
outpng <- paste0(outprefix, ".png")
outpdf <- paste0(outprefix, ".pdf")

#读取OTU表信息
OtuData <- read.table(otufile,
                      row.names = 1,
                      header = T,
                      sep = "\t",
                      skip = 1,
                      comment.char = "",
                      check.names = F)

#转置
OtuData <- t(OtuData)

#读取分组信息表
MetaData <- read.table(metafile,
                       row.names = 1,
                       header = T,
                       sep = "",
                       comment = "",
                       check.names = F)

MetaData$Groups <- MetaData[[groupby]]

bool_idx <- rownames(OtuData) %in% rownames(MetaData)
OtuData <- OtuData[bool_idx,]
#head(OtuData)

#按照OtuData样本顺序获取分组
grouplist <- NULL
for (spname in rownames(OtuData)) {
  #groupname <- MetaData[spname, groupby]
  groupname <- MetaData[spname, "Groups"]
  grouplist <- c(grouplist, groupname)
}
GroupData <- data.frame(group = grouplist)
rownames(GroupData) <- rownames(OtuData)
#head(GroupData)

########################################### PLS-DA ################################################
#PLS-DA分析，这里也是选取2个主成分
plsdaRes <- plsda(OtuData, GroupData$group, ncomp = 2)
#消除object对象
plsdaRes <- unclass(plsdaRes)
print(names(plsdaRes))
# 提取绘图主成分score值数据
plotData <- plsdaRes$variates$X
colnames(plotData) <- paste0("comp", seq_len(ncol(plotData))) #重命名表头列名
plotData <- data.frame(plotData, group = GroupData$group)
head(plotData)
#写入文件
write.table(
  data.frame(SampleID = rownames(plotData), plotData),
  outfile,
  sep = "\t", quote = F, row.names = F, col.names = T)
# 提取解释度数据
eig <- plsdaRes$explained_variance$X
head(eig)

# 绘图
p <- ggplot(plotData,
            aes(x = comp1, y = comp2, group = group, color = group)) +
  geom_point(size = 2) +
  stat_ellipse(
    level = 0.95,
    linetype = "solid",
    linewidth = 0.5) +
  #geom_encircle(s_shape = 1, expand=0) +
  labs(
    title = paste0("OPLS-DA(", ptitle, ")"),
    x = paste0("X-variate 1(", format(100 * eig[1]), "%)"),
    y = paste0("X-variate 2(", format(100 * eig[2]), "%)")) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top",
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10) #unit="cm"
  ) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted")

ggsave(outpdf, p)
ggsave(outpng, p)
