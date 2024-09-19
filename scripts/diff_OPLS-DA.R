### 绘图异常需要修改
#PLS_DA分析
suppressMessages(library(mixOmics))
suppressMessages(library(ggplot2))
suppressMessages(library(ropls))
suppressMessages(library(dplyr)) #管道符筛选 %>%, %in%, mutate
suppressMessages(library(patchwork)) #图片排版

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


#按照OtuData样本顺序获取分组
grouplist <- NULL
for (spname in rownames(OtuData)) {
  groupname <- MetaData[spname, groupby]
  #groupname <- MetaData[spname, "Groups"]
  grouplist <- c(grouplist, groupname)
}
GroupData <- data.frame(group = grouplist)
rownames(GroupData) <- rownames(OtuData)
head(GroupData)

############################################## OPLS-DA#########################################
#OPLS-DA分析
oplsdaRes <- opls(OtuData, GroupData$group, predI = 0, orthoI = NA)
#print(oplsdaRes)
#print(names(oplsdaRes))
#提取score值会图
oplsdaRes@scoreMN
oplsScore <- oplsdaRes@scoreMN %>%
  as.data.frame() %>%
  mutate(group = GroupData$group, o1 = oplsdaRes@orthoScoreMN[, 1])
head(oplsScore)

#写入文件
write.table(
  data.frame(SampleID = rownames(oplsScore),
             scoreNM = oplsScore$p1,
             orthoScoreMN = oplsScore$o1,
             group = oplsScore$group),
  outfile,
  sep = "\t", quote = F, row.names = F, col.names = T)

p <- ggplot(oplsScore,
            aes(x = p1, y = o1, color = group, group = group)) +
  geom_point(size = 2) +
  stat_ellipse(
    level = 0.95,
    linetype = "solid",
    linewidth = 0.5) +
  labs(
    title = paste0("OPLS-DA(", ptitle, ")"),
    x = "scoreNM",
    y = "orthoScoreMN") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top",
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  ) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted")

ggsave(outpdf, p)
ggsave(outpng, p)
