#韦恩图
library(ggVennDiagram)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
#library(extrafont) #显示中文
#font_import(paths="/lustre_new3.0/rawData/analysis/CNV-seq/workflow-16S/fonts", pattern="simkai", prompt=F) #载入楷体, 只需一次, 使用变量名为：KaiTi
#font_import(paths="/lustre_new3.0/rawData/analysis/CNV-seq/workflow-16S/fonts", pattern="simfang", prompt=F) #载入仿宋体, 只需一次, 使用变量名为：FangSong
#载入字体库(包括上述安装的中文)
#loadfonts()
#设置全局系统字体变量
#Sys.setenv(CAIRO_PDF_FONT="KaiTi")

#修改默认中文编码
options(encoding = "UTF-8")


args <- commandArgs(T)

otufile <- args[1] # OUTs/ASVs特征表
metafile <- args[2] # 样本分组信息文件
groupby <- args[3] # 分组列名
outprefix <- args[4] # 输出图片文件前缀
outpng <- paste0(outprefix, ".png")
outpdf <- paste0(outprefix, ".pdf")

#读取特征表feature-table.tsv
OtuData <- read.table(otufile,
                      header = T,
                      sep = "\t",
                      row.names = 1,
                      comment.char = "",
                      skip = 1,
                      check.names = F)
#转置
OtuData <- t(OtuData)

#读取样本分组信息文件
MetaData <- read.table(metafile,
                       header = T,
                       sep = "\t",
                       row.names = 1,
                       comment.char = "",
                       check.names = F)

MetaData$Groups <- MetaData[[groupby]]

# 获取分组列 去重后的组内名称
#uniqgroup <- unique(MetaData[[groupby]])
uniqgroup <- unique(MetaData$Groups)

print(uniqgroup)

#传教韦恩图数据列表
plotdata <- list()
for (groupname in uniqgroup) {
  # 获取每个分组名称的样本信息
  #GroupMetaData <-  MetaData[MetaData[[groupby]] == groupname, ]
  GroupMetaData <- MetaData[MetaData$Groups == groupname,]
  # 获取每个分组名称的样本的OTU信息(行号rownames都是样本名称)
  GroupOtuData <- OtuData[rownames(OtuData) %in% rownames(GroupMetaData),]
  #print(rownames(GroupOtuData))
  # 对每列ASVs/OTUs求和
  colsum <- colSums(GroupOtuData)
  # 获取非0的OTUs/ASVs名称
  OtuIDs <- names(colsum[colsum > 0])
  # 把数据添加到plotdata
  plotdata[[groupname]] <- OtuIDs
}
print(names(plotdata))

p <- ggVennDiagram(plotdata,
                   category.names = names(plotdata),
                   label = "count") +
  scale_x_continuous(expand = expansion(mult = .3)) +
  scale_fill_distiller(palette = "Blues", direction = -1) +
  scale_color_brewer(palette = "Set3") +
  labs(title = paste("OTUs/ASVs Venn Diagram by", groupby, " Sets", sep = " ")) +
  theme(
    plot.title = element_text(hjust = 0.5), #标题居中
    text = element_text(family = "serif")
  )

ggsave(outpng, p, dpi = 150)
ggsave(outpdf, p, height = 10, width = 10)

