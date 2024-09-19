#PCA主成分分析
library(ggplot2)
library(dplyr)

args <- commandArgs(T)

otufile <- args[1] #OTU特征表文件
metafile <- args[2] #样本分组信息metadata
groupby <- args[3] #分组列名
outprefix <- args[4] #输出文件前缀
ptitle <- args[5] #图片中的标题

outfile <- paste0(outprefix, ".txt")  #PCA结果
outpng <- paste0(outprefix, ".png")      #png绘图结果
outpdf <- paste0(outprefix, ".pdf")      #pdf绘图结果

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
#转换为dataframe
OtuData <- data.frame(OtuData)

#读取分组信息表
MetaData <- read.table(metafile,
                       row.names = 1,
                       header = T,
                       sep = "\t",
                       comment = "",
                       check.names = F)

MetaData$Groups <- MetaData[[groupby]]

#按照OtuData样本顺序获取分组
grouplist <- NULL
for (spname in rownames(OtuData)) {
  #groupname <- MetaData[spname, groupby]
  groupname <- MetaData[spname, "Groups"]
  grouplist <- c(grouplist, groupname)
}
#dataframe才能这样添加group列信息
OtuData$group <- grouplist

################# PCA 计算 ################################

#1. 数据标准化，对原数据进行z-score归一化
dt <- as.matrix(scale(OtuData[, 1:(ncol(OtuData) - 1)])) #不含group列
head(dt)

#2. 计算协方差
rm1 <- cor(dt)
#head(rm1)

#3. 计算特征值和相应的特征向量
#### 特征分解
rs1 <- eigen(rm1)
val <- rs1$values #转化为标准差standard deviation
Standard_deviation <- sqrt(val)
#Standard_deviation
#### 计算方差贡献率和累积贡献率
Proportion_of_Variance <- val / sum(val)
#Proportion_of_Variance
Cumulative_Proportion <- cumsum(Proportion_of_Variance)
#Cumulative_Proportion

#4. 计算主成分得分
#### 提取结果中的特征向量(也称为Loadins 载荷矩阵)
U <- as.matrix(rs1$vectors)

#### 进行矩阵乘法, 获得PC score
PC <- dt %*% U
colnames(PC) <- paste0("PC", seq_len(ncol(PC))) #重命名表头
#print(colnames(PC))
### 合并分组列
plotdata <- data.frame(PC, group = OtuData$group)

write.table(
  data.frame(SampleID = rownames(plotdata), plotdata[, c("PC1", "PC2", "PC3", "group")]),
  outfile,
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

head(plotdata$PC1)
head(plotdata$PC2)
unique(plotdata$group)

#### 开始绘图
p <- ggplot(plotdata, aes(x = PC1, y = PC2)) +
  geom_point(
    aes(fill = group),
    size = 2,
    pch = 21) +
  stat_ellipse(
    type = "norm",
    aes(fill = group),
    linewidth = 0.5,
    geom = "polygon",
    alpha = 0.2
  ) +
  labs(
    x = paste0("PC1(", round(Proportion_of_Variance[1] * 100, 2), "%)"),
    y = paste0("PC2(", round(Proportion_of_Variance[2] * 100, 2), "%)")
  ) +
  geom_vline(aes(xintercept = 0), linetype = "dotted") +
  geom_hline(aes(yintercept = 0), linetype = "dotted") +
  ggtitle(ptitle) +
  #theme_classic(base_line_size = 1) +
  theme_bw() +
  theme(
    #panel.background = element_rect(fill = "white", color = "black"),
    #panel.grid = element_blank(),
    axis.title.x = element_text(size = 15, color = "black", face = "bold"),
    axis.title.y = element_text(size = 15, color = "black", face = "bold", angle = 90),
    axis.text.x = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5),
    axis.text.y = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5),
    legend.title = element_text(size = 11, color = "black", face = "bold"),
    legend.text = element_text(size = 10, color = "black", face = "bold"),
    legend.position = "bottom"
  )
ggsave(outpng, p)
ggsave(outpdf, p)
