suppressMessages(library(vegan))
suppressMessages(library(ape))
suppressMessages(library(ggplot2))
suppressMessages(library(grid))
suppressMessages(library(dplyr))
suppressMessages(library(multcomp))
suppressMessages(library(patchwork))

args <- commandArgs(T)

metafile = args[1] #输入分组信息metadata
otufile = args[2] #输入OTU特征表文件
matrixfile = args[3] #距离矩阵文件, 例如：Diversity/core-metrics/bray_curtis_distance_matrix/distance-matrix.tsv, 可以为空(none), 此时matrix由otufile计算
groupby = args[4] #样本按照metadata文件中哪一列分组(特殊字符都需要转换成点.)
dismethod = args[5] #距离矩阵计算方法和多元差异分析方法, 比如：bray, jaccard, euclidean, manhattan, mahalanobis, binomial等
outprefix = args[6] #输出的pcoa结果前缀, 将输出 {outprefix}.pcoa.txt, {outprefix}.png, {outprefix}.pdf

outfile <- paste0(outprefix, ".pcoa.txt") #pcoa结果
outpng <- paste0(outprefix, ".png")      #png绘图结果
outpdf <- paste0(outprefix, ".pdf")      #pdf绘图结果

#读取分组信息
MetaData <- read.table(metafile, header = T, sep = "\t")
#metafile示例：
#sample-id   groupby1      groupby2     ...    groupbyN
#sample1     groupA        groupX       ...    groupN1
#sample2     groupB        groupY       ...    groupN2
#...         ...           ...          ...    ...
#NC2024      Control_NC    Control_NC   ...    Control_NC
#PC2024      Control_PC    Control_PC   ...    Control_PC
MetaData$Groups <- MetaData[[groupby]]

#获取质控样本
ControlData <- MetaData[grepl("Control_[NP]C", MetaData[[groupby]]),]

# 如果不存在matrixfile, 就用otufile计算matrixfile
# 因为计算距离矩阵函数vegdist的算法中没有考虑系统发育树的算法(Weighted-Unifrac/Unweighted-Unifrac), 因此需要用到qiime中计算的这个距离矩阵
# 但是后续permanova多元差异分析中还是考虑其他算法计算差异(比如欧几里得算法euclidean, 马氏距离算法mahalanobis, 二项分布算法binomial)
OtuData <- read.table(otufile,
                      row.names = 1,  #首列样本名称作为rowname)
                      header = T,     #因为样本名可能是数值, 如果作为表头会导致和Group中数据无法匹配上
                      sep = "\t",
                      skip = 1,       #跳过首行
                      comment.char = "%", #默认是#, 会导致第二行表头被忽略, 所以选择用%代替
                      check.names = F, #因为样本名可能是数值, 这个参数防止在数值前自动加X
)
#筛选样本
bool_idx <- !(colnames(OtuData) %in% ControlData$sample.id) &  #矩阵样本不属于metafile中标记的质控样本分组
  !grepl("^[NP]C", colnames(OtuData)) &              #矩阵样本名称不是NC/PC前缀, 因为这些是质控样本
  colnames(OtuData) %in% MetaData$sample.id         #矩阵样本名称需要在metafile中有分组信息
OtuData <- OtuData[, bool_idx]
#转置
OtuData <- t(OtuData)

if (matrixfile != "none") {
  #读取矩阵
  MatrixData <- read.table(matrixfile, header = T, row.names = 1, sep = "\t")
  #matrixfile示例
  #          sample1   sample2   sample3   ...   sampleN
  #sample1   0         0.1111    0.2222    ...   0.4545
  #sample2   0.1111    0         0.1212    ...   0.1234
  #sample3   0.2222    0.1212    0         ...   0.3333
  #...       ...       ...       ...       ...   ...
  #sampleN   0.4545    0.1234    0.3333    ...   0
  #把index替换colname, 因为colname中可能是数值, 变成表头时会自动添加首字母X, 导致不一致
  #并且默认行列的样本名称排列顺序是一致的
  colnames(MatrixData) <- rownames(MatrixData)
}else {
  MatrixData <- vegdist(
    OtuData,
    method = dismethod,
  )
  #dist数据类型转换成matrix数据类型
  MatrixData <- as.matrix(MatrixData)
}

#样本筛选
bool_idx <- !(rownames(MatrixData) %in% ControlData$sample.id) & #矩阵样本不属于metafile中标记的质控样本分组
  !grepl("^[NP]C", rownames(MatrixData)) &                  #矩阵样本名称不是NC/PC前缀, 因为这些是质控样本
  rownames(MatrixData) %in% MetaData$sample.id              #矩阵样本名称需要在metafile中有分组信息
MatrixData <- MatrixData[bool_idx, bool_idx]

#把分组信息顺序与矩阵样本顺序一一对应
grouplist <- c()
for (i in 1:nrow(MatrixData)) {
  groupname <- MetaData$Groups[MetaData$sample.id == rownames(MatrixData)[i]]
  grouplist <- append(grouplist, groupname)
}
sorted_group <- unique(grouplist)
#分组表
groups <- data.frame(rownames(MatrixData), grouplist)
colnames(groups) <- c("V1", "V2")

#计算主坐标
MatrixData <- as.dist(MatrixData) #需要重新把数据转换成dist类型
pcoa <- pcoa(MatrixData, correction = "none", rn = NULL)
PC1 <- pcoa$vectors[, 1] #主坐标1
PC2 <- pcoa$vectors[, 2] #主坐标2

#计算主坐标贡献
pc1 <- round(pcoa$values$Relative_eig[1] * 100, 2)
pc2 <- round(pcoa$values$Relative_eig[2] * 100, 2)

#主图数据整理
plotdata <- data.frame(rownames(pcoa$vectors), PC1, PC2, grouplist)
colnames(plotdata) <- c("sample", "PC1", "PC2", "Group")
#可以自定义levels顺序, 关系到图中分组显示的先后顺序
plotdata$Group <- factor(plotdata$Group, levels = sorted_group)

#输出主坐标绘图结果
write.table(plotdata, outfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, append = FALSE)

##################################################################################
####      PC1和PC2的显著性检验(Tukey 分析)
##################################################################################
# Tukey检验，也称为Tukey的事后多重比较方法，是方差分析（ANOVA）的后续分析中经常采用的一种统计方法
# 它的目标是通过比较各组均值之间的差异来揭示群体之间的显著性差异
# 对没有显著差异的组, 标记相同的重要字母; 存在显著差异的组, 标记不同的重要字母

yf <- plotdata
#获取主坐标每组分类的最大值
yd1 <- yf %>%
  group_by(Group) %>%
  summarize(Max = max(PC1))
yd2 <- yf %>%
  group_by(Group) %>%
  summarize(Max = max(PC2))

# 稍稍增加最大值范围, 使得重要字母标签和绘图箱线有所间隔
yd1$Max <- yd1$Max + max(yd1$Max) * 0.02
yd2$Max <- yd2$Max + max(yd2$Max) * 0.02

# 对主坐标PC1进行Tukey检验
fit1 <- aov(PC1 ~ Group, data = plotdata)
tuk1 <- glht(fit1, linfct = mcp(Group = "Tukey"))
res1 <- cld(tuk1, alpah = 0.05)
#print(res1)

# 对主坐标PC2进行Tukey检验
fit2 <- aov(PC2 ~ Group, data = plotdata)
tuk2 <- glht(fit2, linfct = mcp(Group = "Tukey"))
res2 <- cld(tuk2, alpah = 0.05)

#提取检验的重要字母组成绘图矩阵表
test <- data.frame(
  PC1 = res1$mcletters$Letters,
  PC2 = res2$mcletters$Letters,
  yd1 = yd1$Max,
  yd2 = yd2$Max,
  Group = yd1$Group
)

#分组顺序保持和plotdata一致
test$Group <- factor(test$Group, levels = sorted_group)
# print(test)

##################################################################################
##########  相须图(箱线/盒图) 绘制
##################################################################################
# 绘制的是箱线图, 每个箱子代表数据的4分位范围, 箱子中的水平线表示中位数
# 并且分组两两进行了Tukey显著性检验, 对每组数据标记重要字母(a~z)
# 图中的重要字母(a~z)表示每一组数据, 相同字母的组之间没有显著差异, 不同字母的组存在显著差异
p1 <- ggplot(plotdata, aes(Group, PC1)) +
  geom_boxplot(aes(fill = Group)) +
  geom_text(data = test, aes(x = Group, y = yd1, label = PC1),
            size = 7, color = "black", fontface = "bold") +
  coord_flip() +
  theme_bw() + #ylim(1.5*min(plotdata$PC1), 1.5*max(plotdata$PC1)) +
  theme(axis.ticks.length = unit(0.4, "lines"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(color = "black", size = 20, face = "bold"),
        axis.text.x = element_blank(),
        legend.position = "none")

p3 <- ggplot(plotdata, aes(Group, PC2)) +
  geom_boxplot(aes(fill = Group)) +
  geom_text(data = test, aes(x = Group, y = yd2, label = PC2),
            size = 7, color = "black", fontface = "bold") +
  theme_bw() + #ylim(1.5*min(plotdata$PC2), 1.5*max(plotdata$PC2)) +
  theme(axis.ticks.length = unit(0.4, "lines"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "black", size = 20, face = "bold",
                                   vjust = 1, hjust = 1, angle = 45),
        axis.text.y = element_blank(),
        legend.position = "none")

########################### PCOA 图绘制  #########################
p2 <- ggplot(plotdata, aes(PC1, PC2)) +
  geom_point(aes(fill = Group), size = 8, pch = 21) +
  xlab(paste("PCo1 (", pc1, "%", ")", sep = "")) +
  ylab(paste("PCo2 (", pc2, "%", ")", sep = "")) +
  theme(text = element_text(size = 30)) +
  geom_vline(aes(xintercept = 0), linetype = "dotted") +
  geom_hline(aes(yintercept = 0), linetype = "dotted") +
  theme(panel.background = element_rect(fill = "white", color = "black"),
        panel.grid = element_blank(),
        axis.title = element_text(color = "black", size = 34),
        axis.ticks.length = unit(0.4, "lines"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"),
        axis.title.x = element_text(color = "black", size = 34, vjust = 0),
        axis.title.y = element_text(color = "black", size = 34, vjust = 0),
        axis.text = element_text(color = "black", size = 24),
        legend.title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 20),
        legend.key = element_blank(),
        legend.position = c(1, 1), #图例的坐标位置
        legend.justification = c(1, 1),
        legend.background = element_rect(color = 'black'),
        legend.key.height = unit(1, "cm")) +
  guides(fill = guide_legend(ncol = 1)) +
  stat_ellipse(data = plotdata,
               geom = "polygon",
               linewidth = 0.5,
               aes(fill = Group),
               alpha = 0.2,
  )
#使得相须盒图与Pcoa图的对齐坐标(反过来对齐可能会导致pcoa的置信椭圆绘制不完整)
p1 <- p1 + ylim(ggplot_build(p2)$layout$panel_scales_x[[1]]$
                  range$
                  range)
p3 <- p3 + ylim(ggplot_build(p2)$layout$panel_scales_y[[1]]$
                  range$
                  range)

##################### permanova 分析结果绘制 ###############
# PERMANOVA是多元方差分析的非参数变体。它用来比较多组观测样本的统计指标值的异同。
# 它利用距离矩阵(如欧式距离、Bray-Curtis距离)对总方差进行分解, 分析不同分组因素或不同环境因子对样品差异的解释度, 并使用置换检验对各个变量解释的统计学意义进行显著性分析。
# 目的是检测不同分组的响应变量如菌群构成是否有显著差异。因主要用函数adonis进行分析，有时也称为adonis检验。
# 需要注意的是, 分组group列不能是数值类型, 因为R脚本中会把数值只分为0和非0两组, 自由度df就一直等于1, 会与实际分组不符合
# Diffs：差异比较组
# Df：自由度
# SumsOfSqs：总方差，又称离差平方和
# MeanSqs ：均方（差），即Sums Of Sqs/Df
# Fvalue：F 检验值
# R2：不同分组对样品差异的解释度,(0~1)范围, 即分组方差与总方差的比值，R2 越大表示分组对差异的解释度越高, 差异越大# Pvalue：P值，值越小可信度越高

otu.adonis = adonis(OtuData ~ V2, data = groups, method = dismethod)

p4 <- ggplot(plotdata, aes(PC1, PC2)) +
  geom_text(aes(x = -0.5,
                y = 0.6,
                label = paste("PERMANOVA:",
                              "\nMethod = ",
                              dismethod,
                              "\ndf = ",
                              otu.adonis$aov.tab$Df[1],
                              "\nR2 = ",
                              round(otu.adonis$aov.tab$R2[1], 4),
                              "\np-value = ",
                              otu.adonis$aov.tab$"Pr(>F)"[1],
                              sep = "")),
            size = 6,
  ) +
  theme_bw() +
  xlab("") +
  ylab("") +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()
  )

# 拼接图像
p5 <- p1 +
  p4 +
  p2 +
  p3 +
  plot_layout(heights = c(1, 4), widths = c(4, 1), ncol = 2, nrow = 2)

# 保存图片, 可以是pdf,png
ggsave(outpng, p5, height = 12, width = 15, dpi = 300)
ggsave(outpdf, p5, height = 12, width = 15, dpi = 300)
