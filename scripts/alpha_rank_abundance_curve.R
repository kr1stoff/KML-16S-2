library(BiodiversityR)  #BiodiversityR 包 rankabundance() 实现 OTU 排序
library(ggplot2)  #作图
library(dplyr)

args <- commandArgs(T)

inputfile <- args[1]  #OTU/ASV表 feature-table
outprefix <- args[2]  #输出图片文件前缀
outpng  <- paste(outprefix, ".png", sep="")
outpdf  <- paste(outprefix, ".pdf", sep="")

# 读取OTU表
data <- read.table(inputfile, 
    row.names = 1,
    header = T,
    sep = "\t",
    check.names = F,
    skip =1,
    comment.char = "%")
# 转置
data <- t(data)
#print(head(data))

# 计算相对丰度
rdata <- data / rowSums(data)

# 每个样本的OTU?ASV 丰度进行从大到小排序, 丰度等级rank为排序后的序号
# 先初始化数据列
rank_ranks      <- c()
rank_abundances <- c()
rank_samples    <- c()
samples <- rownames(rdata) 
for (sample in samples){
    #print(sample)
    #筛选样本数据
    sampledata <- subset(rdata, rownames(rdata)==sample)
    #丰度排序划分Rank, 数据保留6位小数
    sampleRankData <- data.frame(rankabundance(sampledata, digits=6))
    #rankdata <- rbind(rankdata, sampleRankData)
    #rankdata <- bind_rows(rankdata, sampleRankData)
    rank_ranks      <- c(rank_ranks, sampleRankData$rank)
    rank_abundances <- c(rank_abundances, sampleRankData$abundance)
    rank_samples    <- c(rank_samples, rep(sample, nrow(sampleRankData)))
}

# 合并成绘图的丰度等级数据
rankdata <- data.frame(
    rank = rank_ranks, 
    abundance = rank_abundances, 
    sample = rank_samples)
# 去除0值
rankdata <- subset(rankdata, abundance !=0)
print(head(rankdata))

#获取样本最后的丰度等级的丰度, 在此位注释样本信息
plotAnno <- rankdata %>% group_by(sample) %>% slice(which.max(rank))
print(head(plotAnno))

#print(nrow(rankdata))
#绘制等级丰度曲线
p <- ggplot(rankdata, aes(x = rank, y = log(abundance, 10), color = sample)) + 
    geom_step(direction = 'vh') + #绘制梯度线, vh表示先水平后垂直
    geom_text(dat = plotAnno,
        aes(x = rank, y = log(abundance, 10), label = sample, color = sample),
        size = 2,
        angle = 90,
        ) + 
    xlab("OTUs/ASVs Rank") + 
    ylab("Log(abundance(%), 10)") +
    theme_classic(base_line_size = 1) +
    theme(
        #panel.grid = element_blank(), #无网格
        #panel.background = element_rect(fill = "white", color = "black"), #背景纯白黑框
        legend.position = "none", #不添加图例, 因为样本太多会导致挤压绘图区空间
    )

#ggsave(outpdf, p, height=10, width=10, dpi=300)
#ggsave(outpng, p, height=10, width=10, dpi=300)
ggsave(outpdf, p)
ggsave(outpng, p)

