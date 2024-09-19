library(ggplot2)
library(stringr) #字符串匹配。 str_extract, str_extract_all
library(dplyr)
library(patchwork)

args <- commandArgs(T)

inputfile = args[1] #qiime2 rarefraction的抽样结果表
groupby   = args[2] #分组表头
outprefix = args[3] #输出文件前缀, 最终输出png和pdf
outpng  = paste(outprefix, ".png", sep="")
outpdf  = paste(outprefix, ".pdf", sep="")

#读取rarefraction的结果表, 否则得自己输出otu去做抽样计算
data <- read.table(inputfile, 
    header = T, 
    sep=",",
    row.names = 1,
    check.names = F,)
#inputfile 示例
#sample-id,depth-1_iter-1,depth-1_iter-2,.....
#sample1,1,1...
#sample2,1,1...
#sampleN,1,1...

#获取抽样深度
depths <- c()
depthheaders <- c()
for(colheader in colnames(data)){
    if(grepl("depth-\\d+_iter-\\d+", colheader)){ #匹配深度抽样和迭代次数得表头
        dpvalues <- str_extract_all(colheader, "\\d+") #获取深度和迭代次数
        depth <- dpvalues[[1]][1]  #深度
        iter  <- dpvalues[[1]][2]  #迭代次数
        depths <- c(depths, depth)
        depthheaders <- c(depthheaders, colheader)
    }
}
#去重
depths <- unique(depths)
print(depthheaders)

# 初始值depth和指数值都赋值为0
all_depths <- c() #rep(0, nrow(data))
all_means  <- c() #rep(0, nrow(data))
all_samples<- c() #rownames(data)
all_groups <- c() #rownames(data)

for(depth in depths){
    target_depth_col_idx  <- grepl(paste("depth-", depth, "_iter-", sep=""), colnames(data))
    target_depth_data <- data[, target_depth_col_idx] #深度depth所有迭代值
    target_depth_mean <- rowMeans(target_depth_data)  #迭代值平均
    #print(depth)
    #print(target_depth_data)
    #print(target_depth_mean)
    all_depths  <- c(all_depths,  rep(as.numeric(depth), nrow(data)))
    all_means   <- c(all_means,   target_depth_mean)
    all_samples <- c(all_samples, rownames(data))
    #all_groups  <- c(all_groups,  paste(data[[groupby]]))
    all_groups  <- c(all_groups, paste(data[[groupby]], data$visitPoints, sep="-"))
}

plotdata <- data.frame(all_depths, all_means, all_samples, all_groups)
colnames(plotdata) <- c("Depth", "Value", "Sample", "Group")
# 筛选
plotdata <- plotdata[!grepl("提前终止", plotdata$Group),]
print(head(plotdata))

#去除空值NA
plotdata <- na.omit(plotdata)

#获取Group(样本)分组的Depth最大值数据，用于在此位置标注Group(样本)分组标签
#plotdata %>% group_by(Group) %>% summarize(Depth_max=max(Depth, na.rm = TRUE))
#plotdata %>% group_by(Group) %>% filter(Depth == max(Depth, na.rm = TRUE))
plotAnno <- plotdata %>% group_by(Sample) %>% slice(which.max(Depth))
#plotAnno <- aggregate(plotdata, by=plotdata["Sample"], FUN=max)

#绘制样本稀疏曲线
p1 <- ggplot(plotdata) + 
    geom_line(
        dat = plotdata, 
        aes(x = Depth, y = Value, group = Sample, color = Sample),
        #color = as.numeric(plotdata$Group),
        #size = 2
        ) + 
    ggtitle("Observed") + 
    geom_text(dat = plotAnno, #在最大深度位置, 添加样本分组注释标签
        aes(x=Depth, y=Value, label = Sample, color = Sample),
        size = 3) +
    theme_classic(base_line_size =1) +
    theme(
        #legend.key.size = unit(0.5, "cm"),
        #legend.position = "right",
        #legend.text = element_text(size = 8),
        legend.position = "none", #无需图例
    )

################## 分组统计绘图  ###################
plotgroupdata <- plotdata %>% group_by(Depth, Group) %>% summarize(meanvalue = mean(Value), se=sd(Value))
p2 <- ggplot(plotgroupdata) + 
    geom_line(
        dat = plotgroupdata,
        aes(x = Depth, y = meanvalue, group = Group, color = Group)
    ) +
    geom_errorbar(aes(x = Depth, ymax = meanvalue + se, ymin = meanvalue - se, color = Group), size=0.5) +
    theme_classic(base_line_size =1) +
    theme(
        legend.position = c(1,1),
        legend.justification = c(1,1),
        axis.title.y = element_blank(),
    )
#使p2和p1的y轴坐标一致
p2 <- p2 + ylim(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range)

#图片布局
p <- p1 + p2 + plot_layout(heights = c(1,1), widths = c(1, 1), ncol=2, nrow=1)

ggsave(outpdf, p, height=8, width=16, dpi=300)
ggsave(outpng, p, height=8, width=16, dpi=300)
#ggsave(outpdf, p)
#ggsave(outpng, p)

