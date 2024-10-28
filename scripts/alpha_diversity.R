library(ggplot2)
library(ggpubr)
library(dplyr)

args <- commandArgs(T)

metafile <- args[1]
diversityfile <- args[2]
DivTag <- args[3]
outprefix <- args[4]

outpng <- paste0(outprefix, ".", DivTag, ".png")
outpdf <- paste0(outprefix, ".", DivTag, ".pdf")
print(outpng)
print(outpdf)

# 读取样本指数文件
AlphaData <- read.table(diversityfile,
    header = T,
    sep = "\t",
    row.names = 1,
    check.names = F
)

# 读取样本分组信息文件
MetaData <- read.table(metafile,
    header = T,
    sep = "\t",
    check.names = F,
    row.names = 1,
)

# 按照左侧AlphaData行号合并MetaFata
data <- merge(AlphaData, MetaData, by = "row.names", all.x = T)
# 筛选异常值
# data <- data %>% as.data.frame() %>% filter(visitPoints != c("提前终止"))
# data$Group <- paste(data$group, data$visitPoints, sep="-")
data$Group <- data$group
data$DrawValue <- data[[DivTag]]
head(data)

compare_list <- list(unique(data$Group))
# compare_list <- list(
#    c("Case-Ⅱ-V6", "Case-Ⅱ-V1"),
#    c("Control-Ⅱ-V6", "Control-Ⅱ-V1"),
#    c("Case-Ⅱ-V6", "Control-Ⅱ-V6"),
#    c("Case-Ⅱ-V1", "Control-Ⅱ-V1")
# )
compare_list <- combn(unique(data$Group), 2, simplify = F)
# print(compare_list)

p <- ggplot(data, aes(x = Group, y = DrawValue, fill = Group, color = Group)) +
    geom_boxplot(
        width = 0.5,
        alpha = 0.6,
        lwd = 1.15,
        outlier.shape = NA
    ) +
    geom_jitter(
        width = 0.3,
        size = 3,
        alpha = 0.75
    ) +
    labs(y = DivTag) +
    # theme_minimal() +
    # theme_classic() +
    theme_bw() +
    theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.1),
    ) +
    stat_compare_means(
        comparisons = compare_list,
        method = "wilcox.test",
        label = "p.signif",
        hide.ns = F
    )

ggsave(outpng, p)
ggsave(outpdf, p)
