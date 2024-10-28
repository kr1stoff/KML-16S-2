#adonis差异分析
suppressMessages(library(vegan))
suppressMessages(library(dplyr))

args <- commandArgs(T)

featurefile <- args[1] #OTUs/ASVs 特征表 feature-table.tsv
metafile <- args[2] #样本分组信息文件
outfile <- args[3]
groupby <- args[4]


#读取OTUs/ASVs特征表
OtuData <- read.table(featurefile,
                      header = T,
                      sep = "\t",
                      check.names = F,
                      row.names = 1,
                      skip = 1,
                      comment.char = "")
OtuData <- t(OtuData)
#print(head(OtuData))

#读取样本信息(默认样本都在MetaData中)
MetaData <- read.table(metafile,
                       header = T,
                       sep = "\t",
                       check.names = F,
                       row.names = 1,
                       comment.char = "")

# MetaData <- MetaData %>% filter(visitPoints != c("提前终止"))
# MetaData$Groups <- paste(MetaData[[groupby]], MetaData$visitPoints, sep="-")
# MetaData <- MetaData %>% filter(Groups %in% c("Case-Ⅱ-V1", "Control-Ⅱ-V1"))
# MetaData <- MetaData %>% filter(Groups %in% c("Control-Ⅱ-V6", "Control-Ⅱ-V1"))
# MetaData <- MetaData %>% filter(Groups %in% c("Case-Ⅱ-V6", "Case-Ⅱ-V1"))
# MetaData <- MetaData %>% filter(Groups %in% c("Case-Ⅱ-V6", "Control-Ⅱ-V6"))

bool_idx <- rownames(OtuData) %in% rownames(MetaData)
OtuData <- OtuData[bool_idx,]

# 筛选分组列, 滨崎重命名为固定名称Group, 方便后续调用
#GroupData <- subset(MetaData, select = groupby)
GroupData <- subset(MetaData, select = "group")
colnames(GroupData) <- "Group"
head(GroupData)

# 对分组筛选
# GroupData <- GroupData %>% filter(Group != c("提前终止"))

#保证分组和OtuData数据行一致(可能删掉部分不在group中的样本)
OtuData <- OtuData[rownames(OtuData) %in% rownames(GroupData),]
head(OtuData)

# 对分组做差异分析
ado <- adonis2(OtuData ~ Group, GroupData)
print(ado)

# 判定显著性等级
sig <- "NS" #默认NS( no significant), 不显著
if (ado$`Pr(>F)`[1] <= 0.0001) {
  sig <- "****"
}else if (ado$`Pr(>F)`[1] <= 0.001) {
  sig <- "***"
}else if (ado$`Pr(>F)`[1] <= 0.01) {
  sig <- "**"
}else if (ado$`Pr(>F)`[1] <= 0.05) {
  sig <- "*"
}else { }

outputdata <- data.frame(
  Group = paste(unique(GroupData$Group), sep = ",", collapse = ","),
  Df = ado$Df[1],
  SumOfSqs = ado$SumOfSqs[1],
  R2 = ado$R2[1],
  Fvalue = ado$F[1],
  Pvalue = ado$`Pr(>F)`[1], #p value
  Significant = sig
)
#保存结果
write.table(outputdata,
            outfile,
            row.names = F,
            col.names = T,
            quote = F,
            sep = "\t"
)
