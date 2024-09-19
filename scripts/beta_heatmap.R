#热图绘制
suppressMessages(library(pheatmap))
suppressMessages(library(tools))
suppressMessages(library(dplyr))

args <- commandArgs(T)

matrixfile <- args[1] #距离矩阵文件
metafile <- args[2] #分组信息文件
groupby <- args[3] #分组表头列
outimg <- args[4] #输出文件


#读取距离矩阵文件
MatrixData <- read.table(matrixfile,
                         header = T,
                         row.names = 1,
                         sep = "\t",
                         check.names = F)
#MatrixData <- MatrixData[1:50, 1:50]

#读取分组信息文件
MetaData <- read.table(metafile,
                       header = T,
                       row.names = 1,
                       sep = "\t",
                       check.names = F)
#head(MetaData)
MetaData$Groups <- MetaData[[groupby]]

# matrix样本筛选
bool_idx <- rownames(MatrixData) %in% rownames(MetaData)
MatrixData <- MatrixData[bool_idx, bool_idx]

#按照样本顺序合并距离矩阵和分组信息表
grouplist <- c()
for (spname in rownames(MatrixData)) {
  #print(spname)
  #groupname <- MetaData[spname, groupby]
  groupname <- MetaData[spname, "Groups"]
  #print(groupname)
  grouplist <- c(grouplist, groupname)
}

annotation <- data.frame(Group = grouplist)
colnames(annotation) <- c(groupby)
rownames(annotation) <- rownames(MatrixData)
#head(annotation)

# 绘图
if (file_ext(outimg) == "pdf") {
  pdf(outimg)
}else if (file_ext(outimg) == "png") {
  #png(out_img, width=10, height=10, units="in", bg="white", res=300)
  png(outimg)
}else {
  stop(cat("Error: 暂不支持输出除了png/pdf之外的文件格式\n"))
}
pheatmap(MatrixData,
         annotation_col = annotation,
         fontsize_col = min(c(15, 350 / ncol(MatrixData))), #设置数据列(横轴)字体, 最大15
         fontsize_row = min(c(15, 400 / nrow(MatrixData)))  #设置数据行(纵轴)字体，最大15
)
dev.off()
