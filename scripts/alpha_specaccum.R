#vegan 包是进行群落数据分析最常用的R包，其中的 specaccum 函数用来计算物种的累计曲线
library(vegan)
args <- commandArgs(T)

otu_table <- args[1]
out_img   <- args[2]

rare = read.table(otu_table, header=T, row.names=1, sep="\t", comment.char="",skip=1) #feature-table

############################################################物种累积曲线
sp1=specaccum(t(rare),method = "random")
#pdf(args[3])
library(tools)
if (file_ext(out_img) == "pdf"){
     pdf(out_img)
}else if (file_ext(out_img) == "png"){
    #png(out_img, width=10, height=10, units="in", bg="white", res=300)
    png(out_img)
}else{
    stop(cat("Error: 暂不支持输出除了png/pdf之外的文件格式\n"))    
}

plot(sp1, ci.type="polygon", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
boxplot(sp1, col="yellow", add=TRUE, pch="+")

dev.off()
