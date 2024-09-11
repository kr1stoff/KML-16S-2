 #调用suppressMessages载入包, 忽略包载入时的messages信息
suppressMessages(library(ggplot2))
suppressMessages(library(ggtree))
suppressMessages(library(treeio))
library(tools)
library(optparse)

#定义脚本参数
option_list <- list(
    make_option(c("-i", "--input"), type="character", action="store", default="tree.newick",
        help = "输入进化树newick格式文件, default:%default"),
    make_option(c("-o", "--outfile"), type="character", action="store", default="tree.png",
        help = "输出文件, 支持png和pdf两种后缀格式, default:%default"),
    make_option(c("-l", "--layout"), type="character", action="store", default="circular",
        help = "进化树绘图模式, default: circular(环形树), defaullt:%default。其余可选:rectangular, dendrogram, slanted, ellipse, roundrect, fan, inward_circular, radial, equal_angle, daylight, ape"),    
    make_option(c("-b", "--branchlength"), type="character", action="store", default="branch.length",
        help = "根据进化距离显示树枝长度, default:%default, 其余值: branch.length(根据距离绘制树枝)")
)

# 解析参数, 参数存储在字典哈希对象 opt中
opt = parse_args(
    OptionParser(
        usage = "usage: %prog [options]",
        option_list = option_list,
        add_help_option = TRUE,
        prog = NULL,
        description = "\v用ggtree绘制系统发育树图\n"
    )
)
#print(opt$input)
# 参数空值/不存在异常报错提醒
if(!file.exists(opt$input)){
    stop({cat("Error: the argument -i[--input] newick文件路径不存在\n")})
}

#用treeio的read.newick函数读取newick树文件
tree <- read.newick(opt$input, node.label="label")
#print(tree$tip) #标签
#print(tree$node) #节点距离

#print(tree$tip[10])
tree$tip <- gsub("'", "", tree$tip) #去除标签中的单引号
tree$tip.label <- gsub("'", "", tree$tip.label) #这个才是有效去除单引号
#print(tree$tip[10])
#print(tree$tip.label)
leaf_counts <- length(tree$tip)  #标签计数


#用tools中的file_ext函数获取文件后缀
outfile_ext = file_ext(opt$outfile)
width_d = 10
height_d = 10
#if(is.null(grep('rectangular', opt$layout))){
#    width_d  = 10
#    height_d = 10
#}else{
#    width_d  = 10
#    height_d = 20
#}

if(outfile_ext == "png"){    
    png(opt$outfile, width=width_d, height=height_d, units="in", bg="white", res=300)    
}else if (outfile_ext == "pdf"){
    pdf(opt$outfile, width=width_d, height=height_d)
}else{
    stop(cat("Error: 暂不支持输出除了png/pdf之外的文件格式\n"))
}

#根据tip/leaf数目自动设置字体大小
leaf_txt_size <- 2*360/leaf_counts

if(is.null(opt$branchlength)){opt$branchlength="none"}
ggtree(
    tree, 
    branch.length = opt$branchlength, 
    layout = opt$layout, #one of 'rectangular', 'dendrogram', 'slanted', 'ellipse','roundrect', 'fan', 'circular', 'inward_circular', 'radial','equal_angle', 'daylight' or 'ape'
    ) + 
geom_tiplab(size=leaf_txt_size, color="seagreen")

dev.off()
