#关联网络
library(tidyverse)
library(igraph)
library(psych)
args <- commandArgs(T)
pvalue=read.table(args[2],head=F,sep="\t",row.names = 1) 
group=read.table(args[1],header=T, sep="\t", comment.char="", stringsAsFactors = F)  #分组信息
#group=group[-1,]
ptemp=read.table(args[2],head=F,sep="\t") #fastspar生成的pvalue结果
rtemp=read.table(args[3],head=T,comment.char = "",sep="\t") #fastspar生成的相关性结果
colnames(ptemp)=c("ASV",rownames(pvalue))
ptemp$node1 = rownames(pvalue)
rtemp$node1 = rownames(pvalue) 
r = rtemp %>% 
  gather(key = "node2", value = "r", -node1) %>%
  data.frame()
df = ptemp %>% 
  gather(key = "node2", value = "padj", -node1) %>%
  data.frame()
cor.data <- merge(r,df,by=c("node1","node2"))
temp=cor.data
cor.data$linetype=c()
cor.data$linesize=c()
a=1
for(i in 1:nrow(temp)){
  if(temp$node1[i] != temp$node2[i]){
    if(abs(as.numeric(temp$r[i])) >= 0.6){
      if(as.numeric(temp$padj[i]) <= 0.01){
        cor.data[a,1:4]=temp[i,]
        cor.data$linesize[a]=abs(as.numeric(temp$r[i]))
        if(as.numeric(temp$r[i])>0){
          cor.data$linetype[a]="positive"
        }
        else{
          cor.data$linetype[a]="negative"
        }
        a=a+1
      }
    }
  }
}
cor.data=cor.data[1:(a-1),]
vertices=c(as.character(cor.data$node1),as.character(cor.data$node2)) %>%
  as_tibble() %>% group_by(value) %>%  summarize(n=n()) 
vertices=vertices %>% select(-n)   # 因为此处的n不准确，所以先删除。
graph <- graph_from_data_frame(cor.data, vertices = vertices, directed = FALSE )
#vcount(graph) # 节点数目：47
#ecount(graph) # 链接数:124
g1 <- igraph::simplify(graph,
                      remove.multiple = TRUE,
                      remove.loops = TRUE,
                      edge.attr.comb = "first")
layout1 <- layout_in_circle(g1) # 径向布局适合节点较少的数据。
layout2 <- layout_with_fr(g1) # fr布局。
layout3 <- layout_on_grid(g1) # grid布局。
V(g1)$degree <- degree(g1)
pdf(args[4])
plot.igraph(g1,layout=layout2,vertex.label.cex=0.8,vertex.label.dist=0,edge.arrow.size=0.5,shape = 1,vertex.size=V(g1)$degree*2) #layout可选
dev.off()
