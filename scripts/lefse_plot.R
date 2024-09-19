pacman::p_load(tidyverse,microeco,magrittr)
args <- commandArgs(T)
feature_table <- read.table(args[1],sep="\t", row.names = 1,skip = 1,comment.char = "%",header=T) #feature-table
sample_table <- read.table(args[2],sep="\t", row.names = 4,header=T) #sample-group
tax_table <- read.table(args[3],sep="\t",header=T, row.names = 1) #taxonomy_v/metadata.tsv
split <- strsplit(tax_table$V2, ";")
tax_table$V2=rep("d__",nrow(tax_table))
tax_table$V3=rep("p__",nrow(tax_table))
tax_table$Class=rep("c__",nrow(tax_table))
tax_table$Order=rep("o__",nrow(tax_table))
tax_table$Family=rep("f__",nrow(tax_table))
tax_table$Genus=rep("g__",nrow(tax_table))
colnames(tax_table)[1]="Domain"
colnames(tax_table)[2]="Phylum"

for(i in 1:length(split)){
    temp=unlist(split[[i]])
    for(j in 1:length(temp)){
        if(length(grep("d__",temp[j]))>0){
            tax_table[i,1]=temp[j]
        }
        if(length(grep("p__",temp[j]))>0){
            tax_table[i,2]=temp[j]
        }
        if(length(grep("c__",temp[j]))>0){
            tax_table[i,3]=temp[j]
        }
        if(length(grep("o__",temp[j]))>0){
            tax_table[i,4]=temp[j]
        }
        if(length(grep("f__",temp[j]))>0){
            tax_table[i,5]=temp[j]
        }
        if(length(grep("g__",temp[j]))>0){
            tax_table[i,6]=temp[j]
        }
    }
}

dataset <- microtable$new(sample_table = sample_table, otu_table = feature_table, tax_table = tax_table)
lefse <- trans_diff$new(dataset = dataset,  method = "lefse", group = "groups", alpha = 0.5,lefse_subgroup = NULL)
#head(lefse$res_diff)
bar=lefse$plot_diff_bar(use_number = 1:30, width = 0.8, group_order = c("Case", "Control")) +
  ggsci::scale_color_npg() + ggsci::scale_fill_npg()
  ggsave(args[4],bar)
# 需要调用ggtree包
  p=lefse$plot_diff_cladogram(use_taxa_num = 150, use_feature_num = 50, clade_label_level = 5, group_order = c("Case", "Control"))
  ggsave(args[5],p)
#图片美化,展示前200个分类单元和前50个特征
#use_labels <- c("c__Deltaproteobacteria", "c__Actinobacteria", "o__Rhizobiales", "p__Proteobacteria", "p__Bacteroidetes", "o__Micrococcales", "p__Acidobacteria", "p__Verrucomicrobia", "p__Firmicutes", "p__Chloroflexi", "c__Acidobacteria", "c__Gammaproteobacteria", "c__Betaproteobacteria", "c__KD4-96","c__Bacilli", "o__Gemmatimonadales","f__Gemmatimonadaceae", "o__Bacillales", "o__Rhodobacterales")
# then use parameter select_show_labels to show
#lefse$plot_diff_cladogram(use_taxa_num = 200, use_feature_num = 50, select_show_labels = use_labels)

