#载入分析包
pacman::p_load(tidyverse, microeco, magrittr, data.table, aplot, ggtree)
library(patchwork)
# 载入数据
args <- commandArgs(T)

otufile <- args[1]
taxfile <- args[2]
metafile <- args[3]
groupby <- args[4]
outprefix <- args[5]
outpng_lda <- paste0(outprefix, ".LDAbar.png")
outpdf_lda <- paste0(outprefix, ".LDAbar.pdf")
outpng_clad <- paste0(outprefix, ".cladogram.png")
outpdf_clad <- paste0(outprefix, ".cladogram.pdf")

feature_table <- read.table(otufile,
                            sep = "\t",
                            row.names = 1,
                            skip = 1,
                            comment.char = "",
                            header = T,
                            check.names = F)

sample_table <- read.table(metafile,
                           sep = "\t",
                           row.names = 1,
                           #skip = 1,
                           comment.char = "",
                           header = T,
                           check.names = F)

tax_table <- read.table(taxfile,
                        sep = "\t",
                        row.names = 1,
                        #skip = 1,
                        comment.char = "",
                        header = T,
                        check.names = F)

sample_table$Groups <- sample_table[[groupby]]
feature_table <- feature_table[, colnames(feature_table) %in% rownames(sample_table)]

head(feature_table[, 1:4])
head(sample_table)
head(tax_table)

#stop("xxx")

# 对Taxon列分割
taxdf <- data.frame(
  domain = rep(NA, nrow(tax_table)),
  phyml = rep(NA, nrow(tax_table)),
  class = rep(NA, nrow(tax_table)),
  order = rep(NA, nrow(tax_table)),
  family = rep(NA, nrow(tax_table)),
  genus = rep(NA, nrow(tax_table)),
  species = rep(NA, nrow(tax_table)))
taxSplit <- strsplit(tax_table$Taxon, ";")
for (i in seq_along(taxSplit)) {
  rowtax <- taxSplit[[i]]
  rowtax <- unlist(rowtax)  #变成一维列表
  rowtax <- trimws(rowtax)  #去除首位空白符

  domain <- rowtax[grepl("d__", rowtax)]
  if (length(domain) > 0) { taxdf[i, "domain"] <- domain }
  phyml <- rowtax[grepl("p__", rowtax)]
  if (length(phyml) > 0) { taxdf[i, "phyml"] <- phyml }
  class <- rowtax[grepl("c__", rowtax)]
  if (length(class) > 0) { taxdf[i, "class"] <- class }
  order <- rowtax[grepl("o__", rowtax)]
  if (length(order) > 0) { taxdf[i, "order"] <- order }
  family <- rowtax[grepl("f__", rowtax)]
  if (length(family) > 0) { taxdf[i, "family"] <- family }
  genus <- rowtax[grepl("g__", rowtax)]
  if (length(genus) > 0) { taxdf[i, "genus"] <- genus }
  species <- rowtax[grepl("s__", rowtax)]
  if (length(species) > 0) { taxdf[i, "species"] <- species }
}
rownames(taxdf) <- rownames(tax_table)
#tax_table <- data.frame(tax_table, taxdf)
head(taxdf)

# 创建microtable对象
dataset <- microtable$new(
  sample_table = sample_table,
  otu_table = feature_table,
  tax_table = taxdf
)

pvalue <- 0.05

#执行lefse分析
lefse <- trans_diff$new(
  dataset = dataset,
  method = "lefse",
  #group = "group",
  group = "Groups",
  taxa_level = "all",
  alpha = pvalue,
  lefse_subgroup = NULL,
  p_adjust_method = "none"
)

print(lefse)

p1 <- lefse$plot_diff_bar(
  use_number = 1:30,
  width = 0.8,
  #group_order = unique(sample_table$group)
  group_order = unique(sample_table$Groups)
) +
  labs(title = paste0("Lefse LDA Score (p-value <", pvalue, ")")) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(outpng_lda, p1, height = 7, width = 14)
ggsave(outpdf_lda, p1, height = 7, width = 14)

p2 <- lefse$plot_diff_cladogram(
  use_taxa_num = 200,
  use_feature_num = 30,
  clade_label_level = 7,
  group_order = unique(sample_table$Groups)
) +
  labs(title = paste0("Lefse Cladogram(p-value <", pvalue, ")")) +
  theme(plot.title = element_text(hjust = 0.5))


#p <- p1 + p2 + plot_layout(widths = c(1,1), ncol=2, nrow=1)

ggsave(outpng_clad, p2, height = 7, width = 14)
ggsave(outpdf_clad, p2, height = 7, width = 14)

