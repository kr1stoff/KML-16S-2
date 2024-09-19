#载入分析包
pacman::p_load(tidyverse, microeco, magrittr, data.table, aplot)
# 载入数据
args <- commandArgs(T)

otufile <- args[1]
taxfile <- args[2]
metafile <- args[3]
groupby <- args[4]
outprefix <- args[5]
outpng <- paste0(outprefix, ".png")
outpdf <- paste0(outprefix, ".pdf")

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
#tax_table <- data.frame(tax_table, taxdf)
rownames(taxdf) <- rownames(tax_table)
head(taxdf)

# 创建microtable对象
dataset <- microtable$new(
  sample_table = sample_table,
  otu_table = feature_table,
  #tax_table = tax_table
  tax_table = taxdf
)

pvalue <- 0.05
p_adjust_m <- "none"

#执行随机森林分析
rf <- trans_diff$new(
  dataset = dataset,
  method = "rf",
  #group = "group",
  group = "Groups",
  taxa_level = "genus",
  alpha = pvalue,
  p_adjust_method = p_adjust_m
)

g1 <- rf$plot_diff_bar(
  use_number = 1:20,
  #group_order = unique(sample_table$group),
  group_order = unique(sample_table$Groups)
)
g2 <- rf$plot_diff_abund(
  use_number = 1:20,
  #group_order = unique(sample_table$group),
  group_order = unique(sample_table$Groups),
  select_taxa = rf$plot_diff_bar_taxa
)
p <- g1 %>% insert_right(g2)

ggsave(outpng, p, height = 7, width = 14)
ggsave(outpdf, p, height = 7, width = 14)
