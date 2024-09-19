suppressMessages(library(readr))
suppressMessages(library(ggpicrust2))
suppressMessages(library(tibble))
suppressMessages(library(tidyverse))
suppressMessages(library(GGally))
suppressMessages(library(ggh4x))
suppressMessages(library(ggprism))
suppressMessages(library(patchwork))
suppressMessages(library(dplyr))

args <- commandArgs(T)

meta_file <- args[1]
abundance_file <- args[2]
groupby <- args[3]
outprefix <- args[4]
#outpng <- paste(outprefix, ".png", sep="")
#outpdf <- paste(outprefix, ".pdf", sep="")

metadata <- read_delim(
  meta_file,
  delim = "\t",
  comment = "#",
  escape_double = FALSE,
  trim_ws = TRUE
)

metadata$group <- metadata[[groupby]]


#head(metadata)

metacyc_abundance <- read_delim(
  abundance_file,
  delim = "\t"
)
head(metacyc_abundance)

#colnames(metacyc_abundance)
#rownames(metadata)

bool_idx <- colnames(metacyc_abundance) %in% c('pathway', metadata$`sample-id`)
print(bool_idx)
metacyc_abundance <- metacyc_abundance[, bool_idx]
head(metacyc_abundance)


#rowSums(metacyc_abundance[, 3:ncol(metacyc_abundance)-1])
#stop("xxx")

metacyc_daa_results_df <- pathway_daa(
  abundance = metacyc_abundance %>% column_to_rownames("pathway"),
  metadata = metadata,
  group = "group",
  daa_method = "LinDA"
)

metacyc_daa_annotated_results_df <- pathway_annotation(
  pathway = "MetaCyc",
  daa_results_df = metacyc_daa_results_df,
  ko_to_kegg = F)

head(metacyc_daa_annotated_results_df)
tail(metacyc_daa_annotated_results_df)

#select_pathway <- metacyc_daa_annotated_results_df[metacyc_daa_annotated_results_df$p_adjust < 0.95, ]$feature
#print(select_pathway[1:3])

#p1 <- pathway_errorbar(
#    abundance = metacyc_abundance %>% column_to_rownames("pathway"),
#    daa_results_df = metacyc_daa_annotated_results_df,
#    Group = groupby,
#    p_values_threshold = 0.95,
#    p_value_bar = T,
#    select = select_pathway[1:30],
#    x_lab = "description"
#    )

#ggsave(outpng, p1)
#ggsave(outpdf, p1)

print("############################################")

################## heatmap ################
feature_with_p <- metacyc_daa_results_df %>% filter(p_values < 0.9)
head(feature_with_p$feature)

pfilter_metacyc_abundance <- metacyc_abundance[metacyc_abundance$pathway %in% feature_with_p$feature,]
head(pfilter_metacyc_abundance)
print(pfilter_metacyc_abundance$pathway[1:20])

#stop("xx")
xxx <- metacyc_abundance %>%
  right_join(metacyc_daa_annotated_results_df %>% select(all_of(c("feature", "description"))), by = c("pathway" = "feature")) %>%
  filter(pathway %in% feature_with_p$feature) %>%
  select(-"pathway") %>%
  column_to_rownames("description")
#print(head(xxx[, 1:4]))
#print(rowSums(xxx))

p_heatmap <- pathway_heatmap(
  #abundance = pfilter_metacyc_abundance %>% column_to_rownames("pathway"),
  abundance = metacyc_abundance %>%
    right_join(
      metacyc_daa_annotated_results_df %>%
        select(all_of(c("feature", "description"))), by = c("pathway" = "feature")) %>%
    filter(pathway %in% feature_with_p$feature) %>%
    select(-"pathway") %>%
    column_to_rownames("description"),
  metadata = metadata,
  group = "group"
)

ggsave(paste(outprefix, ".heatmap.pdf", sep = ""), p_heatmap, height = 7, width = 14)
ggsave(paste(outprefix, ".heatmap.png", sep = ""), p_heatmap, height = 7, width = 14)


#p_pca <- pathway_pca(
#    abundance = metacyc_abundance %>% column_to_rownames("pathway"),
#    metadata = metadata,
#    group = "group"
#)
#ggsave(paste(outprefix, ".pca.pdf", sep=""), p_pca)
#ggsave(paste(outprefix, ".pca.png", sep=""), p_pca)

p_errorbar <- pathway_errorbar(
  abundance = metacyc_abundance %>% column_to_rownames("pathway"),
  daa_results_df = metacyc_daa_annotated_results_df,
  Group = metadata$group,
  ko_to_kegg = FALSE,
  p_values_threshold = 1,
  order = "group",
  select = pfilter_metacyc_abundance$pathway[1:30],
  p_value_bar = T,
  colors = NULL,
  x_lab = "description"
)

ggsave(paste(outprefix, ".errorbar.pdf", sep = ""), p_errorbar, height = 7, width = 14)
ggsave(paste(outprefix, ".errorbar.png", sep = ""), p_errorbar, height = 7, width = 14)
