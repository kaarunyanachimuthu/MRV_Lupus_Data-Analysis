setwd("C:/Users/kaaru/Box/Bigley Lab/Kaarunya Nachimuthu/Nachimuthu Projects/Nachimuthu MRV R848 Spleen scRNAseq vdj/VDJ/tcr")

library(dplyr)
library(ggplot2)
library(scRepertoire)
library(Seurat)
library(tidyverse)
library(dplyr)
library(RColorBrewer)
library(ggrepel)

Mock_Cont <- read.csv("./Mock_Cont_TCR/outs/filtered_contig_annotations.csv")

MRV_Cont <- read.csv("./MRV_Cont_TCR/outs/filtered_contig_annotations.csv")

Mock_R848 <- read.csv("./Mock_R848_TCR/outs/filtered_contig_annotations.csv")

MRV_R848 <- read.csv("./MRV_R848_TCR/outs/filtered_contig_annotations.csv")

contig_list <- list(Mock_Cont, MRV_Cont, Mock_R848, MRV_R848)

combined.TCR <- combineTCR(
  contig_list,
  samples = c("mock_Cont", "MRV_Cont", "mock_R848", "MRV_R848"),
  removeNA = FALSE,
  removeMulti = FALSE,
  filterMulti = FALSE
)


# barplot of cdr3 aa distribution
tab <- clonalLength(combined.TCR, cloneCall="aa", chain = "both", group.by = "sample", exportTable = TRUE)

tab_counts <- tab %>%
  group_by(sample, length) %>%
  summarise(count = n(), .groups = "drop")

tab_percent <- tab_counts %>%
  group_by(sample) %>%
  mutate(percent = count / sum(count) * 100)

ggplot(tab_percent, aes(x = factor(length), y = percent, fill = sample)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "CDR3 length", y = "CDR3 (%)", fill = "Sample", title = ' Percent of CDR3(AA) Length Distribution (TRB chain)') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5)) + scale_fill_manual(
    values = c("red", "blue", "green", "yellow")
  ) 

# Graph of aa length
positionalEntropy(
  combined.TCR,
  chain = "TRB",
  aa.length = 20,
  method = "gini.simpson"
) +
  scale_color_manual(
    values = c(
      "mock_Cont" = "black",
      "MRV_Cont"  = "red",
      "mock_R848" = "blue",
      "MRV_R848"  = "green"
    )
  )

# gini simpson score
# grouped by sample
dp1 <- clonalDiversity(
  combined.TCR,
  cloneCall = "gene",
  group.by = "sample",
  metric = "gini.simpson",
  skip.boots = TRUE
)

dp1$layers[[2]]$aes_params$size <- 5
dp1 +
  scale_fill_manual(
    values = c(
      "mock_Cont" = "black",
      "MRV_Cont"  = "red",
      "mock_R848" = "blue",
      "MRV_R848"  = "green"
    )
  )

# split by sample
dp2 <- clonalDiversity(
  combined.TCR,
  cloneCall = "gene",
  x.axis = "sample",
  metric = "gini.simpson",
  skip.boots = TRUE
)

dp2$layers[[2]]$aes_params$size <- 5
dp2 +
  scale_fill_manual(
    values = c(
      "mock_Cont" = "black",
      "MRV_Cont"  = "red",
      "mock_R848" = "blue",
      "MRV_R848"  = "green"
    )
  )

# Combining with single cell data
tcell <- readRDS('./tcells_updated.rds')

tcell$sample <- tcell$Sample
tcell <- combineExpression(combined.TCR, 
                                   tcell, 
                                   cloneCall="gene", 
                                   group.by = "sample", 
                                   proportion = FALSE, 
                                   cloneSize=c(Single=1, Small=5, Medium=10, Large=14, Hyperexpanded=20))

tcell$cloneSize <- dplyr::recode(tcell$cloneSize,
                                 "Single (0 < X <= 1)"        = "1",
                                 "Small (1 < X <= 5)"         = "2 - 5",
                                 "Medium (5 < X <= 10)"       = "6 - 10",
                                 "Large (10 < X <= 14)"       = "11 - 15",
                                 "Hyperexpanded (14 < X <= 20)" = "16 - 20",
                                 "None ( < X <= 0)"           = NA_character_
)

# Set factor order for proper stacking
tcell$cloneSize <- factor(tcell$cloneSize,
                          levels = c("1", "2 - 5", "6 - 10", "11-20","21 - 50"))


tcell <- subset(tcell, idents = c('CD4','CD8'))
dfa <- tcell@meta.data

dfa_tcr <- dfa %>%
  filter(!is.na(cloneSize))

plot_df <- dfa_tcr %>%
  group_by(Sample, annotations, cloneSize) %>%
  summarise(cell_count = n(), .groups = "drop")

plot_df_pct <- plot_df %>%
  group_by(Sample, annotations) %>%
  mutate(percent = cell_count / sum(cell_count) * 100) %>%
  ungroup()


# count plot
ggplot(plot_df, aes(x = Sample, y = cell_count, fill = cloneSize)) +
  geom_bar(stat = "identity", width = 0.7) +
  facet_grid(
    rows = vars(annotations),
    scales = "free_y",
    switch = "y"
  ) +
  scale_fill_manual(
    values = c(
      "1" = "#B388FF",
      "2 - 5" = "#2EC4C6",
      "6 - 10" = "#6DBE45",
      "11 - 15" = "#1F7A1F"
    ),
    name = "Cells per clonotype"
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_text(angle = 0, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank()
  ) +
  ylab("Cell count")

# percent plot
ggplot(plot_df_pct, aes(x = Sample, y = percent, fill = cloneSize)) +
  geom_bar(stat = "identity", width = 0.7) +
  facet_grid(
    rows = vars(annotations),
    scales = "free_y",
    switch = "y"
  ) +
  scale_fill_manual(
    values = c(
      "1" = "#B388FF",
      "2 - 5" = "#2EC4C6",
      "6 - 10" = "#6DBE45",
      "11 - 15" = "#1F7A1F"
    ),
    name = "Percentage of Cells per clonotype"
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_text(angle = 0, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank()
  ) +
  ylab("% of Cell count")

# Cell count per clonosize per cell type
clonalOccupy(
  tcell,
  x.axis = "sample",
  facet.by = "annotations", label = FALSE
) +
  scale_fill_manual(
    values = c("#B388FF", "#2EC4C6", "#6DBE45", "#1F7A1F")
  )


clonalOccupy(
  tcell,
  x.axis = "sample",
  facet.by = "annotations", label = FALSE, proportion = TRUE
) +
  scale_fill_manual(
    values = c("#B388FF", "#2EC4C6", "#6DBE45", "#1F7A1F")
  )

#### Plotting clonal expansion in CD8 and CD4 T cells
# CD8
cd8_expanded <- WhichCells(tcell, 
                           idents = "CD8", 
                           expression = cloneSize != "Single")

DimPlot(tcell, 
        group.by = "annotations",  # or "seurat_clusters" or "cloneSize"
        cells.highlight = list("CD8 Expanded" = cd8_expanded),
        cols.highlight = "red",    # color for CD8 expanded clonotypes
        cols = "grey",             # all other cells grey
        pt.size = 0.6,
        order = TRUE,
        sizes.highlight = 0.6) +            # draw highlighted cells on top
  labs(title = "Expanded CD8 Clones") + NoLegend()

# CD4
cd4_expanded <- WhichCells(tcell, 
                           idents = "CD4", 
                           expression = cloneSize != "Single")

DimPlot(tcell, 
        group.by = "annotations",  # or "seurat_clusters" or "cloneSize"
        cells.highlight = list("CD4 Expanded" = cd4_expanded),
        cols.highlight = "red",    # color for CD8 expanded clonotypes
        cols = "grey",             # all other cells grey
        pt.size = 0.6,
        order = TRUE,
        sizes.highlight = 0.6) +            # draw highlighted cells on top
  labs(title = "Expanded CD4 Clones") + NoLegend()

# Treg
treg_expanded <- WhichCells(tcell, 
                            idents = "Treg", 
                            expression = cloneSize != "Single")

DimPlot(tcell, 
        group.by = "annotations",  # or "seurat_clusters" or "cloneSize"
        cells.highlight = list("Treg Expanded" = treg_expanded),
        cols.highlight = "red",    # color for CD8 expanded clonotypes
        cols = "grey",             # all other cells grey
        pt.size = 0.6,
        order = TRUE,
        sizes.highlight = 0.6) +            # draw highlighted cells on top
  labs(title = "Expanded Treg Clones") + NoLegend()


#### clonal expansion tables
cd8_meta <- dfa_tcr %>%
  filter(annotations == "CD8")

cd8_summary <- cd8_meta %>%
  group_by(Sample) %>%
  summarise(
    total_cells = n(),
    single_cells = sum(cloneSize == "Single (0 < X <= 1)", na.rm = TRUE),
    expanded_cells = sum(cloneSize != "Single (0 < X <= 1)", na.rm = TRUE),
    expanded_pct = round(expanded_cells / total_cells * 100, 1),
    .groups = "drop"
  )

cd8_summary <- cd8_summary %>%
  mutate(
    Sample = factor(
      Sample,
      levels = c("mock_Cont", "MRV_Cont", "mock_R848", "MRV_R848")
    )
  ) %>%
  arrange(Sample)

#cd4
# dfa <- tcell@meta.data
# dfa_tcr <- dfa %% filter(!is.na(cloneSize))

cd4_meta <- dfa_tcr %>%
  filter(annotations == "CD4")

cd4_summary <- cd4_meta %>%
  group_by(Sample) %>%
  summarise(
    total_cells = n(),
    single_cells = sum(cloneSize == "Single (0 < X <= 1)", na.rm = TRUE),
    expanded_cells = sum(cloneSize != "Single (0 < X <= 1)", na.rm = TRUE),
    expanded_pct = round(expanded_cells / total_cells * 100, 1),
    .groups = "drop"
  )

cd4_summary <- cd4_summary %>%
  mutate(
    Sample = factor(
      Sample,
      levels = c("mock_Cont", "MRV_Cont", "mock_R848", "MRV_R848")
    )
  ) %>%
  arrange(Sample)


