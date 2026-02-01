---
layout: post
title:  "Homework 2"
author: Aarna Sanghai
jhed: asangha1
categories: [ HW2 ]
image: homework/hw2/HW2_asangha1.png
featured: false
---

### Question: How do tSNE coordinates change as you increase or decrease the perplexity?

### I attempted to design a visualization to enhance salience of how perplexity affects cluster structure and cell-type separability. The scatterplots use points as geometric primitives and position on the x axis (tSNE1) and y axis (tSNE2) as primary visual channels. According to the chart from Lesson 0, position is a highly perceptible encoding, allowing viewers to accurately judge distances and cluster boundaries. Color hue is used as a visual channel to encode categorical  data of cell type for rapid class identification.
The use of small multiples enables direct visual comparison across different perplexity settings. Because all panels share the same scales, alignment and consistency are preserved, reducing cognitive load and preventing misleading comparisons. This layout supports the Gestalt principle of similarity, since each panel has the same structure, and the principle of proximity, since related plots are placed close together. Within each panel as well, clusters emerge through the Gestalt principles of similarity and proximity: points with the same color and spatial closeness are perceived as belonging to the same group. 
I discovered that at moderate perplexities, these principles are strongest, producing visually salient and well-separated clusters. At extreme perplexities, these grouping cues weaken due to overlap and dispersion.
Specifically, I saw that: 
At low perplexity values (5), the algorithm emphasizes very small local neighborhoods. As a result, points that are nearby in high-dimensional space are mapped close together, but broader structure is lost. This leads to fragmented clusters and increased overlap between cell types, reducing perceptual separability. 
At moderate perplexity values (15-30), local and global structure are better balanced, producing more stable and coherent clusters. Cell types occupy more distinct regions and proximity more accurately reflects biological similarity. 
At high perplexity values (>50), the embedding prioritizes global relationships, which causes clusters to become elongated and compressed together. There is to reduced cluster separation.
Overall, this visualization effectively enhances saliency by combining strong visual channels, consistent geometric primitives, and Gestalt grouping principles.


### Code

```r

#Xenium-IRI-ShamR_matrix.csv.gz 

data <- read.csv("~/Desktop/Xenium-IRI-ShamR_matrix.csv.gz")


dim(data)
class(data)
head(data)
data[1:5, 1:5]
pos <- data[,c('x','y')]
rownames(pos) <- data[,1]
head(pos)
gexp <- data[, 4:ncol(data)]
rownames(gexp) <- data[,1]
gexp[1:5,1:5]

install.packages("patchwork")
install.packages("Rtsne")

# Load libraries
library(Rtsne)
library(ggplot2)
library(patchwork)

gexp <- data[, 4:ncol(data)]
rownames(gexp) <- data[,1]

#using help from CHATGPT:
#prompts:
# write a script that performs dimensionality reduction and visualization on a kidney
# gene expression matrix (gexp) using PCA and t-SNE to study the effect of
# perplexity on cluster separation.
#
# Step to follow:
# 1. Normalize gene expression using CPM-style normalization and log transform.
# 2. Run PCA on the normalized data and extract the top 10 principal components.
# 3. Run t-SNE on the top PCs using multiple perplexity values (5, 15, 30, 50, 75, 100).
# 4. Classify cells using marker genes:
#    - Aqp1     → Proximal Tubule
#    - Slc12a1  → Loop of Henle
#    - Aqp2     → Collecting Duct
#    - Havcr1   → Injured
#    - Other    → Remaining cells
# 5. Create six ggplot2 scatterplots (one per perplexity) showing tSNE1 vs tSNE2,
#    colored by cell type, with consistent axes and custom colors.
# 6. Arrange the six t-SNE plots using patchwork.
# 7. Plot this quality metric as a line chart with points, and highlight
#    perplexity = 30 as the "optimal" value using a dashed vertical line
#    and annotation.
# 8. Create a legend at the bottom for all cell types.
# 9. Combine all plots into one final figure

# Normalize and PCA
totgexp <- rowSums(gexp)
mat <- log10(gexp/totgexp*1e6+1)
pcs <- prcomp(mat, center = TRUE, scale = FALSE)
toppcs <- pcs$x[, 1:10]

# =============================================================================
# Run t-SNE with 6 perplexities
# =============================================================================

perplexities <- c(5, 15, 30, 50, 75, 100)
all_results <- list()

set.seed(42)
for (perp in perplexities) {
  cat("Running perplexity =", perp, "\n")
  tsne <- Rtsne(toppcs, dims=2, perplexity=perp)
  
  all_results[[as.character(perp)]] <- data.frame(
    tSNE1 = tsne$Y[, 1],
    tSNE2 = tsne$Y[, 2],
    perplexity = perp,
    Aqp1 = mat[, 'Aqp1'],
    Slc12a1 = mat[, 'Slc12a1'],
    Aqp2 = mat[, 'Aqp2'],
    Havcr1 = mat[, 'Havcr1']
  )
}

df_all <- do.call(rbind, all_results)

# Classify cells
df_all$cell_type <- "Other"
df_all$cell_type[df_all$Aqp1 > 0.5] <- "Proximal Tubule"
df_all$cell_type[df_all$Slc12a1 > 0.5] <- "Loop of Henle"
df_all$cell_type[df_all$Aqp2 > 0.5] <- "Collecting Duct"
df_all$cell_type[df_all$Havcr1 > 0.8] <- "Injured"

df_all$cell_type <- factor(df_all$cell_type,
                           levels = c("Proximal Tubule", "Loop of Henle",
                                      "Collecting Duct", "Injured", "Other"))

# Define colors once
cell_colors <- c(
  "Proximal Tubule" = "#FF7F0E",  # Orange
  "Loop of Henle" = "#2CA02C",    # Green
  "Collecting Duct" = "#1F77B4",  # Blue
  "Injured" = "#D62728",          # Red
  "Other" = "gray85"
)

# =============================================================================
# Create 6 t-SNE plots (NO individual legends)
# =============================================================================

plot_list <- list()
for (i in 1:length(perplexities)) {
  perp <- perplexities[i]
  df_sub <- df_all[df_all$perplexity == perp, ]
  
  plot_list[[i]] <- ggplot(df_sub, aes(x=tSNE1, y=tSNE2, color=cell_type)) +
    geom_point(size=0.3, alpha=0.6) +
    scale_color_manual(values = cell_colors, name = "Cell Type") +
    labs(title = paste("Perplexity =", perp),
         x = "tSNE1", y = "tSNE2") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      legend.position = "none",  # Remove individual legends
      panel.grid.minor = element_blank()
    )
}

# =============================================================================
# Create variance quality plot
# =============================================================================

perp_range <- c(5, 10, 15, 20, 30, 40, 50, 75, 100)
quality_scores <- data.frame(perplexity = perp_range)

set.seed(42)
quality_scores$variance <- sapply(perp_range, function(p) {
  tsne_temp <- Rtsne(toppcs, dims=2, perplexity=p, verbose=FALSE)
  var(tsne_temp$Y[,1]) + var(tsne_temp$Y[,2])
})

quality_plot <- ggplot(quality_scores, aes(x=perplexity, y=variance)) +
  geom_line(linewidth=1.5, color="steelblue") +
  geom_point(size=4, color="steelblue") +
  geom_vline(xintercept=30, linetype="dashed", color="red", linewidth=1) +
  annotate("text", x=35, y=max(quality_scores$variance)*0.95,
           label="← Optimal", color="red", size=5, fontface="bold", hjust=0) +
  labs(
    title = "Cluster Separation Quality",
    subtitle = "Higher variance = better separation",
    x = "Perplexity",
    y = "Total Variance"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank()
  )

# =============================================================================
# Combine with SINGLE LEGEND at the bottom
# =============================================================================

# Combine the 6 t-SNE plots
tsne_grid <- wrap_plots(plot_list, ncol = 3)

# Create a dummy plot JUST for the legend
legend_plot <- ggplot(df_all[df_all$perplexity == 5, ], 
                      aes(x=tSNE1, y=tSNE2, color=cell_type)) +
  geom_point() +
  scale_color_manual(values = cell_colors, name = "Cell Type") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold")
  )

# Extract just the legend
library(ggplot2)
get_legend <- function(p) {
  tmp <- ggplot_gtable(ggplot_build(p))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  tmp$grobs[[leg]]
}

shared_legend <- get_legend(legend_plot)

# Final layout: 6 plots + quality plot + shared legend
final_figure <- (tsne_grid / quality_plot / shared_legend) +
  plot_layout(heights = c(2, 0.8, 0.15)) +
  plot_annotation(
    title = "Effect of Perplexity on t-SNE Dimensionality Reduction",
    subtitle = "Kidney cell type separation across different perplexity values",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 13, hjust = 0.5)
    )
  )

print(final_figure)

```
