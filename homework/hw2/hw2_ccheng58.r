set.seed(123)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(viridis)
library(RColorBrewer)

#load data
data <- read.csv('~/Desktop/Github/genomic-data-visualization-2026/data/Xenium-IRI-ShamR_matrix.csv.gz')
id <- data[[1]]
pos <- data[, c("x","y")]
rownames(pos) <- id
gexp <- data[, 4:ncol(data), drop = FALSE]
gexp <- gexp[, sapply(gexp, is.numeric), drop = FALSE]
rownames(gexp) <- id

#normalize
lib <- rowSums(gexp)
mat <- log10(sweep(gexp, 1, lib, "/") * 1e6 + 1)
gene_mean <- colMeans(mat)
gene_var <- apply(mat, 2, var)

#PCA
pcs <- prcomp(mat, center = TRUE, scale. = FALSE)
loading_pc1 <- pcs$rotation[, 1]

dfg <- data.frame(
  gene = names(loading_pc1),
  loading_pc1 = as.numeric(loading_pc1),
  mean_expr = as.numeric(gene_mean[names(loading_pc1)]),
  var_expr = as.numeric(gene_var[names(loading_pc1)])
)

dfg$loading_strength <- cut(abs(dfg$loading_pc1), 
                            breaks = c(0, 0.05, 0.10, 0.15, Inf),
                            labels = c("Weak", "Moderate", "Strong", "Very Strong"))

#identify top genes for labeling
top_lab <- dfg[order(abs(dfg$loading_pc1), decreasing = TRUE), ][1:15, ]

#using a gradient from weak (blue) to strong (red) for loading strength
strength_colors <- c("Weak" = "#3182bd", 
                     "Moderate" = "#9ecae1",
                     "Strong" = "#fc9272", 
                     "Very Strong" = "#de2d26")

#PLOT 1: PC1 loadings vs mean expression
p1 <- ggplot(dfg, aes(x = mean_expr, y = loading_pc1)) +
  geom_point(aes(color = loading_strength), 
             alpha = 0.6, 
             size = 2) +
  #reference line
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  geom_text_repel(data = top_lab, 
                  aes(label = gene), 
                  size = 3.5, 
                  fontface = "bold",
                  max.overlaps = 50,
                  box.padding = 0.5,
                  point.padding = 0.3,
                  segment.color = "gray50",
                  segment.size = 0.3) +
  scale_color_manual(values = strength_colors, name = "Loading Strength") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 10),
    panel.grid.major = element_line(color = "gray95", linewidth = 0.3),
    axis.title = element_text(face = "bold", size = 11)
  ) +
  labs(
    title = "PC1 loadings vs mean expression",
    x = "Mean(log-normalized expression)",
    y = "PC1 loading"
  )

#PLOT 2: PC1 loadings vs variance
p2 <- ggplot(dfg, aes(x = var_expr, y = loading_pc1)) +
  #add points with color by loading strength
  geom_point(aes(color = loading_strength), 
             alpha = 0.6, 
             size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  #labels for top genes
  geom_text_repel(data = top_lab, 
                  aes(label = gene), 
                  size = 3.5, 
                  fontface = "bold",
                  max.overlaps = 50,
                  box.padding = 0.5,
                  point.padding = 0.3,
                  segment.color = "gray50",
                  segment.size = 0.3) +
  scale_color_manual(values = strength_colors, name = "Loading Strength") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 10),
    panel.grid.major = element_line(color = "gray95", linewidth = 0.3),
    axis.title = element_text(face = "bold", size = 11)
  ) +
  labs(
    title = "PC1 loadings vs variance",
    x = "Variance(log-normalized expression)",
    y = "PC1 loading"
  )

fig <- p1 | p2
fig

ggsave("hw2_pc1_loading_mean_var_improved.png", fig, width = 14, height = 6, dpi = 300)
       
       
#spatial visualization of PC1 scores

pc1_scores <- pcs$x[, 1]

spatial_df <- data.frame(
  x = pos$x,
  y = pos$y,
  PC1_score = pc1_scores
)

p_spatial <- ggplot(spatial_df, aes(x = x, y = y, color = PC1_score)) +
  geom_point(size = 0.8, alpha = 0.8) +
  scale_color_gradient2(
    low = "#2166ac",      #blue for negative scores
    mid = "#f7f7f7",      #white/light gray for zero
    high = "#b2182b",     #red for positive scores
    midpoint = 0,
    name = "PC1 Score"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    axis.title = element_text(face = "bold", size = 11),
    aspect.ratio = 1,  #keep square aspect ratio for spatial data
    panel.background = element_rect(fill = "gray98")
  ) +
  labs(
    title = "Spatial Distribution of PC1 Scores",
    x = "X-Coordinate",
    y = "Y-Coordinate"
  ) +
  coord_fixed()  #ensure equal scaling on both axes

p_spatial

ggsave("hw2_pc1_spatial_distribution.png", p_spatial, width = 8, height = 7, dpi = 300)
