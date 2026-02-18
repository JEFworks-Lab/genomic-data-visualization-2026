---
author: Alex Kim
jhed: ykim276
categories: [ HW2 ]
image: homework/hw2/hw2_ykim276.png
featured: false
---


### Write a description explaining what you are trying to make salient and why you believe your data visualization is effective, using vocabulary terms from Lesson 1. (How do the genes with high versus low loadings relate to each other? How are they patterned relative to each other in the tissue?)


I analyzed how genes with high vs. low PC1 loadings relate to each other and spatial tissue patterns. The visualization includes quantitative data (PC1 loadings, gene expression, correlations) and spatial data (x–y tissue coordinates). It shows both individual gene properties and pairwise gene relationships.


Panel 1 uses bar area and color hue (red/blue/gray) to show PC1 loading distributions and highlight selected genes. Panel 2 uses colored tiles (red–white–blue) to encode correlation strength. Panels 3–4 use point position for spatial coordinates and color intensity to show average gene expression.


The visualization shows that genes with opposite PC1 loadings are anti-correlated and spatially separated. Panel 1 highlights a few extreme-loading genes driving major variation. Panel 2 reveals strong negative correlations between high-positive and high-negative genes. Panels 3–4 show these groups occupy distinct tissue regions, meaning PC1 reflects spatial zonation.


Similarity groups genes by color (red vs blue). Proximity forms clear correlation blocks in the heatmap. Continuity appears in smooth spatial color gradients. Position encodes space (most accurate), while color shows expression and correlation. The consistent red = positive, blue = negative scheme keeps the biological meaning clear across panels.










```r


library(dplyr)
library(ggplot2)
library(patchwork)


cat("Loading data...\n")
data <- read.csv("~/genomic-data-visualization-2026/data/Xenium-IRI-ShamR_matrix.csv.gz", row.names = 1)
pos <- data %>% select(x, y)
gexp <- data %>% select(-x, -y)


#Normalize
totals <- rowSums(gexp)
mat <- log10((gexp / totals * 10000) + 1)


cat("Running PCA...\n")
pca <- prcomp(mat, center = TRUE, scale = FALSE)
pc1_loadings <- pca$rotation[, 1]


#Top 3 positive and top 3 negative loading genes
top_pos <- names(sort(pc1_loadings, decreasing = TRUE)[1:3])
top_neg <- names(sort(pc1_loadings, decreasing = FALSE)[1:3])
selected_genes <- c(top_pos, top_neg)


cat("\nTop positive loadings:", top_pos, "\n")
cat("Top negative loadings:", top_neg, "\n")


#PANEL 1: LOADING DISTRIBUTION
loading_df <- data.frame(
 gene = names(pc1_loadings),
 loading = pc1_loadings,
 category = ifelse(names(pc1_loadings) %in% top_pos, "High Positive",
                   ifelse(names(pc1_loadings) %in% top_neg, "High Negative", "Other"))
)


p1 <- ggplot(loading_df, aes(x = loading, fill = category)) +
 geom_histogram(bins = 40, alpha = 0.7) +
 scale_fill_manual(
   values = c("High Positive" = "#D7191C", "High Negative" = "#2C7BB6", "Other" = "gray80"),
   name = "Gene Category"
 ) +
 theme_minimal() +
 labs(
   title = "PC1 Loading Distribution",
   x = "PC1 Loading Value",
   y = "Number of Genes"
 ) +
 theme(
   plot.title = element_text(size = 12, face = "bold"),
   legend.position = "top"
 )


#PANEL 2: GENE-GENE CORRELATION
selected_expr <- mat[, selected_genes]
gene_cor <- cor(selected_expr)


cor_long <- reshape2::melt(gene_cor)
colnames(cor_long) <- c("Gene1", "Gene2", "Correlation")


cor_long$Gene1_cat <- ifelse(cor_long$Gene1 %in% top_pos, "Positive", "Negative")
cor_long$Gene2_cat <- ifelse(cor_long$Gene2 %in% top_pos, "Positive", "Negative")


p2 <- ggplot(cor_long, aes(x = Gene1, y = Gene2, fill = Correlation)) +
 geom_tile(color = "white") +
 scale_fill_gradient2(
   low = "#2C7BB6", mid = "white", high = "#D7191C",
   midpoint = 0, limits = c(-1, 1)
 ) +
 theme_minimal() +
 theme(
   axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
   axis.text.y = element_text(size = 8),
   plot.title = element_text(size = 12, face = "bold"),
   legend.position = "right"
 ) +
 labs(
   title = "Correlation Between High/Low Loading Genes",
   x = "", y = ""
 ) +
 coord_fixed()


#PANELS 3-4: SPATIAL PATTERNS
plot_data <- data.frame(
 x = pos$x,
 y = pos$y,
 pos_score = rowMeans(mat[, top_pos]),
 neg_score = rowMeans(mat[, top_neg])
)


p3 <- ggplot(plot_data, aes(x = x, y = y, color = pos_score)) +
 geom_point(size = 0.3, alpha = 0.6) +
 scale_color_gradient(low = "gray90", high = "#D7191C", name = "Expression") +
 theme_minimal() +
 labs(
   title = "High Positive Loading Genes",
   subtitle = paste0("Average: ", paste(top_pos, collapse = ", ")),
   x = "X (μm)", y = "Y (μm)"
 ) +
 theme(
   plot.title = element_text(size = 10, face = "bold"),
   plot.subtitle = element_text(size = 7),
   legend.position = "right"
 ) +
 coord_fixed()


p4 <- ggplot(plot_data, aes(x = x, y = y, color = neg_score)) +
 geom_point(size = 0.3, alpha = 0.6) +
 scale_color_gradient(low = "gray90", high = "#2C7BB6", name = "Expression") +
 theme_minimal() +
 labs(
   title = "High Negative Loading Genes",
   subtitle = paste0("Average: ", paste(top_neg, collapse = ", ")),
   x = "X (μm)", y = "Y (μm)"
 ) +
 theme(
   plot.title = element_text(size = 10, face = "bold"),
   plot.subtitle = element_text(size = 7),
   legend.position = "right"
 ) +
 coord_fixed()


#COMBINE PANELS
final_plot <- (p1 | p2) / (p3 | p4) +
 plot_annotation(
   title = "Relationship Between High and Low PC1 Loading Genes",
   theme = theme(plot.title = element_text(size = 14, face = "bold"))
 )


#SAVE
ggsave("hw2_tzhan104.png", final_plot, width = 14, height = 10, dpi = 300, bg = "white")
cat("\nSaved as hw2_tzhan104.png\n")
print(final_plot)


```


### AI Prompts
Can you make me a data visualization for my homework with like 2+ panels that does something with PCA or t-SNE on this spatial gene data? I need to pick one of these questions about how genes relate to each other or how changing parameters affects the results, and the code needs to be simple and not copy the examples I showed you.
