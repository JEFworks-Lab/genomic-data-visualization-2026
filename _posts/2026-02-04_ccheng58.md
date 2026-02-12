---
author: Catherine Cheng
jhed: ccheng58
categories: [ HW2 ]
image: homework/hw2/hw2_PC1loading.png
featured: false
---

### Write a description explaining what you are trying to make salient and why you believe your data visualization is effective, using vocabulary terms from Lesson 1. (Question: How do the gene loadings on the first PC relate to features of the genes such as its mean or variance?)

Genes with the largest absolute PC1 loadings generally have moderate to high mean expression and high variance, suggesting that PC1 is mostly driven by genes that vary a lot across cells. In the plots, genes with strong loadings (such as Spink1, Hnf4a, and Cyp4b1) tend to appear toward the upper right of both panels, showing a positive relationship between PC1 loading magnitude and both mean expression and variance. In contrast, genes with low mean expression and low variance cluster near zero loading and contribute very little to PC1. Overall, this indicates that PC1 separates cells mainly based on genes that are highly variable rather than genes that are consistently low or uniform across the tissue.




My visualization is designed to highlight a few main relationships: how gene loading strength on PC1 relates to basic expression features like mean and variance, how cells are spatially organized based on PC1 scores, and which genes contribute most to PC1. I use a blue-to-red color scale to represent PC1 loading strength because color is processed quickly and makes strong positive and negative loadings easy to notice right away. Genes with similar loading values appear similar in color, which follows the Gestalt principle of similarity, while genes with very high loadings stand out more clearly against the background, helping them “pop out” from the rest. Position is used to show quantitative relationships in the scatterplots since it allows for more accurate comparisons, while color helps with grouping and pattern recognition. Instead of labeling every gene, I only label the top 15 genes with the strongest loadings to reduce clutter and keep the visualization readable while still drawing attention to the most important genes.






```r
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






```


### AI PROMPTS:
given the specifications in the hw description below and the example student code as well as lecture code, help me write an R script that aligns with both




i have pasted my code and what my visualization looks like right now in R (black and white) and also a students who i think has a good visualization. right now, i think i can color code the genes instead of putting labels. i also want it to increase salience and look colorful while using gestault principles and visual channels to look better and communicate code more effectively
