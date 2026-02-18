---
layout: post
title:  "Alex Kim_HW1"
author: Alex Kim
jhed: ykim276
categories: [ HW1 ]
image: homework/hw1/hw1_ykim276.png
featured: false
---


### 1. What data types are you visualizing?
All variables I visualized are quantitative. x and y are continuous coordinate variables, and gene expression is a continuous numeric variable. A log transform is used to better show differences between the expression.


### 2. What data encodings (geometric primitives and visual channels) are you using to visualize these data types?
Points are the geometric primitive. X and Y Position encodes spatial location, and color level encodes gene expression level.


### 3. What about the data are you trying to make salient through this data visualization?
This visualization I created makes the spatial pattern of Trpv5 gene expression salient. Bright colors highlight where the gene is highly expressed, while darker areas show low expression and provide tissue context.


### 4. What Gestalt principles or knowledge about perceptiveness of visual encodings are you using to accomplish this?
The plot uses figure–ground separation, where low-expression parts form the background and bright parts stand out. Color contrast helps high expression regions more noticable, and spatial proximity shows the clusters.


### 5. Code (paste your code in between the ``` symbols)


data <- read.csv("~/genomic-data-visualization-2026/data/Xenium-IRI-ShamR_matrix.csv.gz")


pos <- data[,c('x', 'y')]
rownames(pos) <- data[,1]
> head(pos)
gexp <- data[,4:ncol(data)]
> rownames(gexp) <- data[,1]
library(ggplot2)
ggplot(df) + geom_point(aes(x = x, y= y, col = gene), size = 0.001)
> ggplot(df) + geom_point(aes(x = x, y= y, col = log1p(gene)), size = 0.001)
> ggplot(df) + geom_point(aes(x = x, y= y, col = log1p(gene)), size = 0.001) + scale_color_viridis_c(name = "Trpv5\nexpression")
> ggplot(df) + geom_point(aes(x = x, y= y, col = log1p(gene)), size = 0.001) + scale_color_viridis_c(name = "Trpv5\nexpression") + coord_equal() +
+     theme_void()