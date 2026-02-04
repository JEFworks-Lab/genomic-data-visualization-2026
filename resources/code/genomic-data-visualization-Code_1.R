# Install R
# Install R-studio
# Feel free to use VS-Code, etc
?read.csv
data <- read.csv('~/Desktop/genomic-data-visualization-2026/data/Visium-IRI-ShamR_matrix.csv.gz')
dim(data) 
## ~80k cells, ~300 genes for imaging
## ~1000 spots, ~20k genes for sequencing

class(data)
head(data)
data[1:5,1:5] ## first 5 rows, first 5 cols
pos <- data[,c('x', 'y')]
rownames(pos) <- data[,1]
head(pos)
gexp <- data[, 4:ncol(data)]
rownames(gexp) <- data[,1]
gexp[1:5,1:5]
dim(gexp)

#install.packages('ggplot2')
library(ggplot2)
'Mki67' %in% colnames(gexp) ## Protein Ki67 
'Nphs2' %in% colnames(gexp) ## podocyte
df <- data.frame(pos, gene=gexp[,'Cyp4b1'])
head(df)
ggplot(df) + geom_point(aes(x = x, y = y, col=log10(gene+1)), size=2)

## help from Claude
## prompt: "Given a data frame `df` where the first two columns are spatial positional information (x, y) and the third column is quantitative data of a gene's expression (gene) across a population of cells, create a data visualization using ggplot2 where each cell is represented as a point where the gene expression magnitude is the x axis, the spatial y position is the y axis, and the x-position is the color hue"
ggplot(df, aes(x = gene, y = y, color = x)) +
  geom_point() +
  scale_color_viridis_c() +  
  labs(
    x = "Gene Expression",
    y = "Spatial Y Position",
    color = "X Position"
  ) +
  theme_minimal()


