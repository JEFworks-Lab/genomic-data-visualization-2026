data <- read.csv('~/Desktop/genomic-data-visualization-2026/data/Visium-IRI-ShamR_matrix.csv.gz')
pos <- data[,c('x', 'y')]
rownames(pos) <- data[,1]
gexp <- data[, 4:ncol(data)]
rownames(gexp) <- data[,1]
head(pos)
gexp[1:5,1:5]
dim(gexp)

# How many total genes are detected per cell/spot?
# How many unique genes are detected per cell/spot?
# What is the total expression of each gene across all cells/spots?
totgexp = rowSums(gexp)
head(totgexp)
library(ggplot2)
df <- data.frame(
  name = rownames(gexp),
  totgexp = totgexp
  )
head(df)
ggplot(df, aes(x = name, y = totgexp)) + geom_col()

## very skewed
ggplot(df, aes(x = totgexp)) + geom_histogram()
## after log transforming for biological interpretability
ggplot(df, aes(x = log10(totgexp))) + geom_histogram()

## Apply this to your data
# What is the average expression of each gene across all cells/spots?
# How many cells/spots is a particular gene detected in?

####### normalize
## counts per million with a pseudocount of 1
log10(0)
log10(1)
mat <- log10(gexp/totgexp * 1e6 + 1)
#rowSums(mat)

######## PCA
?prcomp
pcs <- prcomp(mat, center=TRUE, scale=FALSE)
names(pcs)
head(pcs$sdev)
pcs$rotation[1:5,1:5]

## What happens if you don't normalize? Don't log transform?
## What gene has the highest loading value on PC1?
## Visualize gene expression of this gene with respect to the PCs
library(ggplot2)
?sort
head(sort(pcs$rotation[,1], decreasing=TRUE))
head(sort(pcs$rotation[,1], decreasing=FALSE))
df <- data.frame(pos, pcs$x[,1:2], gene1 = mat[, 'Cda'], gene2 = mat[, 'Nccrp1'])
head(df)
ggplot(df, aes(x=gene1, y=gene2, col=PC1)) + geom_point()
ggplot(df, aes(x=PC1, y=PC2, col=gene)) + geom_point()
ggplot(df, aes(x=x, y=y, col=PC1)) + geom_point()

       