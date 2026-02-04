# Read in data
data <- read.csv('~/Desktop/genomic-data-visualization-2026/data/Xenium-IRI-ShamR_matrix.csv.gz')
data[1:5,1:5]
pos <- data[,c('x', 'y')]
rownames(pos) <- data[,1]
gexp <- data[, 4:ncol(data)]
rownames(gexp) <- data[,1]
head(pos)
gexp[1:5,1:5]
dim(pos)
dim(gexp)

# Normalize (what happens if you don't? -> how does that affect your PCA for example?)
totgexp <- rowSums(gexp)
head(totgexp)
head(sort(totgexp, decreasing=TRUE))
mat <- log10(gexp/totgexp * 1e6 + 1)
dim(mat)

# PCA
## what happens if you scale? (you may need to remove certain genes)
pcs <- prcomp(mat, center=TRUE, scale=FALSE) 
names(pcs)
head(pcs$sdev)
length(pcs$sdev)
plot(pcs$sdev[1:50])

# tSNE
toppcs <- pcs$x[, 1:10] 
## what happens if we use more PCs?
tsne <- Rtsne::Rtsne(toppcs, dims=3, perplexity = 30) 
## what happens if we change perplexity?

## What is the relationship between genes and tSNE coordinates? PCs? Loadings of genes?
## What can we read off from this tSNE plot about relationships between cells? Between genes?
emb <- tsne$Y
rownames(emb) <- rownames(mat)
colnames(emb) <- c('tSNE1', 'tSNE2')
head(emb)

library(ggplot2)
ratio = log2((mat[, 'Cd24a']+1)/(mat[, 'Spink1']+1))
names(ratio) <- rownames(mat)
head(ratio)
hist(ratio)
#table(is.na(ratio))
#goodcells <- names(ratio)[!is.na(ratio)]
#ratio[goodcells]
df <- data.frame(pos, emb, pcs$x[,1:10], ratio=ratio)
ggplot(df, aes(x= tSNE1, y=tSNE2, col=PC1)) + geom_point(size=0.01, alpha=0.2)
head(sort(pcs$rotation[,1], decreasing=TRUE))
ggplot(df, aes(x= tSNE1, y=tSNE2, col=ratio)) + geom_point(size=0.01, alpha=0.2)
ggplot(df, aes(x= PC1, y=ratio, col=tSNE1)) + geom_point(size=0.01, alpha=0.2)

ggplot(df, aes(x= PC1, y=PC2, col=gene)) + geom_point(size=0.01, alpha=0.2)
ggplot(df, aes(x= x, y=y, col=tSNE1)) + geom_point(size=0.01, alpha=0.2)

ggplot(df, aes(x = x, y = y, col = tSNE1, size = tSNE2)) + 
  geom_point(alpha = 0.2) +
  scale_color_gradient2(low = "red", mid = "white", high = "blue") +
  scale_size(range = c(0.01, 0.5)) +  # Adjust size range as needed
  theme_minimal()  # Optional: to use a cleaner theme

