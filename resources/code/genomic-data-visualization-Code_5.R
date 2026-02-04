# read in data
data <- read.csv('~/Desktop/genomic-data-visualization-2026/data/Visium-IRI-ShamR_matrix.csv.gz')
pos <- data[,c('x', 'y')]
rownames(pos) <- data[,1]
gexp <- data[, 4:ncol(data)]
rownames(gexp) <- data[,1]
head(pos)
gexp[1:5,1:5]
dim(gexp)

# normalize
totgexp = rowSums(gexp)
mat <- log10(gexp/totgexp * 1e6 + 1)
dim(mat)

library(ggplot2)
df <- data.frame(pos, totgexp)
ggplot(df, aes(x=x, y=y, col=log10(totgexp+1))) + geom_point() +
  coord_fixed()
?coord_fixed

# limit to highly variable genes 
?apply
?var
vg <- apply(mat, 2, var)
head(vg)
var(mat[, 'Xkr4']) ## equivalent
vargenes <- names(sort(vg, decreasing=TRUE)[1:300])
matsub <- mat[, vargenes]
dim(matsub) ## much smaller, hopefully won't crash

head(matsub)
ggplot(matsub, aes(x=Serpina1f, y=Ppp1r1b, col=Slc8a1)) + geom_point()

# PCA
pcs <- prcomp(matsub)
df <- data.frame(pcs$x, totgexp)
ggplot(df, aes(x=PC1, y=PC2, col=totgexp)) + geom_point()
#pcs$rotation
#pcs$x
# tSNE
tn <- Rtsne::Rtsne(pcs$x[, 1:10], dim=2) # what happens if we use less PCs?
emb <- tn$Y
df <- data.frame(emb, totgexp)
ggplot(df, aes(x=X1, y=X2, col=totgexp)) + geom_point()
# kmeans
?kmeans
km <- kmeans(matsub, centers=10) # what happens if we change k? on PCs?
names(km)
km$tot.withinss
cluster <- as.factor(km$cluster) # convert from numeric to categorical
df <- data.frame(emb, cluster, pcs$x, pos)
ggplot(df, aes(x=X1, y=X2, col=cluster)) + geom_point()
ggplot(df, aes(x=PC1, y=PC2, col=cluster)) + geom_point()
ggplot(df, aes(x=PC2, y=PC3, col=cluster)) + geom_point()
ggplot(df, aes(x=x, y=y, col=cluster)) + geom_point()

#save.image()
