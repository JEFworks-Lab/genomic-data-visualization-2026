# given a mystery dataset
data <- read.csv('~/Desktop/genomic-data-visualization-2026/data/GDV_practicequiz.csv.gz')
dim(data)
head(data)
data[1:5,1:10]
pos1 <- data[, 2:3] ## grab 2nd and 3rd column
head(pos1)
pos2 <- data[, c('x', 'y')] ## grab columns named x and y
head(pos2)
pos <- pos2

genes <- colnames(data)[4:ncol(data)]
gexp <- data[, genes]
gexp <- data[, 4:ncol(data)]
head(gexp)

rownames(gexp) <- rownames(pos) <- data[,1] ## keep track of cell names
dim(pos)
dim(gexp)

library(ggplot2)
df <- data.frame(pos)
ggplot(df) + geom_point(aes(x=x,y=y), size=0.1)

totgexp <- rowSums(gexp)
hist(totgexp)
df <- data.frame(pos, totgexp)
ggplot(df) + geom_point(aes(x=x,y=y,col=totgexp), size=0.5)

## normalize
#mat <- log10(gexp/totgexp*1e6+1) ## counts per million, log transform
mat <- gexp/totgexp*1e6 ## counts per million
mat[1:5,1:5]
rowSums(mat)

hist(mat[,1])
hist(mat[,10])

## looking at genes
i = 10
df <- data.frame(pos, gene = mat[,i])
ggplot(df) + geom_point(aes(x=x,y=y,col=gene), size=0.5)

## dimensionality via tSNE
seed(0) 
emb <- Rtsne::Rtsne(mat)
embedding <- emb$Y
colnames(embedding) <- c('tSNE1', 'tSNE2')
rownames(embedding) <- rownames(mat)

## tip: write functions to help you
runTSNE <- function(mat) {
  emb <- Rtsne::Rtsne(mat)
  embedding <- emb$Y
  colnames(embedding) <- c('tSNE1', 'tSNE2')
  rownames(embedding) <- rownames(mat)
  return(embedding)
}
embedding <- runTSNE(mat)

df <- data.frame(embedding)
ggplot(df) + geom_point(aes(x=tSNE1, y=tSNE2))

## kmeans
set.seed(10)
set.seed(1)
clusters <- kmeans(mat, centers=10)
com <- as.factor(clusters$cluster) ## set as categorical variable (as opposed to quantitative)

df <- data.frame(embedding, com)
ggplot(df) + geom_point(aes(x=tSNE1, y=tSNE2, col=com))

## differential expression analysis
cluster2cells <- names(com[com==1])
cluster3cells <- names(com[com==9])
results <- sapply(colnames(mat), function(g) {
  print(g)
  wilcox.test(gexp[cluster2cells, g], gexp[cluster3cells, g])$p.value
})
table(p.adjust(results) < 0.05)
siggenes <- names(which(p.adjust(results) < 0.05))

## average fold change
log2(colMeans(mat[cluster2cells, siggenes])/colMeans(mat[cluster3cells, siggenes]))
log2(colMeans(mat[cluster2cells, 1:10])/colMeans(mat[cluster3cells, 1:10]))



## PCA first? then tSNE
pcs <- prcomp(mat)
plot(pcs$sdev[1:1000]) ## scree plot
plot(pcs$x[,1], totgexp)
cor.test(pcs$x[,1], totgexp)
embedding2 <- runTSNE(pcs$x[,1:10])
head(embedding2)
df <- data.frame(embedding2, com, pcs$x[,1:10])
ggplot(df) + geom_point(aes(x=PC1, y=PC2, col=com))
ggplot(df) + geom_point(aes(x=tSNE1, y=tSNE2, col=com))

## visualize genes
df <- data.frame(embedding2, com, pcs$x[,1:10], gene=mat[,"Haus6"])
ggplot(df) + geom_point(aes(x=tSNE1, y=tSNE2, col=gene))


## visualize in tissue
df <- data.frame(pos, embedding2, com, pcs$x[,1:10], gene=mat[,"Haus6"])
ggplot(df) + geom_point(aes(x=tSNE1, y=tSNE2, col=com))

## merge clusters
comnew <- com
comnew[comnew == 2] <- 3
comnew[comnew %in% c(5,6)] <- 3
df <- data.frame(pos, embedding2, comnew, pcs$x[,1:10], gene=mat[,"Haus6"])
ggplot(df) + geom_point(aes(x=tSNE1, y=tSNE2, col=comnew))
ggplot(df) + geom_point(aes(x=x, y=y, col=comnew))

library(patchwork)
g1 <- ggplot(df) + geom_point(aes(x=tSNE1, y=tSNE2, col=comnew))
g2 <- ggplot(df) + geom_point(aes(x=x, y=y, col=comnew))
g1 + g2

g1 <- ggplot(df) + geom_point(aes(x=tSNE1, y=tSNE2, col=gene))
g2 <- ggplot(df) + geom_point(aes(x=x, y=y, col=gene))
g1 + g2


## redo kmeans
clusters <- kmeans(pcs$x[,1:10], centers=10)
com <- as.factor(clusters$cluster) ## set as categorical variable (as opposed to quantitative)
df <- data.frame(pos, embedding2, com, pcs$x[,1:10], gene=mat[,"Haus6"])
ggplot(df) + geom_point(aes(x=tSNE1, y=tSNE2, col=com))
ggplot(df) + geom_point(aes(x=x, y=y, col=com))
## if we didn't normalize, you may be misled into thinking there's no pattern

cluster2cells <- names(com[com==1])
cluster3cells <- names(com[com==9])
results <- sapply(colnames(mat), function(g) {
  print(g)
  wilcox.test(mat[cluster2cells, g], mat[cluster3cells, g])$p.value
})
table(p.adjust(results) < 0.05)
siggenes <- names(which(p.adjust(results) < 0.05))
log2(colMeans(mat[cluster2cells, siggenes])/colMeans(mat[cluster3cells, siggenes]))
log2(colMeans(mat[cluster2cells, 1:10])/colMeans(mat[cluster3cells, 1:10]))

df <- data.frame(pos, embedding2, com, pcs$x[,1:10], gene=mat[,"Ndufa8"])
ggplot(df) + geom_point(aes(x=tSNE1, y=tSNE2, col=gene))
ggplot(df) + geom_point(aes(x=x, y=y, col=gene))

## these are all different ways to check yourself
## you could always be misled
## double dipping problem 

## differences?
clusters <- kmeans(embedding2, centers=10)
clusters <- kmeans(pcs$x[,1:10], centers=10) ## preferred
clusters <- kmeans(mat, centers=10)

## feature select first (looking for variable genes)? then PCA? then tSNE
mat <- log10(gexp/totgexp*1e6+1) ## good idea if there are lots of genes
vg <- apply(mat, 2, var)
hist(vg)
table(vg > 0.001)
goodgenes <- names(which(vg > 0.001))
pcs <- prcomp(mat[, goodgenes]) ## good idea if there are lots of genes
plot(pcs$sdev[1:10]) ## scree plot
embedding3 <- runTSNE(pcs$x[,1:10])
head(embedding3)
df <- data.frame(embedding3, pcs$x[,1:10])
ggplot(df) + geom_point(aes(x=tSNE1, y=tSNE2)) 


