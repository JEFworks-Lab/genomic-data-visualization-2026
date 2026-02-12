# read in data
data <- read.csv('~/Desktop/genomic-data-visualization-2026/data/Visium-IRI-ShamR_matrix.csv.gz')
data[1:5,1:5]
pos <- data[,c('x', 'y')]
rownames(pos) <- data[,1]
gexp <- data[, 4:ncol(data)]
rownames(gexp) <- data[,1]
head(pos)
gexp[1:5,1:5]
dim(gexp)

library(MERINGUE)
dim(pos)
W <- getSpatialNeighbors(pos, filterDist = 100)
dim(W)
sum(W)
plotNetwork(pos, W, cex=0.5)

?MERINGUE::moranTest

x <- gexp[, 'Aqp2'] ## just use gene expression counts for now
names(x) <- rownames(gexp)
x
## maybe you would want to normalize?
moranTest(x, W)
moranPermutationTest(x, W) ## compare with permutation testing (very slow)

## use a loop to test more genes
## for i in 1 to 1000, do something with i
results <- do.call(rbind, lapply(1:1000, function(i) {
  print(i)
  x <- gexp[, i] 
  names(x) <- rownames(gexp)
  moranTest(x, W)
  #moranTest(x, W, alternative = "less")
}))
rownames(results) <- colnames(gexp)[1:1000]
head(results)

df <- data.frame(pos, gene=gexp[, 'Gsta3'])
library(ggplot2)
ggplot(df, aes(x=x, y=y, col=gene)) + geom_point(cex=2) + theme_void() +
  scale_color_gradient2(low="blue", mid="grey", high="red", midpoint=mean(df$gene))

head(results[order(results[, 'p.value']),])
head(results[order(results[, 'observed'], decreasing=TRUE),])
tail(results[order(results[, 'observed'], decreasing=TRUE),])


