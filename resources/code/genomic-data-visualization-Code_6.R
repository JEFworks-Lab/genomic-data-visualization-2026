data <- read.csv('~/Desktop/genomic-data-visualization-2026/data/Xenium-IRI-ShamR_matrix.csv.gz')
dim(data)
## for students struggling with too many cells
data <- data[sample(1:nrow(data), 10000),]
dim(data)
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

# PCA
pcs <- prcomp(mat)
df <- data.frame(pcs$x, pos)
ggplot(df, aes(x=x, y=y, col=PC1)) + geom_point(cex=0.5)

# kmeans
clusters <- as.factor(kmeans(pcs$x[,1:2], centers=2)$cluster)
df <- data.frame(pcs$x, pos, clusters)
ggplot(df, aes(x=PC1, y=PC2, col=clusters)) + geom_point(cex=0.5) 
ggplot(df, aes(x=x, y=y, col=clusters)) + geom_point(cex=0.5)

# differential expression
?t.test
?wilcox.test
'Nphs2' %in% colnames(mat)
df <- data.frame(pcs$x, pos, clusters, gene=mat[,'Nphs2'])
ggplot(df, aes(x=PC5, y=PC6, col=gene)) + geom_point(cex=0.5, alpha=0.5)

# between cells of cluster 2 vs all others
clusterofinterest <- names(clusters)[clusters == 2]
othercells <- names(clusters)[clusters != 2]
x1 <- mat[clusterofinterest, 'Nphs2']
x2 <- mat[othercells, 'Nphs2']
wilcox.test(x1, x2, alternative='two.sided') # differentially expressed
wilcox.test(x1, x2, alternative='greater') # differentially upregulated in cluster of interest
wilcox.test(x2, x1, alternative='greater') # not differentially upregulated in other cells
wilcox.test(x2, x1, alternative='greater')$p.value # grab the p-value

# loops
?apply
?sapply
?lapply
out <- sapply(colnames(mat), function(gene) {
  x1 <- mat[clusterofinterest, gene]
  x2 <- mat[othercells, gene]
  wilcox.test(x1, x2, alternative='greater')$p.value
})
head(sort(out))
df <- data.frame(pcs$x, pos, clusters, gene=mat[,'Cd24a'])
ggplot(df, aes(x=PC1, y=PC2, col=gene)) + geom_point(cex=0.5, alpha=0.5)
ggplot(df, aes(x=x, y=y, col=gene)) + geom_point(cex=0.5, alpha=0.5)
head(sort(pcs$rotation[,1])) 

# do t test and wilcox test disagree?
out2 <- sapply(colnames(mat), function(gene) {
  x1 <- mat[clusterofinterest, gene]
  x2 <- mat[othercells, gene]
  t.test(x1, x2, alternative='greater')$p.value
})
df <- data.frame(wilcox = out, ttest=out2)
ggplot(df, aes(x=wilcox, y=ttest)) + geom_point()
df <- data.frame(wilcox = -log10(out), ttest=-log10(out2))
ggplot(df, aes(x=wilcox, y=ttest)) + geom_point()

which(out < 0.05 & out2 > 0.05) # signficant disagreement
out[c('Calm1','Calm2')]
out2[c('Calm1','Calm2')]
df <- data.frame(pcs$x, pos, clusters, gene=mat[,'Calm1'])
ggplot(df, aes(x=PC1, y=PC2, col=gene)) + geom_point(cex=0.5, alpha=0.5)
ggplot(df, aes(x=x, y=y, col=gene)) + geom_point(cex=0.5, alpha=0.5)
