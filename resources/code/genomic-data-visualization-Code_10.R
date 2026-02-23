# read in Visium data
data <- read.csv('~/Desktop/genomic-data-visualization-2026/data/Visium-IRI-ShamR_matrix.csv.gz')
pos <- data[,c('x', 'y')]
rownames(pos) <- data[,1]
gexp <- data[, 4:ncol(data)]
rownames(gexp) <- data[,1]
head(pos)
gexp[1:5,1:5]
dim(gexp)

# read in xenium data to get gene list
genes <- read.csv('~/Desktop/genomic-data-visualization-2026/data/Xenium-IRI-ShamR_matrix.csv.gz')
genes <- colnames(genes)[4:ncol(genes)]
genes.have <- intersect(genes, colnames(gexp))
length(genes.have)

# subset my Visium to just Xenium genes
gexp.sub <- gexp[, genes.have]
range(rowSums(gexp.sub))

# apply deconvolution to Visium data
# compare with Xenium data
library(STdeconvolve)

# deconvolution
ldas <- fitLDA(gexp.sub, Ks = c(3,4,5,6,7))
optLDA <- optimalModel(models = ldas, opt = "3")
results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
# theta is cell-type proportions
# beta is cell-type-specific gene expression
head(results$theta)
head(results$beta)
head(sort(results$beta[1,], decreasing=TRUE))

deconProp <- results$theta
deconGexp <- results$beta
g1 <- vizAllTopics(deconProp, pos, r = 40)
g1

## what are the genes associated with cell-type 7?
## if you visualize those genes in the tissue?
## HW EC2 - compare with clustering both for Visium and Xenium
## what happens if you pick a different K?

pcs <- prcomp(gexp.sub)
df <- data.frame(pcs$x[,1:2], deconProp)
emb <- Rtsne::Rtsne(gexp.sub)$Y
colnames(emb) <- c('tSNE1', 'tSNE2')
df <- data.frame(emb, deconProp)
head(df)
library(ggplot2)
g2 <- ggplot(df, aes(x=tSNE1, y=tSNE2, col=X3)) + geom_point()
g2

head(sort(deconGexp[3,], decreasing=TRUE))
g <- 'Umod'
deconGexp[, g]
df <- data.frame(emb, pos, gene=log10(gexp.sub[rownames(pos),g]+1))
g3 <- ggplot(df, aes(x=x, y=y, col=gene)) + geom_point()
g3

g4 <- ggplot(df, aes(x=tSNE1, y=tSNE2, col=gene)) + geom_point()
g4

library(patchwork)
(g1 + g2) / (g3 + g4)


## cell density impacts clustering but not deconvolution
df <- data.frame(emb, pos, totgexp = log10(rowSums(gexp)[rownames(pos)]))
head(df)
foo <- ggplot(df, aes(x=x, y=y, col=totgexp)) + geom_point()
bar <- ggplot(df, aes(x=tSNE1, y=tSNE2, col=totgexp)) + geom_point()
foo + bar



