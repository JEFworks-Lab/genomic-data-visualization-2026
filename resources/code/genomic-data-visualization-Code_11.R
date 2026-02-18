# gzipped to save space
data <- read.csv('~/Desktop/genomic-data-visualization-2026/data/codex_spleen2.csv.gz')
data[1:8,1:8]

dim(data)
pos <- data[, 2:3]
head(pos)
pexp <- data[, 5:ncol(data)]
head(pexp)
rownames(pos) <- rownames(pexp) <- data[,1]
area <- data[, 4]
names(area) <- data[,1]
head(area)

## QC -> should I remove cells? should I normalize?
# look at the area
library(ggplot2)
df = data.frame(area)
head(df)
ggplot(df, aes(x = 1, y = area)) + geom_violin() 
# correspondence with total protein expression
df = data.frame(area, totpexp = rowSums(pexp))
ggplot(df, aes(x=log10(area), y=log10(totpexp))) + geom_point()
ggplot(df, aes(x = 1, y = log10(totpexp))) + geom_violin() 
# what does this look like in the tissue?
df = data.frame(area, totpexp = rowSums(pexp), pos)
ggplot(df, aes(x=x, y=y, col=log10(totpexp))) + geom_point(size=0.5)

## look at some of the proteins
dim(pexp)
sort(apply(pexp, 2, var), decreasing=TRUE)
sort(apply(pexp, 2, sum), decreasing=TRUE)
df = data.frame(pos, pexp)
ggplot(df, aes(y=1, x=log10(CD15))) + geom_violin()
ggplot(df, aes(x=x, y=y, col=log10(CD15))) + geom_point(size=0.5)

lapply(colnames(pexp)[1:2], function(p) {
  print(p)
  df = data.frame(pos, protein = pexp[, p])
  g1 <- ggplot(df, aes(y=1, x=log10(protein))) + geom_violin()
  g2 <- ggplot(df, aes(x=x, y=y, col=log10(protein))) + geom_point(size=0.5)
  print(g1 + g2)
})

p = "CD20"
df = data.frame(pos, protein = pexp[, p])
g1 <- ggplot(df, aes(y=1, x=log10(protein))) + geom_violin() + ggtitle(p)
g2 <- ggplot(df, aes(x=x, y=y, col=log10(protein))) + geom_point(size=0.5)
print(g1 + g2)

df = data.frame(tcell = pexp[, 'CD8'], bcell = pexp[, 'CD20'])
ggplot(df, aes(x=tcell, y=bcell)) + geom_point()

table(pexp[, 'CD8'] > 100 & pexp[, 'CD20'] > 100)


