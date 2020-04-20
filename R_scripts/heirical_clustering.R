

clusdata <- read.csv("clustering.csv", row.names = 1)

clusdata <- t(clusdata)
d <- dist(clusdata, method = "euclidean")
hc1 <- hclust(d, method = "complete" )
plot(hc1, cex = 0.6, hang = -1, las= .025, cex = .1 )
