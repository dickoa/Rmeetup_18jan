tempsen <- read.csv("tempsen.csv", row.names = 1)
mat_dist <- dist(tempsen)

tempclust <- hclust(mat_dist, method = "average")
tempclust
plot(tempclust, hang = -1)

as.matrix(mat_dist)[1:5, 1:5]

set.seed(1)
tempkm <- kmeans(tempsen, centers = 2, nstart = 100)

library(FactoMineR)
temp_pca <- PCA(tempsen, graph = FALSE)
temp_hcpc <- HCPC(temp_pca, graph = FALSE)
plot.HCPC
plot(temp_hcpc, choice = "map", title = "", xaxt = "n", yaxt = "n")
axis(1, at = c(-10, -5, 0, 5), c(-10, -5, 0, 5), lwd = 0.2)
axis(2, at = c(-10, -5, 0, 5), c(-10, -5, 0, 5), lwd = 0.2)

par()







