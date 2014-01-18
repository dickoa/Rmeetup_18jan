### Package nécessaire pour l'analyse
library(Cairo)
library(MASS)
library(ggplot2)
library(FactoMineR)
library(reshape2)
library(ade4)
source("panel.lm.R")


### R comme une calculatrice
x <- c(pi, 4*pi / 3)

log(-cos(x))

### R collection (1/2)
c(1, 2, 10)

c("Ali", "Modou", "Marie")

### R collection (1/2)
c(1, 2, "a")

list("Ali", 10, "Marie")

### R table de donnees
df <- data.frame(nom = c("Ali", "Modou", "Marie"),
                 taille = c(170, 185, 165))

df$nom

df[df$taille > 165, ]


#### ACP
## lire base de donnees
auto <- read.csv("data/auto.csv", row.names = 1)
summary(auto)

## regarder les 6 premieres lignes
head(auto)

### stars plot
stars(auto, len = 0.65,
      key.loc = c(4.5, 12),
      main = "", draw.segments = FALSE)


### Pairs plot
pairs(auto, col = c(rep("black", 18), "red", rep("black", 5)),
      pch = c(rep(3, 18), 15, rep(3, 5)),
      gap=0, panel = panel.lm)


### ACP
auto_acp <- PCA(auto, ncp = 2, graph = FALSE)

### Scree plot
eig_val <- auto_acp$eig[,1]
barplot(rev(eig_val),
        col = "steelblue", horiz = TRUE,
        names.arg = paste0("PC", rev(seq_along(eig_val))),
        las = 1, xaxt = "n")

### Eigenvalue + cumulative
auto_acp$eig


####
tab1 <- cbind(
    dist = auto_acp$ind$dist^2,
    coord =  auto_acp$ind$coord,
    cos2 = auto_acp$ind$cos2
    )

colnames(tab1) <- c("dist", paste0("coord_dim", 1:2), paste0("cos2_dim", 1:2))
head(tab1)

#### coord des variables
auto_acp$var$coord
auto_acp$var$cor


### Graphique
plot(auto_acp, cex = 0.8, title = "")
graph.var(auto_acp, title = "", new.plot = FALSE)

pca <- prcomp(auto, scale. = TRUE)
biplot(pca, las = 1, cex = 0.8)
abline(v = 0, h = 0, lty = "dashed")

### AFC
nice <- read.csv("data/nice.csv")
names(nice) <- tolower(names(nice))
str(nice)

### reshape to long format
nice <- melt(nice, id = "csp", variable.name = "filiere")
str(nice)

### turn the data to contigency table
cont_table <- xtabs(value ~ csp + filiere, data = nice)
cont_table

### test de khi2
chisq.test(cont_table)

### AFC
nice_afc <- CA(cont_table, graph = FALSE)

plot(nice_afc, title = "", cex = 0.8)



### Contrib
rbind(
    nice_afc$row$contrib,
    nice_afc$col$contrib
    )

### Cos2
rbind(
    nice_afc$row$cos2,
    nice_afc$col$cos2
    )


### ACM
chien <- read.csv("data/chien.csv", row.names = 1, colClasses = "factor")
names(chien) <- tolower(names(chien))
str(chien)

acm.disjonctif(chien)

chien_acm <- MCA(chien,
                 quali.sup = match("fonction", names(chien)),
                 graph = FALSE)

barplot(rev(chien_acm$eig[,1]),
        col = "steelblue", horiz = TRUE,
        names.arg = paste0("axes ", rev(seq_along(chien_acm$eig[,1]))),
        las = 1, xaxt = "n")

plot(chien_acm, title = "", cex = 0.8)


### CAH
tempsen <- read.csv("data/tempsen.csv", row.names = 1)
tempsen

mat_dist <- dist(tempsen)
as.matrix(mat_dist)[1:5, 1:5]

### hclust pour cah
tempclust <- hclust(mat_dist, method = "average")

plot(tempclust, hang = -1, main = "", cex = 0.8, las = 1)


barplot(sort(tempclust$height, decreasing = TRUE),
        col = c("steelblue", "red")[(sort(tempclust$height, decreasing = TRUE) > 10) + 1], las = 1)
abline(h = 10, lty = "dashed", col = "red")

### Decoupage de l'arbre
temphcut <- cutree(tempclust, k = 3)

plot(tempclust, hang = -1, main = "")
abline(h = 10, lty = "dashed", col = "red")

### Avec HCPC
temp_pca <- PCA(tempsen, graph = FALSE)
temp_hcpc <- HCPC(temp_pca, graph = FALSE)

plot(temp_hcpc, choice = "map", title = "", xaxt = "n", yaxt = "n")


#### Example classification de ESPS 2005
household <- readRDS("data/household.rds")

### Nettoyage des donnees
is.na(household$cons_tot) <- household$cons_tot == 0
household$cons_tot[is.na(household$cons_tot)] <- median(household$cons_tot, na.rm = TRUE)
household$cons_tot <- log(household$cons_tot)
is.na(household$deptotjr) <- household$deptotjr < 50
household$deptotjr[is.na(household$deptotjr)] <- median(household$deptotjr, na.rm = TRUE)
household$deptotjr <- log(household$deptotjr)


household$hh_size <- cut(household$hh_size,
                         breaks = quantile(household$hh_size,
                         probs = seq(0, 1, 1/3)),
                         right = FALSE,
                         ordered_result = TRUE,
                         include.lowest = TRUE,
                         labels = c("low", "medium", "high"))

household$cons_tot <- cut(household$cons_tot,
                          breaks = quantile(household$cons_tot,
                          probs = seq(0, 1, 1/3)),
                          right = FALSE,
                          ordered_result = TRUE,
                          include.lowest = TRUE,
                          labels = c("low", "medium", "high"))

household$deptotjr <- cut(household$deptotjr,
                          breaks = quantile(household$deptotjr,
                          probs = seq(0, 1, 1/3)),
                          right = FALSE,
                          ordered_result = TRUE,
                          include.lowest = TRUE,
                          labels = c("low", "medium", "high"))

summary(household)

naxes_max  <- sum(sapply(household[ ,-1], nlevels)) - ncol(household[ ,-1])
mca <- dudi.acm(household[,-1], nf = naxes_max, scannf = FALSE)

poids <- list()
for(i in 2:ncol(household)) {
    poids[[i]] <- table(household[,i])
}


poids <- cbind(unlist(poids))
str(mca)
variable  <- mca$co[ ,1:3]
variable$poids <- poids
names(variable) <- c("dim.1", "dim.2", "dim.3", "poids")
str(variable)

ggplot() +
    geom_point(aes(x = dim.1, y = dim.2, size = poids), pch = 21,
               colour = "steelblue", fill = "white", data = variable) +
    geom_text(aes(x = dim.1 - 0.05, y = dim.2 + 0.07, label = rownames(variable)), size = 1.5, data = variable) +
    scale_x_continuous(limits = c(-1.6, 1.2)) +
    geom_vline(aes(xintercept = 0), linetype = "dotted") +
    geom_hline(aes(yintercept = 0), linetype = "dotted")+
    labs(x = "axe 1", y = "axe 2") +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.ticks = element_blank())


ggplot(data.frame(axes = seq_along(mca$eig), eig = mca$eig)) +
    geom_bar(stat = "identity", size = 0.2,
             aes(x = axes, y = eig, fill = eig > 0.1)) +
    scale_fill_manual(values = c("grey", "steelblue")) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.ticks = element_blank())


mcadata <- mca$l1[ ,1:4]
sapply(mcadata, function(x) list(mean = round(mean(x), 2), std = sd(x)))


set.seed(1234)
hh_part <- kmeans(mcadata, centers = 100,
                  nstart = 20, iter.max = 500)

hh_hclust <- hclust(dist(hh_part$centers), method = "ward")


plot(hh_hclust, hang = -1)
abline(h = 20, lty = "dashed", col = "red")


indice <- hh_hclust$height
barplot(sort(indice, decreasing = TRUE),
        col = c("steelblue", "red")[(sort(indice, decreasing = TRUE) > 20) + 1])
abline(h = 20, lty = "dashed", col = "red")

### Decoupage
hh_tree <- cutree(hh_hclust, k = 3)
table(hh_tree)

###
dataclust <- data.frame(kclust = seq_along(hh_tree), hclust = hh_tree)
household$kclust <- hh_part$cluster
household <- merge(household, dataclust, by = "kclust")

### Groupe
table(household$hclust)

centers <- by(mcadata, household$hclust, function(x) apply(x, 2, median))
centers <- matrix(unlist(centers), ncol = ncol(mcadata), byrow = TRUE)
consol <- kmeans(mcadata, centers = centers, iter.max = 50, nstart = 10)

indiv  <- mca$li[ ,1:3]
names(indiv) <- c("dim.1", "dim.2", "dim.3")
indiv$hclust <- factor(household$hclust)
indiv$hclust2 <- factor(consol$cluster)

names(indiv)[5] <- "cluster"
ggplot() +
    geom_point(aes(x = dim.1, y = dim.2, colour = cluster), pch = 19, data = indiv, size = 0.8) +
    geom_vline(aes(xintercept = 0), linetype = "dotted", size = 0.5) +
    geom_hline(aes(yintercept = 0), linetype = "dotted", size = 0.5) +
    labs(x = "axe 1", y = "axe 2") +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.ticks = element_blank())
