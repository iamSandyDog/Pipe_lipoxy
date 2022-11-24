library(tidyverse)
library(RColorBrewer)
library(magrittr)
library(ggfortify)

f.table <- read.table("https://raw.githubusercontent.com/iamSandyDog/Pipe_lipoxy/main/.github/workflows/full_biochem.txt", header = TRUE)

f2 <- na.omit(f.table)
#Доделать
fun1 <- function(x) {
  c.mean <- f2 %>%
    filter(gr == x)
  c <- mean(c.mean$X1lipoxy.pr)
}
f.table$X1lipoxy.pr[9] <- fun1('c')
f.table$X1lipoxy.pr[11] <- fun1('a')
f.table$X1lipoxy.pr[10] <- fun1('cg')
f.table$X1lipoxy.pr[12] <- fun1('ag')

fun2 <- function(x) {
  c2.mean <- f2 %>%
    filter(gr == x)
  c2 <- mean(c2.mean$X2lipoxy.pr)
}
f.table$X2lipoxy.pr[9] <- fun2('c')
f.table$X2lipoxy.pr[11] <- fun2('a')
f.table$X2lipoxy.pr[10] <- fun2('cg')
f.table$X2lipoxy.pr[12] <- fun2('ag')

fun3 <- function(x) {
  c3.mean <- f2 %>%
    filter(gr == x)
  c3 <- mean(c3.mean$X3lipoxy.pr)
}
f.table$X3lipoxy.pr[9] <- fun3('c')
f.table$X3lipoxy.pr[11] <- fun3('a')
f.table$X3lipoxy.pr[10] <- fun3('cg')
f.table$X3lipoxy.pr[12] <- fun3('ag')

#claster
biochem.scale <-scale(f.table[, 4:21], center = TRUE, scale = TRUE)

dist.biochem <- dist(biochem.scale, method = "euclidean")
clust.biochem <- hclust(dist.biochem, "ward.D")
plot(1:23, clust.biochem$height, type = "b")
plot(clust.biochem, labels = f.table$gr)

dist.biochem <- dist(biochem.scale, method = "euclidean")
clust.biochem <- hclust(dist.biochem, "mcquitty")
plot(1:23, clust.biochem$height, type = "b")
plot(clust.biochem, labels = f.table$gr)

#PCA
pca <- prcomp(f.table[, 4:21], center = TRUE, scale = TRUE)
summary(pca)
plot(pca)
autoplot(pca, data = f.table, colour = "gr", loadings = TRUE, loadings.label = TRUE, size = 5) +
  theme_bw()
  