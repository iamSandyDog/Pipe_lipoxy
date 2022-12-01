library(tidyverse)
library(RColorBrewer)
library(magrittr)
library(ggfortify)
library(corrplot)
library(rgl)


f.table <- read.table("https://raw.githubusercontent.com/iamSandyDog/Pipe_lipoxy/main/.github/workflows/full_biochem.txt", header = TRUE)

f2 <- na.omit(f.table)

#заполнение NA по среднему в группе
fun1 <- function(x) {
  c.mean <- f2 %>%
    filter(gr == x)
  c <- mean(c.mean$f1_lipoxy)
}
f.table$f1_lipoxy[9] <- fun1('c')
f.table$f1_lipoxy[11] <- fun1('a')
f.table$f1_lipoxy[10] <- fun1('cg')
f.table$f1_lipoxy[12] <- fun1('ag')

fun2 <- function(x) {
  c2.mean <- f2 %>%
    filter(gr == x)
  c2 <- mean(c2.mean$f2_lipoxy)
}
f.table$f2_lipoxy[9] <- fun2('c')
f.table$f2_lipoxy[11] <- fun2('a')
f.table$f2_lipoxy[10] <- fun2('cg')
f.table$f2_lipoxy[12] <- fun2('ag')

fun3 <- function(x) {
  c3.mean <- f2 %>%
    filter(gr == x)
  c3 <- mean(c3.mean$f3_lipoxy)
}
f.table$f3_lipoxy[9] <- fun3('c')
f.table$f3_lipoxy[11] <- fun3('a')
f.table$f3_lipoxy[10] <- fun3('cg')
f.table$f3_lipoxy[12] <- fun3('ag')
colorrr = c("#4C9900", "#7F00FF", "#FF0000", "#FF00FF", 
            "#FF00FF", "#FF0000", "#7F00FF", "#4C9900", 
            "#4C9900", "#7F00FF", "#FF0000", "#FF00FF",
            "#FF00FF", "#FF0000", "#7F00FF", "#4C9900", 
            "#4C9900", "#7F00FF", "#FF0000", "#FF00FF",
            "#FF00FF", "#FF0000", "#7F00FF", "#4C9900"
)
f.table <-  f.table %>% mutate(colorrr)

#claster
biochem.scale <-scale(f.table[, 4:21], center = TRUE, scale = TRUE)

dist.biochem <- dist(biochem.scale, method = "euclidean")
clust.biochem <- hclust(dist.biochem, "ward.D")
plot(1:23, clust.biochem$height, type = "b")
plot_clust_eu <- plot(clust.biochem, labels = f.table$gr)

dist.biochem <- dist(biochem.scale, method = "euclidean")
clust.biochem <- hclust(dist.biochem, "mcquitty")
plot(1:23, clust.biochem$height, type = "b")
plot(clust.biochem, labels = f.table$gr)

#PCA
pca <- prcomp(f.table[, 4:21], center = TRUE, scale = TRUE)
summary(pca)
plot(pca)
pca_plot <- autoplot(pca, data = f.table, colour = "gr", loadings = TRUE, loadings.label = TRUE, size = 5) +
  theme_bw()
pca_plot <-  pca_plot +
  scale_colour_manual(values = c("#FF0000", "#FF00FF", "#4C9900", "#7F00FF"))
pca_plot
cor_data <- cor(f.table[, 4:21])
corrplot(cor_data)  
pca_plot2 <- plot(f.table[, 4:21], col = rep(1:4))

#3d plot of pca
pca2 <- princomp(f.table[, 4:21])

scores <- as.data.frame(pca$x)
plot_3d_pca <- plot3d(scores[,1:3], col = f.table$colorrr, size = .75, type = "s") 
text3d(scores[,1:3], 
       texts=f.table$gr, 
       cex= 0.8, pos=4)
play3d(spin3d(axis = c(1,0,0), rpm = 10), duration = 40)

plot_3d_pca2 <- plot3d(scores[,4:6], col = f.table$colorrr, size = .75, type = "s") 
text3d(scores[,4:6], 
       texts=f.table$gr, 
       cex= 0.8, pos=4)
play3d(spin3d(axis = c(1,0,0), rpm = 1), duration = 20)

plot_3d_pca3 <- plot3d(scores[,7:9], col = f.table$colorrr, size = .75, type = "s") 
text3d(scores[,7:9], 
       texts=f.table$gr, 
       cex= 0.8, pos=4)

plot_3d_pca4 <- plot3d(scores[,10:12], col = f.table$colorrr, size = .75, type = "s") 
text3d(scores[,10:12], 
       texts=f.table$gr, 
       cex= 0.8, pos=4)

plot_3d_pca5 <- plot3d(scores[,13:15], col = f.table$colorrr, size = .75, type = "s") 
text3d(scores[,13:15], 
       texts=f.table$gr, 
       cex= 0.8, pos=4)

plot_3d_pca6 <- plot3d(scores[,16:18], col = f.table$colorrr, size = .75, type = "s") 
text3d(scores[,16:18], 
       texts=f.table$gr, 
       cex= 0.8, pos=4)
