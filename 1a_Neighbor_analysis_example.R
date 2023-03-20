#This file contains data and code to reproduce the graphical explanation of the 
#neighbour smoothing method
library(DepecheR)
series11 <- rnorm(10)
series12 <- rnorm(10)
series21 <- rnorm(10)
series22 <- rnorm(10)
series31 <- rnorm(10)
series32 <- rnorm(10)

series11 <- series11-2.5
series12 <- series12+2.5
series31 <- series31+2.5
series32 <- series32+2.5

fullX <- c(series11, series21, series31)
fullY <- c(series12, series22, series32)
plot(fullX, fullY)

numVec <- c(12, 12, 12, 12, 7, 12, 7, 12, 12, 12, 7, 7,7,12,1,7,7,12,7,1,1,1,1,1,1,7,1,1,7,1)

colVec <- dColorVector(1:12, colorScale = c("blue", "white", "red"))

plot(fullX, fullY, pch = 21, col = colVec)

fullDf <- data.frame(fullX, fullY, facVec = factor(numVec, levels = c(1,7,12)))


dir.create("Neighbor_analysis_example")
library(ggplot2)
ggplot(fullDf, aes(x=fullX, y=fullY, fill = facVec)) +
  geom_point(size=4, shape=21) + 
  scale_fill_manual(values = c("#0000FF", "#FFE7E7", "#FF0000")) +
  theme_bw() + coord_fixed()
ggsave("Neighbor_analysis_example/Raw_values.pdf")
#Here we create the neighbor analysis

smoothVec <- round(neighSmooth(numVec, kNeighK = 10, fullDf[,1:2]), digits = 0)

fullDf$smooth <- smoothVec

#And for convenience, we add two extra dots that will be removed afterwards
fullDf[31:32,] <- NA
fullDf[31,1] <- -3
fullDf[31,2] <- -2
fullDf[32,1] <- -2.5
fullDf[32,2] <- -2
fullDf$smooth[31] <- 1
fullDf$smooth[32] <- 12

ggplot(fullDf, aes(x=fullX, y=fullY, fill = smooth)) +
  geom_point(size=4, shape=21) + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 6.5) +
  theme_bw()
ggsave("Neighbor_analysis_example/Smooth_values.pdf")

ggplot(fullDf, aes(x=fullX, y=fullY, fill = facVec)) +
  geom_point(size=4, shape=21) + 
  scale_fill_manual(values = c("#0000FF", "#FFE7E7", "#FF0000")) +
  theme_bw() +theme(legend.position="none")
ggsave("Neighbor_analysis_example/Raw_values_no_legend.pdf", width = 5, height = 5)

ggplot(fullDf, aes(x=fullX, y=fullY, fill = smooth)) +
  geom_point(size=4, shape=21) + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 6.5) +
  theme_bw() +theme(legend.position="none")
ggsave("Neighbor_analysis_example/Smooth_values_no_legend.pdf", width = 5, height = 5)

