#Data sets:
#
# 1) Isotopic Composition Plutonium Batches
#https://rdrr.io/cran/cluster/man/pluton.html
#https://wis.kuleuven.be/stat/robust/papers/1996/rousseeuwkaufmantrauwaert-fuzzyclusteringusingscat.pdf

# library(cluster)
# data(pluton)
# plot(pluton)
# 
# #2) Diday and Govaert 
# set.seed(123)
# mu1 <- c(0,3)
# mu2 <- c(0,0)
# mu3 <- c(4,3)
# sigma1 <- matrix(c(0.25, 0, 0, 0.25), nrow = 2)
# sigma2 <- matrix(c(4, 1.7, 1.7, 1), nrow = 2)
# sigma3 <- matrix(c(4, -1.7, -1.7, 1), nrow = 2)
# library(MASS)
# data1 <- mvrnorm(n = 50, mu = mu1, Sigma = sigma1)
# data2 <- mvrnorm(n = 50, mu = mu2, Sigma = sigma2)
# data3 <- mvrnorm(n = 50, mu = mu3, Sigma = sigma3)
# data <- rbind(data1, data2, data3)
# plot(data)
# 
# #Example 3 (Tibshirani)
# 
# tions, with mean (0; 0) and (m; 0);
# 
# where m = 1; : : : ; 7:
