library(mvtnorm)
P <- matrix(c(1, 0.4, 0.7, 0.3, 0.4, 1, -0.1, -0.2, 0.7, -0.1, 1, 0.5, 0.3, -0.2, 0.5, 1), nrow = 4)
eigen(P)
data <- rmvnorm(1000, sigma = P)

res <- quadrants(sign(data), diagonal = TRUE)
table(res)/length(res)

sgns <- signsofpairs(data)
res <- quadrants(sgns, diagonal = TRUE)
weights <- table(res)/length(res)
weights
kappa <- Amatrix(4) %*% weights
kappa
2*kappa -1
cor(data, method = "kendall")

res <- quadrants(sgns)
weights <- table(res)/length(res)
sum(weights[c(1,length(weights))])
