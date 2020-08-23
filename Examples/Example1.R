library(extremalcopula)

codingset(3,6)

data <- rextremal(k = rep(2,100), d = 4)
data
cor(data)
cor(data, method = "kendall")


w <- c(1,2,3,4,2,3,3,2)/20
data <- rextremalmix(w, 1000)
cor(data)
cor(data, method ="kendall")
cor(data, method ="spearman")
# although the population matrices are identical, the estimates typically differ

# Equicorrelation example
# Sample extremal mixture with same concordance structure as an elliptical
rho <- 0.5
kappa = 1-acos(rho)/pi

P <- rho* matrix(1, nrow = 4, ncol = 4) + (1-rho)*diag(4)
P
kappa4 <- Gaussorthant(P)
w1 <- (3*kappa-1)/2 - kappa4
w2 <- 1 - 2*kappa + kappa4
w <- c(kappa4, w1, w1, w2, w1, w2, w2, w1)
data <- rextremalmix(w, 5000)
cor(data)
cor(data, method ="kendall")
cor(data, method ="spearman")

# True correlation of model
2*asin(P/pi)


# Examples
Amatrix(2)
Amatrix(3)
Amatrix(4)
Amatrix(5)

d <- 6
corvals <- (1:15) / 16
P <- copula::p2P(corvals)

consignature(P)
findextremal(consignature(P))
tmp <- extremaltable(consignature(P))
tmp

print(xtable::xtable(tmp,digits=4),sanitize.text.function=function(x){x})

d <- 4
rho <- -1/4
corvals <- rep(rho,choose(d,2))

P <- copula::p2P(corvals)
P
B <- consignature(P)
B
A <- Amatrix(nrow(P))
A
weights <- solve(A, B)
sum(weights)
U <- rextremalmix(weights, 1000)
cor(U)
2*asin(P)/pi
