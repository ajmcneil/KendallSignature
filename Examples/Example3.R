library(qrmdata)
library(extremalcopula)
library(rgl)
library(gMOIP)
library(xts)

# Scenario: 4 variables 1,2,3,4
# In the absence of any further information what is the set of attainable ranl
# correlations for {1,2}, {1,3} and {1,4}?
partkappa <- 1
names(partkappa) <- "empty"
vert <- findpolytope(partkappa, 4)
kappa <- Amatrix(4) %*% vert
kappawith4 <- t(kappa[c("{1,2}", "{1,3}", "{1,4}"),])
taus <- 2*kappawith4 -1
plot3d(taus, plot=FALSE)
plotHull3D(taus, drawPoints = TRUE,drawLines=FALSE, drawPolygons=TRUE)

# A real 4-dimensional dataset
data(crypto)
?crypto
data <- crypto['2016-12-31/2017-12-31']
X <- (diff(log(data))[-1]) * 100
plot(X, multi.panel=TRUE)
X <- as.matrix(X)

# Suppose we only knew the rank correlations for first 3 cryptocurrencies
# What could we infer about the rank correlations with fourth?
d <- 4
fullsig <- estsignature(X)
fullsig
sig_3missing <- fullsig[c("empty", "{1,2}", "{1,3}", "{2,3}")]

#############################################

sig_3missing
vert <- findpolytope(sig_3missing, d = d)
vert
nvert <- dim(vert)[1]
kappas <- Amatrix(d) %*% t(vert)

# A check that the vertices give the correct rank correlations
weights <- apply(kappas, 2, findextremal)
weights[ weights < 10e-10] <- 0
output <- vector("list", length = nvert)
for (i in 1:nvert)
  output[[i]] <- extmixcorr(weights[, i])
output
2 * sig_3missing - 1

kappas_3missing <- t(kappas[c("{1,4}", "{2,4}", "{3,4}"), ])
taus3 <- 2*kappas_3missing -1
plot3d(taus3, plot=FALSE)

taus3_hull <- taus3[convexHull(taus3)$pts[,"vtx"],]
taus3_hull <- round(taus3_hull, 8)
plotHull3D(taus3_hull, drawPoints = TRUE, drawLines=FALSE, drawPolygons=TRUE)

####################

# Now suppose we are told the rank correlation {1,4}
sig_2missing <- fullsig[c("empty", "{1,2}", "{1,3}", "{1,4}", "{2,3}")]
sig_2missing

vert <- findpolytope(sig_2missing, d = d)
vert
nvert <- dim(vert)[1]
kappas <- Amatrix(d) %*% t(vert)

# A check that the vertices give the correct rank correlations
weights <- apply(kappas, 2, findextremal)
output <- vector("list", length = nvert)
for (i in 1:nvert)
  output[[i]] <- extmixcorr(weights[, i])
output
2 * sig_2missing - 1

kappas_2missing <- t(kappas[c("{2,4}", "{3,4}"), ])
taus2 <- 2*kappas_2missing -1

plot(taus2, col = 2, xlim = c(-1,1), ylim = c(-1,1))
taus2_hull <- taus2[convexHull(taus2)$pts[,"vtx"],]
lines(rbind(taus2_hull, taus2_hull[1,]))
points(2*fullsig["{2,4}"]-1, 2*fullsig["{3,4}"]-1)

#################################################

plane_pts <- cbind(rep(as.numeric(2*fullsig["{1,4}"]-1), dim(taus2_hull)[1]), taus2_hull)
dimnames(plane_pts)[[2]] <- c("{1,4}", "{2,4}", "{3,4}")
plotHull3D(plane_pts)
