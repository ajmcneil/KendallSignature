library(qrmdata)
library(extremalcopula)
library(rgl)
library(gMOIP)
library(xts)

data(crypto)
data <- crypto['2016-12-31/2017-12-31']
X <- (diff(log(data))[-1]) * 100
plot(X, multi.panel=TRUE)
X <- as.matrix(X)

d <- 4
fullsig <- estsignature(X)
fullsig
sig_3missing <- fullsig[c("empty", "{1,2}", "{1,3}", "{2,3}")]
missing3 <- fullsig[c("{1,4}", "{2,4}", "{3,4}")]
sig_2missing <- fullsig[c("empty", "{1,2}", "{1,3}", "{1,4}", "{2,3}")]
missing2 <- fullsig[c("{2,4}", "{3,4}")]

####################

vert <- findpolytope(sig_2missing, d = d)
vert
nvert <- dim(vert)[1]
kappas <- Amatrix(d) %*% t(vert)

weights <- apply(kappas, 2, findextremal)

taumat <- function(v) {
  Reduce('+', Map('*', v, extcorr(1:2 ^ (d - 1), d)))
}
output <- vector("list", length = nvert)
for (i in 1:nvert)
  output[[i]] <- taumat(weights[, i])
output
2 * sig_2missing - 1


kappas_2missing <- t(kappas[names(missing2), ])

plot(kappas_2missing, col = 2, xlim = c(0,1), ylim = c(0,1))
kappas_hull2 <- kappas_2missing[convexHull(kappas_2missing)$pts[,"vtx"],]
lines(rbind(kappas_hull2,kappas_hull2[1,]))
points(missing2[1], missing2[2])

#############################################

vert <- findpolytope(sig_3missing, d = d)
vert
nvert <- dim(vert)[1]
kappas <- Amatrix(d) %*% t(vert)

weights <- apply(kappas, 2, findextremal)
weights[ weights < 10e-10] <- 0

output <- vector("list", length = nvert)
for (i in 1:d ^ 2)
  output[[i]] <- taumat(weights[, i])
output
2 * sig_3missing - 1

kappas_3missing <- t(kappas[names(missing3), ])
kappas_hull3 <- kappas_3missing[convexHull(kappas_3missing)$pts[,"vtx"],]
kappas_hull3 <- round(kappas_hull3, 8)

plot3d(kappas_hull3, plot=FALSE)
plotHull3D(kappas_hull3, drawPoints = TRUE,drawLines=FALSE, drawPolygons=TRUE)



plane_pts <- cbind(rep(as.numeric(missing3[1]), dim(kappas_hull2)[1]), kappas_hull2)
dimnames(plane_pts)[[2]] <- names(missing3)
plotHull3D(plane_pts)
