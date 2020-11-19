library(qrmdata)
library(extremalcopula)
library(rgl)
library(gMOIP)
library(xts)


# A real 4-dimensional dataset
data(crypto)
?crypto
data <- crypto['2016-12-31/2017-12-31']
X <- (diff(log(data))[-1]) * 100
plot(X, multi.panel=TRUE)
X <- as.matrix(X)

# First estimate the signature and find a matching extremal mixture copula

fullsig <- estsignature(X)
fullsig
findextremal(fullsig)

# Suppose we only knew the rank correlations for first 3 cryptocurrencies
# What could we infer about the rank correlations with fourth?

sig_3missing <- fullsig[c("empty", "{1,2}", "{1,3}", "{2,3}")]
sig_3missing
d <- 4
vert <- findpolytope(sig_3missing, d = d)
vert
nvert <- dim(vert)[1]

# A check that the vertices give the correct rank correlations
output <- vector("list", length = nvert)
for (i in 1:nvert)
  output[[i]] <- extmixcorr(vert[i,])
output[[nvert]] # for example for final vertex
2 * sig_3missing - 1

# now we find the corresponding even concordance signatures (kappas)
# and the corresponding values for the missing pairwise concordance probabilities
kappas <- Amatrix(d) %*% t(vert)
(kappas_3missing <- t(kappas[c("{1,4}", "{2,4}", "{3,4}"), ]))

# Now we show the set of attainable Kendall's tau values
taus3 <- 2*kappas_3missing -1
taus3_hull <- taus3[convexHull(taus3)$pts[,"vtx"],]
taus3_hull <- round(taus3_hull, 10)
par3d(cex=2.0)
plot3d(taus3_hull, plot=FALSE)
plotHull3D(taus3_hull, drawPoints = TRUE, drawLines=FALSE, drawPolygons=TRUE)

####################

# Now suppose we are told the rank correlation {1,4}
sig_2missing <- fullsig[c("empty", "{1,2}", "{1,3}", "{1,4}", "{2,3}")]
sig_2missing

# we again find the corresponding even concordance signatures (kappas)
# and the corresponding values for the missing pairwise concordance probabi
vert <- findpolytope(sig_2missing, d = d)
vert
nvert <- dim(vert)[1]
kappas <- Amatrix(d) %*% t(vert)
(kappas_2missing <- t(kappas[c("{2,4}", "{3,4}"), ]))
taus2 <- 2*kappas_2missing -1

# The plot is now in the plane
taus2_hull <- taus2[convexHull(taus2)$pts[,"vtx"],]

pdf("crypto2D.pdf")
par(cex=1.2, cex.axis=1.2)
plot(taus2_hull, col = 2, xlim = c(-1,1), ylim = c(-1,1))
lines(rbind(taus2_hull, taus2_hull[1,]))


# Of course the true values must be in this set
points(2*fullsig["{2,4}"]-1, 2*fullsig["{3,4}"]-1)
dev.off()
#################################################

plane_pts <- cbind(rep(as.numeric(2*fullsig["{1,4}"]-1), dim(taus2_hull)[1]), taus2_hull)
dimnames(plane_pts)[[2]] <- c("{1,4}", "{2,4}", "{3,4}")


plotHull3D(plane_pts,drawPolygons = FALSE)

rgl.postscript("crypto3D.pdf",fmt="pdf")
