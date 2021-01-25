#' Identify quadrants of multivariate data points
#'
#' @param signs a matrix containing signs of the data points
#' @param diagonal whether quadrants on same diagonal should be amalgamated
#' @param pretty whether sets should be printed in pretty fashion
#'
#' @return a vector containing the quadrant of each point
#' @export
#'
#' @examples
#' data <- mvtnorm::rmvnorm(n = 1000, sigma = diag(3))
#' table(quadrants(sign(data)))
#' table(quadrants(sign(data), diagonal = TRUE))
quadrants <- function(signs, diagonal = FALSE, pretty = FALSE) {
  if (diagonal)
    signs <- -signs * signs[, 1]
  bdata <- (signs + 1) / 2
  if (pretty)
    output <- apply(bdata, 1, setpaste)
  else
    output <- apply(bdata, 1, paste, collapse="")
  output
}

#' Compute signs of differences of pairs
#'
#' @param data a matrix of data of length n
#'
#' @return a matrix of length n*(n-1)/2
#' @export
#'
#' @examples
#' data <- mvtnorm::rmvnorm(n = 10, sigma = diag(3))
#' signsofpairs(data)
signsofpairs <- function(data){
  d <- ncol(data)
  n <- nrow(data)
  output <- matrix(NA, nrow =n*(n-1)/2, ncol = d)
  k <- 1
  for (i in 1:(n-1))
    for (j in (i+1):n){
      output[k, ] <- sign(data[i,]- data[j,])
      k <- k+1
    }
  output
}

#' Estimate a concordance signature
#'
#' @param data a matrix of data observations
#'
#' @return a vector giving the estimated concordance signature
#' @export
#'
#' @examples
#' data <- as.matrix(iris[,1:4])
#' estsignature(data)
estsignature <- function(data){
  sgns <- signsofpairs(data)
  n <- dim(data)[1]
  d <- dim(data)[2]
  vals <- apply(binarycoding((1:2^(d-1)), d), 1, paste, collapse = "")
  ties <- apply(sgns, 1, function(v){0 %in% v})
  t1 <- quadrants(sgns[!ties, ], diagonal = TRUE)
  t1 <- table(factor(t1, levels = vals))
  if (sum(ties) > 0){
    v2 <- sgns[ties,]
    v3 <- v2
    v2[v2==0] <- 1
    v3[v3==0] <- -1
    t2 <- table(factor(quadrants(v2, diagonal = TRUE), levels = vals))
    t3 <- table(factor(quadrants(v3, diagonal = TRUE), levels = vals))
    t1 <- t1 + (t2+t3)/2
  }
  weights <- 2*t1/(n*(n-1))
  kappa <- Amatrix(d, named = TRUE) %*% weights
  drop(kappa)
}
