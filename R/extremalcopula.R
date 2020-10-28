#' Binary coding set for extremal copulas
#'
#' @param k a vector of integers between one and 2^(d-1) specifying extremal
#' copulas in dimension d
#' @param d the dimension of the extremal copula
#'
#' @return a matrix with as many rows as the length of k and d columns containing zeros and ones
#' @export
#'
#' @examples
#' codingset(c(1:8), 4)
codingset <- function(k, d){
  if (sum(!(k %in% (1:2^(d-1)))) > 0)
    stop("Illegal value for k")
  scode <- binaryLogic::as.binary(k - 1, n = d)
  t(sapply(scode, function(listitem){
    as.integer(as.logical(listitem))
  }))
}

#' Sample from extremal copulas
#'
#' @param k a vector of integers between one and 2^(d-1) specifying extremal
#' copulas in dimension d
#' @param d the dimension of the extremal copula
#'
#' @return a matrix with as many rows as the length of k
#' such that each row is a realisation from an extremal copula in dimension d
#' @export
#'
#' @examples
#' rextremal(rep(2, 10), 4)
rextremal <- function(k, d){
  V <- codingset(k, d)
  U <- matrix(runif(length(k)), nrow = length(k), ncol = d)
  (U^V)*((1-U)^(1-V))
}

#' Sample from extremal mixture copula
#'
#' @param w vector of weights specifing the extremal mixture copula
#' @param n number of realisations to generate
#'
#' @return a matrix with n rows where each row is a a realisation from the
#' extremal mixture copula
#' @export
#'
#' @examples
#' rextremalmix(rep(1, 8)/8, 100)
rextremalmix <-function(w, n = 1){
  if (abs(sum(w)-1) > 1e-5)
    stop("Infeasible weight vector w (must sum to 1)")
  d <- as.integer(log(length(w))/log(2)+1)
  if (!(identical(as.integer(2^(d-1)), length(w))))
    stop("Infeasible weight vector w (wrong length)")
  mmat <- rmultinom(n, 1, w)
  Nvec <- apply(mmat, 2, function(v){
    (1:length(v))[as.logical(v)]})
  rextremal(Nvec, d)
}

#' Orthant probability for a multivariate Gaussian vector
#'
#' @param P correlation matrix of the normal distribution
#'
#' @return the lower orthant probability
#' @export
#'
#' @examples
#' P <- 0.7 * matrix(1, nrow = 3, ncol = 3) + 0.3 * diag(3)
#' Gaussorthant(P)
Gaussorthant<- function(P){
  d <- length(diag(P))
  2*as.numeric(mvtnorm::pmvnorm(lower = rep(0, d), sigma = P))
}

#' Construct A matrix in dimension d
#'
#' @param d the dimension
#' @param named logical value to determine whether rows and columns are named
#'
#' @return a square matrix with 2^(d-1) rows and columns
#' @export
#'
#' @examples
#' Amatrix(4)
Amatrix <- function(d, named = TRUE){
  kmax <- 2 ^ (d - 1)
  A <- matrix(1, nrow = kmax, ncol  = kmax)
  if (named){
    paste2 <- function(v){paste("{",paste(v,collapse=","),"}",sep="")}
    dimnames(A)[[1]] <- rep("empty", kmax)
  }
  coding <- codingset(1:kmax, d)
  layers <- floor(d / 2)
  n_layer <- sapply(1:layers, function(x, d){choose(d, 2 * x)}, d = d)
  pos <- 1
  for (j in 1:layers) {
    nelements <- 2 * j
    combs <- gtools::combinations(d, nelements, repeats.allowed = FALSE)
    for (i in 1:n_layer[j]) {
      thecomb <- combs[i, ]
      A[pos + 1,] <-
        apply(
          coding[, thecomb],
          1,
          FUN = function(v) {
            (length(unique(v)) == 1)
          }
        )
      if (named)
        dimnames(A)[[1]][pos+1] <- paste2(thecomb)
      pos <- pos + 1
    }
  }
  if (named)
    dimnames(A)[[2]] <- apply(codingset(1:kmax, d), 1, paste2)
  A
}

#' Elliptical concordance signature
#'
#' @param P correlation matrix of elliptical distribution
#' @param limit whether a limit should be applied to the size of the correlation matrix
#' for accuracy in the calculation of Gaussian orthant probabilities
#' @param named logical value to determine whether signature is named
#'
#' @return a vector containing the concordance signature
#' @export
#'
#' @examples
#' P <- 0.7 * matrix(1, nrow = 3, ncol = 3) + 0.3 * diag(3)
#' consignature(P)
consignature <- function(P, limit = TRUE, named = TRUE){
  if (!is.matrix(P))
    stop("P is not a matrix")
  if (nrow(P) != ncol(P))
    stop("P is not a square matrix")
  evalues <- eigen(P)$values
  if (min(evalues) < 0)
    stop("P is not a psd matrix")
  d <- nrow(P)
  if ((sum(diag(P)) != d) | (min(P) < -1) | (max(P) > 1))
    stop("P is not a correlation matrix")
  layers <- floor(d / 2)
  if (limit & (layers > 3))
    stop("Maximum of 3 layers imposed for accuracy")
  kmax <- 2 ^ (d - 1)
  B <- rep(1, kmax)
  if (named){
    paste2 <- function(v){paste("{",paste(v,collapse=","),"}",sep="")}
    names(B)[1] <- "empty"
  }
  coding <- codingset(1:kmax, d)
  layers <- floor(d / 2)
  n_layer <- sapply(1:layers, function(x, d){choose(d, 2 * x)}, d = d)
  pos <- 1
  for (j in 1:layers) {
    nelements <- 2 * j
    combs <- gtools::combinations(d, nelements, repeats.allowed = FALSE)
    for (i in 1:n_layer[j]) {
      thecomb <- combs[i, ]
      B[pos+1] <- Gaussorthant(P[thecomb, thecomb])
       if (named)
         names(B)[pos+1] <- paste2(thecomb)
      pos <- pos + 1
    }
  }
  B
}

#' Determine extremal mixture representation from full concordance signature
#'
#' @param B full concordance signature
#' @param named logical value to determine whether weights are named
#'
#' @return vector of weights of extremal copula
#' @export
#'
#' @examples
#' P <- 0.7 * matrix(1, nrow = 3, ncol = 3) + 0.3 * diag(3)
#' B <- consignature(P)
#' findextremal(B)
findextremal <- function(B, named = TRUE) {
  d <- as.integer(log(length(B))/log(2)+1)
  if (!(identical(as.integer(2^(d-1)), length(B))))
    stop("Infeasible full signature (wrong length)")
  A <- Amatrix(d, named = named)
  solve(A, B)
}

#' Make table of extremal representation
#'
#' @param B full concordance signature
#'
#' @return a dataframe that can be readily turned into latex table or knitr table
#' @export
#'
#' @examples
#' P <- 0.7 * matrix(1, nrow = 3, ncol = 3) + 0.3 * diag(3)
#' B <- consignature(P)
#' extremaltable(B)
extremaltable <- function(B){
  weights <- findextremal(B)
  d <- as.integer(log(length(B))/log(2)+1)
  kmax <- length(B)
  names(B)[1] <- paste("$\\emptyset$")
  names(B)[2:kmax] <- paste("$\\{", names(B)[2:kmax], "\\}$", sep ="")
  paste2 <- function(v){paste("$\\{",paste(v,collapse=","),"\\}$",sep="")}
  names(weights) <- apply(codingset(1:kmax, d), 1, paste2)
  example <- data.frame(S=names(weights),
                        w=weights,
                        I=names(B),
                        kappa=B,
                        row.names = (1:kmax))
  names(example) <- c("S", "w", "I", "$\\mathcal{K}$")
  example
}

#' Extremal correlation matrices
#'
#' @param k a vector of integers between one and 2^(d-1) specifying extremal
#' copulas in dimension d
#' @param d the dimension of the extremal correlation matrix
#'
#' @return a list of matrices if k is a vector or a single matrix if k is scalar
#' @export
#'
#' @examples
#' extcorr(2, 4)
#' extcorr(1:8, 4)
extcorr <- function(k, d){
  if (sum(!(k %in% (1:2^(d-1)))) > 0)
    stop("Illegal value for k")
  output <- lapply(k, function(k,d){
    S <- as.vector(codingset(k, d))
    V <- (2*S -1)
    outer(V, V)
  },
  d = d)
  if (length(output) == 1)
    output <- output[[1]]
  output
}

#' Determine  mixture representation from partial concordance signature
#'
#' @param sig partial concordance signature
#' @param maximize a logical variable stating whether unspecified concordance
#' probabilities should be jointly maximized rather than minimized which is the
#' default
#'
#' @return vector of weights of a possible extremal copula
#' @export
#'
#' @examples
#' kappa2 <- 0.7
#' sig <- c(1, rep(kappa2, 6))
#' attainable(sig)
#' attainable(sig, maximize = TRUE)
attainable <- function(sig, maximize = FALSE) {
  d <- ceiling(log(length(sig))/log(2)+1)
  n <- length(sig)
  if (identical(as.integer(2^(d-1)), n))
    findextremal(sig)
  A <- Amatrix(d)
  A1 <- A[1:n, ]
  A2 <- A[-(1:n), , drop = FALSE]
  B <- rep(0, nrow(A2))
  if (maximize)
    B <- rep(1, nrow(A2))
  tmp <- limSolve::lsei(A = A2, B = B, E = A1, F = sig, G = diag(2^(d-1)), H = rep(0,2^(d-1)))
  if (tmp$IsError)
    stop("Solution on simplex not found")
  return(tmp$X)
}

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
#' data <- rmvnorm(n = 1000, sigma = diag(3))
#' table(quadrants(sign(data)))
#' table(quadrants(sign(data), diagonal = TRUE))
quadrants <- function(signs, diagonal = FALSE, pretty = FALSE) {
  if (diagonal)
    signs <- -signs * signs[, 1]
  bdata <- (signs + 1) / 2
  if (pretty)
    paste2 <- function(v){paste("{",paste(v,collapse=","),"}",sep="")}
  else
    paste2 <- function(v){paste(v,collapse="")}
  apply(bdata, 1, paste2)
}


#' Compute signs of differences of pairs
#'
#' @param data a matrix of data of length n
#'
#' @return a matrix of length n*(n-1)/2
#' @export
#'
#' @examples
#' data <- rmvnorm(n = 10, sigma = diag(3))
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
  vals <- apply(codingset((1:2^(d-1)), d), 1, paste, collapse = "")
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
  kappa <- Amatrix(d, name = TRUE) %*% weights
  drop(kappa)
}

#' Find polytope of attainable weights from partial concordance signature
#'
#' @param sig partial concordance signature
#' @param d dimension of copula
#'
#' @return a matrix where each row is a vertex of the polytope
#' @export
#'
#' @examples
#' tau <- c(-0.19, -0.29, 0.49, -0.34, 0.3, -0.79)
#' kappa <- (1+tau)/2
#' findpolytope(c(1, kappa), d = 4)
findpolytope <- function(sig, d){
  if (sig[1] != 1)
    stop("Signature must have a 1 in first position for sum constraint")
  m <- 2^(d-1)
  n <- length(sig)
  if (identical(m, n))
    findextremal(sig)
  A <- Amatrix(d)[(1:n),]
  qux <- rcdd::makeH(- diag(m), rep(0, m), A, sig)
  results <- scdd(qux)
  if (dim(results$output)[1] == 0)
    stop("Signature not attainable")
  results$output[,-c(1,2)]
}
