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
#' kappa <- c(1, (1 + tau)/2)
#' findpolytope(kappa, d = 4)
findpolytope <- function(sig, d){
  if (sig[1] != 1)
    stop("Signature must have a 1 in first position for sum constraint")
  if (is.null(names(sig)))
    names(sig) <- evenpowerset(d)[1:length(sig)]
  m <- 2^(d-1)
  n <- length(sig)
  if (identical(m, n))
    findextremal(sig)
  A <- Amatrix(d)[names(sig),]
  qux <- rcdd::makeH(- diag(m), rep(0, m), A, sig)
  results <- rcdd::scdd(qux)
  if (dim(results$output)[1] == 0)
    stop("Signature not attainable")
  output <- results$output[,-c(1,2)]
  rnames <- paste("v", (1:dim(output)[1]), sep="")
  cnames <- apply(binarycoding((1:2^(d-1)), d), 1, setpaste)
  dimnames(output) <- list(rnames, cnames)
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
