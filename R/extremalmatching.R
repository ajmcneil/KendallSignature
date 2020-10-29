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
  coding <- binarycoding(1:kmax, d)
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
      pos <- pos + 1
    }
  }
  if (named){
    dimnames(A)[[1]] <- evenpowerset(d)
    dimnames(A)[[2]] <- apply(coding, 1, setpaste)
  }
  A
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
  names(weights) <- apply(binarycoding(1:kmax, d), 1, paste2)
  example <- data.frame(S=names(weights),
                        w=weights,
                        I=names(B),
                        kappa=B,
                        row.names = (1:kmax))
  names(example) <- c("S", "w", "I", "$\\mathcal{K}$")
  example
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
  coding <- binarycoding(1:kmax, d)
  layers <- floor(d / 2)
  n_layer <- sapply(1:layers, function(x, d){choose(d, 2 * x)}, d = d)
  pos <- 1
  for (j in 1:layers) {
    nelements <- 2 * j
    combs <- gtools::combinations(d, nelements, repeats.allowed = FALSE)
    for (i in 1:n_layer[j]) {
      thecomb <- combs[i, ]
      B[pos+1] <- Gaussorthant(P[thecomb, thecomb])
      pos <- pos + 1
    }
  }
  if (named)
    names(B) <- evenpowerset(d)
  B
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





