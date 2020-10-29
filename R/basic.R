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
#' binarycoding(c(1:8), 4)
binarycoding <- function(k, d){
  if (sum(!(k %in% (1:2^(d-1)))) > 0)
    stop("Illegal value for k")
  scode <- binaryLogic::as.binary(k - 1, n = d)
  t(sapply(scode, function(listitem){
    as.integer(as.logical(listitem))
  }))
}

#' Lexicographically ordered even power set
#'
#' @param d the dimension of the set
#'
#' @return a character vector with the elements of the even power set
#' @export
#'
#' @examples
#' evenpowerset(4)
evenpowerset <- function(d){
  output <- rep("empty",2^(d-1))
  pointer <- 1
  for (i in 1:floor(d/2))
  {
    combs <- gtools::combinations(d, 2*i, repeats.allowed = FALSE)
    combsaschars <- apply(combs, 1, setpaste)
    output[(pointer + 1): (pointer + choose(d, 2*i))] <- combsaschars
    pointer <- pointer + choose(d, 2*i)
  }
  output
}

#' Paste vector in set notation
#'
#' @param v a numeric vector
#'
#' @return a character string gving vector as set
#' @export
#'
setpaste <- function(v){
  paste("{",paste(v,collapse=","),"}",sep="")
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
    S <- as.vector(binarycoding(k, d))
    V <- (2*S -1)
    outer(V, V)
  },
  d = d)
  if (length(output) == 1)
    output <- output[[1]]
  output
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
  V <- binarycoding(k, d)
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
