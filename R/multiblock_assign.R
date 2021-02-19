#' @title Random matrix
#'
#' Create a random data matrix
#' @param nrow Number of rows
#' @param ncol Number of colums
#' @return A matrix of dimension (nrow, ncol)
#' @importFrom stats runif
#' @export
randomData <- function(nrow, ncol) {
    return (matrix(runif(nrow*ncol, min=0, max=1), nrow, ncol))
}


#' @title Trace
#'
#' Trace of a matrix
#' @param x A symmetric numeric matrix or data frame.
#' @return tr(x)
tr <- function(x) {
    if (nrow(x) != ncol(x)) {
    	stop("x should be symmetric.")
    }
    y <- apply(x, c(1,2), as.numeric)
    return (sum(diag(y)))
}


#' @title Matrix scaling
#'
#' Scale a matrix to trace of 1
#' @param x A numeric matrix or data frame.
#' @return A scaled matrix with tr(xx') = 1
#' @export
trscale <- function(x, na.rm = TRUE) {
    y <- apply(x, c(1,2), as.numeric)
    if (na.rm) {
        y <- y[!is.na(rowSums(y)), , drop=F]
    } else {
        y[is.na(y)] <- 0
    }
    y <- y[order(rownames(y)), ]
    return (y/sqrt(tr(crossprod(y))))
}


#' @title Dot product
#'
#' Dot product of two matrices
#' @param x A numeric matrix or data frame.
#' @param y A numeric matrix or data frame.
#' @return xy
#' @export
dotprod <- function(x, y) {
    #stop if y is trivial
    return ()
}
