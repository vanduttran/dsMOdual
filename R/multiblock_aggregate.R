#' @title Row means
#'
#' Row means of linear transformations of a matrix 
#' @param x A numeric matrix
#' @param y A list of symmetric matrices of dimension (ncol(x), ncol(x))
#' @return Row means of x if y = NULL, or row means x %*% yy for each matrix yy in y otherwise
#' @export
rowmeans <- function(x, y = NULL) {
    ## if (is.null(dim(x)) || min(dim(x)) < 10) {
    ##     stop("x should be a matrix with two dimensions higher than 10.")
    ## }
    if (is.null(y)) {
        return (matrix(rowMeans(x), ncol=1))
    } else if (!all(sapply(y, isSymmetric))) {
        stop("y is not all symmetric.")
    } else {
        return (lapply(y, function(yy) matrix(rowMeans(tcrossprod(x, yy)), ncol=1)))
    }
}

#' @title Matrix cross product
#' 
#' Calculates the cross product t(x) \%*\% x
#' @param x A numeric matrix
#' @return t(x) \%*\% x
#' @export
crossProd <- function(x) {
    ## if (is.null(dim(x)) || min(dim(x)) < 10) {
    ##     stop("x should be a matrix with two dimensions higher than 10.")
    ## }
    return (crossprod(x))
}

#' @title Matrix cross product
#' 
#' Calculates the cross product x \%*\% t(x)
#' @param x A numeric matrix
#' @return x \%*\% t(x)
#' @export
tcrossProd <- function(x) {
    ## if (is.null(dim(x)) || min(dim(x)) < 10) {
    ##     stop("x should be a matrix with two dimensions higher than 10.")
    ## }
    return (tcrossprod(x))
}

#' @title Matrix triple product
#'
#' Calculate the triple product
#' @param x A numeric matrix
#' @param y A list of symmatric numeric matrices of dimension (ncol(x), ncol(x))
#' @return List of x %*% y %*% t(x)
#' @export
tripleProd <- function(x, y) {
    if (!all(sapply(y, isSymmetric))) {
        stop("y is not all symmetric.")
    } else {
        return (lapply(y, function(yy) tcrossprod(x, tcrossprod(x, yy))))
    }
}
