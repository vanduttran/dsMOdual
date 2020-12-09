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
        return (matrix(rowMeans(x), ncol=1, dimnames=list(rownames(x), "mean")))
    } else if (!all(sapply(y, isSymmetric))) {
        stop("y is not all symmetric.")
    } else {
        return (lapply(y, function(yy) matrix(rowMeans(tcrossprod(x, yy)), ncol=1, dimnames=list(rownames(x), "mean"))))
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


#' @title Cross login
#'
#' Cross login # rather on client side
#' @param logins An encoded dataframe with server, url, user, password, and driver fields.
#' @export
crossLogin <- function(logins) {
    loginfo <- dsSwissKnife:::.decode.arg(logins)
    myDf <- data.frame(server=loginfo$server,
                       url=loginfo$url,
                       user=loginfo$user,
                       password=loginfo$password,
                       driver=loginfo$driver)
    x <- tryCatch(DSI::datashield.login(myDf), error=function(e) return (sessionInfo()))
    return (x)
}


#' @title Cross aggregate
#'
#' Cross aggregate # rather on client side
#' @param opal An opal object or list of opal objects.
#' @param expr An encoded expression to evaluate.
#' @param wait See DSI::datashield.aggreate options. Default: FALSE.
#' @param async See DSI::datashield.aggreate options. Default: TRUE.
#' @import DSI
#' @export
crossAggregate <- function(opal, expr, wait = F, async = T) {
    expr <- dsSwissKnife:::.decode.arg(expr)
    DSI::datashield.aggregate(opal=opal, expr=as.symbol(expr), wait=wait, async=async)
}


#' @title Cross assign
#'
#' Cross assign # rather on client side
#' @param opal An opal object or list of opal objects.
#' @param symbol Name of an R symbol.
#' @param value A variable name or an R epxression with allowed assign function calls.
#' @param value.call A logical value, TRUE if value is function call, FALSE if value is a variable name.
#' @param wait See DSI::datashield.aggreate options. Default: FALSE.
#' @param async See DSI::datashield.aggreate options. Default: TRUE.
#' @import DSI
#' @export
crossAssign <- function(opal, symbol, value, value.call, variables = NULL, wait = F, async = T) {
    value <- dsSwissKnife:::.decode.arg(value)
    variables <- dsSwissKnife:::.decode.arg(variables)
    DSI::datashield.assign(opal=opal, symbol=symbol, value=ifelse(value.call, as.symbol(value), value), variables=variables, wait=wait, async=async)
}
