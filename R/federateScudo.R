#' @title Federated ComDim
#' @description Function for ComDim federated analysis on the virtual cohort combining multiple cohorts
#' Finding common dimensions in multitable data (Xk, k=1...K)
#' @usage federateComDim(X,group,H=2,scale="none",option="none",threshold=1e-10)
#'
#' @param loginFD Login information of the FD server
#' @param logins Login information of the servers containing cohort data
#' @param variables Variables
#' @param TOL Tolerance of 0
#' @param XX  :	        list of dataframes XX = X %*% t(X)
#' @param group :       named list of variables for each table
#' @param H :           number of common dimensions
#' @param scale  either value "none" / "sd" indicating the same scaling for all tables or a vector of scaling ("none" / "sd") for each table
#' @param option weighting of te tables \cr
#'        "none" :  no weighting of the tables - (default) \cr
#'      "uniform": weighting to set the table at the same inertia \cr
#' @param threshold if the difference of fit<threshold then break the iterative loop (default 1E-10)
#' @return \item{group}{ input parameter group }
#' @return \item{scale}{ scaling factor applied to the dataset X}
#' @return \item{Q}{common scores (nrow x ndim)}
#' @return \item{saliences}{weights associated to each table for each dimension}
#' @return \item{explained}{retun total variance explained}
#' @return \item{RV}{RV coefficients between each table (Xk) and compromise table}
#' @return \item{Block}{results associated with each table. You will find block component ...
#'         \itemize{
#'                \item {Qk}{: Block component}
#'                \item {Wk}{: Block component}
#'                \item {Pk}{: Block component}
#'        }}
#'
#' @return \item{call}{: call of the method }
#' @importFrom utils setTxtProgressBar
#' @importFrom DSI datashield.aggregate
#' @export
federateComDim <- function(loginFD, logins, queryvar, querytab, size = NA, H = 2, scale = "none", option = "none", threshold = 1e-10, TOL = 1e-10) {
    group <- dsSwissKnife:::.decode.arg(queryvar)
    ## compute SSCP matrix for each centered data table
    XX <- lapply(group, function(variables) {
        federateSSCP(loginFD, logins, .encode.arg(variables), TOL)
    })
    return (XX)
}
