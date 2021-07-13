#' @title Federated Scudo
#' @description Function for ComDim federated analysis on the virtual cohort combining multiple cohorts
#' Finding common dimensions in multitable data (Xk, k=1...K)
#' @usage federateScudo(X,group,H=2,scale="none",option="none",threshold=1e-10)
#'
#' @param loginFD Login information of the FD server
#' @param logins Login information of the servers containing cohort data
#' @param variables Variables
#' @param XX  :	        list of dataframes XX = X %*% t(X)
#' @param TOL tolerance
#' @param group :       named list of variables for each table
#' @return XX
#' @importFrom utils setTxtProgressBar
#' @importFrom DSI datashield.aggregate
#' @export
federateScudo <- function(loginFD, logins, queryvar, querytab, size = NA, TOL = 1e-10) {
    group <- dsSwissKnife:::.decode.arg(queryvar)
    ## compute SSCP matrix for each centered data table

    XX <- lapply(group, function(variables) {
        federateSSCP(loginFD, logins, .encode.arg(variables), TOL)
    })
    return (XX)
}
