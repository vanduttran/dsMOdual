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
    

    ScudoResults <- setClass("ScudoResults",
                           slots = list(
                           distMatrix = "matrix",
                           upSignatures = "data.frame",
                           downSignatures = "data.frame",
                           groupsAnnotation = "factor",
                           consensusUpSignatures = "data.frame",
                           consensusDownSignatures = "data.frame",
                           selectedFeatures = "character",
                           scudoParams = "list"))
    
   # loginFD <-dsSwissKnife:::.decode.arg(loginFD)
   # logins <- dsSwissKnife:::.decode.arg(logins)
   group <- dsSwissKnife:::.decode.arg(queryvar)
    ## compute SSCP matrix for each centered data table

   # opals <- datashield.login(logins=logins)
   # nNode <- length(opals)
   # querytable <- unique(logindata$table)
  
   # datashield.assign(opals, 'rawData', querytable,
   #                 variables=queryvar, async=T)

   # datashield.assign(opals, "center", as.symbol('center(rawData)'), async=T)

    XX <- lapply(group, function(variables) {
        federateSSCP(loginFD, logins, .encode.arg(variables), TOL)
    })
    
    dimensions = list(server1 = c(101,5), server2 = c(101, 5))
    #dimensions = datashield.aggregate(opals, as.symbol('dimDSS(center)'), async=T)
    #print(dimensions)
  
    x1_cov = XX[1:dimensions[[1]][1], 1:dimensions[[1]][1]] /(dimensions[[1]][2]-1)
    print(dim(x1_cov))
    x2_cov = XX[(dimensions[[1]][1]+1): nrow(XX), (dimensions[[1]][1]+1): ncol(XX)] /(dimensions[[2]][2]-1) 
    x1x2_cov = XX[1:dimensions[[1]][1], (dimensions[[1]][1]+1): ncol(XX)] / (dimensions[[1]][2]-1) 
    x2x1_cov = XX[(dimensions[[1]][1]+1): nrow(XX),1:dimensions[[1]][1]] / (dimensions[[1]][2]-1) 
  
    block = rbind(cbind(x1_cov, x1y2_cov), cbind(x2x1_cov,x2_cov))
  
    print(dim(block))
    
    correlation <- function(Cxx){
    
    inv_var_x = diag(1/sqrt(diag(Cxx)), ncol(Cxx), ncol(Cxx))
    
    corr = inv_var_x %*% Cxx %*% inv_var_x
    rownames(corr) = rownames(Cxx)
    colnames(corr) = colnames(Cxx)
    return(corr)
    
  }
  
  
  distances = 1 - correlation(block)
  
  #define output
  pars = list(nTop, nBottom)
  
  ScudoResults(distMatrix = distances, 
               upSignatures = NULL, 
               downSignatures = NULL, 
               groupsAnnotation = groups,
               consensusUpSignatures = NULL, 
               consensusDownSignatures = NULL, 
               selectedFeatures = rownames(expressionData), 
               scudoParams = pars)
  
}



:w2q
:wq

