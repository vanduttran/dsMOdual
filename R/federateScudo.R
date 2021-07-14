#' @title Federated Scudo
#' @description Function for ComDim federated analysis on the virtual cohort combining multiple cohorts
#' Finding common dimensions in multitable data (Xk, k=1...K)
#' @usage federateScudo(X,group,H=2,scale="none",option="none",threshold=1e-10)
#'
#' @param loginFD Login information of the FD server
#' @param logins Login information of the servers containing cohort data
#' @param variables Variables
#' @param XX  :	        list of dataframes XX = X %*% t(X)
#' @param nTop how many of the most expressed features we select
#' @param nBott how many of the least expressed features we select
#' @param labels labels assigned to each group of samples
#' @param TOL tolerance
#' @param group :       named list of variables for each table
#' @return XX
#' @importFrom utils setTxtProgressBar
#' @importFrom DSI datashield.aggregate
#' @import rScudo
#' @export
federateScudo <- function(loginFD, logins, queryvar, querytab, nTop=10, nBott=10, labels = "NA", size = NA, TOL = 1e-10) {
    

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
    queryvariables <- dsSwissKnife:::.decode.arg(queryvar)
    querytable     <- dsSwissKnife:::.decode.arg(querytab)

    ## compute SSCP matrix for each centered data table

   # opals <- datashield.login(logins=logins)
   # nNode <- length(opals)
   # querytable <- unique(logindata$table)
  
   # datashield.assign(opals, 'rawData', querytable,
   #                 variables=queryvar, async=T)

   # datashield.assign(opals, "center", as.symbol('center(rawData)'), async=T)
    XX <- lapply(queryvariables, function(variables) {
        federateSSCP(loginFD, logins, querytable, .encode.arg(variables), TOL)
    })
    
    dimensions = list(server1 = c(101,5), server2 = c(101, 5))
    #dimensions = datashield.aggregate(opals, as.symbol('dimDSS(center)'), async=T)
    #print(dimensions)
  
    XXcov = lapply(XX, function(x) {x/(dimensions[[1]][2]-1)})
    print(lapply(XXcov, dim))
    
    correlation <- function(Cxx){
    
    inv_var_x = diag(1/sqrt(diag(Cxx)), ncol(Cxx), ncol(Cxx))
    
    corr = inv_var_x %*% Cxx %*% inv_var_x
    rownames(corr) = rownames(Cxx)
    colnames(corr) = colnames(Cxx)
    return(corr)
    
  }
  
  
  distances =  lapply(XXcov, function(x) {abs(1- correlation(x))})[[1]]

  
  #define output
  pars = list(nTop, nBottom)
  y <- c(rep(0,101),rep(1,101))
  labels <- factor(y, labels = c("Smoker-tumor","Normal"))

  upSignatures = as.data.frame(matrix(rep("NA", ncol(distances)),nTop, ncol(distances)))
  colnames(upSignatures) = colnames(distances)
  
  downSignatures = as.data.frame(matrix(rep("NA", ncol(distances)),nBott, ncol(distances)))
  colnames(downSignatures) = colnames(distances)
  
  
  consensusUpSignatures = as.data.frame(matrix("NA", nTop, length(unique(labels))))
  colnames(consensusUpSignatures) = unique(labels)
  
  consensusDownSignatures = as.data.frame(matrix("NA", nBott, length(unique(labels))))
  colnames(consensusDownSignatures) = unique(labels)
  
  pars$foldChange = 0
  pars$groupedFoldChange = 0
  
  res = ScudoResults(distMatrix = distances, 
               upSignatures = NULL, 
               downSignatures = NULL, 
               groupsAnnotation = labels,
               consensusUpSignatures = NULL, 
               consensusDownSignatures = NULL, 
               selectedFeatures = rownames(expressionData), 
               scudoParams = pars)


  to_plot = scudoNetwork(res, 0.2)
  scudoPlot(to_plot)

  
}





