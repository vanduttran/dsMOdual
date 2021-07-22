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
    
       
    #loginFD <-dsSwissKnife:::.decode.arg(loginFD)
    #logins <- dsSwissKnife:::.decode.arg(logins)
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
    
    print(lapply(XX, dim))
    
    dimensions = list(server1 = c(101,5), server2 = c(101, 5))
    #dimensions = datashield.aggregate(opals, as.symbol('dimDSS(center)'), async=T)
    print(dimensions)
  
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
    pars = list(nTop, nBott)
    y <- c(rep(0,101),rep(1,101))
    labels <- factor(y, labels = c("Smoker-tumor","Normal"))

    #upSignatures = as.data.frame(matrix(rep("NA", ncol(distances)),nTop, ncol(distances)))
    #colnames(upSignatures) = colnames(distances)
  
    #downSignatures = as.data.frame(matrix(rep("NA", ncol(distances)),nBott, ncol(distances)))
    #colnames(downSignatures) = colnames(distances)
  
  
    #consensusUpSignatures = as.data.frame(matrix("NA", nTop, length(unique(labels))))
    #colnames(consensusUpSignatures) = unique(labels)
  
    #consensusDownSignatures = as.data.frame(matrix("NA", nBott, length(unique(labels))))
    #colnames(consensusDownSignatures) = unique(labels)
  
    pars$foldChange = 0
    pars$groupedFoldChange = 0
  
    res = list(distMatrix = distances, 
               upSignatures = "Na", 
               downSignatures =  "Na", 
               groupsAnnotation = labels,
               consensusUpSignatures =  "Na", 
               consensusDownSignatures =  "Na", 
               selectedFeatures = "NA", 
               scudoParams = pars)

    

    addColors <-function(result, object, colors) {
    if (length(object$groupsAnnotation) == 0) {
      igraph::V(result)$color <- rep("#FFFFFF", dim(object$distMatrix)[1])
    } else {
      igraph::V(result)$group <- as.character(object$groupsAnnotation)
      
      if (length(colors) == 0) {
        pal <- grDevices::rainbow(length(levels(object$groupsAnnotation)))
        pal <- stringr::str_extract(pal, "^#[0-9a-fA-F]{6}")
        igraph::V(result)$color <- pal[as.integer(object$groupsAnnotation)]
      } else {
        igraph::V(result)$color <- stringr::str_extract(colors,
                                                        "^#[0-9a-fA-F]{6}")
      }
    }
    
    result
  }


    scudoNetwork = function(object, N, colors = character()) {
    
       			 # input checks
   			 stopifnot(rScudo:::.isSinglePositiveNumber(N),
             		 N <= 1.0,
             		 is.character(colors),
             		 is.vector(colors) )
    
   		 if (length(colors) != 0) {
     			 if (any(is.na(colors))) stop("colors contains NAs")
     			 if (length(colors) != dim(object$distMatrix)[1]) {
       				 stop(paste("length of colors differs from number of samples",
                  		 "in object"))
     		     }
     		 if (any(is.na(stringr::str_match(colors, "^#[0-9a-fA-F]{6,8}$")))) {
       			 stop(paste("colors contains invalid hexadecimal colors (see",
                  		 "documentation for correct format)"))
     			 }
   				 }
    
   	   # get distance matrix and generate igraph object
   	   result <- rScudo:::.makeNetwork(object$distMatrix, N)
       
   	   # add group and color annotation
    
   	   result <- addColors(result, object, colors)

    	 return(result)
 	 }
   
  to_plot = scudoNetwork(res, 0.2)

  return(list(plot = to_plot, res = res))
}




#' @title Try weighted function
#' @description Function for ComDim federated analysis on the virtual cohort combining multiple cohorts
#' Finding common dimensions in multitable data (Xk, k=1...K)
#' @usage federateScudo(X,group,H=2,scale="none",option="none",threshold=1e-10)
#'
#' @param loginFD Login information of the FD server
#' @param logins Login information of the servers containing cohort data
#' @param variables Variables
#' @param XX  :         list of dataframes XX = X %*% t(X)
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
federateTrial <- function(loginFD, logins, queryvar, querytab, nTop=10, nBott=10, labels = "NA", size = NA, TOL = 1e-10) {


    #loginFD <-dsSwissKnife:::.decode.arg(loginFD)
    #logins <- dsSwissKnife:::.decode.arg(logins)
    queryvariables <- dsSwissKnife:::.decode.arg(queryvar)
    querytable     <- dsSwissKnife:::.decode.arg(querytab)



    XX <- lapply(queryvariables, function(variables) {
        federateSSCPweight(loginFD, logins, querytable, .encode.arg(variables), TOL)
    })

    
    return(XX)

}


#' @title Federate SSCP on weighted data
#' @description Function for computing the federated SSCP matrix
#' @param loginFD Login information of the FD server
#' @param logins Login information of data repositories
#' @param querytab Encoded name of a table reference in data repositories
#' @param queryvar Encoded variables from the table reference
#' @param TOL Tolerance of 0
#' @import DSOpal parallel bigmemory
#' @keywords internal
federateSSCPweight <- function(loginFD, logins, querytab, queryvar, TOL = 1e-10) {
    require(DSOpal)

    loginFDdata    <- dsSwissKnife:::.decode.arg(loginFD)
    logindata      <- dsSwissKnife:::.decode.arg(logins)
    querytable     <- dsSwissKnife:::.decode.arg(querytab)
    queryvariables <- dsSwissKnife:::.decode.arg(queryvar)

    opals <- DSI::datashield.login(logins=logindata)
    nNode <- length(opals)

    datashield.assign(opals, "rawData", querytable, variables=queryvariables, async=T)
    #When I print the matrices, Na appear!! But dssSubset doesn't work
    #dssSubset('filtered', 'rawData', row.filter = 'complete.cases(rawData)', datasources = opals)
    datashield.assign(opals, "indexMatrix", as.symbol('dsRank(rawData)'), async=T)
    datashield.assign(opals, "weightMatrix",as.symbol("computeWeights(rawData, indexMatrix)"), async = T)
    datashield.assign(opals, "centeredData", as.symbol('center(weightMatrix)'), async=T)
    datashield.assign(opals, "crossProdSelf", as.symbol('crossProd(centeredData)'), async=T)
    datashield.assign(opals, "tcrossProdSelf", as.symbol('tcrossProd(centeredData, chunk=50)'), async=T)



     ##- received by each from other nodes ----
    invisible(mclapply(names(opals), mc.cores=1, function(opn) {
        logindata.opn <- logindata[logindata$server != opn, , drop=F]
        logindata.opn$user <- logindata.opn$userserver
        logindata.opn$password <- logindata.opn$passwordserver
        opals.loc <- paste0("crossLogin('", .encode.arg(logindata.opn), "')")
        datashield.assign(opals[opn], 'mates', as.symbol(opals.loc), async=F)

        command.opn <- list(paste0("crossAssign(mates, symbol='rawDataMate', value='",
                                   querytab,
                                   "', value.call=F, variables='",
                                   queryvar,
                                   "', async=F)"),
                            paste0("crossAssign(mates, symbol='centeredDataMate', value='",
                                   .encode.arg("center(rawDataMate)"),
                                   "', value.call=T, async=F)")
        )
        for (command in command.opn) {
            cat("Command: ", command, "\n")
            print(datashield.aggregate(opals[opn], as.symbol(command), async=F))
        }

        command.opn <- paste0("crossAggregate(mates, '", .encode.arg('singularProd(centeredDataMate)'), "', async=F)")
        cat("Command: ", command.opn, "\n")
        print(datashield.assign(opals[opn], "singularProdMate", as.symbol(command.opn), async=F))

        command.opn <- paste0("crossAggregate(mates, '",
                              .encode.arg(paste0("as.call(list(as.symbol('pushValue'), dsSSCP:::.encode.arg(crossProdSelf), dsSSCP:::.encode.arg('", opn, "')))")),
                              "', async=F)")
        cat("Command: ", command.opn, "\n")
        print(datashield.assign(opals[opn], "pidMate", as.symbol(command.opn), async=F))
    }))
    datashield.symbols(opals)
 

      #-----

    ## (X_i) * (X_i)': push this symmetric matrix from each node to FD
    #crossProdSelf     <- datashield.aggregate(opals, as.symbol('tcrossProd(centeredData)'), async=T)
    datashield.assign(opals, 'FD', as.symbol(paste0("crossLogin('", loginFD, "')")), async=T)

    # command <- paste0("crossAggregate(FD, '", 
    #                   .encode.arg(paste0("as.call(list(as.symbol('garbageCollect')", "))")), 
    #                   "', async=T)")
    # cat("Command: ", command, "\n")
    # datashield.assign(opals, "GC", as.symbol(command), async=T)

    command <- paste0("dscPush(FD, '",
                      .encode.arg(paste0("as.call(list(as.symbol('pushSymmMatrix'), dsSSCP:::.encode.arg(tcrossProdSelf)", "))")),
                      "', async=T)")
    cat("Command: ", command, "\n")
    crossProdSelfDSC <- datashield.aggregate(opals, as.symbol(command), async=T)


    crossProdSelf <- mclapply(crossProdSelfDSC, mc.cores=min(length(opals), detectCores()), function(dscblocks) {
        return (as.matrix(attach.big.matrix(dscblocks[[1]])))
        ## retrieve the blocks as matrices: on FD
        matblocks <- lapply(dscblocks[[1]], function(dscblock) {
            lapply(dscblock, function(dsc) {
                as.matrix(attach.big.matrix(dsc))
            })
        })
        uptcp <- lapply(matblocks, function(bl) do.call(cbind, bl))
        ## combine the blocks into one matrix
        if (length(uptcp)>1) {
            ## without the first layer of blocks
            no1tcp <- lapply(2:length(uptcp), function(i) {
                cbind(do.call(cbind, lapply(1:(i-1), function(j) {
                    t(matblocks[[j]][[i-j+1]])
                })), uptcp[[i]])
            })
            ## with the first layer of blocks
            tcp <- rbind(uptcp[[1]], do.call(rbind, no1tcp))
        } else {
            tcp <- uptcp[[1]]
        }
        stopifnot(isSymmetric(tcp))
        return (tcp)
    })
    gc(reset=F)

    ## (X_i) * (X_j)' * ((X_j) * (X_j)')[,1]: push this single-column matrix from each node to FD
    #singularProdCross <- datashield.aggregate(opals, as.symbol('tcrossProd(centeredData, singularProdMate)'), async=T)
    datashield.assign(opals, "singularProdCross", as.symbol('tcrossProd(centeredData, singularProdMate)'), async=T)

     command <- paste0("dscPush(FD, '",
                      .encode.arg(paste0("as.call(list(as.symbol('pushSingMatrix'), dsSSCP:::.encode.arg(singularProdCross)", "))")),
                      "', async=T)")
    cat("Command: ", command, "\n")
    singularProdCrossDSC <- datashield.aggregate(opals, as.symbol(command), async=T)

    singularProdCross <- mclapply(singularProdCrossDSC, mc.cores=length(singularProdCrossDSC), function(dscbigmatrix) {
        dscMatList <- lapply(dscbigmatrix[[1]], function(dsc) {
            dscMat <- matrix(as.matrix(attach.big.matrix(dsc)), ncol=1) #TOCHECK: with more than 2 servers
            stopifnot(ncol(dscMat)==1)
            return (dscMat)
        })
        return (dscMatList)
    })
    gc(reset=F)
    #return(singularProdCross)
    ##  (X_i) * (X_j)' * (X_j) * (X_i)'
    #prodDataCross     <- datashield.aggregate(opals, as.symbol('tripleProd(centeredData, crossProdMate)'), async=F)
    ## N.B. save-load increase numeric imprecision!!!
    prodDataCross     <- datashield.aggregate(opals, as.call(list(as.symbol("tripleProd"),
                                                                  as.symbol("centeredData"),
                                                                  .encode.arg(names(opals)))), async=T)

            ## deduced from received info by federation
    crossProductPair <- lapply(1:(nNode-1), function(opi) {
        crossi <- lapply((opi+1):(nNode), function(opj) {
            opni <- names(opals)[opi]
            opnj <- names(opals)[opj]

            a1 <- solveSSCP(XXt=prodDataCross[[opni]][[opnj]],
                            XtX=prodDataCross[[opnj]][[opni]],
                            r=crossProdSelf[[opnj]][, 1, drop=F],
                            Xr=singularProdCross[[opni]][[opnj]],
                            TOL=TOL)
            a2 <- solveSSCP(XXt=prodDataCross[[opnj]][[opni]],
                            XtX=prodDataCross[[opni]][[opnj]],
                            r=crossProdSelf[[opni]][, 1, drop=F],
                            Xr=singularProdCross[[opnj]][[opni]],
                            TOL=TOL)
            cat("Precision on a1 = t(a2):", max(abs(a1 - t(a2))), "\n")
            return (a1)
        })
        names(crossi) <- names(opals)[(opi+1):(nNode)]
        return (crossi)
    })
    names(crossProductPair) <- names(opals)[1:(nNode-1)]
 

     ## SSCP whole matrix
    XXt <- do.call(rbind, lapply(1:nNode, function(opi) {
        upper.opi <- do.call(cbind, as.list(crossProductPair[[names(opals)[opi]]]))
        lower.opi <- do.call(cbind, lapply(setdiff(1:opi, opi), function(opj) {
            t(crossProductPair[[names(opals)[opj]]][[names(opals)[opi]]])
        }))
        return (cbind(lower.opi, crossProdSelf[[opi]], upper.opi))
    }))
    datashield.logout(opals)


    samples = datashield.aggregate(opals, as.symbol('aggRownames(rawData)'), async=T)
    #print(samples)
    #print(".....")


    return (samples)
}












