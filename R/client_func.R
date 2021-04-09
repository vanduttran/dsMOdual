#' @title Federated ComDim
#' @param logins Login info
#' @param variables Variables
#' @param TOL Tolerance of 0
#' @import DSOpal
#' @export
ComDimFD <- function(logins, variables, TOL = 1e-10) {
    require(DSOpal)
    #opals.cen <- paste0("crossLogin('", dsSSCP:::.encode.arg(logins), "')")
    #datashield.assign(opals[opn], 'mates', as.symbol(opals.cen), async = F)
    logindata <- dsSwissKnife:::.decode.arg(logins)
    vardata <- dsSwissKnife:::.decode.arg(variables)
    
    opals <- DSI::datashield.login(logins=logindata)
    nNode <- length(opals)
    querytable <- unique(logindata$table)

    datashield.assign(opals, "rawData", querytable,
                      variables=vardata, async=T)
    datashield.assign(opals, "centeredData", as.symbol('center(rawData)'), async=T)
    datashield.assign(opals, "crossProdSelf", as.symbol('crossProd(centeredData)'), async=T)

    ##- received by node i from other nodes ----
    invisible(mclapply(names(opals), mc.cores=1, function(opn) {
        logindata.opn <- logindata[logindata$server != opn, , drop=F]
        logindata.opn$user <- logindata.opn$userserver
        logindata.opn$password <- logindata.opn$passwordserver
        opals.loc <- paste0("crossLogin('", dsSwissKnifeClient:::.encode.arg(logindata.opn), "')")
        datashield.assign(opals[opn], 'mates', as.symbol(opals.loc), async = F)
        
        command.opn <- list(paste0("crossAssign(mates, symbol='rawDataMate', value='", 
                                   dsSwissKnifeClient:::.encode.arg(querytable), 
                                   "', value.call=F, variables='",
                                   dsSwissKnifeClient:::.encode.arg(VAR),
                                   "', async=F)"),
                            paste0("crossAssign(mates, symbol='centeredDataMate', value='",
                                   dsSwissKnifeClient:::.encode.arg("center(rawDataMate)"),
                                   "', value.call=T, async=F)")
        )
        for (command in command.opn) {
            cat("Command: ", command, "\n")
            print(datashield.aggregate(opals[opn], as.symbol(command), async=F))
        }
        
        command.opn <- paste0("crossAggregate(mates, '", dsSwissKnifeClient:::.encode.arg('singularProd(centeredDataMate)'), "', async=F)")
        cat("Command: ", command.opn, "\n")
        print(datashield.assign(opals[opn], "singularProdMate", as.symbol(command.opn), async=F))
        
        command.opn <- paste0("crossAggregate(mates, '", 
                              dsSwissKnifeClient:::.encode.arg(paste0("as.call(list(as.symbol('pushValue'), dsSSCP:::.encode.arg(crossProdSelf), dsSSCP:::.encode.arg('", opn, "')))")), 
                              "', async=F)")
        cat("Command: ", command.opn, "\n")
        print(datashield.assign(opals[opn], "pidMate", as.symbol(command.opn), async=F))
    }))
    datashield.symbols(opals)
    
    #-----
    
    ##  (X_i) * (X_i)'
    crossProdSelf     <- datashield.aggregate(opals, as.symbol('tcrossProd(centeredData)'), async=T)
    ##  (X_i) * (X_j)' * ((X_j) * (X_j)')[,1]
    singularProdCross <- datashield.aggregate(opals, as.symbol('tcrossProd(centeredData, singularProdMate)'), async=T)
    ##  (X_i) * (X_j)' * (X_j) * (X_i)'
    #prodDataCross     <- datashield.aggregate(opals, as.symbol('tripleProd(centeredData, crossProdMate)'), async=F)
    ## N.B. save-load increase numeric imprecision!!!
    prodDataCross     <- datashield.aggregate(opals, as.call(list(as.symbol("tripleProd"), 
                                                                  as.symbol("centeredData"), 
                                                                  dsSwissKnifeClient:::.encode.arg(names(opals)))), async=T)
    ## deduced from received info by federation
    crossProductPair <- lapply(1:(nNode-1), function(opi) {
        crossi <- lapply((opi+1):(nNode), function(opj) {
            opni <- names(opals)[opi]
            opnj <- names(opals)[opj]

            a1 <- solveSSCP(XXt=prodDataCross[[opni]][[opnj]],
                            XtX=prodDataCross[[opnj]][[opni]],
                            r=crossProdSelf[[opnj]][, 1, drop=F],
                            Xr=singularProdCross[[opni]][[opnj]])
            a2 <- solveSSCP(XXt=prodDataCross[[opnj]][[opni]],
                            XtX=prodDataCross[[opni]][[opnj]],
                            r=crossProdSelf[[opni]][, 1, drop=F],
                            Xr=singularProdCross[[opnj]][[opni]])
            #cat("Precision on A = a1:", max(abs(As[[opni]][[opnj]] - a1)), "\n")
            #cat("Precision on A = a2:", max(abs(As[[opnj]][[opni]] - a2)), "\n")
            cat("Precision on a1 = a2:", max(a1 - a2)), "\n")
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
    
    return (XXt)
}

