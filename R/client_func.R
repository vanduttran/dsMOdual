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
    opals <- DSI::datashield.login(logins=logindata)
    nNode <- length(opals)
    querytable <- unique(logindata$table)
    print(opals)
    print(querytable)
    
    datashield.assign(opals, "rawData", querytable,
                      variables=variables, async=T)
    print(datashield.symbols(opals))
    print(datashield.erros())
    return (NULL)
    #datashield.assign(opals, "centeredData", as.symbol('center(rawData)'), async=T)
    #print(datashield.symbols(opals))
    #datashield.assign(opals, "crossProdSelf", as.symbol('crossProd(centeredData)'), async=T)
    #print(datashield.symbols(opals))
    #datashield.symbols(opals)
    #ds.summary("centeredData", datasources=opals)
    #ds.summary("crossProdSelf", datasources=opals)
    
    # ##- received by node i from other nodes ----
    # invisible(mclapply(names(opals), mc.cores=1, function(opn) {
    #     logindata.opn <- logins[logins$server != opn, , drop=F]
    #     logindata.opn$user <- logindata.opn$userserver
    #     logindata.opn$password <- logindata.opn$passwordserver
    #     opals.loc <- paste0("crossLogin('", dsSwissKnifeClient:::.encode.arg(logindata.opn), "')")
    #     datashield.assign(opals[opn], 'mates', as.symbol(opals.loc), async = F)
    #     
    #     command.opn <- list(paste0("crossAssign(mates, symbol='rawDataMate', value='", 
    #                                dsSwissKnifeClient:::.encode.arg(querytable), 
    #                                "', value.call=F, variables='",
    #                                dsSwissKnifeClient:::.encode.arg(VAR),
    #                                "', async=F)"),
    #                         paste0("crossAssign(mates, symbol='centeredDataMate', value='",
    #                                dsSwissKnifeClient:::.encode.arg("center(rawDataMate)"),
    #                                "', value.call=T, async=F)")
    #     )
    #     for (command in command.opn) {
    #         cat("Command: ", command, "\n")
    #         print(datashield.aggregate(opals[opn], as.symbol(command), async=F))
    #     }
    #     
    #     command.opn <- paste0("crossAggregate(mates, '", dsSwissKnifeClient:::.encode.arg('singularProd(centeredDataMate)'), "', async=F)")
    #     #command.opn <- paste0("crossAggregate(mates, '", dsSwissKnifeClient:::.encode.arg("as.call(list(as.symbol('singularProd'), as.symbol('centeredDataMate')))"), "', async=F)")
    #     cat("Command: ", command.opn, "\n")
    #     print(datashield.assign(opals[opn], "singularProdMate", as.symbol(command.opn), async=F))
    #     
    #     command.opn <- paste0("crossAggregate(mates, '", 
    #                           dsSwissKnifeClient:::.encode.arg(paste0("as.call(list(as.symbol('pushValue'), dsSSCP:::.encode.arg(crossProdSelf), dsSSCP:::.encode.arg('", opn, "')))")), 
    #                           "', async=F)")
    #     cat("Command: ", command.opn, "\n")
    #     print(datashield.assign(opals[opn], "pidMate", as.symbol(command.opn), async=F))
    # }))
    # datashield.symbols(opals)
    ##-----
    
    ##  (X_i) * (X_i)'
    #crossProdSelf     <- datashield.aggregate(opals, as.symbol('tcrossProd(centeredData)'), async=T)
    ##  (X_i) * (X_j)' * ((X_j) * (X_j)')[,1]
    #singularProdCross <- datashield.aggregate(opals, as.symbol('tcrossProd(centeredData, singularProdMate)'), async=T)
    ##  (X_i) * (X_j)' * (X_j) * (X_i)'
    #prodDataCross     <- datashield.aggregate(opals, as.symbol('tripleProd(centeredData, crossProdMate)'), async=F)
    ## N.B. save-load increase numeric imprecision!!!
    # prodDataCross     <- datashield.aggregate(opals, as.call(list(as.symbol("tripleProd"), 
    #                                                               as.symbol("centeredData"), 
    #                                                               dsSwissKnifeClient:::.encode.arg(names(opals)))), async=T)
    #return (crossProdSelf)
}

