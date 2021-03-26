#' @title Federated ComDim
#' @param logins Login info
#' @param variables Variables
#' @param TOL Tolerance of 0
#' @import DSI
#' @export
ComDimFD <- function(logins, variables, TOL = 1e-10) {
    #opals.cen <- paste0("crossLogin('", dsSSCP:::.encode.arg(logins), "')")
    #datashield.assign(opals[opn], 'mates', as.symbol(opals.cen), async = F)
        
    opals <- datashield.login(logins=logins)
    nNode <- length(opals)
    return (opals)
    # querytable <- unique(logins$table)
    # 
    # datashield.assign(opals, 'rawData', querytable, 
    #                   variables=VAR, async=T)
    # datashield.assign(opals, "centeredData", as.symbol('center(rawData)'), async=T)
    # datashield.assign(opals, "crossProdSelf", as.symbol('crossProd(centeredData)'), async=T)
    # datashield.symbols(opals)
    # ds.summary("centeredData", datasources=opals)
    # ds.summary("crossProdSelf", datasources=opals)
    # 
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
}

