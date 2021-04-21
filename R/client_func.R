#' @title Encode function  arguments
#' @description Serialize to JSON, then encode base64,
#'  then replace '+', '/' and '=' in the result in order to play nicely with the opal sentry.
#'  Used to encode non-scalar function arguments prior to sending to the opal server.
#'  There's a corresponding function in the server package calle .decode_args
#' @param some.object the object to be encoded
#' @return encoded text with offending characters replaced by strings
#' @keywords internal
.encode.arg <- function(some.object){
    encoded <- RCurl::base64Encode(jsonlite::toJSON(some.object, null = 'null'));
    # go fishing for '+', '/' and '=', opal rejects them :
    my.dictionary <- c('\\/' = '-slash-', '\\+' = '-plus-', '\\=' = '-equals-')
    sapply(names(my.dictionary), function(x){
        encoded[1] <<- gsub(x, my.dictionary[x], encoded[1])
    })
    return(paste0(encoded[1],'base64'))
}


#' @title Cross push
#' @description Push matrix data into the federated server
#' @import bigmemory
#' @export
pushValue <- function(value, name) {
    valued <- dsSwissKnife:::.decode.arg(value)
    if (is.list(valued)) valued <- do.call(rbind, valued)
    stopifnot(ncol(valued)!=1)
    dscbigmatrix <- describe(as.big.matrix(valued))
    save(dscbigmatrix, file=paste0("/tmp/", dsSwissKnife:::.decode.arg(name)))
    return (dscbigmatrix)
}


#' @title Find X from XX' and X'X
#' @description Find X from XX' and X'X
#' @param XXt XX'
#' @param XtX X'X
#' @param r A non-null vector of length \code{ncol(t(X)*X)}
#' @param Xr A vector of length \code{nrow(X * t(X))}, equals to the product X %*% r
#' @import parallel
#' @keywords internal
#' @return X
solveSSCP <- function(XXt, XtX, r, Xr, TOL = 1e-10) {
    if (length(r) != ncol(XtX)) {
        stop("r length shoud match ncol(XtX).")
    }
    if (length(Xr) != nrow(XXt)) {
        stop("Xr length shoud match nrow(XXt).")
    }
    if (max(abs(r)) < TOL) {
        stop("Cannot solve with r = 0.")
    }
    
    B1 <- XXt
    B2 <- XtX
    N1 <- nrow(B1)
    N2 <- nrow(B2)
    
    eB1 <- eigen(B1, symmetric=T)
    eB2 <- eigen(B2, symmetric=T)
    vecB1 <- eB1$vectors                    # not unique
    vecB2 <- eB2$vectors                    # not unique
    valB1 <- eB1$values
    valB2 <- eB2$values                     # valB2 == union(valB1, 0)
    vecs <- list("XXt"=vecB1, "XtX"=vecB2)
    vals <- list("XXt"=valB1, "XtX"=valB2)
    if (N2 > N1) {
        tol <- max(abs(valB2[(N1+1):N2]))*10
    } else if (N1 > N2) {
        tol <- max(abs(valB1[(N2+1):N1]))*10
    } else {
        tol <- TOL
    }
    vals <- mclapply(vals, mc.cores=length(vals), function(x) {
        x[abs(x) < tol] <- 0
        return (x)
    })
    eignum <- length(vals[[1]])
    poseignum <- unique(sapply(vals, function(x) {
        max(which(x > 0))
    }))
    cat("Number of strictly positive eigenvalues:", poseignum, "\n")
    stopifnot(length(poseignum)==1)
    ## verify deduced info
    invisible(lapply(1:length(vecs), function(j) {
        vec <- vecs[[j]]
        cat("------", names(vecs)[j], "------\n")
        cat("Determinant:", det(vec), "\n")
        cat("Precision on v' = 1/v:", max(abs(t(vec) - solve(vec))), "\n")
        cat("Precision on Norm_col = 1:", max(abs(apply(vec, 2, function(x) norm(as.matrix(x), "2")) - 1)), "\n")
        cat("Precision on Norm_row = 1:", max(abs(apply(vec, 1, function(x) norm(as.matrix(x), "2")) - 1)), "\n")
        cat("Precision on Orthogonal:", max(sapply(1:(ncol(vec)-1), function(i) {
            max(sum(vec[i,] * vec[i+1,]), sum(vec[, i] * vec[, i+1]))
        })), "\n")
    }))
    
    ## solution S: X * r = vecB1 * E * S * vecB2' * r = Xr
    ## E * S * vecB2' * r = vecB1' * Xr = tmprhs1
    tmprhs1 <- crossprod(vecs[[1]], Xr)
    if (poseignum < N1) cat("Precision on tmprhs1's zero:", max(abs(tmprhs1[(poseignum+1):N1, 1])), "\n")
    ## S * vecB2' * rmX2 = S * lhs1 = 1/E * tmprhs1 = rhs1
    E <- diag(sqrt(vals[[1]][1:poseignum]))
    invE <- diag(1/diag(E))
    rhs1 <- crossprod(t(invE), tmprhs1[1:poseignum, , drop=F])
    lhs1 <- crossprod(vecs[[2]], r)
    signs1 <- rhs1[1:poseignum,]/lhs1[1:poseignum,]
    S <- cbind(diag(signs1), matrix(0, nrow=poseignum, ncol=N2-poseignum)) # S = [signs1 0]
    D <- rbind(crossprod(t(E), S), matrix(0, nrow=N1-poseignum, ncol=N2))  # D = E %*% S
    a1 <- tcrossprod(tcrossprod(vecs[[1]], t(D)), vecs[[2]]) # a = vecs[["A*A'"]] %*% D %*% t(vecs[["A'*A"]])
    
    cat("----------------------\n")
    cat("Precision on XXt = a1*a1':", max(abs(B1 - tcrossprod(a1))), "\n")
    cat("Precision on XtX = a1'*a1:", max(abs(B2 - crossprod(a1))), "\n")
    
    return (a1)
    
    # ## solution S: A' * r = vecB2 * S' * E' * vecB1' * r = Xr
    # ## S' * E' * vecB1' * r = S' * lhs2 = vecB2' * Xr = rhs2
    # rhs2 <- crossprod(vecs[[2]], Xr)
    # if (poseignum < N2) cat("Precision on rhs2's zero:", max(abs(rhs2[(poseignum+1):N2, 1])), "\n")
    # E <- diag(sqrt(vals[[1]][1:poseignum]))
    # 
    # lhs2 <- crossprod(E, crossprod(vecs[[1]][,1:poseignum], r))
    # signs2 <- rhs2[1:poseignum,]/lhs2[1:poseignum,]
    # 
    # ## check signs: signs2 should be identical to signs1
    # cat("Precision on signs double-check:", max(abs(signs1-signs2)), "\n")
    # 
    # S <- cbind(diag(signs2), matrix(0, nrow=poseignum, ncol=N2-poseignum)) # S = [signs1 0]
    # D <- rbind(crossprod(t(E), S), matrix(0, nrow=N1-poseignum, ncol=N2))  # D = E %*% S
    # a2 <- tcrossprod(tcrossprod(vecs[[1]], t(D)), vecs[[2]]) # a = vecs[["A*A'"]] %*% D %*% t(vecs[["A'*A"]])
    # 
    # cat("----------------------\n")
    # cat("Precision on XXt = a2*a2':", max(abs(B1 - tcrossprod(a2))), "\n")
    # cat("Precision on XtX = a2'*a2:", max(abs(B2 - crossprod(a2))), "\n")
    
    # return (a2)
}


#' @title Federated ComDim
#' @param logins Login info
#' @param variables Variables
#' @param TOL Tolerance of 0
#' @import DSOpal parallel bigmemory
#' @export
ComDimFD <- function(loginFD, logins, variables, TOL = 1e-10) {
    require(DSOpal)

    loginFDdata <- dsSwissKnife:::.decode.arg(loginFD)
    logindata <- dsSwissKnife:::.decode.arg(logins)
    vardata <- dsSwissKnife:::.decode.arg(variables)
    
    opals <- DSI::datashield.login(logins=logindata)
    nNode <- length(opals)
    querytable <- unique(logindata$table)

    datashield.assign(opals, "rawData", querytable,
                      variables=vardata, async=T)
    datashield.assign(opals, "centeredData", as.symbol('center(rawData)'), async=T)
    datashield.assign(opals, "crossProdSelf", as.symbol('crossProd(centeredData)'), async=T)
    datashield.assign(opals, "tcrossProdSelf", as.symbol('tcrossProd(centeredData)'), async=T)

    ##- received by node i from other nodes ----
    invisible(mclapply(names(opals), mc.cores=1, function(opn) {
        logindata.opn <- logindata[logindata$server != opn, , drop=F]
        logindata.opn$user <- logindata.opn$userserver
        logindata.opn$password <- logindata.opn$passwordserver
        opals.loc <- paste0("crossLogin('", .encode.arg(logindata.opn), "')")
        datashield.assign(opals[opn], 'mates', as.symbol(opals.loc), async=F)
        
        command.opn <- list(paste0("crossAssign(mates, symbol='rawDataMate', value='", 
                                   .encode.arg(querytable), 
                                   "', value.call=F, variables='",
                                   variables,
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
    
    ##  (X_i) * (X_i)'
    #crossProdSelf     <- datashield.aggregate(opals, as.symbol('tcrossProd(centeredData)'), async=T)
    datashield.assign(opals, 'FD', as.symbol(paste0("crossLogin('", loginFD, "')")), async=T)
    command <- paste0("dscPush(FD, '", 
                      .encode.arg(paste0("as.call(list(as.symbol('pushValue'), dsSSCP:::.encode.arg(crossProdSelf), dsSSCP:::.encode.arg('", names(opals)[1], "')))")), 
                      "', async=T)")
    cat("Command: ", command, "\n")
    crossProdSelfDSC <- datashield.aggregate(opals, as.symbol(command), async=T)
    
    crossProdSelf <- mclapply(crossProdSelfDSC, function(dscbigmatrix) {
        y <- as.matrix(attach.big.matrix(dscbigmatrix[[1]]))
        stopifnot(isSymmetric(y))
        return (y)
    })
        
    ##  (X_i) * (X_j)' * ((X_j) * (X_j)')[,1]
    #singularProdCross <- datashield.aggregate(opals, as.symbol('tcrossProd(centeredData, singularProdMate)'), async=T)
    datashield.assign(opals, "singularProdCross", as.symbol('tcrossProd(centeredData, singularProdMate)'), async=T)
    command <- paste0("dscPush(FD, '", 
                      .encode.arg(paste0("as.call(list(as.symbol('pushValue'), dsSSCP:::.encode.arg(singularProdCross), dsSSCP:::.encode.arg('", names(opals)[1], "')))")), 
                      "', async=T)")
    cat("Command: ", command, "\n")
    singularProdCrossDSC <- datashield.aggregate(opals, as.symbol(command), async=T)
    singularProdCross <- mclapply(singularProdCrossDSC, function(dscbigmatrix) {
        y <- as.matrix(attach.big.matrix(dscbigmatrix[[1]]))
        stopifnot(ncol(y)!=1)
        return (y)
    })
    return(singularProdCross)
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
                            Xr=singularProdCross[[opni]][[opnj]])
            a2 <- solveSSCP(XXt=prodDataCross[[opnj]][[opni]],
                            XtX=prodDataCross[[opni]][[opnj]],
                            r=crossProdSelf[[opni]][, 1, drop=F],
                            Xr=singularProdCross[[opnj]][[opni]])
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
    
    return (XXt)
}

