#' @title Garbage collection
#' @description Call gc on the federated server
#' @keywords internal
garbageCollect <- function() {
    gc(reset=T)
    return (NULL)
}


#' @title Bigmemory description of a matrix
#' @description Bigmemory description of a matrix
#' @param value Encoded value of a matrix
#' @import bigmemory
#' @return Bigmemory description of the given matrix
#' @export
matrix2DscFD <- function(value) {
    valued <- .decode.arg(value)
    tcp <- do.call(rbind, .decode.arg(valued))
    dscbigmatrix <- describe(as.big.matrix(tcp, backingfile = ""))
    rm(list=c("valued", "tcp"))
    return (dscbigmatrix)
}


#' @title Symmetric matrix reconstruction
#' @description Rebuild a matrix from its partition
#' @param matblocks List of lists of matrix blocks, obtained from .partitionMatrix
#' @param mc.cores Number of cores for parallel computing. Default: 1
#' @return The complete symmetric matrix
#' @keywords internal
.rebuildMatrix <- function(matblocks, mc.cores = 1) {
    uptcp <- lapply(matblocks, function(bl) do.call(cbind, bl))
    .printTime("perform .rebuildMatrix")
    print(lapply(uptcp, range))
    .printTime("wrote uptcp")
    ## combine the blocks into one matrix
    if (length(uptcp)>1) {
        if (length(unique(sapply(uptcp, ncol)))==1) {
            tcp <- do.call(rbind, uptcp)
        } else {
            ## without the first layer of blocks
            no1tcp <- mclapply(2:length(uptcp), mc.cores=mc.cores, function(i) {
                cbind(do.call(cbind, lapply(1:(i-1), function(j) {
                    t(matblocks[[j]][[i-j+1]])
                })), uptcp[[i]])
            })
            ## with the first layer of blocks
            tcp <- rbind(uptcp[[1]], do.call(rbind, no1tcp))
            rm(list=c("no1tcp"))
        }
    } else {
        tcp <- uptcp[[1]]
    }
    print(sum(tcp))
    print(tcp[1:3,1:3])
    print(range(tcp))
    print(max(abs(tcp-t(tcp))))
    .printTime("diff tcp - t(tcp)")
    stopifnot(isSymmetric(tcp))
    rm(list=c("uptcp"))
    return (tcp)
}


#' @title Symmetric matrix reconstruction
#' @description Rebuild a symmetric matrix from its partition bigmemory objects
#' @param dscblocks List of lists of bigmemory objects pointed to matrix blocks
#' @param mc.cores Number of cores for parallel computing. Default: 1
#' @return The complete symmetric matrix
#' @keywords internal
.rebuildMatrixDsc <- function(dscblocks, mc.cores = 1) {
    ## access to matrix blocks 
    matblocks <- mclapply(dscblocks, mc.cores=mc.cores, function(y) {
        lapply(y, function(x) {
            return (as.matrix(attach.big.matrix(x)))
        })
    })
    tcp <- .rebuildMatrix(matblocks, mc.cores=mc.cores)
    return (tcp)
}


#' @title Euclidean distance
#' @description Transform XX' matrix into Euclidean between samples (rows) in X
#' @param XXt An SSCP matrix XX'
#' @import parallel
#' @return Euclidean distance
#' @keywords internal
.toEuclidean <- function(XXt) {
    if (!isSymmetric(XXt) || any(rownames(XXt) != colnames(XXt))) stop('Input XXt (an SSCP matrix) should be symmetric.')
    .printTime("before toEuclidean XXt")
    print(class(XXt))
    print(dim(XXt))
    print(sum(XXt))
    print(XXt[1:3,1:3])
    #lowerTri <- cbind(do.call(cbind, mclapply(1:(ncol(XXt)-1), mc.cores=max(2, min(ncol(XXt)-1, detectCores())), function(i) {
    lowerTri <- cbind(do.call(cbind, lapply(1:(ncol(XXt)-1), function(i) {
        res <- sapply((i+1):ncol(XXt), function(j) {
            ## d(Xi, Xj)^2 = Xi'Xi - 2Xi'Xj + Xj'Xj = XXt[i,i] - 2XXt[i,j] + XXt[j,j]
            t(c(1,-1)) %*% XXt[c(i,j),c(i,j)] %*% c(1, -1)
        })
        return (c(rep(0, i), res))
    })), rep(0, ncol(XXt)))
    distmat <- sqrt(lowerTri + t(lowerTri))
    rownames(distmat) <- rownames(XXt)
    colnames(distmat) <- colnames(XXt)
    
    return (distmat)
}


#' @title Find X from XX' and X'X
#' @description Find X from XX' and X'X
#' @param XXt XX'
#' @param XtX X'X
#' @param r A non-null vector of length \code{ncol(X'X)}
#' @param Xr A vector of length \code{nrow(XX'}, equals to the product Xr
#' @param TOL Tolerance of 0
#' @import parallel
#' @importFrom Matrix rankMatrix
#' @keywords internal
#' @return X
.solveSSCP <- function(XXt, XtX, r, Xr, TOL = 1e-10) {
    if (length(r) != ncol(XtX)) {
        stop("r length shoud match ncol(XtX).")
    }
    if (length(Xr) != nrow(XXt)) {
        stop("Xr length shoud match nrow(XXt).")
    }
    if (max(abs(r)) < TOL) {
        stop("Cannot solve with r = 0.")
    }
    ## scale by some exponential of 10 to deal with high values in input
    absval <- setdiff(c(abs(XXt), abs(XtX)), 0)
    scaling <- 10^ceiling(log10(min(absval))/4+1)
    Xrb <- Xr/(scaling^4)
    rb <- r/(scaling^2)
    B1 <- XXt/(scaling^4)
    B2 <- XtX/(scaling^4)
    N1 <- nrow(B1)
    N2 <- nrow(B2)
    Nmin <- min(N1, N2)
    
    eB1 <- eigen(B1, symmetric=T)
    eB2 <- eigen(B2, symmetric=T)
    vecB1 <- eB1$vectors                    # not unique
    vecB2 <- eB2$vectors                    # not unique
    valB1 <- eB1$values
    valB2 <- eB2$values                     # valB2 == [valB1, 0] or reversely
    vecs <- list("XXt"=vecB1, "XtX"=vecB2)
    vals <- list("XXt"=valB1, "XtX"=valB2)
    
    ## NB: numerically imprecise: poseignum = min(Matrix::rankMatrix(B1), Matrix::rankMatrix(B2))
    ## this number of positive eigenvalues can upper-limitted by the number of variables, yet required as argument of the function .solveSSCP
    ## the following formula is empiric, based on the fact that errors of zero-value are random
    poseignum <- min(Nmin+1, which(vals$XXt[1:Nmin]/(vals$XtX[1:Nmin]+.Machine$double.eps) < 0.99 |
                                       vals$XtX[1:Nmin]/(vals$XXt[1:Nmin]+.Machine$double.eps) < 0.99)) - 1

    vals <- mclapply(vals, mc.cores=length(vals), function(x) {
        x[(poseignum+1):length(x)] <- 0
        return (x)
    })

    # if (N2 > N1) {
    #     tol <- max(abs(valB2[(N1+1):N2]))*10
    # } else if (N1 > N2) {
    #     tol <- max(abs(valB1[(N2+1):N1]))*10
    # } else {
    #     tol <- TOL
    # }
    # vals <- mclapply(vals, mc.cores=length(vals), function(x) {
    #     x[abs(x) < tol] <- 0
    #     return (x)
    # })
    # eignum <- length(vals[[1]])
    # poseignum <- unique(sapply(vals, function(x) {
    #     max(which(x > 0))
    # }))
    # cat("Number of strictly positive eigenvalues:", poseignum, "with tolerance of", tol, "\n")
    # stopifnot(length(poseignum)==1)
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
    tmprhs1 <- crossprod(vecs[[1]], Xrb)
    if (poseignum < N1) cat("Precision on tmprhs1's zero:", max(abs(tmprhs1[(poseignum+1):N1, 1])), "\n")
    ## S * vecB2' * rmX2 = S * lhs1 = 1/E * tmprhs1 = rhs1
    E <- diag(sqrt(vals[[1]][1:poseignum]), ncol=poseignum, nrow=poseignum)
    invE <- diag(1/diag(E), ncol=poseignum, nrow=poseignum)
    rhs1 <- crossprod(t(invE), tmprhs1[1:poseignum, , drop=F])
    lhs1 <- crossprod(vecs[[2]], rb)
    signs1 <- rhs1[1:poseignum,]/lhs1[1:poseignum,]
    S <- cbind(diag(signs1, ncol=poseignum, nrow=poseignum), matrix(0, nrow=poseignum, ncol=N2-poseignum)) # S = [signs1 0]
    D <- rbind(crossprod(t(E), S), matrix(0, nrow=N1-poseignum, ncol=N2))  # D = E %*% S
    a1 <- tcrossprod(tcrossprod(vecs[[1]], t(D)), vecs[[2]]) # a = vecs[["A*A'"]] %*% D %*% t(vecs[["A'*A"]])
    
    cat("----------------------\n")
    cat("Precision on XXt = a1*a1':", max(abs(B1 - tcrossprod(a1))), " / (", quantile(abs(B1)), ")\n")
    cat("Precision on XtX = a1'*a1:", max(abs(B2 - crossprod(a1))),  " / (", quantile(abs(B2)), ")\n")
    
    return (a1*(scaling^2))
    
    # ## solution S: A' * r = vecB2 * S' * E' * vecB1' * r = Xr
    # ## S' * E' * vecB1' * r = S' * lhs2 = vecB2' * Xr = rhs2
    # rhs2 <- crossprod(vecs[[2]], Xrb)
    # if (poseignum < N2) cat("Precision on rhs2's zero:", max(abs(rhs2[(poseignum+1):N2, 1])), "\n")
    # E <- diag(sqrt(vals[[1]][1:poseignum]))
    # 
    # lhs2 <- crossprod(E, crossprod(vecs[[1]][,1:poseignum], rb))
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
    
    # return (a2*(scaling^2))
}


#' @title Federated SSCP
#' @description Function for computing the federated SSCP matrix
#' @param loginFD Login information of the FD server
#' @param logins Login information of data repositories
#' @param funcPreProc Definition of a function for preparation of raw data matrices. 
#' Two arguments are required: conns (list of DSConnection-classes), 
#' symbol (name of the R symbol) (see datashield.assign).
#' @param querytables Vector of names of the R symbols to assign in the Datashield R session on each server in \code{logins}.
#' The assigned R variables will be used as the input raw data.
#' Other assigned R variables in \code{funcPreProc} are ignored.
#' @param ind Index of table in querytables to compute SSCP. Default, ind = 1.
#' @param byColumn A logical value indicating whether the input data is centered by column or row.
#' Default, TRUE, centering by column. Constant variables across samples are removed. 
#' If FALSE, centering and scaling by row. Constant samples across variables are removed.
#' @param TOL Tolerance of 0
#' @import DSOpal parallel bigmemory
#' @keywords internal
.federateSSCP <- function(loginFD, logins, funcPreProc, querytables, ind = 1, byColumn = TRUE, scale = FALSE, chunk = 500, mc.cores = 1, TOL = 1e-10) {
    #require(DSOpal)
    #require(dsBaseClient)
    stopifnot((length(querytables) > 0) & (ind %in% 1:length(querytables)))
    
    loginFDdata    <- .decode.arg(loginFD)
    logindata      <- .decode.arg(logins)
    opals <- .login(logins=logindata) #datashield.login(logins=logindata)
    nNode <- length(opals)
    .printTime(".federateSSCP login-ed")
    tryCatch({
        ## take a snapshot of the current session
        safe.objs <- .ls.all()
        safe.objs[['.GlobalEnv']] <- setdiff(safe.objs[['.GlobalEnv']], '.Random.seed')  # leave alone .Random.seed for sample()
        ## lock everything so no objects can be changed
        .lock.unlock(safe.objs, lockBinding)
        
        ## apply funcPreProc for preparation of querytables on opals
        ## TODO: control hacking!
        ## TODO: control identical colnames!
        funcPreProc(conns=opals, symbol=querytables)
        
        ## unlock back everything
        .lock.unlock(safe.objs, unlockBinding)
        ## get rid of any sneaky objects that might have been created in the filters as side effects
        .cleanup(safe.objs)
    }, error=function(e) {
        print(paste0("DATA MAKING PROCESS: ", e))
        return (paste0("DATA MAKING PROCESS: ", e, ' --- ', datashield.symbols(opals), ' --- ', datashield.errors(), ' --- ', datashield.logout(opals)))
    })
    .printTime(".federateSSCP data processed")
    if (nNode==1) {
        tryCatch({
            datashield.assign(opals, "centeredData", as.symbol(paste0("center(", querytables[ind], ", subset=NULL, byColumn=", byColumn, ", scale=", scale, ")")), async=T)
            datashield.assign(opals, "tcrossProdSelf", as.symbol(paste0('tcrossProd(x=centeredData, y=NULL, chunk=', chunk, ')')), async=T)
            .printTime(".federateSSCP intermediate data computed")
            samplenames <- datashield.aggregate(opals, as.symbol("rowNames(centeredData)"), async=T)
            datashield.assign(opals, 'FD', as.symbol(paste0("crossLogin('", loginFD, "')")), async=T)
            tryCatch({
                cat("Command: pushToDscFD(FD, 'tcrossProdSelf')", "\n")
                tcrossProdSelfDSC <- datashield.aggregate(opals, as.symbol("pushToDscFD(FD, 'tcrossProdSelf')"), async=T)
                tcrossProdSelf <- .rebuildMatrixDsc(tcrossProdSelfDSC[[1]])
            },
            error=function(e) print(paste0("FD PROCESS SINGLE: ", e, ' --- ', datashield.symbols(opals), ' --- ', datashield.errors())),
            finally=datashield.assign(opals, 'crossEnd', as.symbol("crossLogout(FD)"), async=T))
            .printTime(".federateSSCP XX' communicated to FD")
            XXt <- tcrossProdSelf
            rownames(XXt) <- colnames(XXt) <- unlist(samplenames, use.names=F)
            gc(reset=F)
        },
        error=function(e) print(paste0("XX' PROCESS SINGLE: ", e, ' --- ', datashield.symbols(opals), ' --- ', datashield.errors())),
        finally=datashield.logout(opals))
    } else {
        tryCatch({
            datashield.assign(opals, "centeredData", as.symbol(paste0("center(", querytables[ind], ", subset=NULL, byColumn=", byColumn, ", scale=", scale, ")")), async=T)
            datashield.assign(opals,  "crossProdSelf",  as.symbol(paste0('crossProd(x=centeredData, y=NULL, chunk=', chunk, ')')), async=T)
            datashield.assign(opals, "tcrossProdSelf", as.symbol(paste0('tcrossProd(x=centeredData, y=NULL, chunk=', chunk, ')')), async=T)
            .printTime(".federateSSCP intermediate data computed")
            samplenames <- datashield.aggregate(opals, as.symbol("rowNames(centeredData)"), async=T)
            
            ##- received by each from other nodes ----
            prodDataCross <- lapply(names(opals), function(opn) {
                ind.opn <- which(logindata$server == opn)
                logindata.opn <- logindata[-ind.opn, , drop=F]
                logindata.opn$user <- logindata.opn$userserver
                logindata.opn$password <- logindata.opn$passwordserver
                ## from each node opn, log in other nodes (mates)
                opals.loc <- paste0("crossLogin('", .encode.arg(logindata.opn), "')")
                datashield.assign(opals[opn], 'mates', as.symbol(opals.loc), async=F)
                tryCatch({
                    ## prepare raw data matrices on mates of opn
                    command.opn <- list(as.symbol("crossAssignFunc"),
                                        as.symbol("mates"),
                                        .encode.arg(funcPreProc, serialize.it=T),
                                        .encode.arg(querytables))
                    cat("Command: crossAssignFunc(mates, funcPreProc, ...", "\n")
                    invisible(datashield.aggregate(opals[opn], as.call(command.opn), async=F))
                    
                    ## center raw data on mates of opn
                    command.opn <- paste0("crossAssign(mates, symbol='centeredDataMate', value='",
                                          .encode.arg(paste0("center(", querytables[ind], ", subset=NULL, byColumn=", byColumn, ", scale=", scale, ")")),
                                          "', value.call=T, async=F)")
                    cat("Command: crossAssign(mates, centeredDataMate, ...", "\n")
                    invisible(datashield.aggregate(opals[opn], as.symbol(command.opn), async=F))
                    
                    ## create singularProd from mates of opn
                    command.opn <- paste0("crossAggregatePrimal(mates, '", .encode.arg('singularProd(centeredDataMate)'), "', async=T)")
                    cat("Command: crossAggregatePrimal(mates, singularProd(centeredDataMate), ...", "\n")
                    datashield.assign(opals[opn], "singularProdMate", as.symbol(command.opn), async=T)
                       
                    ## send X'X from opn to mates of opn
                    command.opn <- list(as.symbol("pushToDscMate"),
                                        as.symbol("mates"),
                                        'crossProdSelf',
                                        opn,
                                        async=T)
                    cat("Command: pushToDscMate(mates, crossProdSelf...", "\n")
                    invisible(datashield.aggregate(opals[opn], as.call(command.opn), async=F))
                    .printTime(paste0(".federateSSCP pairwise X'X communicated: ", opn))
                    
                    ## send (X_i) * (X_j)' * (X_j) * (X_i)' to FD
                    command.opn <- paste0("crossAssign(mates, symbol='prodDataCross_", opn, "', value='",
                                          .encode.arg(paste0("tripleProdChunk(centeredDataMate, mate='", .encode.arg(opn), "', chunk=", chunk, ", mc.cores=", mc.cores, ")")),
                                          "', value.call=T, async=F)")
                    cat("Command: crossAssign(mates, prodDataCross, tripleProdChunk(...", "\n")
                    invisible(datashield.aggregate(opals[opn], as.symbol(command.opn), async=F))
                    
                    # login to FD from mates
                    command.opn <- paste0("crossAssign(mates, symbol='FDmate', value='",
                                          .encode.arg(paste0("crossLogin('", loginFD, "')")),
                                          "', value.call=T, async=F)")
                    cat("Command: crossAssign(mates, FDmate, crossLogin(...", "\n")
                    invisible(datashield.aggregate(opals[opn], as.symbol(command.opn), async=F))
                    tryCatch({
                        # push to FD from mates
                        command.opn <- paste0("crossAggregateDual(mates, '", 
                                              .encode.arg(paste0("pushToDscFD(FDmate, 'prodDataCross_", opn, "')")), 
                                              "', async=T)")
                        cat("Command: crossAggregateDual(mates, pushToDscFD(FDmate, prodDataCross_...", "\n")
                        prodDataCrossDSC <- datashield.aggregate(opals[opn], as.symbol(command.opn), async=T)
                        prodDataCross.opn <- mclapply(prodDataCrossDSC[[opn]], mc.cores=mc.cores, function(dscblocks) {
                            return (.rebuildMatrixDsc(dscblocks[[opn]]))
                        })
                        .printTime(paste0(".federateSSCP XY'YX' tripleProd communicated to FD: ", opn))
                    }, 
                    error=function(e) print(paste0("MATES-FD PROCESS MULTIPLE: ", e, ' --- ', datashield.symbols(opals[opn]), ' --- ', datashield.errors())),
                    finally=datashield.aggregate(opals[opn], 
                                                 as.symbol(paste0("crossAssign(mates, symbol='crossFDEnd', value='",
                                                                  .encode.arg("crossLogout(FDmate)"),
                                                                  "', value.call=T, async=F)")), 
                                                 async=T))
                    return (prodDataCross.opn)
                }, 
                error=function(e) print(paste0("CROSS PROCESS: ", e, ' --- ', datashield.symbols(opals[opn]), ' --- ', datashield.errors())),
                finally=datashield.assign(opals[opn], 'crossEnd', as.symbol("crossLogout(mates)"), async=T))
            })
            names(prodDataCross) <- names(opals)
            #-----
            
            datashield.assign(opals, 'FD', as.symbol(paste0("crossLogin('", loginFD, "')")), async=T)
            tryCatch({
                # command <- paste0("crossAggregate(FD, '", 
                #                   .encode.arg(paste0("as.call(list(as.symbol('garbageCollect')", "))")), 
                #                   "', async=T)")
                # cat("Command: ", command, "\n")
                # datashield.assign(opals, "GC", as.symbol(command), async=T)
                
                ## (X_i) * (X_i)': push this symmetric matrix from each node to FD
                cat("Command: pushToDscFD(FD, tcrossProdSelf)", "\n")
                tcrossProdSelfDSC <- datashield.aggregate(opals, as.symbol("pushToDscFD(FD, 'tcrossProdSelf')"), async=T)
                tcrossProdSelf <- mclapply(tcrossProdSelfDSC, mc.cores=mc.cores, function(dscblocks) {
                    return (.rebuildMatrixDsc(dscblocks))
                })
                .printTime(".federateSSCP XX' communicated to FD")
                
                ## (X_i) * (X_j)' * ((X_j) * (X_j)')[,1]: push this single-column matrix from each node to FD
                datashield.assign(opals, "singularProdCross", as.symbol('tcrossProd(centeredData, singularProdMate)'), async=T)
                command <- list(as.symbol("pushToDscFD"),
                                as.symbol("FD"),
                                'singularProdCross',
                                async=T)
                cat("Command: pushToDscFD(FD, 'singularProdCross')", "\n")
                singularProdCrossDSC <- datashield.aggregate(opals, as.call(command), async=T)
                singularProdCross <- mclapply(singularProdCrossDSC, mc.cores=mc.cores, function(dscbigmatrix) {
                    dscMatList <- lapply(dscbigmatrix, function(dsc) {
                        dscMat <- do.call(rbind, lapply(dsc, function(dsci) {
                            return (matrix(as.matrix(attach.big.matrix(dsci[[1]])), ncol=1))
                        }))
                        stopifnot(ncol(dscMat)==1)
                        return (dscMat)
                    })
                    return (dscMatList)
                })
                .printTime(".federateSSCP Ar communicated to FD")
                gc(reset=F)
            },
            error=function(e) print(paste0("FD PROCESS MULTIPLE: ", e, ' --- ', datashield.symbols(opals), ' --- ', datashield.errors())),
            finally=datashield.assign(opals, 'crossEnd', as.symbol("crossLogout(FD)"), async=T))
        },
        error=function(e) print(paste0("XX' PROCESS MULTIPLE: ", e, ' --- ', datashield.symbols(opals), ' --- ', datashield.errors())),
        finally=datashield.logout(opals))
        
        ## deduced from received info by federation: (X_i) * (X_j)'
        crossProductPair <- mclapply(1:(nNode-1), mc.cores=mc.cores, function(opi) {
            crossi <- lapply((opi+1):(nNode), function(opj) {
                opni <- names(opals)[opi]
                opnj <- names(opals)[opj]
                a1 <- .solveSSCP(XXt=prodDataCross[[opnj]][[opni]],
                                 XtX=prodDataCross[[opni]][[opnj]],
                                 r=tcrossProdSelf[[opnj]][, 1, drop=F],
                                 Xr=singularProdCross[[opni]][[opnj]],
                                 TOL=TOL)
                a2 <- .solveSSCP(XXt=prodDataCross[[opni]][[opnj]],
                                 XtX=prodDataCross[[opnj]][[opni]],
                                 r=tcrossProdSelf[[opni]][, 1, drop=F],
                                 Xr=singularProdCross[[opnj]][[opni]],
                                 TOL=TOL)
                cat("Precision on a1 = t(a2):", max(abs(a1 - t(a2))),  " / (", quantile(abs(a1)), ")\n")
                return (a1)
            })
            names(crossi) <- names(opals)[(opi+1):(nNode)]
            return (crossi)
        })
        names(crossProductPair) <- names(opals)[1:(nNode-1)]
        .printTime(".federateSSCP XY' computed")
        
        ## SSCP whole matrix
        XXt <- do.call(rbind, mclapply(1:nNode, mc.cores=mc.cores, function(opi) {
            upper.opi <- do.call(cbind, as.list(crossProductPair[[names(opals)[opi]]]))
            lower.opi <- do.call(cbind, lapply(setdiff(1:opi, opi), function(opj) {
                t(crossProductPair[[names(opals)[opj]]][[names(opals)[opi]]])
            }))
            return (cbind(lower.opi, tcrossProdSelf[[opi]], upper.opi))
        }))
        rownames(XXt) <- colnames(XXt) <- unlist(samplenames[names(opals)], use.names=F)
        .printTime(".federateSSCP whole XX' computed")
        gc(reset=F)
    }
    
    return (XXt)
}


#' @title Federated ComDim
#' @description Function for ComDim federated analysis on the virtual cohort combining multiple cohorts
#' Finding common dimensions in multitable data (Xk, k=1...K)
#' @usage federateComDim(loginFD, logins, func, symbol, ncomp = 2, scale = "none", option = "uniform", chunk = 500, threshold = 1e-10)
#'
#' @param loginFD Login information of the FD server
#' @param logins Login information of data repositories
#' @param func Encoded definition of a function for preparation of raw data matrices. 
#' Two arguments are required: conns (list of DSConnection-classes), 
#' symbol (name of the R symbol) (see datashield.assign).
#' @param symbol Encoded vector of names of the R symbols to assign in the Datashield R session on each server in \code{logins}.
#' The assigned R variables will be used as the input raw data.
#' Other assigned R variables in \code{func} are ignored.
#' @param H Number of common dimensions
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
#' @import DSI
#' @export
federateComDim <- function(loginFD, logins, func, symbol, ncomp = 2, scale = "none", option = "uniform", chunk = 500, mc.cores = 1, threshold = 1e-10) {
    require(DSOpal)
    .printTime("federateComDim started")
    TOL <- 1e-10
    funcPreProc <- .decode.arg(func)
    querytables <- .decode.arg(symbol)
    ntab <- length(querytables)
    
    ## compute SSCP matrix for each centered data table
    XX_query <- lapply(1:ntab, function(i) {
        .federateSSCP(loginFD=loginFD, logins=logins, funcPreProc=funcPreProc, querytables=querytables, ind=i, byColumn=TRUE, chunk=chunk, mc.cores=mc.cores, TOL=TOL)
    })
    names(XX_query) <- querytables
    XX <- XX_query
    
    ## set up the centered data table on every node
    loginFDdata <- .decode.arg(loginFD)
    logindata <- .decode.arg(logins)
    opals <- .login(logins=logindata) #datashield.login(logins=logindata)
    nNode <- length(opals)
    
    tryCatch({
        ## take a snapshot of the current session
        safe.objs <- .ls.all()
        safe.objs[['.GlobalEnv']] <- setdiff(safe.objs[['.GlobalEnv']], '.Random.seed')  # leave alone .Random.seed for sample()
        ## lock everything so no objects can be changed
        .lock.unlock(safe.objs, lockBinding)
        
        ## apply funcPreProc for preparation of querytables on opals
        ## TODO: control hacking!
        ## TODO: control identical colnames!
        funcPreProc(conns=opals, symbol=querytables)
        
        ## unlock back everything
        .lock.unlock(safe.objs, unlockBinding)
        ## get rid of any sneaky objects that might have been created in the filters as side effects
        .cleanup(safe.objs)
    }, error=function(e) {
        print(paste0("DATA MAKING PROCESS: ", e))
        return (paste0("DATA MAKING PROCESS: ", e, ' --- ', datashield.symbols(opals), ' --- ', datashield.errors(), ' --- ', datashield.logout(opals)))
    })
    
    ## take variables (colnames)
    queryvariables <- lapply(querytables, function(querytable) {
        DSI::datashield.aggregate(opals[1], as.symbol(paste0('colNames(', querytable, ')')), async=F)[[1]]
    })
    names(queryvariables) <- querytables
    
    ## centered cbind-ed centered data matrix
    DSI::datashield.assign(opals, "centeredAllData",
                           as.symbol(paste0('center(list(', paste(querytables, collapse=','), '), byColumn=TRUE, na.rm=FALSE)')),
                           async=T)
    
    ## compute the total variance from a XX' matrix
    inertie <- function(tab) {
        return (sum(diag(tab)))    #Froebenius norm
    }
    ## compute the RV between WX and WY: similarity between two matrices
    coefficientRV <- function(WX, WY) {
        rv <- inertie(WX %*% WY)/(sqrt(inertie(WX %*% WX) * inertie(WY %*% WY)))
        return(rv)
    }
    ## normalize a vector to a unit vector
    normv <- function(x) {
        normx <- norm(x, '2')
        if (normx==0) normx <- 1
        return (x/normx)
    }
    
    ##- 0. Preliminary tests ----
    if (any(sapply(XX, is.na)))
        stop("No NA values are allowed")
    nind <- unique(unlist(apply(sapply(XX, dim), 1, unique))) ## number of samples
    if (length(nind) > 1)
        stop("XX elements should be symmetric of the same dimension")
    samples <- sapply(XX, function(x) union(rownames(x), colnames(x)))
    if (is.list(samples) && max(lengths(samples))==0) {
        XX <- lapply(XX, function(x) {
            rownames(x) <- colnames(x) <- paste("Ind", 1:nind, sep='.')
            return (x)
        })
        samples <- sapply(XX, function(x) union(rownames(x), colnames(x)))
    }
    if (is.list(samples) || is.list(apply(samples, 1, unique)))
        stop("XX elements should have the same rownames and colnames")
    
    ## TOREVIEW
    # if (is.character(scale)) {
    #     if (!scale %in% c("none","sd"))
    #         stop("scale must be either none or sd")
    # }
    # else {
    #     if (!is.numeric(scale) | length(scale)!=ncol(X))
    #         stop("Non convenient scaling parameter")
    # }
    # if (!option %in% c("none","uniform"))
    #     stop("option must be either none or uniform")
    ##-----
    
    ##- 1. Output preparation ----
    nvar <- lengths(queryvariables)
    if (ncomp > 2 && min(sapply(XX, rankMatrix)) <= ncomp) {
        print(paste0("Security issue: maximum ", min(sapply(XX, rankMatrix)) - 1, " components could be inquired. ncomp will be set to 2."))
        ncomp <- 2
    }
    if (ncomp < 1) {
        print("ncomp should be at least 1. ncomp will be set to 2.")
        ncomp <- 2
    }
    compnames <- paste("Dim.", 1:ncomp, sep="")
    
    contrib           <- matrix(0, ntab, ncomp)
    dimnames(contrib) <- list(querytables, compnames)
    
    saliences <- LAMBDA <- NNLAMBDA <- matrix(1, ntab, ncomp) # Specific weights for each dataset and each dimension
    dimnames(saliences) <- dimnames(LAMBDA) <- dimnames(NNLAMBDA) <- list(querytables, compnames)
    
    T <- matrix(0, nrow=nind, ncol=ncomp)          # Global components
    C <- matrix(0, nrow=nind, ncol=ncomp)          # Unnormed global components
    dimnames(Q) <- dimnames(C) <- list(rownames(XX[[1]]), compnames)
    
    T.b <- array(0, dim=c(nind, ncomp, ntab))      # Block components
    dimnames(Q.b) <- list(rownames(XX[[1]]), compnames, querytables)
    
    cor.g.b <- array(0, dim=c(ncomp, ncomp, ntab)) # Correlations between global components and their respective block components
    dimnames(cor.g.b) <- list(compnames, compnames, querytables)
    
    W.b      <- vector("list", length=ntab) # Weights for the block components
    blockcor <- vector("list", length=ntab) # Correlation between the original variables of each block and its block components (loadings)
    for (k in 1:ntab) {
        W.b[[k]] <- blockcor[[k]] <- matrix(0, nrow=nvar[k], ncol=ncomp)
        dimnames(W.b[[k]]) <- dimnames(blockcor[[k]]) <- list(queryvariables[[k]], compnames)
    }
    
    Px           <- matrix(0, nrow=sum(nvar), ncol=ncomp)  # Loadings for the global components
    W            <- Px                                     # Weights for the global components
    Wm           <- Px                                     # Modified Weights to take account of deflation
    dimnames(Px) <- dimnames(W) <- dimnames(Wm) <- list(do.call(c, queryvariables), compnames)
    
    IT.X         <- vector("numeric", length=ntab+1)
    explained.X  <- matrix(0, nrow=ntab+1, ncol=ncomp)     # Percentage of explained inertia of each Xb block and global
    dimnames(explained.X) <- list(c(querytables, 'Global'), compnames)
    cumexplained <- matrix(0, nrow=ncomp, ncol=2)
    dimnames(cumexplained) <- list(compnames, c("%explX", "cum%explX"))

    tb           <- matrix(0, nrow=nind, ncol=ntab)        # Concatenated Xb block components at current dimension comp
    components   <- vector("numeric", length=2)
    
    Xscale       <- NULL
    Block        <- NULL
    res          <- NULL
    ##-----
    
    ##- 2. Required parameters and data preparation ----
    # Xscale$mean   <- apply(X, 2, mean)
    # X             <- scale(X,center=Xscale$mean, scale=FALSE)   # Default centering
    # 
    # if (scale=="none") {
    #     Xscale$scale <-rep(1,times=ncol(X))
    # }   else {
    #     if (scale=="sd") {
    #         sd.tab    <- apply(X, 2, function (x) {return(sqrt(sum(x^2)/length(x)))})   # sd based on biased variance
    #         temp      <- sd.tab < 1e-14
    #         if (any(temp)) {
    #             warning("Variables with null variance not standardized.")
    #             sd.tab[temp] <- 1
    #         }
    #         X         <- sweep(X, 2, sd.tab, "/")
    #         Xscale$scale <-sd.tab
    #     }   else {     # Specific scaling depending on blocks defined as a vector with a scaling parameter for each variable
    #         X         <- sweep(X, 2, scale, "/")
    #         Xscale$scale <-scale
    #     }
    # }
    
    #  Pre-processing: block weighting to set each block inertia equal to 1
    if (option=="uniform") {
        inertia <- sapply(XX, inertie)
        XX <- lapply(1:ntab, function(k) {
            XX[[k]]/inertia[k]
        })
        inertia0.sqrt <- sqrt(inertia)
    }
    XX0 <- XX
    IT.X[1:ntab] <- sapply(XX, inertie) # Inertia of each block of variable
    IT.X[ntab+1] <- sum(IT.X[1:ntab])
    
    # Computation of the cross-product matrices among individuals (or association matrices)
    Trace <- IT.X[1:ntab]
    Itot  <- 0
    ##-----
    
    ##- 3. computation of Q and LAMBDA for the various dimensions ----
    for (comp in 1:ncomp)  { # Iterative computation of the various components
        critt     <- 0
        deltacrit <- 1
        
        ## 3.1 Interatively compute the comp-th common component
        while (deltacrit > threshold) {
            P <- Reduce("+", lapply(1:ntab, function(k) LAMBDA[k, comp]*XX[[k]])) # Weighted sum of XX'
            reseig    <- eigen(P)
            q         <- reseig$vectors[, 1]
            Q[, comp] <- q
            optimalcrit[comp] <- reseig$values[1]
            LAMBDA[, comp]  <- sapply(1:ntab, function(k) {t(q) %*% XX[[k]] %*% q})
            LAMBDA[, comp]  <- normv(LAMBDA[, comp])
            criterion   <- reseig$values[1]
            deltacrit   <- criterion - critt
            critt       <- criterion
        }
        
        ## 3.2 Storage of the results associated with dimension comp
        for (k in 1:ntab) {
            #W.b[[k]][, comp] <- t(X[[k]]) %*% Q[, comp]
            T.b[, comp, k] <- XX[[k]] %*% q
        }
        
        LAMBDA[, comp]   <- sapply(1:ntab, function(k) {t(T[,comp]) %*% T.b[, comp, k]})
        NNLAMBDA[, comp] <- LAMBDA[, comp]          # Non normalized specific weights
        LAMBDA[, comp]   <- normv(LAMBDA[, comp])
        
        ## 3.3 Deflation
        X.exp <- lapply(XX, function(xx) T[, comp] %*% xx %*% t(T[, comp]))
        X0.exp <- lapply(XX0, function(xx) T[, comp] %*% xx %*% t(T[, comp]))
        explained.X[1:ntab, comp] <- sapply(X0.exp, function(x) {sum(x^2)})
        explained.X[ntab+1, comp] <- sum(explained.X[1:ntab, comp])
        proj <- diag(1, nind) - tcrossprod(T[, comp])
        XX <- lapply(XX, function(xx) proj %*% xx %*% t(proj)) # Deflation of XX
    }
    ##-----
    
    ##- 4. loadings ----
    # number of samples on each node
    tryCatch({
        size <- sapply(datashield.aggregate(opals, as.symbol('dsDim(centeredAllData)')), function(x) x[1])
    }, error=function(e) {
        print(paste0("INDIVIDUAL DIMENSION COMPUTATION PROCESS: ", e)); 
        return (paste0("INDIVIDUAL DIMENSION COMPUTATION PROCESS: ", e, ' --- ', datashield.symbols(opals), ' --- ', datashield.errors(), ' --- ', datashield.logout(opals)))
    })
    size <- c(0, size)
    func <- function(x, y) {x %*% y}
    Qlist <- setNames(lapply(2:length(size), function(i) {
        Qi <- Q[(cumsum(size)[i-1]+1):cumsum(size)[i], , drop=F]
        ## As Q is orthonormal, Qi == Qi.iter
        # Qi.iter <- sapply(1:H, function(dimension) {
        #   projs <- lapply(setdiff(1:dimension, dimension), function(dimprev) {
        #     return (diag(1, size[i]) - tcrossprod(Qi[,dimprev]))
        #   })
        #   projs <- c(Id=list(diag(1, size[i])), projs)
        #   return (crossprod(Reduce(func, projs), Qi[,dimension,drop=F]))
        # })
        # return (Qi.iter)
    }), names(opals))
    tryCatch({
        Wbk <- Reduce('+', unlist(mclapply(names(opals), mc.cores=1, function(opn) {
            expr <- list(as.symbol("loadings"),
                         as.symbol("centeredAllData"),
                         .encode.arg(Qlist[[opn]]),
                         "prod")
            loadings <- datashield.aggregate(opals[opn], as.call(expr))
            return (loadings)
        }), recursive = F))
    }, error=function(e) {
        print(paste0("LOADING COMPUTATION PROCESS: ", e, ' --- ', datashield.symbols(opals), ' --- ', datashield.errors()))
    }, 
    finally=datashield.logout(opals))
    colnames(Wbk) <- compnames
    csnvar <- cumsum(nvar)
    W.b <- mclapply(1:ntab, mc.cores=ntab, function(k) {
        if (option=="uniform") return (Wbk[ifelse(k==1, 1, csnvar[k-1]+1):csnvar[k], , drop=F]/inertia0.sqrt[k])
        return (Wbk[ifelse(k==1, 1, csnvar[k-1]+1):csnvar[k], , drop=F])
    })
    
    Px <- do.call(rbind, W.b)
    W <- do.call(rbind, lapply(1:ntab, function(k) tcrossprod(W.b[[k]], diag(LAMBDA[k,]))))
    W <- do.call(cbind, lapply(1:ncomp, function(comp) normv(W[, comp])))
    colnames(W) <- compnames
    
    contrib <- t(t(NNLAMBDA)/colSums(NNLAMBDA))
    
    Wm <- W %*% solve(crossprod(Px, W), tol=1e-150)      # Weights that take into account the deflation procedure
    
    # Unnormed global components
    if (ncomp==1) {
        LambdaMoyen <- apply(NNLAMBDA^2, 2, sum)
        C <- Q * LambdaMoyen
    }
    else {
        LambdaMoyen <- apply(NNLAMBDA^2, 2, sum)
        C <- Q %*% sqrt(diag(LambdaMoyen))
    }
    
    #globalcor <- cor(X00, C)
    
    for (k in 1:ntab) {
        cor.g.b[, , k] <- cor(Q, Q.b[, , k])
        #blockcor[[k]] <- cor(X0[[k]], Q.b[, 1:ncomp, k])
        #if (is.null(rownames(blockcor[[k]]))) rownames(blockcor[[k]]) <- names(group[k])
    }
    
    ## 4.1 Preparation of the results Global
    res$components          <- c(ncomp=ncomp)
    res$optimalcrit         <- optimalcrit[1:ncomp]
    res$saliences           <- round(LAMBDA[, 1:ncomp, drop=FALSE]^2, 2)
    res$Q                   <- Q[, 1:ncomp, drop=FALSE]     # Storage of the normed global components associated with X
    res$C                   <- C[, 1:ncomp, drop=FALSE]     # Storage of the unnormed global components associated with X
    res$explained.X         <- round(100*explained.X[1:ntab, 1:ncomp], 2)
    res$cumexplained        <- round(100*cumexplained[1:ncomp,], 2)
    res$contrib             <- round(100*contrib[1:ntab, 1:ncomp], 2)
    #res$globalcor           <- globalcor[,1:ncomp]
    res$cor.g.b             <- cor.g.b#[1:ncomp, 1:ncomp, ]

    ## 4.2 Preparation of the results Block
    Block$Q.b         <-  Q.b[,1:ncomp,]
    Block$blockcor    <-  blockcor
    res$Block         <-  Block                         # Results for each block

    ##- 5. Return res ----
    res$Xscale <- Xscale
    res$call   <- match.call()
    class(res) <- c("ComDim", "list")
    ##-----
    
    return(Res)
}


#' @title Federated ComDim deprecated
#' @description Function for ComDim federated analysis on the virtual cohort combining multiple cohorts
#' Finding common dimensions in multitable data (Xk, k=1...K)
#' @usage federateComDim(loginFD, logins, func, symbol, H = 2, scale = "none", option = "uniform", threshold = 1e-10, TOL = 1e-10)
#'
#' @param loginFD Login information of the FD server
#' @param logins Login information of data repositories
#' @param func Encoded definition of a function for preparation of raw data matrices. 
#' Two arguments are required: conns (list of DSConnection-classes), 
#' symbol (name of the R symbol) (see datashield.assign).
#' @param symbol Encoded vector of names of the R symbols to assign in the Datashield R session on each server in \code{logins}.
#' The assigned R variables will be used as the input raw data.
#' Other assigned R variables in \code{func} are ignored.
#' @param H Number of common dimensions
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
#' @import DSI
#' @importFrom utils setTxtProgressBar
#' @export
federateComDimRm <- function(loginFD, logins, func, symbol, H = 2, scale = "none", option = "uniform", chunk = 500, mc.cores = 1, threshold = 1e-10) {
    require(DSOpal)
    .printTime("federateComDim started")
    TOL <- 1e-10
    funcPreProc <- .decode.arg(func)
    querytables <- .decode.arg(symbol)
    ntab <- length(querytables)
    
    ## compute SSCP matrix for each centered data table
    XX <- lapply(1:ntab, function(i) {
        .federateSSCP(loginFD=loginFD, logins=logins, funcPreProc=funcPreProc, querytables=querytables, ind=i, byColumn=TRUE, chunk=chunk, mc.cores=mc.cores, TOL=TOL)
    })
    names(XX) <- querytables
    
    ## set up the centered data table on every node
    loginFDdata <- .decode.arg(loginFD)
    logindata <- .decode.arg(logins)
    opals <- .login(logins=logindata) #datashield.login(logins=logindata)
    nNode <- length(opals)
    
    tryCatch({
        ## take a snapshot of the current session
        safe.objs <- .ls.all()
        safe.objs[['.GlobalEnv']] <- setdiff(safe.objs[['.GlobalEnv']], '.Random.seed')  # leave alone .Random.seed for sample()
        ## lock everything so no objects can be changed
        .lock.unlock(safe.objs, lockBinding)
        
        ## apply funcPreProc for preparation of querytables on opals
        ## TODO: control hacking!
        ## TODO: control identical colnames!
        funcPreProc(conns=opals, symbol=querytables)
        
        ## unlock back everything
        .lock.unlock(safe.objs, unlockBinding)
        ## get rid of any sneaky objects that might have been created in the filters as side effects
        .cleanup(safe.objs)
    }, error=function(e) {
        print(paste0("DATA MAKING PROCESS: ", e))
        return (paste0("DATA MAKING PROCESS: ", e, ' --- ', datashield.symbols(opals), ' --- ', datashield.errors(), ' --- ', datashield.logout(opals)))
    })
    
    ## take variables (colnames)
    queryvariables <- lapply(querytables, function(querytable) {
        DSI::datashield.aggregate(opals[1], as.symbol(paste0('colNames(', querytable, ')')), async=F)[[1]]
    })
    names(queryvariables) <- querytables
    
    ## centered cbind-ed centered data matrix
    DSI::datashield.assign(opals, "centeredAllData",
                           as.symbol(paste0('center(list(', paste(querytables, collapse=','), '), byColumn=TRUE, na.rm=FALSE)')),
                           async=T)
    ## compute the total variance of a dataset
    inertie <- function(tab) {
        return (sum(diag(tab)))    #Froebenius norm
    }

    ## compute the RV between WX and WY
    coefficientRV <- function(WX, WY) {
        rv <- sum(diag(WX %*% WY))/((sum(diag(WX %*% WX)) * sum(diag(WY %*% WY)))^0.5)
        return(rv)
    }
    # ---------------------------------------------------------------------------
    # 0. Preliminary tests
    # ---------------------------------------------------------------------------
    if (any(sapply(XX, is.na)))
        stop("No NA values are allowed")
    nsamples <- unique(unlist(apply(sapply(XX, dim), 1, unique))) ## number of samples
    if (length(nsamples) > 1)
        stop("XX elements should be symmetric of the same dimension")
    samples <- sapply(XX, function(x) union(rownames(x), colnames(x)))
    if (is.list(samples) && max(lengths(samples))==0) {
        XX <- lapply(XX, function(x) {
            rownames(x) <- colnames(x) <- paste("X", 1:nsamples, sep='.')
            return (x)
        })
        samples <- sapply(XX, function(x) union(rownames(x), colnames(x)))
    }
    if (is.list(samples) || is.list(apply(samples, 1, unique)))
        stop("XX elements should have the same rownames and colnames")
    
    ## TOREVIEW
    # if (is.character(scale)) {
    #   if (!scale %in% c("none","sd"))
    #     stop("Non convenient scaling parameter")
    # }
    # else {
    #   if (!is.numeric(scale) | length(scale)!=ncol(X))
    #     stop("Non convenient scaling parameter")
    # }
    # if (!option %in% c("none","uniform"))
    #   stop("Non convenient weighting parameter")
    
    
    # ---------------------------------------------------------------------------
    # 1. Output preparation
    # ---------------------------------------------------------------------------
    nvar <- lengths(queryvariables)
    W <- array(0, dim=c(nsamples, nsamples, ntab+1)) # association matrices
    
    LAMBDA <- matrix(0, nrow=ntab, ncol=H)   # will contains the saliences
    Q <- matrix(0, nrow=nsamples, ncol=H)    # will contain common components
    J <- rep(1:ntab, times=nvar)             # indicates which block each variable belongs to
    names.H <- paste("Dim.", 1:H, sep="")
    Q.b <- array(0, dim=c(nsamples, H, ntab))  # components for the block components
    dimnames(Q.b) <- list(rownames(XX[[1]]), names.H, querytables)
    W.b <- vector("list", length=ntab)        # weights for the block components
    P.b <- vector("list", length=ntab)        # loadings for the block components
    for (k in 1:ntab) {
        W.b[[k]] <- matrix(0, nrow=nvar[k], ncol=H)
        P.b[[k]] <- matrix(0, nrow=nvar[k], ncol=H)
        rownames(W.b[[k]]) <- rownames(P.b[[k]]) <- queryvariables[[k]]
        colnames(W.b[[k]]) <- colnames(P.b[[k]]) <- names.H
    }
    We <- Pe <- matrix(0, nrow=sum(nvar), ncol=H)
    Res <- NULL              # Results to be returned
    
    explained.block <- matrix(0, nrow=ntab+1,ncol=H)      # percentage of inertia recovered
    
    # ---------------------------------------------------------------------------
    # 2. Required parameters and data preparation
    # ---------------------------------------------------------------------------
    
    
    #Xscale$mean <- apply(X, 2, mean)
    #X<-scale(X, center=Xscale$mean, scale=FALSE)   #default centering
    
    
    # if (scale=="none") {
    #   Xscale$scale <-rep(1,times=ncol(X))
    # }
    # else {
    #   if (scale=="sd") {
    #     sd.tab <- apply(X, 2, function (x) {return(sqrt(sum(x^2)/length(x)))})   #sd based on biased variance
    #     temp <- sd.tab < 1e-14
    #     if (any(temp)) {
    #       warning("Variables with null variance not standardized.")
    #       sd.tab[temp] <- 1
    #     }
    #     X <- sweep(X, 2, sd.tab, "/")
    #     Xscale$scale <-sd.tab
    #   }
    #   else {     #specific scaling depending on blocks defined as a vector with a scaling parameter for each variable
    #       X <- sweep(X, 2, scale, "/")
    #       Xscale$scale <-scale
    #   }
    # }
    
    # Pre-processing: block weighting
    inertia <- sapply(1:ntab, function(k) inertie(XX[[k]]))
    # if (option=="uniform") {
    #    #set each block inertia equal 1
    #   w.tab <- rep(sqrt(inertia), times=nvar)     # weighting parameter applied to each variable
    #   X <- sweep(X, 2, w.tab, "/")
    #   Xscale$scale<- Xscale$scale*w.tab
    #   inertia <- rep(1, times=ntab)
    # }
    if (option=="uniform") {
        XX <- lapply(1:ntab, function(k) {
            XX[[k]]/inertia[k]
        })
        inertia0.sqrt <- sqrt(inertia)
        inertia <- rep(1, times=ntab)
    }
    
    # Computation of association matrices
    W[,,1:ntab] <- array(as.numeric(unlist(XX)), dim=c(nsamples, nsamples, ntab))
    tvar <- sapply(1:ntab, function(k) sum(as.matrix(W[,,k])^2))
    Itot <- sum(tvar) # Total inertia of all dataset sum(trace(Wj*Wj))
    
    #X0 <- X            #keep initial values with standardization and weighting scheme
    # ---------------------------------------------------------------------------
    # 3. computation of Q and LAMBDA for the various dimensions
    # ---------------------------------------------------------------------------
    explained <- matrix(0, nrow=H, ncol=1)
    
    pb <- txtProgressBar(min=0, max=H, style=3)
    
    for (dimension in 1:H)  {
        Sys.sleep(0.1)
        
        previousfit <- 100000;
        lambda <- rep(1, ntab)
        deltafit <- 1000000;
        while (deltafit > threshold) {
            W[,,ntab+1] <- Reduce("+", lapply(1:ntab, function(k) lambda[k]*W[,,k]))
            Svdw <- svd(as.matrix(W[,,ntab+1]))
            q <- Svdw$u[, 1, drop=F]
            
            fit <- 0
            for (k in 1:ntab) {
                # estimating residuals
                lambda[k]  <- (t(q) %*% as.matrix(W[,,k]) %*% q)
                pred       <- lambda[k]*q %*% t(q)
                aux        <- as.matrix(W[,,k]) - pred
                fit        <- fit + sum(aux^2)
            }
            deltafit <- previousfit - fit
            previousfit <- fit
        } #deltafit>threshold
        
        explained[dimension,1] <- 100*sum(lambda^2)/Itot  ## vca modif
        LAMBDA[,dimension] <- lambda
        Q[,dimension] <- q
        
        # updating association matrices
        proj <- diag(1, nsamples) - tcrossprod(q)
        for (k in 1:ntab)   {
            #W.b[[k]][,dimension] <- t(X[,J==k]) %*% q
            ##TODO
            
            #Q.b[,dimension,k]    <- X[,J==k]%*%matrix(W.b[[k]][,dimension],ncol=1)
            Q.b[,dimension,k] <- W[,,k] %*% q
            
            #P.b[[k]][,dimension] <- t(W.b[[k]][,dimension])#t(q)%*%X[,J==k]
            ##TODO
            
            #Pe[J==k,dimension]   <- P.b[[k]][,dimension]
            ##TODO
            
            #X <- as.matrix(X)
            #X.hat <- q%*%matrix(P.b[[k]][,dimension],nrow=1)
            ##TODO
            
            #explained.block[k,dimension] <- inertie(X.hat)
            explained.block[k,dimension] <- inertie(tcrossprod(q) %*% W[,,k] %*% tcrossprod(q))
            
            #X[,J==k] <- X[,J==k]-X.hat
            ##TODO
            
            #W[,,k] <- X[,J==k]%*%t(X[,J==k]);
            W[,,k] <- proj %*% W[,,k] %*% t(proj)
        }
        explained.block[ntab+1,dimension]<- sum(explained.block[1:ntab,dimension])
        #We[,dimension]   <- unlist(sapply(1:ntab,function(j){lambda[j]*W.b[[j]][,dimension]}))
        #We[,dimension]   <- We[,dimension] / sum(lambda^2)
        setTxtProgressBar(pb, dimension)
    }
    
    ## loadings
    
    # number of samples on each node
    tryCatch({
        size <- sapply(datashield.aggregate(opals, as.symbol('dsDim(centeredAllData)'), async=T), function(x) x[1])
    }, error=function(e) {e; datashield.logout(opals)})
    size <- c(0, size)
    func <- function(x, y) {x %*% y}
    Qlist <- setNames(lapply(2:length(size), function(i) {
        Qi <- Q[(cumsum(size)[i-1]+1):cumsum(size)[i],,drop=F]
        ## As Q is orthonormal, Qi == Qi.iter
        # Qi.iter <- sapply(1:H, function(dimension) {
        #   projs <- lapply(setdiff(1:dimension, dimension), function(dimprev) {
        #     return (diag(1, size[i]) - tcrossprod(Qi[,dimprev]))
        #   })
        #   projs <- c(Id=list(diag(1, size[i])), projs)
        #   return (crossprod(Reduce(func, projs), Qi[,dimension,drop=F]))
        # })
        # return (Qi.iter)
    }), names(opals))
    tryCatch({
        Wbk <- Reduce('+', unlist(mclapply(names(opals), mc.cores=1, function(opn) {
            expr <- list(as.symbol("loadings"),
                         as.symbol("centeredAllData"),
                         .encode.arg(Qlist[[opn]]))
            loadings <- datashield.aggregate(opals[opn], as.call(expr), async=T)
            return (loadings)
        }), recursive = F))
    }, error=function(e) e, finally=datashield.logout(opals))
    colnames(Wbk) <- names.H
    csnvar <- cumsum(nvar)
    W.b <- mclapply(1:length(nvar), mc.cores=length(nvar), function(k) {
        if (option=="uniform") return (Wbk[ifelse(k==1, 1, csnvar[k-1]+1):csnvar[k], , drop=F]/inertia0.sqrt[k])
        return (Wbk[ifelse(k==1, 1, csnvar[k-1]+1):csnvar[k], , drop=F])
    })
    
    # W.b <- lapply(1:ntab, function(k) {
    #     #Wbk <- crossprod(as.matrix(X[,J==k]), Q)
    #     Wbk <- Reduce('+', unlist(mclapply(names(opals), mc.cores=1, function(opn) {
    #         expr <- list(as.symbol("loadings"),
    #                      as.symbol("centeredAllData"),
    #                      .encode.arg(Qlist[[opn]]))
    #         loadings <- datashield.aggregate(opals[opn], as.call(expr), async=F)
    #         return (loadings)
    #     }), recursive = F))
    #     
    #     colnames(Wbk) <- names.H
    #     return (Wbk/inertia0.sqrt)
    # })
    # return (W.b)

    We <- do.call(rbind, lapply(1:ntab, function(k) tcrossprod(W.b[[k]], diag(LAMBDA[k,]))))
    #We <- do.call(rbind, lapply(1:ntab, function(k) W.b[[k]] %*% diag(LAMBDA[k,]))) #crossprod(as.matrix(X[,J==k]), Q)
    
    P.b <- W.b #lapply(W.b, function(x) t(x))
    Pe  <- do.call(rbind, P.b)
    We  <- sapply(1:H, function(dimension) We[,dimension]/sum(LAMBDA[,dimension]^2))
    colnames(We) <- names.H
    
    # ---------------------------------------------------------------------------
    # 4.1 Preparation of the results Global
    # ---------------------------------------------------------------------------
    close(pb)
    
    # Overall agreement
    if (H==1) {
        LambdaMoyen <- apply(LAMBDA, 2, mean)
        C <- Q*LambdaMoyen
    } else {
        LambdaMoyen <- apply(LAMBDA, 2, mean)
        C <- Q %*% sqrt(diag(LambdaMoyen))
    }
    
    rownames(C) <- rownames(XX[[1]])
    colnames(C) <- names.H
    
    RV <- sapply(1:ntab, function(k) {
        coefficientRV(XX[[k]], tcrossprod(C))
    })
    names(RV) <- names(queryvariables)
    
    # global components, saliences and explained variances
    Res$saliences <- LAMBDA
    colnames(Res$saliences) <- names.H
    rownames(Res$saliences) <- names(queryvariables)
    Res$Q <- Q
    rownames(Res$Q) <- rownames(XX[[1]])
    colnames(Res$Q) <- names.H
    Res$C  <- C
    Res$RV <- RV
    Res$W  <- We
    Res$Wm <- We %*% solve(t(Pe) %*% We)
    rownames(Res$Wm) <- rownames(Res$W)
    colnames(Res$Wm) <- colnames(Res$W) <- names.H
    
    fit <- matrix(0, nrow=H, ncol=2)
    fit[,1] <- explained
    fit[,2] <- cumsum(fit[,1])
    Res$fit <- fit
    rownames(Res$fit) <- names.H
    colnames(Res$fit) <- c("%Fit", "%Cumul Fit")
    
    inertia                 <- c(inertia, sum(inertia))
    Res$explained           <- sweep(explained.block, 1, inertia, "/")
    rownames(Res$explained) <- c(rownames(Res$saliences), 'global')
    colnames(Res$explained) <- colnames(Res$saliences)
    
    
    # ---------------------------------------------------------------------------
    # 4.2 Preparation of the results Blocks
    # ---------------------------------------------------------------------------
    Block     <- NULL
    Block$Q.b <- Q.b
    Block$W.b <- W.b
    Block$P.b <- P.b
    Res$Block <- Block
    
    # ---------------------------------------------------------------------------
    # Return Res
    # ---------------------------------------------------------------------------
    Res$call   <- match.call()
    class(Res) <- c("federateComDim")
    
    return(Res)
}


#' @title Federated SNF
#' @description Function for SNF federated analysis on the virtual cohort combining multiple cohorts
#' @usage federateSNF(loginFD, logins, func, symbol, K = 20, sigma = 0.5, t = 20)
#' @param loginFD Login information of the FD server
#' @param logins Login information of data repositories
#' @param func Encoded definition of a function for preparation of raw data matrices. 
#' Two arguments are required: conns (list of DSConnection-classes), 
#' symbol (name of the R symbol) (see datashield.assign).
#' @param symbol Encoded vector of names of the R symbols to assign in the Datashield R session on each server in \code{logins}.
#' The assigned R variables will be used as the input raw data.
#' Other assigned R variables in \code{func} are ignored.
#' @param metric Either \code{euclidean} or \code{correlation} for distance metric between samples. 
#' For Euclidean distance, the data from each cohort will be centered (not scaled) for each variable.
#' For correlation-based distance, the data from each cohort will be centered scaled for each sample.
#' @param K Number of neighbors in K-nearest neighbors part of the algorithm, see \code{SNFtool::SNF}.
#' @param sigma Variance for local model, see \code{SNFtool::affinityMatrix}.
#' @param t Number of iterations for the diffusion process, see \code{SNFtool::SNF}.
#' @return The overall status matrix derived W.
#' @import SNFtool DSI
#' @export
federateSNF <- function(loginFD, logins, func, symbol, metric = 'euclidean', K = 20, sigma = 0.5, t = 20, chunk = 500, mc.cores = 1) {
    require(DSOpal)
    .printTime("federateSNF started")
    TOL <- 1e-10
    funcPreProc <- .decode.arg(func)
    querytables <- .decode.arg(symbol)
    ntab <- length(querytables)
    metric <- match.arg(metric, choices=c('euclidean', 'correlation'))
    logindata <- .decode.arg(logins)
    opals <- .login(logins=logindata) #datashield.login(logins=logindata)

    tryCatch({
        ## take a snapshot of the current session
        safe.objs <- .ls.all()
        safe.objs[['.GlobalEnv']] <- setdiff(safe.objs[['.GlobalEnv']], '.Random.seed')  # leave alone .Random.seed for sample()
        ## lock everything so no objects can be changed
        .lock.unlock(safe.objs, lockBinding)
        
        ## apply funcPreProc for preparation of querytables on opals
        ## TODO: control hacking!
        ## TODO: control identical colnames!
        funcPreProc(conns=opals, symbol=querytables)
        
        ## unlock back everything
        .lock.unlock(safe.objs, unlockBinding)
        ## get rid of any sneaky objects that might have been created in the filters as side effects
        .cleanup(safe.objs)
    }, error=function(e) {
        print(paste0("DATA MAKING PROCESS: ", e))
        return (paste0("DATA MAKING PROCESS: ", e, ' --- ', datashield.symbols(opals), ' --- ', datashield.errors(), ' --- ', datashield.logout(opals)))
    })
    
    ## take variables (colnames)
    queryvariables <- lapply(querytables, function(querytable) {
        DSI::datashield.aggregate(opals[1], as.symbol(paste0('colNames(', querytable, ')')), async=F)[[1]]
    })
    names(queryvariables) <- querytables
    DSI::datashield.logout(opals)
    
    if (metric == "correlation") {
        ## compute (1 - correlation) distance between samples for each data table 
        XX <- lapply(1:ntab, function(i) {
            1 - .federateSSCP(loginFD=loginFD, logins=logins, 
                              funcPreProc=funcPreProc, querytables=querytables, ind=i, 
                              byColumn=FALSE, chunk=chunk, mc.cores=mc.cores, TOL=TOL)/(length(queryvariables[[i]])-1)
        })
    } else if (metric == "euclidean"){
        ## compute Euclidean distance between samples for each data table 
        XX <- lapply(1:ntab, function(i) {
            .toEuclidean(.federateSSCP(loginFD=loginFD, logins=logins,
                                       funcPreProc=funcPreProc, querytables=querytables, ind=i,
                                       byColumn=TRUE, chunk=chunk, mc.cores=mc.cores, TOL=TOL))
        })
    }
    
    ## take common samples
    commons <- Reduce(intersect, lapply(XX, rownames))
    XX <- lapply(XX, function(distmat) {
        distmat[commons,commons]
    })
    ## similarity graphs
    Ws <- lapply(XX, function(distmat) {
        affinityMatrix(distmat, K, sigma)
    })
    
    ## fuse similarity graphs
    W <- SNF(Ws, K, t)
    
    return (W)
}


#' @title Federated UMAP
#' @description Function for UMAP federated analysis on the virtual cohort combining multiple cohorts
#' @usage federateUMAP(loginFD, logins, func, symbol, TOL = 1e-10, metric = 'euclidean', ...)
#' @param loginFD Login information of the FD server
#' @param logins Login information of data repositories
#' @param func Encoded definition of a function for preparation of raw data matrices. 
#' Two arguments are required: conns (list of DSConnection-classes), 
#' symbol (name of the R symbol) (see datashield.assign).
#' @param symbol Encoded vector of names of the R symbols to assign in the Datashield R session on each server in \code{logins}.
#' The assigned R variables will be used as the input raw data.
#' Other assigned R variables in \code{func} are ignored.
#' @param metric Either \code{euclidean} or \code{correlation} for distance metric between samples. 
#' For Euclidean distance, the data from each cohort will be centered (not scaled) for each variable.
#' For correlation-based distance, the data from each cohort will be centered scaled for each sample.
#' @param ... see \code{uwot::umap}
#' @return A matrix of optimized coordinates.
#' @import uwot DSI
#' @importFrom "stats" "as.dist" "cor" "quantile" "setNames"
#' @export
federateUMAP <- function(loginFD, logins, func, symbol, metric = 'euclidean', chunk = 500, mc.cores = 1, ...) {
    require(DSOpal)
    .printTime("federateUMAP started")
    TOL <- 1e-10
    funcPreProc <- .decode.arg(func)
    querytables <- .decode.arg(symbol)
    ntab <- length(querytables)
    metric <- match.arg(metric, choices=c('euclidean', 'correlation'))
    logindata <- .decode.arg(logins)
    opals <- .login(logins=logindata) #datashield.login(logins=logindata)
    
    tryCatch({
        ## take a snapshot of the current session
        safe.objs <- .ls.all()
        safe.objs[['.GlobalEnv']] <- setdiff(safe.objs[['.GlobalEnv']], '.Random.seed')  # leave alone .Random.seed for sample()
        ## lock everything so no objects can be changed
        .lock.unlock(safe.objs, lockBinding)
        
        ## apply funcPreProc for preparation of querytables on opals
        ## TODO: control hacking!
        ## TODO: control identical colnames!
        funcPreProc(conns=opals, symbol=querytables)
        
        ## unlock back everything
        .lock.unlock(safe.objs, unlockBinding)
        ## get rid of any sneaky objects that might have been created in the filters as side effects
        .cleanup(safe.objs)
    }, error=function(e) {
        print(paste0("DATA MAKING PROCESS: ", e))
        return (paste0("DATA MAKING PROCESS: ", e, ' --- ', datashield.symbols(opals), ' --- ', datashield.errors(), ' --- ', datashield.logout(opals)))
    })
    .printTime("federateUMAP data processed")
    
    ## take variables (colnames)
    queryvariables <- lapply(querytables, function(querytable) {
        DSI::datashield.aggregate(opals[1], as.symbol(paste0('colNames(', querytable, ')')), async=F)[[1]]
    })
    names(queryvariables) <- querytables
    DSI::datashield.logout(opals)
    
    if (metric == "correlation") {
        ## compute (1 - correlation) distance between samples for each data table 
        XX <- lapply(1:ntab, function(i) {
            as.dist(1 - .federateSSCP(loginFD=loginFD, logins=logins, 
                                      funcPreProc=funcPreProc, querytables=querytables, ind=i, 
                                      byColumn=FALSE, chunk=chunk, mc.cores=mc.cores, TOL=TOL)/(length(queryvariables[[i]])-1))
        })
    } else if (metric == "euclidean"){
        ## compute Euclidean distance between samples for each data table 
        XX <- lapply(1:ntab, function(i) {
            as.dist(.toEuclidean(.federateSSCP(loginFD=loginFD, logins=logins, 
                                               funcPreProc=funcPreProc, querytables=querytables, ind=i, 
                                               byColumn=TRUE, chunk=chunk, mc.cores=mc.cores, TOL=TOL)))
        })
    }
    
    return (setNames(lapply(1:ntab, function(i) uwot::umap(XX[[i]], ...)), querytables)) #ret_model = FALSE, ret_nn = FALSE, ret_extra = c(),
}


#' @title Federated hdbscan
#' @description Function for hdbscan federated analysis on the virtual cohort combining multiple cohorts
#' @usage federateHdbscan(loginFD, logins, func, symbol, TOL = 1e-10, metric = 'euclidean', ...)
#' @param loginFD Login information of the FD server
#' @param logins Login information of data repositories
#' @param func Encoded definition of a function for preparation of raw data matrices. 
#' Two arguments are required: conns (list of DSConnection-classes), 
#' symbol (name of the R symbol) (see datashield.assign).
#' @param symbol Encoded vector of names of the R symbols to assign in the Datashield R session on each server in \code{logins}.
#' The assigned R variables will be used as the input raw data.
#' Other assigned R variables in \code{func} are ignored.
#' @param metric Either \code{euclidean} or \code{correlation} for distance metric between samples. 
#' For Euclidean distance, the data from each cohort will be centered (not scaled) for each variable.
#' For correlation-based distance, the data from each cohort will be centered scaled for each sample.
#' @param ... see \code{dbscan::hdbscan}
#' @return An object of class \code{hdbscan}.
#' @import dbscan DSI
#' @export
federateHdbscan <- function(loginFD, logins, func, symbol, metric = 'euclidean', minPts = 10, chunk = 500, mc.cores = 1, ...) {
    require(DSOpal)
    .printTime("federateHdbscan started")
    TOL <- 1e-10
    funcPreProc <- .decode.arg(func)
    querytables <- .decode.arg(symbol)
    ntab <- length(querytables)
    metric <- match.arg(metric, choices=c('euclidean', 'correlation'))
    logindata <- .decode.arg(logins)
    opals <- .login(logins=logindata) #datashield.login(logins=logindata)
    
    tryCatch({
        ## take a snapshot of the current session
        safe.objs <- .ls.all()
        safe.objs[['.GlobalEnv']] <- setdiff(safe.objs[['.GlobalEnv']], '.Random.seed')  # leave alone .Random.seed for sample()
        ## lock everything so no objects can be changed
        .lock.unlock(safe.objs, lockBinding)
        
        ## apply funcPreProc for preparation of querytables on opals
        ## TODO: control hacking!
        ## TODO: control identical colnames!
        funcPreProc(conns=opals, symbol=querytables)
        
        ## unlock back everything
        .lock.unlock(safe.objs, unlockBinding)
        ## get rid of any sneaky objects that might have been created in the filters as side effects
        .cleanup(safe.objs)
    }, error=function(e) {
        print(paste0("DATA MAKING PROCESS: ", e))
        return (paste0("DATA MAKING PROCESS: ", e, ' --- ', datashield.symbols(opals), ' --- ', datashield.errors(), ' --- ', datashield.logout(opals)))
    })
    
    ## take variables (colnames)
    queryvariables <- lapply(querytables, function(querytable) {
        DSI::datashield.aggregate(opals[1], as.symbol(paste0('colNames(', querytable, ')')), async=F)[[1]]
    })
    names(queryvariables) <- querytables
    DSI::datashield.logout(opals)
    
    if (metric == "correlation") {
        ## compute (1 - correlation) distance between samples for each data table 
        XX <- lapply(1:ntab, function(i) {
            as.dist(1 - .federateSSCP(loginFD=loginFD, logins=logins, 
                                      funcPreProc=funcPreProc, querytables=querytables, ind=i, 
                                      byColumn=FALSE, chunk=chunk, mc.cores=mc.cores, TOL=TOL)/(length(queryvariables[[i]])-1))
        })
    } else if (metric == "euclidean"){
        ## compute Euclidean distance between samples for each data table 
        XX <- lapply(1:ntab, function(i) {
            as.dist(.toEuclidean(.federateSSCP(loginFD=loginFD, logins=logins, 
                                               funcPreProc=funcPreProc, querytables=querytables, ind=i, 
                                               byColumn=TRUE, chunk=chunk, mc.cores=mc.cores, TOL=TOL)))
        })
    }
    
    return (setNames(lapply(1:ntab, function(i) dbscan::hdbscan(XX[[i]], minPts = minPts, ...)), querytables))
}


#' @title Test federateSSCP
#' @description Deprecated
#' @param loginFD Login information of the FD server
#' @param logins Login information of data repositories
#' @param func Encoded definition of a function for preparation of raw data matrices. 
#' Two arguments are required: conns (list of DSConnection-classes), 
#' symbol (name of the R symbol) (see datashield.assign).
#' @param symbol Encoded vector of names of the R symbols to assign in the Datashield R session on each server in \code{logins}.
#' The assigned R variables will be used as the input raw data.
#' Other assigned R variables in \code{func} are ignored.
#' @param byColumn A logical value indicating whether the input data is centered by column or row.
#' Default, TRUE, centering by column. Constant variables across samples are removed. 
#' If FALSE, centering and scaling by row. Constant samples across variables are removed.
#' @param TOL Tolerance of 0
#' @import DSOpal parallel bigmemory
# ' @export
#' @keywords internal
testSSCP <- function(loginFD, logins, func, symbol, metric = 'euclidean', chunk = 500, mc.cores = 1, TOL = 1e-10) {
    funcPreProc <- .decode.arg(func)
    querytables <- .decode.arg(symbol)
    ntab <- length(querytables)
    
    metric <- match.arg(metric, choices=c('euclidean', 'correlation'))
    logindata <- .decode.arg(logins)
    opals <- .login(logins=logindata) #datashield.login(logins=logindata)
    
    tryCatch({
        ## take a snapshot of the current session
        safe.objs <- .ls.all()
        safe.objs[['.GlobalEnv']] <- setdiff(safe.objs[['.GlobalEnv']], '.Random.seed')  # leave alone .Random.seed for sample()
        ## lock everything so no objects can be changed
        .lock.unlock(safe.objs, lockBinding)
        
        ## apply funcPreProc for preparation of querytables on opals
        ## TODO: control hacking!
        ## TODO: control identical colnames!
        funcPreProc(conns=opals, symbol=querytables)
        
        ## unlock back everything
        .lock.unlock(safe.objs, unlockBinding)
        ## get rid of any sneaky objects that might have been created in the filters as side effects
        .cleanup(safe.objs)
    }, error=function(e) {
        print(paste0("DATA MAKING PROCESS: ", e))
        return (paste0("DATA MAKING PROCESS: ", e, ' --- ', datashield.symbols(opals), ' --- ', datashield.errors(), ' --- ', datashield.logout(opals)))
    })
    
    ## take variables (colnames)
    queryvariables <- lapply(querytables, function(querytable) {
        DSI::datashield.aggregate(opals[1], as.symbol(paste0('colNames(', querytable, ')')), async=F)[[1]]
    })
    names(queryvariables) <- querytables
    DSI::datashield.logout(opals)
    
    if (metric == "correlation") {
        ## compute (1 - correlation) distance between samples for each data table 
        XX <- lapply(1:ntab, function(i) {
            as.dist(1 - .federateSSCP(loginFD=loginFD, logins=logins, 
                                      funcPreProc=funcPreProc, querytables=querytables, ind=i, 
                                      byColumn=FALSE, chunk=chunk, mc.cores=mc.cores, TOL=TOL)/(length(queryvariables[[i]])-1))
        })
    } else if (metric == "euclidean"){
        ## compute Euclidean distance between samples for each data table 
        XX <- lapply(1:ntab, function(i) {
            as.dist(.toEuclidean(.federateSSCP(loginFD=loginFD, logins=logins, 
                                               funcPreProc=funcPreProc, querytables=querytables, ind=i, 
                                               byColumn=TRUE, chunk=chunk, mc.cores=mc.cores, TOL=TOL)))
        })
    }
    
    return (XX)
}
