#' @title Garbage collection
#' @description Call gc on the federated server.
#' @keywords internal
garbageCollect <- function() {
    gc(reset=T)
    return (NULL)
}


#' @title Bigmemory description of a matrix
#' @description Bigmemory description of a matrix.
#' @param value Encoded value of a matrix.
#' @import bigmemory
#' @returns Bigmemory description of the given matrix
#' @keywords internal
matrix2DscFDrm <- function(value) {
    valued <- .decode.arg(value)
    tcp <- do.call(rbind, .decode.arg(valued))
    dscbigmatrix <- describe(as.big.matrix(tcp, backingfile = ""))
    rm(list=c("valued", "tcp"))
    return (dscbigmatrix)
}


#' @title Bigmemory description of a matrix
#' @description Bigmemory description of a matrix.
#' @param value Encoded value of a matrix.
#' @returns Bigmemory description of the given matrix
#' @importFrom arrow read_ipc_stream
#' @importFrom bigmemory as.big.matrix describe
#' @export
matrix2DscFD <- function(value) {
    tcp <- as.matrix(read_ipc_stream(.decode.arg(value)))
    dscbigmatrix <- describe(as.big.matrix(tcp, backingfile = ""))
    rm(list=c("tcp"))
    return (dscbigmatrix)
}


#' @title Matrix partition
#' @description Partition a matrix into blocks.
#' @param x A matrix.
#' @param seprow A numeric vectors indicating sizes of blocks in rows.
#' @param sepcol A numeric vectors indicating sizes of blocks in columns.
#' @returns List of blocks
#' @importFrom arrow arrow_table
#' @keywords internal
.partitionMatrix <- function(x, seprow, sepcol = seprow) {
    stopifnot(sum(seprow)==nrow(x) && sum(sepcol)==ncol(x))
    csseprow <- cumsum(seprow)
    indrow <- lapply(1:length(seprow), function(i) {
        return (c(ifelse(i==1, 0, csseprow[i-1])+1, csseprow[i]))
    })
    cssepcol <- cumsum(sepcol)
    indcol <- lapply(1:length(sepcol), function(i) {
        return (c(ifelse(i==1, 0, cssepcol[i-1])+1, cssepcol[i]))
    })
    parMat <- lapply(1:length(indrow), function(i) {
        lapply(ifelse(isSymmetric(x), i, 1):length(indcol), function(j) {
            xij <- x[indrow[[i]][1]:indrow[[i]][2],
                     indcol[[j]][1]:indcol[[j]][2], drop=F]
            return (arrow_table(as.data.frame(xij)))
        })
    })
    
    return (parMat)
}


#' @title Matrix reconstruction
#' @description Rebuild a matrix from its partition.
#' @param matblocks List of lists of matrix blocks, obtained from
#' .partitionMatrix.
#' @param mc.cores Number of cores for parallel computing. Default, 1.
#' @returns The reconstructed matrix.
#' @importFrom parallel mclapply
#' @keywords internal
.rebuildMatrix <- function(matblocks, mc.cores = 1) {
    uptcp <- lapply(matblocks, function(bl) do.call(cbind, bl))
    ## combine the blocks into one matrix
    if (length(uptcp) > 1) {
        if (length(unique(sapply(uptcp, ncol)))==1) {
            tcp <- do.call(rbind, uptcp)
        } else {
            ## without the first layer of blocks
            no1tcp <- mclapply(
                2:length(uptcp),
                mc.cores=mc.cores,
                function(i) {
                    cbind(do.call(cbind, lapply(1:(i-1), function(j) {
                        t(matblocks[[j]][[i-j+1]])
                    })), uptcp[[i]])
                })
            ## with the first layer of blocks
            tcp <- rbind(uptcp[[1]], do.call(rbind, no1tcp))
            rm(list=c("no1tcp"))
            rownames(tcp) <- colnames(tcp)
            stopifnot(isSymmetric(tcp))
        }
    } else {
        ## for either asymmetric matrix or column-wise blocks
        tcp <- uptcp[[1]]
    }
    rm(list=c("uptcp"))
    
    return (tcp)
}


#' @title Symmetric matrix reconstruction
#' @description Rebuild a matrix from its partition.
#' @param matblocks List of lists of matrix blocks, obtained from
#' .partitionMatrix.
#' @param mc.cores Number of cores for parallel computing. Default, 1.
#' @returns The complete symmetric matrix
#' @keywords internal
.rebuildMatrixrm <- function(matblocks, mc.cores = 1) {
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
#' @returns The complete symmetric matrix
#' @keywords internal
.rebuildMatrixDscrm <- function(dscblocks, mc.cores = 1) {
    ## access to matrix blocks 
    matblocks <- mclapply(dscblocks, mc.cores=mc.cores, function(y) {
        lapply(y, function(x) {
            return (as.matrix(attach.big.matrix(x)))
        })
    })
    tcp <- .rebuildMatrix(matblocks, mc.cores=mc.cores)
    return (tcp)
}


#' @title Symmetric matrix reconstruction
#' @description Rebuild a symmetric matrix from its partition bigmemory objects
#' @param dscblocks List of lists of bigmemory objects pointed to matrix blocks
#' @param mc.cores Number of cores for parallel computing. Default: 1
#' @returns The complete symmetric matrix
#' @importFrom parallel mclapply
#' @importFrom bigmemory attach.big.matrix
#' @keywords internal
.rebuildMatrixDsc <- function(dscblocks, mc.cores = 1) {
    ## access to matrix blocks 
    matblocks <- mclapply(dscblocks, mc.cores=mc.cores, function(y) {
        lapply(y, function(x) {
            return ((attach.big.matrix(x))[,,drop=F])
        })
    })
    tcp <- .rebuildMatrix(matblocks, mc.cores=mc.cores)
    return (tcp)
}


#' @title Euclidean distance
#' @description Transform XX' matrix into Euclidean distance between samples
#' (rows) in X.
#' @param XXt An SSCP matrix XX'.
#' @returns Euclidean distance
#' @keywords internal
.toEuclidean <- function(XXt) {
    if (!isSymmetric(XXt) || any(rownames(XXt) != colnames(XXt)))
        stop('Input XXt (an SSCP matrix) should be symmetric.')
    .printTime("toEuclidean XXt started")
    lowerTri <- cbind(do.call(cbind, lapply(1:(ncol(XXt)-1), function(i) {
        res <- sapply((i+1):ncol(XXt), function(j) {
            ## d(Xi, Xj)^2 = Xi'Xi - 2Xi'Xj + Xj'Xj
            ## = XXt[i,i] - 2XXt[i,j] + XXt[j,j]
            t(c(1,-1)) %*% XXt[c(i,j), c(i,j)] %*% c(1, -1)
        })
        return (c(rep(0, i), res))
    })), rep(0, ncol(XXt)))
    distmat <- sqrt(lowerTri + t(lowerTri))
    rownames(distmat) <- rownames(XXt)
    colnames(distmat) <- colnames(XXt)
    
    return (distmat)
}


#' @title Find X from XX' and X'X
#' @description Find X from XX' and X'X.
#' @param XXt XX'
#' @param XtX X'X
#' @param r A non-null vector of length \code{ncol(X'X)}
#' @param Xr A vector of length \code{nrow(XX'}, equals to the product Xr.
#' @param TOL Tolerance of 0
#' @importFrom parallel mclapply
#' @importFrom Matrix rankMatrix
#' @keywords internal
#' @returns X
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
    poseignum <- min(
        Nmin+1,
        which(
            vals$XXt[1:Nmin]/(vals$XtX[1:Nmin]+.Machine$double.eps) < 0.99 |
                vals$XtX[1:Nmin]/(vals$XXt[1:Nmin]+.Machine$double.eps) < 0.99)
        ) - 1
    
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
        cat("Determinant:",
            det(vec), "\n")
        cat("Precision on v' = 1/v:",
            max(abs(t(vec) - solve(vec))), "\n")
        cat("Precision on Norm_col = 1:",
            max(abs(apply(vec, 2, function(x) norm(as.matrix(x), "2")) - 1)),
            "\n")
        cat("Precision on Norm_row = 1:",
            max(abs(apply(vec, 1, function(x) norm(as.matrix(x), "2")) - 1)),
            "\n")
        cat("Precision on Orthogonal:",
            max(sapply(1:(ncol(vec)-1), function(i) {
                max(sum(vec[i,] * vec[i+1,]), sum(vec[, i] * vec[, i+1]))
            })), "\n")
    }))
    
    ## solution S: X * r = vecB1 * E * S * vecB2' * r = Xr
    ## E * S * vecB2' * r = vecB1' * Xr = tmprhs1
    tmprhs1 <- crossprod(vecs[[1]], Xrb)
    if (poseignum < N1)
        cat("Precision on tmprhs1's zero:",
            max(abs(tmprhs1[(poseignum+1):N1, 1])), "\n")
    ## S * vecB2' * rmX2 = S * lhs1 = 1/E * tmprhs1 = rhs1
    E <- diag(sqrt(vals[[1]][1:poseignum]), ncol=poseignum, nrow=poseignum)
    invE <- diag(1/diag(E), ncol=poseignum, nrow=poseignum)
    rhs1 <- crossprod(t(invE), tmprhs1[1:poseignum, , drop=F])
    lhs1 <- crossprod(vecs[[2]], rb)
    signs1 <- rhs1[1:poseignum,]/lhs1[1:poseignum,]
    ## S = [signs1 0]
    S <- cbind(diag(signs1, ncol=poseignum, nrow=poseignum),
               matrix(0, nrow=poseignum, ncol=N2-poseignum))
    ## D = E %*% S
    D <- rbind(crossprod(t(E), S), matrix(0, nrow=N1-poseignum, ncol=N2))
    ## a = vecs[["A*A'"]] %*% D %*% t(vecs[["A'*A"]])
    a1 <- tcrossprod(tcrossprod(vecs[[1]], t(D)), vecs[[2]])
    
    cat("----------------------\n")
    cat("Precision on XXt = a1*a1':", max(abs(B1 - tcrossprod(a1))),
        " / (", quantile(abs(B1)), ")\n")
    cat("Precision on XtX = a1'*a1:", max(abs(B2 - crossprod(a1))),
        " / (", quantile(abs(B2)), ")\n")
    
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
#' @description Function for computing the federated SSCP matrix.
#' @param loginFD Login information of the FD server
#' @param logins Login information of data repositories
#' @param funcPreProc Definition of a function for preparation of raw data
#' matrices. Two arguments are required: conns (list of DSConnection-classes), 
#' symbol (name of the R symbol) (see datashield.assign).
#' @param querytables Vector of names of the R symbols to assign in the
#' DataSHIELD R session on each server in \code{logins}.
#' The assigned R variables will be used as the input raw data.
#' Other assigned R variables in \code{funcPreProc} are ignored.
#' @param byColumn A logical value indicating whether the input data is
#' centered by column or row. Default, TRUE, centering by column, with constant
#' variables across samples are removed. If FALSE, centering and scaling by
#' row, with constant samples across variables are removed.
#' @param width.cutoff Default, 500. See \code{deparse1}.
#' @param TOL Tolerance of 0.
#' @param connRes A logical value indicating if the connection to \code{logins}
#' is returned. Default, FALSE, connections are closed.
#' @importFrom parallel mclapply
#' @importFrom DSI datashield.aggregate datashield.assign datashield.logout
#' datashield.errors datashield.symbols
#' @keywords internal
.federateSSCP <- function(loginFD, logins, funcPreProc, querytables,
                          byColumn = TRUE, scale = FALSE,
                          chunk = 500, mc.cores = 1,
                          width.cutoff = 500L,
                          TOL = 1e-10,
                          connRes = FALSE) {
    loginFDdata <- .decode.arg(loginFD)
    logindata   <- .decode.arg(logins)
    opals <- .login(logins=logindata)
    nnode <- length(opals)
    .printTime(".federateSSCP Login-ed")
    
    tryCatch({
        ## take a snapshot of the current session
        safe.objs <- .ls.all()
        ## leave alone .Random.seed for sample()
        safe.objs[['.GlobalEnv']] <- setdiff(safe.objs[['.GlobalEnv']],
                                             '.Random.seed')
        ## lock everything so no objects can be changed
        .lock.unlock(safe.objs, lockBinding)
        
        ## apply funcPreProc for preparation of querytables on opals
        ## TODO: control hacking!
        ## TODO: control identical colnames!
        funcPreProc(conns=opals, symbol=querytables)
        
        ## unlock back everything
        .lock.unlock(safe.objs, unlockBinding)
        ## get rid of any sneaky objects that might have been created in the
        ## filters as side effects
        .cleanup(safe.objs)
    }, error=function(e) {
        .printTime(paste0("DATA MAKING PROCESS: ", e))
        return (paste0("DATA MAKING PROCESS: ", e,
                       ' --- ', datashield.symbols(opals),
                       ' --- ', datashield.errors(),
                       ' --- ', datashield.logout(opals)))
    })
    .printTime(".federateSSCP Input data processed")
    
    ## compute XX' on opals
    tryCatch({
        ## center data
        datashield.assign(opals,
                          "centeredData",
                          as.call(c(as.symbol("center"),
                                    x=as.call(c(as.symbol("list"),
                                                setNames(
                                                    lapply(querytables,
                                                           as.symbol),
                                                    querytables))),
                                    subset=NULL,
                                    byColumn=byColumn,
                                    scale=scale
                                    )),
                          async=T)
        
        ## samples
        samples <- lapply(
            datashield.aggregate(opals,
                                 as.symbol("rowNames(centeredData)"),
                                 async=T),
            function(x) x[[1]])
        
        ## variables
        variables <- datashield.aggregate(
            opals[1],
            as.symbol('colNames(centeredData)'),
            async=T)[[1]]
        
        ## compute XX'
        datashield.assign(opals, "tcrossProdSelf", 
                          as.call(list(as.symbol("tcrossProd"),
                                       x=as.symbol("centeredData"),
                                       y=NULL,
                                       chunk=chunk)),
                          async=T)
            
        ## connection from non-FD servers to FD-assigned server:
        ## user and password for login between servers are required
        loginFDdata$user     <- loginFDdata$userserver
        loginFDdata$password <- loginFDdata$passwordserver
        datashield.assign(opals, 'FD',
                          as.symbol(paste0("crossLogin('",
                                           .encode.arg(loginFDdata), "')")),
                          async=T)
    }, error=function(e) {
        .printTime(paste0("SSCP MAKING PROCESS: ", e))
        return (paste0("SSCP MAKING PROCESS: ", e,
                       ' --- ', datashield.symbols(opals),
                       ' --- ', datashield.errors(),
                       ' --- ', datashield.logout(opals)))
    })
    .printTime(".federateSSCP XX' data computed")
    
    ## send XX' from opals to FD
    tryCatch({
        ## send decomposed XX' to FD memory
        command <- list(as.symbol("pushToDscFD"),
                        as.symbol("FD"),
                        as.symbol("tcrossProdSelf"),
                        async=T)
        cat("Command: pushToDscFD(FD, tcrossProdSelf)", "\n")
        tcrossProdSelfDSC <- datashield.aggregate(opals,
                                                  as.call(command),
                                                  async=T)
        
        ## rebuild XX'
        tcrossProdSelf <- lapply(tcrossProdSelfDSC, function(dscblocks) {
            tcps <- lapply(dscblocks, function(dscblocki) {
                return (.rebuildMatrixDsc(dscblocki, mc.cores=mc.cores))
            })
            if (is.null(names(tcps))) names(tcps) <- querytables
            return (tcps)
        })
        gc(reset=F)
    }, error = function(e) {
        .printTime(paste0("SSCP PUSH PROCESS: ", e))
        datashield.assign(opals, 'crossEnd',
                          as.symbol("crossLogout(FD)"), async=T)
        return (paste0("SSCP PUSH PROCESS: ", e,
                       ' --- ', datashield.symbols(opals),
                       ' --- ', datashield.errors(),
                       ' --- ', datashield.logout(opals)))
    })
    .printTime(".federateSSCP Single XX' communicated to FD")
    
    if (nnode==1) {
        XXt <- tcrossProdSelf[[1]]
        for (tab in querytables) {
            rownames(XXt[[tab]]) <- colnames(XXt[[tab]]) <- samples[[1]]
        }
    } else {
        tryCatch({
            ## compute X'X
            datashield.assign(opals, "crossProdSelf", 
                              as.call(list(as.symbol("crossProd"),
                                           x=as.symbol("centeredData"),
                                           pair=F,
                                           chunk=chunk)),
                              async=T)
        
            ##- received by each from other nodes ----
            prodDataCross <- lapply(names(opals), function(opn) {
                ind.opn <- which(logindata$server == opn)
                logindata.opn <- logindata[-ind.opn, , drop=F]
                logindata.opn$user <- logindata.opn$userserver
                logindata.opn$password <- logindata.opn$passwordserver
                ## from each node opn, log in other nodes (mates)
                opals.loc <- paste0("crossLogin('",
                                    .encode.arg(logindata.opn),
                                    "')")
                datashield.assign(opals[opn], 'mates',
                                  as.symbol(opals.loc), async=F)
                tryCatch({
                    ## prepare raw data matrices on mates of opn
                    command.opn <- list(as.symbol("crossAssignFunc"),
                                        as.symbol("mates"),
                                        .encode.arg(funcPreProc),
                                        .encode.arg(querytables))
                    cat("Command: crossAssignFunc(mates, funcPreProc, ...)",
                        "\n")
                    invisible(datashield.aggregate(opals[opn],
                                                   as.call(command.opn),
                                                   async=F))
                    
                    ## center raw data on mates of opn
                    command.opn <- list(
                        as.symbol("crossAssign"),
                        as.symbol("mates"), 
                        symbol="centeredDataMate",
                        value=.encode.arg(deparse1(
                            as.call(c(as.symbol("center"),
                                      x=as.call(
                                          c(as.symbol("list"),
                                            setNames(
                                                lapply(querytables,
                                                       as.symbol),
                                                querytables))),
                                      byColumn=byColumn,
                                      scale=scale)),
                            width.cutoff=width.cutoff)),
                        value.call=T,
                        async=F
                    )
                    cat("Command: crossAssign(mates, centeredDataMate)", "\n")
                    invisible(datashield.aggregate(opals[opn],
                                                   as.call(command.opn),
                                                   async=F))
                    
                    ## create singularProd from mates of opn
                    command.opn <- paste0(
                        "crossAggregatePrimal(mates, '",
                        .encode.arg('singularProd(centeredDataMate)'),
                        "', async=T)")
                    cat("Command: crossAggregatePrimal(mates, 
                        singularProd(centeredDataMate), ...", "\n")
                    datashield.assign(opals[opn],
                                      "singularProdMate",
                                      as.symbol(command.opn), async=F)
                    
                    ## send X'X from opn to mates of opn
                    command.opn <- list(as.symbol("pushToDscMate"),
                                        as.symbol("mates"),
                                        as.symbol("crossProdSelf"),
                                        opn,
                                        async=F)
                    cat("Command: pushToDscMate(mates, crossProdSelf...", "\n")
                    invisible(datashield.aggregate(opals[opn],
                                                   as.call(command.opn),
                                                   async=F))
                    .printTime(paste0(
                        ".federateSSCP Pairwise X'X communicated: ",
                        opn))
                    
                    ## compute (X_i) * (X_j)' * (X_j) * (X_i)' on mates
                    command.opn <- list(
                        as.symbol("crossAssign"),
                        as.symbol("mates"),
                        symbol=paste0("prodDataCross_",opn),
                        value=.encode.arg(deparse1(
                            as.call(c(as.symbol("tripleProdChunk"),
                                      x=as.symbol("centeredDataMate"),
                                      mate=opn,
                                      chunk=chunk,
                                      mc.cores=mc.cores)),
                            width.cutoff=width.cutoff)),
                        value.call=T,
                        async=F
                    )
                    cat("Command: crossAssign(mates, prodDataCross,
                        tripleProdChunk(...", "\n")
                    invisible(datashield.aggregate(opals[opn],
                                                   as.call(command.opn),
                                                   async=F))
                    
                    ## login to FD from mates
                    command.opn <- paste0(
                        "crossAssign(mates, symbol='FDmate', value='",
                        .encode.arg(paste0("crossLogin('", loginFD, "')")),
                        "', value.call=T, async=F)")
                    cat("Command: crossAssign(mates, FDmate,
                        crossLogin(...", "\n")
                    invisible(datashield.aggregate(opals[opn],
                                                   as.symbol(command.opn),
                                                   async=F))
                    
                    ## push X_i * X_j * X_j' * X_i' to FD from mates
                    tryCatch({
                        command.opn <- paste0(
                            "crossAggregateDual(mates, '", 
                            .encode.arg(
                                paste0("pushToDscFD(FDmate, prodDataCross_",
                                       opn, ")")), 
                            "', async=T)")
                        cat("Command: crossAggregateDual(mates,
                            pushToDscFD(FDmate, prodDataCross_...", "\n")
                        prodDataCrossDSC <- datashield.aggregate(
                            opals[opn],
                            as.symbol(command.opn),
                            async=T)
                        
                        prodDataCross.opn <- lapply(
                            prodDataCrossDSC[[opn]],
                            function(dscblocks) {
                                pdcs <- lapply(dscblocks, function(dscblocki) {
                                    return (.rebuildMatrixDsc(
                                        dscblocki, mc.cores=mc.cores))
                                })
                                if (is.null(names(pdcs)))
                                    names(pdcs) <- querytables
                                return (pdcs)
                            })
                        .printTime(paste0(".federateSSCP XY'YX' tripleProd
                                          communicated to FD: ", opn))
                    }, error = function(e) {
                        .printTime(paste0("CROSS SSCP PUSH PROCESS: ", e))
                        return (paste0("CROSS SSCP PUSH PROCESS: ", e,
                                       ' --- ', datashield.symbols(opals[opn]),
                                       ' --- ', datashield.errors()))
                    }, finally = {
                        datashield.aggregate(
                            opals[opn], 
                            as.symbol(paste0(
                                "crossAssign(mates, symbol='crossFDEnd', ",
                                "value='",
                                .encode.arg("crossLogout(FDmate)"),
                                "', value.call=T, async=F)")),
                            async=T)
                    })
                    return (prodDataCross.opn)
                }, error = function(e) {
                    .printTime(paste0("CROSS PROCESS: ", e))
                    return (paste0("CROSS PROCESS: ", e,
                                   ' --- ', datashield.symbols(opals[opn]),
                                   ' --- ', datashield.errors()))
                }, finally = {
                    datashield.assign(opals[opn], 'crossEnd',
                                      as.symbol("crossLogout(mates)"),
                                      async=T)
                })
            })
            names(prodDataCross) <- names(opals)
            ##-----
            
            tryCatch({
                # command <- paste0("crossAggregate(FD, '",
                #                   .encode.arg(paste0("as.call(list(as.symbol('garbageCollect')", "))")),
                #                   "', async=T)")
                # cat("Command: ", command, "\n")
                # datashield.assign(opals, "GC", as.symbol(command), async=T)

                ## send the single-column matrix from opals to FD
                ## (X_i) * (X_j)' * ((X_j) * (X_j)')[,1]
                datashield.assign(opals, "singularProdCross",
                                  as.symbol('tcrossProd(centeredData,
                                            singularProdMate)'),
                                  async=T)
                command <- list(as.symbol("pushToDscFD"),
                                as.symbol("FD"),
                                as.symbol("singularProdCross"),
                                async=T)
                cat("Command: pushToDscFD(FD, singularProdCross)", "\n")
                singularProdCrossDSC <- datashield.aggregate(opals,
                                                             as.call(command),
                                                             async=T)
                singularProdCross <- lapply(singularProdCrossDSC,
                                            function(dscspc) {
                    lapply(dscspc, function(dscblocks) {
                        spcs <- lapply(dscblocks, function(dscblocki) {
                            return (.rebuildMatrixDsc(dscblocki,
                                                      mc.cores=mc.cores))
                        })
                        if (is.null(names(spcs))) names(spcs) <- querytables
                        return (spcs)
                    })
                })
                gc(reset=F)
                .printTime(".federateSSCP Ar communicated to FD")
            }, error = function(e) {
                .printTime(paste0("FD PROCESS MULTIPLE: ", e))
                datashield.assign(opals, 'crossEnd',
                                  as.symbol("crossLogout(FD)"), async=T)
                return (paste0("FD PROCESS MULTIPLE: ", e,
                               ' --- ', datashield.symbols(opals),
                               ' --- ', datashield.errors(),
                               ' --- ', datashield.logout(opals)))
            })
        }, error = function(e) {
            .printTime(paste0("XX' PROCESS MULTIPLE: ", e))
            return (paste0("XX' PROCESS MULTIPLE: ", e,
                           ' --- ', datashield.symbols(opals),
                           ' --- ', datashield.errors(),
                           ' --- ', datashield.logout(opals)))
        })
        
        ## deduced from received info by federation: (X_i) * (X_j)'
        crossProductPair <- lapply(querytables, function(tab) {
            cptab <- mclapply(
                1:(nnode-1),
                mc.cores=mc.cores,
                function(opi) {
                    crossi <- lapply((opi+1):(nnode), function(opj) {
                        opni <- names(opals)[opi]
                        opnj <- names(opals)[opj]
                        a1 <- .solveSSCP(
                            XXt=prodDataCross[[opnj]][[opni]][[tab]],
                            XtX=prodDataCross[[opni]][[opnj]][[tab]],
                            r=tcrossProdSelf[[opnj]][[tab]][, 1, drop=F],
                            Xr=singularProdCross[[opni]][[opnj]][[tab]],
                            TOL=TOL)
                        a2 <- .solveSSCP(
                            XXt=prodDataCross[[opni]][[opnj]][[tab]],
                            XtX=prodDataCross[[opnj]][[opni]][[tab]],
                            r=tcrossProdSelf[[opni]][[tab]][, 1, drop=F],
                            Xr=singularProdCross[[opnj]][[opni]][[tab]],
                            TOL=TOL)
                        cat("Precision on a1 = t(a2):", max(abs(a1 - t(a2))),
                            " / (", quantile(abs(a1)), ")\n")
                        return (a1)
                    })
                    names(crossi) <- names(opals)[(opi+1):(nnode)]
                    return (crossi)
                })
            names(cptab) <- names(opals)[1:(nnode-1)]
            return (cptab)
        })
        names(crossProductPair) <- querytables
        .printTime(".federateSSCP XY' computed")
        
        ## SSCP whole matrix
        XXt <- lapply(querytables, function(tab) {
            XXt.tab <- do.call(rbind, mclapply(
                1:nnode,
                mc.cores=mc.cores,
                function(opi) {
                    upper.opi <- do.call(cbind, as.list(
                        crossProductPair[[tab]][[names(opals)[opi]]]))
                    lower.opi <- do.call(cbind, lapply(
                        setdiff(1:opi, opi),
                        function(opj) {
                            t(crossProductPair[[tab]][[names(opals)[opj]]][[
                                names(opals)[opi]]])
                        }))
                    return (cbind(lower.opi,
                                  tcrossProdSelf[[opi]][[tab]],
                                  upper.opi))
                }))
            rownames(XXt.tab) <- colnames(XXt.tab) <- 
                unlist(samples, use.names=F)
    
            return (XXt.tab)
        })
        names(XXt) <- querytables
        .printTime(".federateSSCP Whole XX' computed")
        gc(reset=F)
    }
    
    if (!connRes) {
        datashield.assign(opals, 'crossEnd',
                          as.symbol("crossLogout(FD)"), async=T)
        datashield.logout(opals)
    }
    return (list(sscp=XXt,
                 samples=samples,
                 variables=variables,
                 conns=opals))
}


#' @title Federated ComDim
#' @description Function for ComDim federated analysis on the virtual cohort
#' combining multiple cohorts: finding common dimensions in multiblock data.
#' @usage federateComDim(loginFD,
#'                       logins,
#'                       func,
#'                       symbol,
#'                       ncomp = 2,
#'                       scale = FALSE,
#'                       scale.block = TRUE,
#'                       threshold = 1e-8,
#'                       chunk = 500,
#'                       mc.cores = 1,
#'                       width.cutoff = 500L
#'                       )
#' @param loginFD Login information of the FD server
#' @param logins Login information of data repositories
#' @param func Encoded definition of a function for preparation of raw data
#' matrices. Two arguments are required: conns (list of Opal connections), 
#' symbol (name of the R symbol) (see datashield.assign).
#' @param symbol Encoded vector of names of the R symbols to assign in the
#' DataSHIELD R session on each server in \code{logins}.
#' The assigned R variables will be used as the input raw data.
#' Other assigned R variables in \code{func} are ignored.
#' @param ncomp Number of common dimensions
#' @param scale A logical value indicating if variables are scaled to unit
#' variance. Default, FALSE. See \code{MBAnalysis::ComDim}.
#' @param scale.block A logical value indicating if each block of variables is
#' divided by the square root of its inertia. Default, TRUE.
#' See \code{MBAnalysis::ComDim}.
#' @param threshold if the difference of fit<threshold then break the iterative
#' loop (default 1e-8)
#' @param width.cutoff Default, 500. See \code{deparse1}.
#' @param chunk Size of chunks into what the resulting matrix is partitioned.
#' Default, 500.
#' @param mc.cores Number of cores for parallel computing. Default, 1.
#' @returns A \code{ComDim} object. See \code{MBAnalysis::ComDim}.
#' @importFrom DSI datashield.logout datashield.errors datashield.symbols
#' datashield.assign
#' @importFrom arrow write_to_raw
#' @export
federateComDim <- function(loginFD, logins, func, symbol,
                           ncomp = 2,
                           scale = FALSE,
                           scale.block = TRUE,
                           threshold = 1e-8,
                           chunk = 500, mc.cores = 1,
                           width.cutoff = 500L) {
    #require(DSOpal)
    .printTime("federateComDim started")
    TOL <- 1e-10
    funcPreProc <- .decode.arg(func)
    querytables <- .decode.arg(symbol)
    ntab <- length(querytables)
    
    if (!is.logical(scale)) stop("scale must be logical.")
    if (!is.logical(scale.block)) stop("scale.block must be logical.")
    
    ## compute SSCP matrix for each centered data table
    XX_query <- .federateSSCP(loginFD=loginFD, logins=logins, 
                              funcPreProc=funcPreProc, querytables=querytables,
                              byColumn=TRUE,
                              chunk=chunk, mc.cores=mc.cores,
                              width.cutoff=width.cutoff, TOL=TOL, 
                              connRes=T)
    XX <- XX_query$sscp
    querysamples <- XX_query$samples
    queryvariables <- XX_query$variables
    opals <- XX_query$conns
    nnode <- length(opals)
    
    ## function: compute the total variance from a XX' matrix
    inertie <- function(tab) {
        return (sum(diag(tab)))
    }
    
    ## function: compute the RV between WX and WY: 
    ## similarity between two matrices
    # coefficientRV <- function(WX, WY) {
    #     rv <- inertie(WX %*% WY)/(sqrt(inertie(WX %*% WX) *
    #                                        inertie(WY %*% WY)))
    #     return (rv)
    # }
    
    ## function: normalize a vector to a unit vector
    normv <- function(x) {
        normx <- norm(x, type='2')
        if (normx==0) normx <- 1
        return (x/normx)
    }
    
    ##### ComDim algorithm inherited from MBAnalysis #####
    
    ##- 0. Preliminary tests ----
    if (any(sapply(XX, is.na)))
        stop("No NA values are allowed.")
    if (!all(sapply(XX, isSymmetric, check.attributes=T)))
        stop("XX elements should be symmetric.")
    if (length(unique(sapply(XX, nrow))) != 1)
        stop("XX elements should be the same dimension.")
    ## number of samples
    # nind <- unique(unlist(apply(sapply(XX, dim), 1, unique)))
    # if (length(nind) > 1)
    #     stop("XX elements should be symmetric of the same dimension")
    # samples <- sapply(XX, function(x) intersect(rownames(x), colnames(x)))
    # if (is.list(samples) && max(lengths(samples))==0) {
    #     XX <- lapply(XX, function(x) {
    #         rownames(x) <- colnames(x) <- paste("Ind", 1:nind, sep='.')
    #         return (x)
    #     })
    #     samples <- sapply(XX, function(x) union(rownames(x), colnames(x)))
    # }
    # if (is.list(samples) || is.list(apply(samples, 1, unique)))
    #     stop("XX elements should have the same rownames and colnames")
    # 
    ##-----
    
    ##- 1. Output preparation ----
    nind <- nrow(XX[[1]])
    samples <- unlist(querysamples)
    if (any(samples!=rownames(XX[[1]])))
        stop("Wrong samples order.")
    nvar <- lengths(queryvariables)
    if (ncomp > 2 && min(sapply(XX, rankMatrix)) <= ncomp) {
        print(paste0("Security issue: maximum ",
                     min(sapply(XX, rankMatrix)) - 1,
                     " components could be inquired. ncomp will be set to 2."))
        ncomp <- 2
    }
    if (ncomp < 1) {
        print("ncomp should be at least 1. ncomp will be set to 2.")
        ncomp <- 2
    }
    compnames <- paste("Dim.", 1:ncomp, sep="")
    
    ## optimal value of the criterion
    optimalcrit  <- vector("numeric", ncomp)
    names(optimalcrit) <- compnames
    
    #contrib           <- matrix(0, nrow=ntab, ncol=ncomp)
    #dimnames(contrib) <- list(querytables, compnames)
    
    ## Specific weights for each dataset and each dimension
    LAMBDA <- #NNLAMBDA <-
        matrix(1, nrow=ntab, ncol=ncomp,
               dimnames=list(querytables, compnames))
    
    ## Normed global components
    Q <- matrix(0, nrow=nind, ncol=ncomp,
                dimnames=list(samples, compnames))

    ## Block components
    Q.b <- array(0, dim=c(nind, ncomp, ntab),
                 dimnames=list(samples, compnames, querytables))

    ## Correlations between global components and respective block components
    # cor.g.b <- array(0, dim=c(ncomp, ncomp, ntab),
    #                  dimnames=list(compnames, compnames, querytables))

    ## Weights for the block components
    # W.b      <- vector("list", length=ntab)
    # ## Correlation between the original variables of each block and its block
    # ## components (loadings)
    # blockcor <- vector("list", length=ntab)
    # for (k in 1:ntab) {
    #     W.b[[k]] <- blockcor[[k]] <- matrix(0, nrow=nvar[k], ncol=ncomp)
    #     dimnames(W.b[[k]]) <- dimnames(blockcor[[k]]) <-
    #         list(queryvariables[[k]], compnames)
    # }
    
    ## Weights for the block components
    ## Correlation between the original variables of each block and its block
    ## components (loadings)
    # W.b <- blockcor <- lapply(querytables, function(tab) {
    #     matrix(0, nrow=nvar[tab], ncol=ncomp,
    #            dimnames=list(queryvariables[[tab]], compnames))
    # })
    # names(W.b) <- querytables
    
    ## Weights for the global components
    ## Modified Weights to take account of deflation
    # W <- Wm <- matrix(0, nrow=sum(nvar), ncol=ncomp,
    #                         dimnames=list(do.call(c, queryvariables),
    #                                       compnames))
    
    ## Percentage of explained inertia of each Xb block and global
    explained.X  <- matrix(0, nrow=ntab+1, ncol=ncomp,
                           dimnames=list(c(querytables, 'Global'), compnames))
    cumexplained <- matrix(0, nrow=ncomp, ncol=2,
                           dimnames=list(compnames, c("%explX", "cum%explX")))
    
    ## Concatenated Xb block components at current dimension comp
    #tb           <- matrix(0, nrow=nind, ncol=ntab)
    #components   <- vector("numeric", length=2)
    
    ## Block results
    Block <- NULL
    ## Output
    res   <- NULL
    ##-----
    
    ##- 2. Required parameters and data preparation ----
    ## Pre-processing: block weighting to set each block inertia equal to 1
    if (scale.block) {
        inertia <- sapply(XX, inertie)
        XX <- lapply(querytables, function(tab) {
            XX[[tab]]/inertia[tab]
        })
        names(XX) <- querytables
        inertia0.sqrt <- sqrt(inertia)
    }
    XX0 <- XX
    IT.X <- sapply(XX, inertie) # Inertia of each block of variable
    IT.X <- c(IT.X, sum(IT.X[1:ntab]))
    
    ## Computation of the cross-product matrices among individuals
    ## (or association matrices)
    Trace <- IT.X[1:ntab]
    Itot  <- 0
    ##-----
    
    ##- 3. Computation of Q, LAMBDA and scores for the various dimensions ----
    ## Iterative computation of the various components
    for (comp in 1:ncomp)  {
        previousfit <- 0
        deltafit <- threshold + 1
        
        ## 3.1 Interatively compute the comp-th common component
        while (deltafit > threshold) {
            ## Weighted sum of XX'
            P <- Reduce("+", lapply(1:ntab, function(k)
                LAMBDA[k, comp] * XX[[k]]
            ))
            ## NB. svd is sensitive to noise (~1e-15) in P
            #reseig    <- eigen(P)
            #Q[, comp] <- reseig$vectors[, 1]
            ressvd    <- svd(P, nu=1, nv=1)
            Q[, comp] <- ressvd$u
            
            #optimalcrit[comp] <- reseig$values[1]
            LAMBDA[, comp] <- sapply(1:ntab, function(k) {
                crossprod(Q[, comp, drop=F],
                          crossprod(XX[[k]],
                                    Q[, comp, drop=F]))
            })
            fit <- sum(LAMBDA[, comp]^2)
            
            if (previousfit==0)
                deltafit <- threshold + 1
            else
                deltafit <- (fit - previousfit)/previousfit
            previousfit <- fit
        }
        optimalcrit[comp] <- sum(LAMBDA[, comp]^2)
        #LAMBDA[, comp]    <- normv(LAMBDA[, comp])
        
        ## 3.2 Storage of the results associated with dimension comp
        for (k in 1:ntab) {
            explained.X[k, comp] <- inertie(tcrossprod(Q[, comp]) %*% XX[[k]]
                                            %*% tcrossprod(Q[, comp]))
            Q.b[, comp, k] <- XX[[k]] %*% Q[, comp]
        }
        explained.X[ntab+1, comp] <- sum(explained.X[1:ntab, comp])
        
        # LAMBDA[, comp]   <- sapply(1:ntab, function(k) {
        #     t(Q[, comp]) %*% Q.b[, comp, k]
        # })
        
        ## Non normalized specific weights
        #NNLAMBDA[, comp] <- LAMBDA[, comp]
        #LAMBDA[, comp]   <- normv(LAMBDA[, comp])
        
        # 3.3 Deflation
        ## Deflation of XX
        proj <- diag(1, nind) - tcrossprod(Q[, comp])
        XX <- lapply(XX, function(xx) proj %*% xx %*% t(proj))
    }
    
    ## Explained inertia for X
    explained.X       <- sweep(explained.X, 1, IT.X, "/")
    cumexplained[, 1] <- explained.X[ntab+1, 1:ncomp]
    cumexplained[, 2] <- cumsum(cumexplained[, 1])
    
    ## Global components (scores of individuals)
    Scor.g <- tcrossprod(Q, sqrt(diag(colSums(LAMBDA), ncol = ncol(LAMBDA))))
    # Unnormed global components
    #LambdaMoyen <- apply(NNLAMBDA^2, 2, sum)
    #C <- Q %*% sqrt(diag(LambdaMoyen, ))
    ##-----
    
    ##- 4. Loadings ----
    size <- c(0, lengths(querysamples))
    Qlist <- setNames(lapply(2:length(size), function(i) {
        Qi <- Q[(cumsum(size)[i-1]+1):cumsum(size)[i], , drop=F]
        rownames(Qi) <- querysamples[[i-1]]
        return (Qi)
    }), names(opals))
    chunkList <- mclapply(Qlist, mc.cores=mc.cores, function(xx) {
        nblocksrow <- ceiling(nrow(xx)/chunk)
        sepblocksrow <- rep(ceiling(nrow(xx)/nblocksrow), nblocksrow-1)
        sepblocksrow <- c(sepblocksrow, nrow(xx) - sum(sepblocksrow))
        tcpblocks <- .partitionMatrix(xx, seprow=sepblocksrow, sepcol=ncomp)
        return (lapply(tcpblocks, function(tcpb) {
            return (lapply(tcpb, function(tcp) {
                return (.encode.arg(write_to_raw(tcp)))
            }))
        }))
    })
    tryCatch({
        ## send Qlist from FD to opals
        # TOCHECK: security on pushed data
        invisible(sapply(names(opals), function(opn) {
            lapply(1:length(chunkList[opn]), function(i) {
                lapply(1:length(chunkList[[i]]), function(j) {
                    lapply(1:length(chunkList[[i]][[j]]), function(k) {
                        datashield.assign(
                            opals[opn], paste(c('FD', i, j, k),
                                              collapse="__"),
                            as.call(list(as.symbol("matrix2DscMate"),
                                         chunkList[opn][[i]][[j]][[k]])),
                            async=T)
                    })
                })
            })
            datashield.assign(
                opals[opn],
                paste("pushed", 'FD', sep="_"),
                as.call(list(as.symbol("rebuildMatrixVar"),
                             symbol='FD',
                             len1=length(chunkList[opn]),
                             len2=lengths(chunkList[opn]),
                             len3=lapply(chunkList[opn], lengths),
                             querytables="common")),
                async=T)
            # gc on opals
            # command <- paste0("crossAggregate(FD, '",
            #                   .encode.arg(paste0("as.call(list(as.symbol('garbageCollect')", "))")),
            #                   "', async=T)")
            # cat("Command: ", command, "\n")
            # datashield.assign(opals, "GC", as.symbol(command), async=T)
        }))
        ## compute loadings X'*Qlist
        datashield.assign(
            opals,
            "loadings", 
            as.call(list(as.symbol("crossProd"),
                         x=as.symbol("centeredData"),
                         y=as.symbol("pushed_FD"),
                         chunk=chunk)),
            async=T)
        ## send loadings from opals to FD
        command <- list(as.symbol("pushToDscFD"),
                        as.symbol("FD"),
                        as.symbol("loadings"),
                        async=T)
        cat("Command: pushToDscFD(FD, loadings)", "\n")
        loadingsDSC <- datashield.aggregate(opals, as.call(command), async=T)
        loadingsLoc <- lapply(loadingsDSC, function(dscblocks) {
            cps <- lapply(names(dscblocks), function(dscn) {
                cpsi <- .rebuildMatrixDsc(dscblocks[[dscn]], mc.cores=mc.cores)
                colnames(cpsi) <- paste0("Comp.", 1:ncomp)
                rownames(cpsi) <- queryvariables[[gsub("__common" ,"", dscn)]]
                return (cpsi)
            })
            names(cps) <- gsub("__common" ,"", names(dscblocks))
            return (cps)
        })
        gc(reset=F)
        .printTime("federatedComDim Loadings communicated to FD")
    }, error = function(e) {
        .printTime(paste0("LOADING COMPUTATION PROCESS: ", e))
        return (paste0("LOADING COMPUTATION PROCESS: ", e,
                       ' --- ', datashield.symbols(opals),
                       ' --- ', datashield.errors()))
    }, finally = {
        datashield.assign(opals, 'crossEnd',
                          as.symbol("crossLogout(FD)"), async=T)
        datashield.logout(opals)
    })

    ## Weights for the block components
    W.b <- lapply(querytables, function(tab) {
        Reduce('+', lapply(loadingsLoc, function(ll) ll[[tab]]))
    })
    names(W.b) <- querytables
    if (scale.block) {
        W.b <- lapply(querytables, function(tab) {
            W.b[[tab]]/inertia0.sqrt[tab]
        })
    }
    names(W.b) <- querytables
    
    W.g[, comp] <- unlist(sapply(1:ntab, function(j) {
        LAMBDA[j, comp] * W.b[[j]][, comp]
    }))
    W.g[, comp] <- normv(W.g[, comp])
    # Wbk <- Reduce('+', unlist(mclapply(names(opals),
    #                                    mc.cores=1, function(opn) {
    #     expr <- list(as.symbol("loadings"),
    #                  as.symbol("centeredAllData"),
    #                  .encode.arg(Qlist[[opn]]))
    #     loadings <- datashield.aggregate(opals[opn], as.call(expr))
    #     return (loadings)
    # }), recursive = F))
    
    #colnames(Wbk) <- compnames
    #csnvar <- cumsum(nvar)
    # W.b <- mclapply(1:ntab, mc.cores=ntab, function(k) {
    #     if (option=="uniform") return (Wbk[ifelse(k==1, 1, csnvar[k-1]+1):csnvar[k], , drop=F]/inertia0.sqrt[k])
    #     return (Wbk[ifelse(k==1, 1, csnvar[k-1]+1):csnvar[k], , drop=F])
    # })
    
    ## Loadings for the global components: normed
    Load.g <- tcrossprod(do.call(rbind, W.b),
                         sqrt(diag(colSums(LAMBDA), ncol = ncol(LAMBDA))))
    Load.g <- t(t(Load.g)/colSums(Scor.g^2)) 
    
    ## Global weights: normed
    W.g <- do.call(rbind, mclapply(1:ntab, mc.cores=mc.cores, function(k)
        tcrossprod(W.b[[k]], diag(LAMBDA[k,], nrow=ncomp, ncol=ncomp))))
    W.g <- do.call(cbind, mclapply(1:ncomp, mc.cores=mc.cores, function(comp)
        normv(W.g[, comp])))
    colnames(W.g) <- compnames
    
    ## Global projection
    Proj.g <- W.g %*% solve(crossprod(Load.g, W.g),
                            tol=9.99999999999999e-301)
    ##-----
    
    #contrib <- t(t(NNLAMBDA)/colSums(NNLAMBDA))
    
    ## Weights that take into account the deflation procedure
    #Wm <- W %*% solve(crossprod(Px, W), tol=1e-150)
    
    # globalcor <- NA #cor(X00, C)
    # 
    # for (k in 1:ntab) {
    #     cor.g.b[, , k] <- cor(Q, Q.b[, , k])
    #     blockcor[[k]] <- NA #cor(X0[[k]], Q.b[, 1:ncomp, k])
    #     #if (is.null(rownames(blockcor[[k]]))) rownames(blockcor[[k]]) <- names(group[k])
    # }
    
    ##- 5. Results ----
    ## Global
    res$components   <- c(ncomp=ncomp)
    res$optimalcrit  <- optimalcrit[1:ncomp]
    res$saliences    <- LAMBDA[, 1:ncomp, drop=FALSE]
    res$T.g          <- Q
    res$Scor.g       <- Scor.g
    res$W.g          <- W.g
    res$Load.g       <- Load.g
    res$Proj.g       <- Proj.g
    res$explained.X  <- explained.X
    res$cumexplained <- cumexplained*100
    
    ## Blocks
    Block$T.b <- Q.b
    Block$W.b <- W.b
    res$Block <- Block
    res$call  <- match.call()
    class(res) <- c("ComDim", "list")
    ##-----
    
    return (res)
}


#' @title Federated ComDim
#' @description Function for ComDim federated analysis on the virtual cohort
#' combining multiple cohorts: finding common dimensions in multiblock data.
#' @usage federateComDim(loginFD,
#'                       logins,
#'                       func,
#'                       symbol,
#'                       ncomp = 2,
#'                       scale = "none",
#'                       option = "uniform",
#'                       chunk = 500,
#'                       mc.cores = 1,
#'                       width.cutoff = 500L,
#'                       threshold = 1e-10)
#' @param loginFD Login information of the FD server
#' @param logins Login information of data repositories
#' @param func Encoded definition of a function for preparation of raw data
#' matrices. Two arguments are required: conns (list of Opal connections), 
#' symbol (name of the R symbol) (see datashield.assign).
#' @param symbol Encoded vector of names of the R symbols to assign in the
#' DataSHIELD R session on each server in \code{logins}.
#' The assigned R variables will be used as the input raw data.
#' Other assigned R variables in \code{func} are ignored.
#' @param ncomp Number of common dimensions
#' @param scale  either value "none" / "sd" indicating the same scaling for
#' all tables or a vector of scaling ("none" / "sd") for each table. 
#' Only "none" is considered in this version.
#' @param option weighting of te tables \cr
#'        "none" : no weighting of the tables - (default) \cr
#'      "uniform": weighting to set the table at the same inertia \cr
#' @param chunk Size of chunks into what the resulting matrix is partitioned.
#' Default, 500.
#' @param mc.cores Number of cores for parallel computing. Default, 1.
#' @param threshold if the difference of fit<threshold then break the iterative
#' loop (default 1E-10)
#' @returns A \code{ComDim} object. See \code{MBAnalysis::ComDim}.
#' @importFrom DSI datashield.logout datashield.errors datashield.symbols
#' datashield.assign
#' @importFrom arrow write_to_raw
#' @keywords internal
federateComDim1 <- function(loginFD, logins, func, symbol,
                           ncomp = 2,
                           scale = "none", option = "uniform",
                           chunk = 500, mc.cores = 1,
                           width.cutoff = 500L, threshold = 1e-10) {
    #require(DSOpal)
    .printTime("federateComDim started")
    TOL <- 1e-10
    funcPreProc <- .decode.arg(func)
    querytables <- .decode.arg(symbol)
    ntab <- length(querytables)
    
    ## compute SSCP matrix for each centered data table
    XX_query <- .federateSSCP(loginFD=loginFD, logins=logins, 
                              funcPreProc=funcPreProc, querytables=querytables,
                              byColumn=TRUE,
                              chunk=chunk, mc.cores=mc.cores,
                              width.cutoff=width.cutoff, TOL=TOL, 
                              connRes=T)
    XX <- XX_query$sscp
    querysamples <- XX_query$samples
    queryvariables <- XX_query$variables
    opals <- XX_query$conns
    nnode <- length(opals)
    
    ## function: compute the total variance from a XX' matrix
    inertie <- function(tab) {
        return (sum(diag(tab,)))    #Froebenius norm
    }
    
    ## function: compute the RV between WX and WY: 
    ## similarity between two matrices
    coefficientRV <- function(WX, WY) {
        rv <- inertie(WX %*% WY)/(sqrt(inertie(WX %*% WX) *
                                           inertie(WY %*% WY)))
        return(rv)
    }
    
    ## function: normalize a vector to a unit vector
    normv <- function(x) {
        normx <- norm(x, '2')
        if (normx==0) normx <- 1
        return (x/normx)
    }
    
    ##### ComDim algorithm inherited from MBAnalysis #####
    
    ##- 0. Preliminary tests ----
    if (any(sapply(XX, is.na)))
        stop("No NA values are allowed.")
    if (!all(sapply(XX, isSymmetric, check.attributes=T)))
        stop("XX elements should be symmetric.")
    if (length(unique(sapply(XX, nrow))) != 1)
        stop("XX elements should be the same dimension.")
    ## number of samples
    # nind <- unique(unlist(apply(sapply(XX, dim), 1, unique)))
    # if (length(nind) > 1)
    #     stop("XX elements should be symmetric of the same dimension")
    # samples <- sapply(XX, function(x) intersect(rownames(x), colnames(x)))
    # if (is.list(samples) && max(lengths(samples))==0) {
    #     XX <- lapply(XX, function(x) {
    #         rownames(x) <- colnames(x) <- paste("Ind", 1:nind, sep='.')
    #         return (x)
    #     })
    #     samples <- sapply(XX, function(x) union(rownames(x), colnames(x)))
    # }
    # if (is.list(samples) || is.list(apply(samples, 1, unique)))
    #     stop("XX elements should have the same rownames and colnames")
    # 
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
    nind <- nrow(XX[[1]])
    nvar <- lengths(queryvariables)
    if (ncomp > 2 && min(sapply(XX, rankMatrix)) <= ncomp) {
        print(paste0("Security issue: maximum ",
                     min(sapply(XX, rankMatrix)) - 1,
                     " components could be inquired. ncomp will be set to 2."))
        ncomp <- 2
    }
    if (ncomp < 1) {
        print("ncomp should be at least 1. ncomp will be set to 2.")
        ncomp <- 2
    }
    compnames <- paste("Dim.", 1:ncomp, sep="")

    ## optimal value of the criterion
    optimalcrit  <- vector("numeric", ncomp)
    names(optimalcrit) <- compnames
    
    contrib           <- matrix(0, ntab, ncomp)
    dimnames(contrib) <- list(querytables, compnames)
    
    ## Specific weights for each dataset and each dimension
    saliences <- LAMBDA <- NNLAMBDA <- matrix(1, ntab, ncomp)
    dimnames(saliences) <- dimnames(LAMBDA) <-
        dimnames(NNLAMBDA) <- list(querytables, compnames)
    
    ## Global components
    Q <- matrix(0, nrow=nind, ncol=ncomp)
    ## Unnormed global components
    C <- matrix(0, nrow=nind, ncol=ncomp)
    dimnames(Q) <- dimnames(C) <- list(rownames(XX[[1]]), compnames)
    
    ## Block components
    Q.b <- array(0, dim=c(nind, ncomp, ntab))
    dimnames(Q.b) <- list(rownames(XX[[1]]), compnames, querytables)
    
    ## Correlations between global components and respective block components
    cor.g.b <- array(0, dim=c(ncomp, ncomp, ntab))
    dimnames(cor.g.b) <- list(compnames, compnames, querytables)
    
    ## Weights for the block components
    W.b      <- vector("list", length=ntab)
    ## Correlation between the original variables of each block and its block
    ## components (loadings)
    blockcor <- vector("list", length=ntab)
    for (k in 1:ntab) {
        W.b[[k]] <- blockcor[[k]] <- matrix(0, nrow=nvar[k], ncol=ncomp)
        dimnames(W.b[[k]]) <- dimnames(blockcor[[k]]) <-
            list(queryvariables[[k]], compnames)
    }
    
    ## Loadings for the global components
    Px           <- matrix(0, nrow=sum(nvar), ncol=ncomp)
    ## Weights for the global components
    W            <- Px
    ## Modified Weights to take account of deflation
    Wm           <- Px                                     
    dimnames(Px) <- dimnames(W) <- dimnames(Wm) <-
        list(do.call(c, queryvariables), compnames)
    
    IT.X         <- vector("numeric", length=ntab+1)
    ## Percentage of explained inertia of each Xb block and global
    explained.X  <- matrix(0, nrow=ntab+1, ncol=ncomp)
    dimnames(explained.X) <- list(c(querytables, 'Global'), compnames)
    cumexplained <- matrix(0, nrow=ncomp, ncol=2)
    dimnames(cumexplained) <- list(compnames, c("%explX", "cum%explX"))

    ## Concatenated Xb block components at current dimension comp
    tb           <- matrix(0, nrow=nind, ncol=ntab)
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
    
    ## Pre-processing: block weighting to set each block inertia equal to 1
    if (option=="uniform") {
        inertia <- sapply(XX, inertie)
        XX <- lapply(querytables, function(tab) {
            XX[[tab]]/inertia[tab]
        })
        names(XX) <- querytables
        inertia0.sqrt <- sqrt(inertia)
    }
    XX0 <- XX
    IT.X[1:ntab] <- sapply(XX, inertie) # Inertia of each block of variable
    IT.X[ntab+1] <- sum(IT.X[1:ntab])
    
    ## Computation of the cross-product matrices among individuals
    ## (or association matrices)
    Trace <- IT.X[1:ntab]
    Itot  <- 0
    ##-----
    
    ##- 3. computation of Q and LAMBDA for the various dimensions ----
    ## Iterative computation of the various components
    for (comp in 1:ncomp)  {
        critt     <- 0
        deltacrit <- 1
        
        ## 3.1 Interatively compute the comp-th common component
        while (deltacrit > threshold) {
            ## Weighted sum of XX'
            P <- Reduce("+", lapply(1:ntab, function(k)
                LAMBDA[k, comp] * XX[[k]])
            )
            reseig    <- eigen(P)
            Q[, comp] <- reseig$vectors[, 1]
            optimalcrit[comp] <- reseig$values[1]
            LAMBDA[, comp]    <- sapply(1:ntab, function(k) {
                t(Q[, comp]) %*% XX[[k]] %*% Q[, comp]
            })
            LAMBDA[, comp]    <- normv(LAMBDA[, comp])
            criterion <- reseig$values[1]
            deltacrit <- criterion - critt
            critt     <- criterion
        }
        
        ## 3.2 Storage of the results associated with dimension comp
        for (k in 1:ntab) {
            explained.X[k, comp] <- inertie(tcrossprod(Q[, comp]) %*% XX[[k]]
                                            %*% tcrossprod(Q[, comp]))
            #Q.b[, comp, k] <- XX[[k]] %*% Q[, comp]
        }
        
        Q.b[, comp,] <- sapply(1:ntab, function(k) {
            XX[[k]] %*% Q[, comp]
        })
        
        LAMBDA[, comp]   <- sapply(1:ntab, function(k) {
            t(Q[, comp]) %*% Q.b[, comp, k]
        })
        ## Non normalized specific weights
        NNLAMBDA[, comp] <- LAMBDA[, comp]
        LAMBDA[, comp]   <- normv(LAMBDA[, comp])
        
        ## 3.3 Deflation
        explained.X[ntab+1, comp] <- sum(explained.X[1:ntab, comp])
        proj <- diag(1, nind) - tcrossprod(Q[, comp])
        ## Deflation of XX
        XX <- lapply(XX, function(xx) proj %*% xx %*% t(proj))
    }
    
    ## Explained inertia for X
    explained.X       <- sweep(explained.X, 1, IT.X, "/")
    cumexplained[, 1] <- explained.X[ntab+1, 1:ncomp]
    cumexplained[, 2] <- cumsum(cumexplained[, 1])
    ##-----
    
    ##- 4. loadings ----
    size <- c(0, lengths(querysamples))
    Qlist <- setNames(lapply(2:length(size), function(i) {
        Qi <- Q[(cumsum(size)[i-1]+1):cumsum(size)[i], , drop=F]
        rownames(Qi) <- querysamples[[i-1]]
        return (Qi)
    }), names(opals))
    chunkList <- mclapply(Qlist, mc.cores=mc.cores, function(xx) {
        nblocksrow <- ceiling(nrow(xx)/chunk)
        sepblocksrow <- rep(ceiling(nrow(xx)/nblocksrow), nblocksrow-1)
        sepblocksrow <- c(sepblocksrow, nrow(xx) - sum(sepblocksrow))
        tcpblocks <- .partitionMatrix(xx, seprow=sepblocksrow, sepcol=ncomp)
        return (lapply(tcpblocks, function(tcpb) {
            return (lapply(tcpb, function(tcp) {
                return (.encode.arg(write_to_raw(tcp)))
            }))
        }))
    })
    tryCatch({
        ## send Q from FD to opals
        # TOCHECK: security on pushed data
        invisible(sapply(names(opals), function(opn) {
            lapply(1:length(chunkList[opn]), function(i) {
                lapply(1:length(chunkList[[i]]), function(j) {
                    lapply(1:length(chunkList[[i]][[j]]), function(k) {
                        datashield.assign(
                            opals[opn], paste(c('FD', i, j, k),
                                            collapse="__"),
                            as.call(list(as.symbol("matrix2DscMate"),
                                         chunkList[opn][[i]][[j]][[k]])),
                            async=T)
                    })
                })
            })
            datashield.assign(
                opals[opn],
                paste("pushed", 'FD', sep="_"),
                as.call(list(as.symbol("rebuildMatrixVar"),
                             symbol='FD',
                             len1=length(chunkList[opn]),
                             len2=lengths(chunkList[opn]),
                             len3=lapply(chunkList[opn], lengths),
                             querytables="common")),
                async=T)
            # gc on opals
            # command <- paste0("crossAggregate(FD, '",
            #                   .encode.arg(paste0("as.call(list(as.symbol('garbageCollect')", "))")),
            #                   "', async=T)")
            # cat("Command: ", command, "\n")
            # datashield.assign(opals, "GC", as.symbol(command), async=T)
        }))
        ## compute loadings X'*Qlist
        datashield.assign(
            opals,
            "loadings", 
            as.call(list(as.symbol("crossProd"),
                         x=as.symbol("centeredData"),
                         y=as.symbol("pushed_FD"),
                         chunk=chunk)),
            async=T)
        ## send loadings from opals to FD
        command <- list(as.symbol("pushToDscFD"),
                        as.symbol("FD"),
                        as.symbol("loadings"),
                        async=T)
        cat("Command: pushToDscFD(FD, loadings)", "\n")
        loadingsDSC <- datashield.aggregate(opals, as.call(command), async=T)
        loadingsLoc <- lapply(loadingsDSC, function(dscblocks) {
            cps <- lapply(names(dscblocks), function(dscn) {
                cpsi <- .rebuildMatrixDsc(dscblocks[[dscn]], mc.cores=mc.cores)
                colnames(cpsi) <- paste0("Comp.", 1:ncomp)
                rownames(cpsi) <- queryvariables[[gsub("__common" ,"", dscn)]]
                return (cpsi)
            })
            names(cps) <- gsub("__common" ,"", names(dscblocks))
            return (cps)
        })
        gc(reset=F)
        .printTime("federatedComDim Loadings communicated to FD")
    }, error = function(e) {
        .printTime(paste0("LOADING COMPUTATION PROCESS: ", e))
        return (paste0("LOADING COMPUTATION PROCESS: ", e,
                       ' --- ', datashield.symbols(opals),
                       ' --- ', datashield.errors()))
    }, finally = {
        datashield.assign(opals, 'crossEnd',
                          as.symbol("crossLogout(FD)"), async=T)
        datashield.logout(opals)
    })
    
    W.b <- lapply(querytables, function(tab) {
        Reduce('+', lapply(loadingsLoc, function(ll) ll[[tab]]))
    })
    names(W.b) <- querytables
    if (option=="uniform") {
        W.b <- lapply(querytables, function(tab) {
            W.b[[tab]]/inertia0.sqrt[tab]
        })
    }
    names(W.b) <- querytables
    
    # Wbk <- Reduce('+', unlist(mclapply(names(opals),
        #                                    mc.cores=1, function(opn) {
        #     expr <- list(as.symbol("loadings"),
        #                  as.symbol("centeredAllData"),
        #                  .encode.arg(Qlist[[opn]]))
        #     loadings <- datashield.aggregate(opals[opn], as.call(expr))
        #     return (loadings)
        # }), recursive = F))
    
    #colnames(Wbk) <- compnames
    #csnvar <- cumsum(nvar)
    # W.b <- mclapply(1:ntab, mc.cores=ntab, function(k) {
    #     if (option=="uniform") return (Wbk[ifelse(k==1, 1, csnvar[k-1]+1):csnvar[k], , drop=F]/inertia0.sqrt[k])
    #     return (Wbk[ifelse(k==1, 1, csnvar[k-1]+1):csnvar[k], , drop=F])
    # })
    
    Px <- do.call(rbind, W.b)
    W  <- do.call(rbind, lapply(1:ntab, function(k)
        tcrossprod(W.b[[k]], diag(LAMBDA[k,], nrow=ncomp, ncol=ncomp))))
    W  <- do.call(cbind, lapply(1:ncomp, function(comp)
        normv(W[, comp])))
    colnames(W) <- compnames
    
    contrib <- t(t(NNLAMBDA)/colSums(NNLAMBDA))
    
    ## Weights that take into account the deflation procedure
    Wm <- W %*% solve(crossprod(Px, W), tol=1e-150)
    
    
    # Unnormed global components
    LambdaMoyen <- apply(NNLAMBDA^2, 2, sum)
    C <- Q %*% sqrt(diag(LambdaMoyen, ))
    
    globalcor <- NA #cor(X00, C)
    
    for (k in 1:ntab) {
        cor.g.b[, , k] <- cor(Q, Q.b[, , k])
        blockcor[[k]] <- NA #cor(X0[[k]], Q.b[, 1:ncomp, k])
        #if (is.null(rownames(blockcor[[k]]))) rownames(blockcor[[k]]) <- names(group[k])
    }
    
    ## 4.1 Preparation of the results Global
    res$components          <- c(ncomp=ncomp)
    res$optimalcrit         <- optimalcrit[1:ncomp]
    res$saliences           <- round(LAMBDA[, 1:ncomp, drop=FALSE]^2, 2)
    ## Storage of the normed global components associated with X
    res$T                   <- Q[, 1:ncomp, drop=FALSE]
    ## Storage of the unnormed global components associated with X
    res$C                   <- C[, 1:ncomp, drop=FALSE]
    res$explained.X         <- round(100*explained.X[1:ntab, 1:ncomp], 2)
    res$cumexplained        <- round(100*cumexplained[1:ncomp,], 2)
    res$contrib             <- round(100*contrib[1:ntab, 1:ncomp], 2)
    res$globalcor           <- globalcor
    res$cor.g.b             <- cor.g.b

    ## 4.2 Preparation of the results Block
    Block$T.b         <-  Q.b[, 1:ncomp, ]
    Block$blockcor    <-  blockcor
    ## Results for each block
    res$Block         <-  Block
    ##-----

    ##- 5. Return res ----
    res$call   <- NA
    class(res) <- c("ComDim", "list")
    ##-----
    
    return (res)
}


#' @title Federated SNF
#' @description Function for SNF federated analysis on the virtual cohort
#' combining multiple cohorts.
#' @usage federateSNF(loginFD,
#'                    logins,
#'                    func,
#'                    symbol,
#'                    metric = 'euclidean',
#'                    K = 20,
#'                    sigma = 0.5,
#'                    t = 20,
#'                    chunk = 500,
#'                    mc.cores = 1,
#'                    width.cutoff = 500L)
#' @param loginFD Login information of the FD server
#' @param logins Login information of data repositories
#' @param func Encoded definition of a function for preparation of raw data
#' matrices. Two arguments are required: conns (list of DSConnection-classes), 
#' symbol (name of the R symbol) (see datashield.assign).
#' @param symbol Encoded vector of names of the R symbols to assign in the
#' DataSHIELD R session on each server in \code{logins}.
#' The assigned R variables will be used as the input raw data.
#' Other assigned R variables in \code{func} are ignored.
#' @param metric Either \code{euclidean} or \code{correlation} for distance
#' metric between samples. 
#' For Euclidean distance, the data from each cohort will be centered (not
#' scaled) for each variable.
#' For correlation-based distance, the data from each cohort will be centered
#' scaled for each sample.
#' @param K Number of neighbors in K-nearest neighbors part of the algorithm,
#' see \code{SNFtool::SNF}.
#' @param sigma Variance for local model, see \code{SNFtool::affinityMatrix}.
#' @param t Number of iterations for the diffusion process,
#' see \code{SNFtool::SNF}.
#' @param chunk Size of chunks into what the resulting matrix is partitioned.
#' Default, 500.
#' @param mc.cores Number of cores for parallel computing. Default, 1.
#' @param width.cutoff Default, 500. See \code{deparse1}.
#' @returns The overall status matrix derived W.
#' @importFrom SNFtool SNF affinityMatrix
#' @export
federateSNF <- function(loginFD, logins, func, symbol,
                        metric = "euclidean",
                        K = 20, sigma = 0.5, t = 20,
                        chunk = 500, mc.cores = 1,
                        width.cutoff = 500L) {
    #require(DSOpal)
    .printTime("federateSNF started")
    TOL <- 1e-10
    funcPreProc <- .decode.arg(func)
    querytables <- .decode.arg(symbol)
    ntab <- length(querytables)
    metric <- match.arg(metric, choices=c("euclidean", "correlation"))
    
    ## compute SSCP matrix for each centered data table
    XX_query <- .federateSSCP(loginFD=loginFD, logins=logins, 
                              funcPreProc=funcPreProc,
                              querytables=querytables,
                              byColumn=(metric=="euclidean"),
                              chunk=chunk, mc.cores=mc.cores,
                              width.cutoff=width.cutoff, TOL=TOL, 
                              connRes=F)
    XX <- lapply(1:ntab, function(i) {
        if (metric == "correlation") {
            ## compute (1 - correlation) distance between samples
            1 - XX_query$sscp[[i]]/(length(XX_query$variables[[i]])-1)
        } else if (metric == "euclidean") {
            ## compute Euclidean distance between samples
            .toEuclidean(XX_query$sscp[[i]])
        }
    })
    names(XX) <- querytables
    
    # if (metric == "correlation") {
    #     ## compute (1 - correlation) distance between samples for each data table 
    #     XX <- lapply(1:ntab, function(i) {
    #         xxi <- .federateSSCP(loginFD=loginFD, logins=logins, 
    #                              funcPreProc=funcPreProc, querytables=querytables, ind=i, 
    #                              byColumn=FALSE, chunk=chunk, mc.cores=mc.cores, TOL=TOL)
    #         return (1 - xxi$sscp/(length(xxi$var)-1))
    #     })
    # } else if (metric == "euclidean"){
    #     ## compute Euclidean distance between samples for each data table 
    #     XX <- lapply(1:ntab, function(i) {
    #         xxi <- .federateSSCP(loginFD=loginFD, logins=logins,
    #                              funcPreProc=funcPreProc, querytables=querytables, ind=i,
    #                              byColumn=TRUE, chunk=chunk, mc.cores=mc.cores, TOL=TOL)
    #         return (.toEuclidean(xxi$sscp))
    #     })
    # }
    
    ## take common samples
    commons <- Reduce(intersect, lapply(XX, rownames))
    XX <- lapply(XX, function(distmat) {
        distmat[commons, commons]
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
#' @description Function for UMAP federated analysis on the virtual cohort
#' combining multiple cohorts.
#' @usage federateUMAP(loginFD,
#'                     logins,
#'                     func,
#'                     symbol,
#'                     metric = 'euclidean',
#'                     chunk = 500,
#'                     mc.cores = 1,
#'                     width.cutoff = 500L,
#'                     ...)
#' @param loginFD Login information of the FD server
#' @param logins Login information of data repositories
#' @param func Encoded definition of a function for preparation of raw data
#' matrices. Two arguments are required: conns (list of DSConnection-classes), 
#' symbol (name of the R symbol) (see datashield.assign).
#' @param symbol Encoded vector of names of the R symbols to assign in the
#' DataSHIELD R session on each server in \code{logins}.
#' The assigned R variables will be used as the input raw data.
#' Other assigned R variables in \code{func} are ignored.
#' @param metric Either \code{euclidean} or \code{correlation} for distance
#' metric between samples. 
#' For Euclidean distance, the data from each cohort will be centered (not
#' scaled) for each variable.
#' For correlation-based distance, the data from each cohort will be centered
#' scaled for each sample.
#' @param chunk Size of chunks into what the resulting matrix is partitioned.
#' Default, 500.
#' @param mc.cores Number of cores for parallel computing. Default, 1.
#' @param width.cutoff Default, 500. See \code{deparse1}.
#' @param ... see \code{uwot::umap}
#' @returns A matrix of optimized coordinates.
#' @importFrom uwot umap
#' @importFrom stats as.dist
#' @export
federateUMAP <- function(loginFD, logins, func, symbol,
                         metric = "euclidean",
                         chunk = 500, mc.cores = 1,
                         width.cutoff = 500L, ...) {
    #require(DSOpal)
    .printTime("federateUMAP started")
    TOL <- 1e-10
    funcPreProc <- .decode.arg(func)
    querytables <- .decode.arg(symbol)
    ntab <- length(querytables)
    metric <- match.arg(metric, choices=c("euclidean", "correlation"))

    ## compute SSCP matrix for each centered data table
    XX_query <- .federateSSCP(loginFD=loginFD, logins=logins, 
                              funcPreProc=funcPreProc,
                              querytables=querytables,
                              byColumn=(metric=="euclidean"),
                              chunk=chunk, mc.cores=mc.cores,
                              width.cutoff=width.cutoff, TOL=TOL, 
                              connRes=F)
    XX <- lapply(1:ntab, function(i) {
        if (metric == "correlation") {
            ## compute (1 - correlation) distance between samples
            return (as.dist(
                1 - XX_query$sscp[[i]]/(length(XX_query$variables[[i]])-1)
                ))
        } else if (metric == "euclidean") {
            ## compute Euclidean distance between samples
            return (as.dist(
                .toEuclidean(XX_query$sscp[[i]])
                ))
        }
    })
    names(XX) <- querytables
    
    return (setNames(
        lapply(1:ntab, function(i) {
            uwot::umap(XX[[i]], ...)
            #ret_model = FALSE, ret_nn = FALSE, ret_extra = c(),
        }),
        querytables))
    
}


#' @title Federated hdbscan
#' @description Function for hdbscan federated analysis on the virtual cohort
#' combining multiple cohorts.
#' @usage federateHdbscan(loginFD,
#'                        logins,
#'                        func,
#'                        symbol,
#'                        metric = 'euclidean',
#'                        minPts = 10,
#'                        chunk = 500,
#'                        mc.cores = 1,
#'                        width.cutoff = 500L,
#'                        ...)
#' @param loginFD Login information of the FD server
#' @param logins Login information of data repositories
#' @param func Encoded definition of a function for preparation of raw data
#' matrices. Two arguments are required: conns (list of DSConnection-classes), 
#' symbol (name of the R symbol) (see datashield.assign).
#' @param symbol Encoded vector of names of the R symbols to assign in the
#' DataSHIELD R session on each server in \code{logins}.
#' The assigned R variables will be used as the input raw data.
#' Other assigned R variables in \code{func} are ignored.
#' @param metric Either \code{euclidean} or \code{correlation} for distance
#' metric between samples. 
#' For Euclidean distance, the data from each cohort will be centered (not
#' scaled) for each variable.
#' For correlation-based distance, the data from each cohort will be centered
#' scaled for each sample.
#' @param minPts Minimum size of clusters, see \code{dbscan::hdbscan}.
#' Default, 10.
#' @param chunk Size of chunks into what the resulting matrix is partitioned.
#' Default, 500.
#' @param mc.cores Number of cores for parallel computing. Default, 1.
#' @param width.cutoff Default, 500. See \code{deparse1}.
#' @param ... see \code{dbscan::hdbscan}
#' @returns An object of class \code{hdbscan}.
#' @importFrom dbscan hdbscan
#' @export
federateHdbscan <- function(loginFD, logins, func, symbol,
                            metric = "euclidean", minPts = 10, 
                            chunk = 500, mc.cores = 1,
                            width.cutoff = 500L, ...) {
    #require(DSOpal)
    .printTime("federateHdbscan started")
    TOL <- 1e-10
    funcPreProc <- .decode.arg(func)
    querytables <- .decode.arg(symbol)
    ntab <- length(querytables)
    metric <- match.arg(metric, choices=c("euclidean", "correlation"))
    
    ## compute SSCP matrix for each centered data table
    XX_query <- .federateSSCP(loginFD=loginFD, logins=logins, 
                              funcPreProc=funcPreProc,
                              querytables=querytables,
                              byColumn=(metric=="euclidean"),
                              chunk=chunk, mc.cores=mc.cores,
                              width.cutoff=width.cutoff, TOL=TOL, 
                              connRes=F)
    XX <- lapply(1:ntab, function(i) {
        if (metric == "correlation") {
            ## compute (1 - correlation) distance between samples
            return (as.dist(
                1 - XX_query$sscp[[i]]/(length(XX_query$variables[[i]])-1)
            ))
        } else if (metric == "euclidean") {
            ## compute Euclidean distance between samples
            return (as.dist(
                .toEuclidean(XX_query$sscp[[i]])
            ))
        }
    })
    names(XX) <- querytables
    
    return (setNames(
        lapply(1:ntab, function(i) {
            dbscan::hdbscan(XX[[i]], minPts = minPts, ...)
        }),
        querytables))
}


#' @title Test federateSSCP
#' @description Deprecated
#' @param loginFD Login information of the FD server
#' @param logins Login information of data repositories
#' @param func Encoded definition of a function for preparation of raw data
#' matrices. 
#' Two arguments are required: conns (list of DSConnection-classes), 
#' symbol (name of the R symbol) (see datashield.assign).
#' @param symbol Encoded vector of names of the R symbols to assign in the
#' DataSHIELD R session on each server in \code{logins}.
#' The assigned R variables will be used as the input raw data.
#' Other assigned R variables in \code{func} are ignored.
#' @param byColumn A logical value indicating whether the input data is centered by column or row.
#' Default, TRUE, centering by column. Constant variables across samples are removed. 
#' If FALSE, centering and scaling by row. Constant samples across variables are removed.
#' @param chunk Size of chunks into what the resulting matrix is partitioned. Default: 500
#' @param mc.cores Number of cores for parallel computing. Default: 1.
#' @param width.cutoff Default, 500. See \code{deparse1}.
#' @param TOL Tolerance of 0
#' @export
# #' @keywords internal
testSSCP <- function(loginFD, logins, func, symbol,
                     metric = "euclidean",
                     chunk = 500, mc.cores = 1, TOL = 1e-10,
                     width.cutoff = 500L, ...) {
    #require(DSOpal)
    .printTime("testSSCP started")
    TOL <- 1e-10
    funcPreProc <- .decode.arg(func)
    querytables <- .decode.arg(symbol)
    ntab <- length(querytables)
    metric <- match.arg(metric, choices=c("euclidean", "correlation"))
    
    ## compute SSCP matrix for each centered data table
    XX_query <- .federateSSCP(loginFD=loginFD, logins=logins, 
                              funcPreProc=funcPreProc,
                              querytables=querytables,
                              byColumn=(metric=="euclidean"),
                              chunk=chunk, mc.cores=mc.cores,
                              width.cutoff=width.cutoff, TOL=TOL, 
                              connRes=F)
    XX <- lapply(1:ntab, function(i) {
        if (metric == "correlation") {
            ## compute (1 - correlation) distance between samples
            return (as.dist(
                1 - XX_query$sscp[[i]]/(length(XX_query$variables[[i]])-1)
            ))
        } else if (metric == "euclidean") {
            ## compute Euclidean distance between samples
            return (as.dist(
                .toEuclidean(XX_query$sscp[[i]])
            ))
        }
    })
    names(XX) <- querytables
    
    # if (metric == "correlation") {
    #     ## compute (1 - correlation) distance between samples for each data table 
    #     XX <- lapply(1:ntab, function(i) {
    #         xxi <- .federateSSCP(loginFD=loginFD, logins=logins, 
    #                              funcPreProc=funcPreProc, querytables=querytables, ind=i, 
    #                              byColumn=FALSE, chunk=chunk, mc.cores=mc.cores, TOL=TOL)
    #         return (as.dist(1 - xxi$sscp/(length(xxi$var)-1)))
    #     })
    # } else if (metric == "euclidean"){
    #     ## compute Euclidean distance between samples for each data table 
    #     XX <- lapply(1:ntab, function(i) {
    #         xxi <- .federateSSCP(loginFD=loginFD, logins=logins, 
    #                              funcPreProc=funcPreProc, querytables=querytables, ind=i, 
    #                              byColumn=TRUE, chunk=chunk, mc.cores=mc.cores, TOL=TOL)
    #         return (as.dist(.toEuclidean(xxi$sscp)))
    #     })
    # }
    
    return (XX)
}
