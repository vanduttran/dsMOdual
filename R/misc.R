#' @title Encode function argument
#' @description Serialize to JSON, then encode base64,
#'  then replace '+', '/' and '=' in the result in order to play nicely with the opal entry.
#'  Used to encode non-scalar function arguments prior to sending to the opal server.
#'  There's a corresponding function in the server package called .decode.arg.
#'  See \code{dsSwissKnifeClient:::.encode.arg}.
#' @param some.object the object to be encoded
#' @return encoded text with offending characters replaced by strings
#' @keywords internal
.encode.arg <- function(some.object, serialize.it = FALSE){
    if(serialize.it){
        encoded <- paste0(RCurl::base64Encode(jsonlite::serializeJSON(some.object)), 'serialized')
    } else {
        encoded <- RCurl::base64Encode(jsonlite::toJSON(some.object, null = 'null'))
    }
    # go fishing for '+', '/' and '=', opal rejects them :
    my.dictionary <- c('\\/' = '-slash-', '\\+' = '-plus-', '\\=' = '-equals-')
    sapply(names(my.dictionary), function(x){
        encoded[1] <<- gsub(x, my.dictionary[x], encoded[1])
    })
    return(paste0(encoded[1],'base64'))
}


#' @title Decode from base64 and deserialize from json if necessary
#' @description Work around the restrictions imposed by the Opal server on function arguments
#' The Opal server is very picky as regards function arguments. The workaround is
#' to serialize and encode them on the client and strip the right padding.
#' See \code{dsSwissKnife:::.decode.arg}.
#' @details It looks for the string 'base64' in the argument to determine if it's encoded
#' @param some.thing the thing to be decoded and deserialized from json if necessary
#' @return the decoded and deserialized argument
#' @keywords internal
.decode.arg <- function(some.thing, simplifyMatrix = FALSE){
    if(length(some.thing) == 1 && grepl('base64', some.thing, ignore.case = TRUE)){
        some.thing <- gsub('base64', '', some.thing, ignore.case =TRUE)
        serialized <- FALSE
        if(grepl('serialized', some.thing, ignore.case = TRUE)){
            serialized <- TRUE
            some.thing <- gsub('serialized', '', some.thing, ignore.case =TRUE)
        }
        my.dictionary = c('-plus-' = '+', '-slash-' = '/', '-equals-' = '=')
        sapply(names(my.dictionary), function(x){
            some.thing <<- gsub(x, my.dictionary[x], some.thing)
        })
        #
        if(serialized){
            some.thing <- jsonlite::unserializeJSON(RCurl::base64Decode(some.thing))
        } else {
            some.thing <- jsonlite::fromJSON(RCurl::base64Decode(some.thing), simplifyMatrix = simplifyMatrix)
        }
    }
    return(some.thing)
}


#' @title List all objects in all environments
#' @description \code{ls.all} returns all objects in all environments
#' See \code{dsSwissKnife:::.ls.all}.
#' @param start a character the environment name where to start (default .GlobalEnv)
#' @return  a list of environment names and the respective objects defined in each environment
#' @keywords internal
.ls.all <- function(start = '.GlobalEnv'){
    envir <- get(start)
    objs <- ls(envir, all.names = TRUE)
    ret <- list()
    ret[[start]] <- objs
    more.envs <- names(which(sapply(objs, function(x) is.environment(get(x)))==TRUE))
    c(ret,sapply(more.envs,function(x) ls(get(x), all.names = TRUE), USE.NAMES = TRUE, simplify = FALSE))
    
}


#' @title Locks or unlocks bindings in environments
#' @description Locks or unlocks bindings in environments
#' See \code{dsSwissKnife:::.lock.unlock}.
#' @param what a list of  environments and their respective objects - the output of ls.all above
#' @param func a function, either lockBinding or unlockBinding
#' @keywords internal
.lock.unlock <- function(what, func){
    stopifnot(deparse(substitute(func)) %in% c('lockBinding', 'unlockBinding'))
    invisible(lapply(names(what), function(x){
        lapply(what[[x]], function(y){
            func(y,get(x))
        })
    }))
}


#' @title Remove objects from the current workspace
#' @description This function removes objects from the current workspace
#' See \code{dsSwissKnife:::.cleanup}.
#' @param what a list of  environments and their respective objects - the output of a previous call to ls.all
#' @param start a character the environment name where to start (default .GlobalEnv)
#' @keywords internal
.cleanup <- function(initial, start = '.GlobalEnv'){
    objs <- .ls.all(start)
    new.envs <- setdiff(names(objs), names(initial))
    Map(function(x){
        rm(get(x))
        objs[x] <- NULL
    }, new.envs)
    invisible(Map(function(x){
        new.objs <- setdiff(objs[[x]], initial[[x]])
        rm(list = new.objs, pos = get(x))
    }, names(objs)))
}
