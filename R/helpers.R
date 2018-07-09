#' @title Get pooled-member local statistics
#' @description Retrieve statistics for each grid point by pooling together all members
#' @param grid A climate4R CDM structure
#' @param fun Character string. Currently in use \code{"mean"} or \code{"sd"}, but any other valid function may work.
#' @param dim.along Dimension to pool along. Default to \code{"latlon"}. Otherwise \code{"time"}
#' @param ... Further args passed to \code{fun}
#' @return A vector, in the appropriate coordinates ordering (consistent with \code{\link[transformeR]{array3Dto2Dmat}} in the case of
#' spatial aggregations or \code{\link[transformeR]{getRefDates}} in the case of time aggregations).
#' @keywords internal
#' @importFrom abind abind
#' @importFrom magrittr %>% 
#' @importFrom transformeR getDim array3Dto2Dmat

getPooledMemberStat <- function(grid, fun, dim.along = "latlon", ...) {
    stopifnot(is.character(fun))
    stopifnot(is.function(eval(mean)))
    dim.along <- match.arg(dim.along, choices = c("latlon", "time"))
    if (dim.along == "latlon") {
        mar <- grep(c("lat|lon"), getDim(grid))
        s <- apply(grid$Data, MARGIN = mar, FUN = fun, ...) %>% abind(along = -1L) %>% unname()
        attr(s, "dimensions") <- c("time", "lat", "lon")
        array3Dto2Dmat(s)[1,]    
    } else {
        mar <- grep("time", getDim(grid))
        apply(grid$Data, MARGIN = mar, FUN = fun, ...)
    }
}


#' @importFrom magrittr %<>% 
#' @import transformeR
#' @keywords internal

prepareCalInputs <- function(x, y, newdata, crossval) {
    stopifnot(isGrid(y) | isGrid(x))
    stopifnot(is.logical(crossval))
    x %<>% redim(x, member = TRUE)
    y %<>% redim(y, member = TRUE)
    if (typeofGrid(y) == "station") {
        checkDim(x, y, dimensions = "time")    
        x <- suppressMessages(interpGrid(x, y))
    } else {
        checkDim(x, y, dimensions = c("time", "lat", "lon"))    
    }
    checkSeason(x, y)
    if (!is.null(newdata)) {
        crossval <- FALSE
        stopifnot(isGrid(newdata))
        checkDim(x, newdata, dimensions = c("member", "lat", "Lon"))
    } else {
        newdata <- x 
    }
    y %<>% redim(drop = TRUE) %>% redim(member = FALSE)
    return(list(x = x, y = y, newdata = newdata, crossval = crossval))
}


#' @importFrom stats cor.test
#' @keywords internal

corfun <- function(ens, obs, ...) {
    ctest <- cor.test(x = obs, y = rowMeans(ens), ...)
    return(list("pval" = ctest$p.value, "coef" = ctest$estimate))
}

