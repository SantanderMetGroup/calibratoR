#' @title Get pooled-member local statistics
#' @description Retrieve statistics for each grid point by pooling together all members
#' @param grid A climate4R CDM structure
#' @param fun Character string. Currently in use \code{"mean"} or \code{"sd"}, but any other valid function may work.
#' @param ... Further args passed to \code{fun}
#' @return A vector, in the appropriate coordinates ordering (consistent with \code{\link{array3Dto2Dmat}})
#' @keywords internal
#' @importFrom abind abind
#' @importFrom magrittr %>% 
#' @importFrom transformeR getDim array3Dto2Dmat

getPooledMemberStat <- function(grid, fun, ...) {
    stopifnot(is.character(fun))
    stopifnot(is.function(eval(mean)))
    mar <- grep(c("lat|lon"), getDim(grid))
    s <- apply(grid$Data, MARGIN = mar, FUN = fun, ...) %>% abind(along = -1L) %>% unname()
    attr(s, "dimensions") <- c("time", "lat", "lon")
    array3Dto2Dmat(s)[1,]    
}

