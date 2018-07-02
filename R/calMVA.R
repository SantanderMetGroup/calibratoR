##     calMVA.R Calibration of monthly/seasonal forecasts
##
##     Copyright (C) 2018 Santander Meteorology Group (http://www.meteo.unican.es)
##
##     This program is free software: you can redistribute it and/or modify
##     it under the terms of the GNU General Public License as published by
##     the Free Software Foundation, either version 3 of the License, or
##     (at your option) any later version.
## 
##     This program is distributed in the hope that it will be useful,
##     but WITHOUT ANY WARRANTY; without even the implied warranty of
##     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##     GNU General Public License for more details.
## 
##     You should have received a copy of the GNU General Public License
##     along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' @title Calibration of seasonal climate forecasts.
#' @description This function implements the method 1 described in Torralba et al. 2017, which has been previously applied 
#' in the context of seasonal forecasting to correct temperature and precipitation in Leung et al. 1999.
#' The ensemble mean produced by this method has the same mean and standard deviation as the observations. 
#' It is based on the assumption that both the reference and predicted distributions are approximated well by a Gaussian (normal) distribution.
#' @return climate4R grid. Calibrated forecasts.
#' @import transformeR 
#' @importFrom magrittr %>% extract2 %<>% 
#' @importFrom stats sd
#' @export
#' @template templateInputPars
#' @template templateParallelParams
#' @template templateParallel
#' @references
#' \itemize{
#' \item Torralba, V., F.J. Doblas-Reyes, D. MacLeod, I. Christel, and M. Davis, 2017: Seasonal Climate Prediction: A New Source of Information for the Management of Wind Energy Resources. J. Appl. Meteor. Climatol., 56, 1231–1247, https://doi.org/10.1175/JAMC-D-16-0204.1 
#' \item Leung, L.R., A.F. Hamlet, D.P. Lettenmaier, and A. Kumar, 1999: Simulations of the ENSO hydroclimate signals in the Pacific Northwest Columbia River basin. Bull. Amer. Meteor. Soc., 80, 2313–2329, doi:10.1175/1520-0477(1999)080,2313: SOTEHS.2.0.CO;2
#' }
#' @family calibration
#' @author R. Manzanas and V. Torralba.
#' @examples \dontrun{
#' ## loading seasonal forecasts (CFS) and observations (NCEP) of boreal winter temperature over Iberia
#' data("CFS_Iberia_tas"); fcst = CFS_Iberia_tas
#' data("NCEP_Iberia_tas"); obs = NCEP_Iberia_tas
#' ## passing from daily data to seasonal averages
#' fcst = aggregateGrid(fcst, aggr.y = list(FUN = "mean", na.rm = TRUE))
#' obs = aggregateGrid(obs, aggr.y = list(FUN = "mean", na.rm = TRUE))
#' ## interpolating forecasts to the observations' resolution
#' fcst = interpGrid(fcst, new.coordinates = getGrid(obs))
#' ## applying calibration
#' fcst.cal = calMVA(fcst, obs, crossval = TRUE)
#' ## plotting climatologies
#' library(visualizeR)
#' spatialPlot(makeMultiGrid(climatology(obs),
#'                          climatology(fcst, by.member = FALSE),
#'                          climatology(fcst.cal, by.member = FALSE)),
#'            backdrop.theme = "coastline",
#'            layout = c(3, 1),
#'            names.attr = c("NCEP", "CFS (raw)", "CFS (calibrated)"))
#' }

# data("cfs.forecast")
# data("cfs.hindcast")
# data("ncep.ref")
# newdata <- cfs.forecast
# x <- cfs.hindcast     
# y <- ncep.ref
# a <- calMVA(y = ncep.ref, x = cfs.hindcast, newdata = cfs.forecast)



calMVA <- function(y, x, newdata = NULL, crossval = TRUE, parallel = FALSE, max.ncores = 16, ncores = NULL) {
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
    ppars <- parallelCheck(parallel, max.ncores, ncores)
    lapplyfun <- selectPar.pplyFun(parallel.pars = ppars, .pplyFUN = "lapply")
    if (ppars$hasparallel) on.exit(parallel::stopCluster(ppars$cl))
    n.mem <- getShape(x, dimension = "member")
    coords <- getCoordinates(x)
    x.coords <- coords$x
    y.coords <- coords$y
    if (crossval) {
        yrs <- getYearsAsINDEX(x) %>% unique()
        l1 <- lapplyfun(yrs, function(j) {
            yr.ind <- setdiff(yrs, j)
            a <- subsetGrid(x, years = yr.ind)
            clim.fcst <- getPooledMemberStat(a, fun = "mean", na.rm = TRUE)
            sigma.e <- getPooledMemberStat(a, fun = "sd", na.rm = TRUE)
            a <- subsetGrid(y, years = yr.ind)
            clim.obs <- getPooledMemberStat(a, fun = "mean", na.rm = TRUE)
            sigma.ref <- getPooledMemberStat(a, fun = "sd", na.rm = TRUE)
            l <- lapply(1:n.mem, function(i) {
                aux <- subsetGrid(x, members = i) %>% redim(member = FALSE)
                x.train <- subsetGrid(aux, years = j) %>% redim(member = FALSE)
                aux2 <- x.train %>% extract2("Data") %>% array3Dto2Dmat() %>% t()
                x.train$Data <- ((aux2 - clim.fcst) * (sigma.ref/sigma.e) + clim.obs) %>% t() %>% mat2Dto3Darray(x = x.coords, y = y.coords)
                return(x.train)
            })
            do.call("bindGrid", c(l, dimension = "member")) %>% return()
        })
        cal <- do.call("bindGrid", c(l1, dimension = "time")) %>% redim(drop = TRUE)
    } else {
        clim.fcst <- getPooledMemberStat(x, fun = "mean", na.rm = TRUE)
        sigma.e <- getPooledMemberStat(x, fun = "sd", na.rm = TRUE)
        clim.obs <- getPooledMemberStat(y, fun = "mean", na.rm = TRUE)
        sigma.ref <- getPooledMemberStat(y, fun = "sd", na.rm = TRUE)
        l <- lapplyfun(1:n.mem, function(i) {
            x.train <- subsetGrid(newdata, members = i) %>% redim(member = FALSE) 
            aux2 <- x.train %>% extract2("Data") %>% array3Dto2Dmat() %>% t()
            x.train$Data <- ((aux2 - clim.fcst) * (sigma.ref/sigma.e) + clim.obs) %>% t() %>% mat2Dto3Darray(x = x.coords, y = y.coords)
            return(x.train)
        })
        cal <- do.call("bindGrid", c(l, dimension = "member")) %>% return()
    }
    attr(cal$Variable, "correction") <- "MVA"
    attr(cal$Variable, "correction_descr") <- "Mean-Variance Adjustment"
    invisible(cal)
}
    

# 
# library(visualizeR)
# 
# library(microbenchmark)
# out <- microbenchmark(
# par = calMVA(y = ncep.ref, x = cfs.hindcast),
# ser = calDebias(obs.grid = ncep.ref, fcst.grid = cfs.hindcast, crossval = FALSE), ntimes = 5
# )
# 
# boxplot(out)
# 
# a <- calMVA(y = ncep.ref, x = cfs.hindcast)
# 
# a %>% climatology(by.member = FALSE)  %>% spatialPlot()
# b %>% climatology(by.member = FALSE)  %>% spatialPlot()
# cfs.hindcast %>% climatology(by.member = FALSE)  %>% spatialPlot()
 # b <- calDebias(cfs.hindcast, obs = ncep.obs, crossval = T)
# 
 # b %>% climatology() %>% spatialPlot()
# a %>% climatology() %>% spatialPlot()
# err <- gridArithmetics(a, b, operator = "-")
# 
# 


calDebias <- function(fcst.grid, obs.grid, crossval = TRUE) {
    .Deprecated(new = "calMVA", old = "calDebias")
    ## Method 1 in Torralba et al. 2017: http://www.bsc.es/ESS/sites/default/files/imce/amspaper_final.pdf
    
    fcst = fcst.grid$Data
    obs = obs.grid$Data
    
    stopifnot(identical(dim(fcst)[-1], dim(obs)))  # check for equality of dimensions between fcst and obs
    
    nmemb = getShape(fcst.grid, "member")
    ntimes = getShape(fcst.grid, "time")
    nlat = getShape(fcst.grid, "lat")
    nlon = getShape(fcst.grid, "lon")
    
    fcst.cal = NA*fcst
    for (ilat in 1:nlat) {
        if (!(ilat/10) - trunc(ilat/10)) {
            message(sprintf("... lat %d of %d ...", ilat, nlat))
        }
        for (ilon in 1:nlon) {
            if (crossval) {
                ## leave-one-out cross-validation
                aux = sapply(1:ntimes, function(x) {
                    fcst.train = fcst[,-x,ilat,ilon]
                    fcst.test = fcst[,x,ilat,ilon]
                    obs.train = obs[-x,ilat,ilon]
                    clim.obs = mean(obs.train, na.rm = T)
                    clim.fcst = mean(fcst.train, na.rm = T)
                    sigma.e = sd(fcst.train, na.rm = T)
                    sigma.ref = sd(obs.train, na.rm = T)
                    
                    ((fcst.test - clim.fcst) * (sigma.ref/sigma.e)) + clim.obs
                })
                fcst.cal[,,ilat,ilon] = aux; rm(aux)
            } else {
                fcst.train = fcst[,,ilat,ilon]
                fcst.test = fcst.train
                obs.train = obs[,ilat,ilon]
                clim.obs = mean(obs.train, na.rm = T)
                clim.fcst = mean(fcst.train, na.rm = T)
                sigma.e = sd(fcst.train , na.rm = T)
                sigma.ref = sd(obs.train, na.rm = T)
                
                fcst.cal[,,ilat,ilon] = ((fcst.test - clim.fcst)*(sigma.ref/sigma.e)) + clim.obs
            }
        }
    }
    fcst.out = fcst.grid
    fcst.out$Data = fcst.cal
    attributes(fcst.out$Data)$dimensions = attributes(fcst.grid$Data)$dimensions
    return(fcst.out)
}

