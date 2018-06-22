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
#' @param fcst.grid climate4R grid. Forecasts to be calibrated (typically on a monthly/seasonal basis). At the moment, only gridded data are supported.
#' @param obs.grid climate4R grid. Reference observations the forecasts are calibrated towards (typically on a monthly/seasonal basis).
#' @param crossval Logical. TRUE (default) for leave-one-out cross-validation. FALSE for not cross-validation.
#' @return climate4R grid. Calibrated forecasts.
#' @importFrom transformeR getShape
#' @importFrom stats sd
#' @export
#' @references
#' \itemize{
#' \item Torralba, V., F.J. Doblas-Reyes, D. MacLeod, I. Christel, and M. Davis, 2017: Seasonal Climate Prediction: A New Source of Information for the Management of Wind Energy Resources. J. Appl. Meteor. Climatol., 56, 1231–1247, https://doi.org/10.1175/JAMC-D-16-0204.1 
#' \item Leung, L.R., A.F. Hamlet, D.P. Lettenmaier, and A. Kumar, 1999: Simulations of the ENSO hydroclimate signals in the Pacific Northwest Columbia River basin. Bull. Amer. Meteor. Soc., 80, 2313–2329, doi:10.1175/1520-0477(1999)080,2313: SOTEHS.2.0.CO;2
#' }
#' @seealso \code{\link{calInflation}}, \code{\link{calLM}}, \code{\link{calRPC}}, \code{\link{calNGR}}
#' @author R. Manzanas and V. Torralba.
#' @examples{
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

calMVA <- function(fcst.grid, obs.grid, crossval = TRUE) {
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

