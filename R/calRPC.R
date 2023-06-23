##     calRPC.R Calibration of monthly/seasonal forecasts
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

#' @title Calibration using the Ratio of Predictable Components for seasonal climate forecasts.
#' @description This function implements the EMOS method described in Eade et. 2014.
#' It uses the ensemble to reduce noise and adjust the forecast variance so that the ratio of predictable components (RPC) in the model and in the observations is the same. 
#' In Eade et al. 2014, this method was used to adjust seasonal forecasts of the North Atlantic Oscillation (NAO), temperature and pressure in the North Atlantic region.
#' @note Ensemble Model Output Statistics (EMOS) methods use the correspondence between the ensemble mean and the observations in the calibration process.
#' @param fcst.grid climate4R grid. Forecasts to be calibrated (typically on a monthly/seasonal basis). At the moment, only gridded data are supported.
#' @param obs.grid climate4R grid. Reference observations the forecasts are calibrated towards (typically on a monthly/seasonal basis).
#' @param crossval Logical. TRUE (default) for leave-one-out-out cross-validation. FALSE for not cross-validation.		
#' @param apply.to Character. If \code{"all"} is selected, all forecasts are calibrated. Alternatively, if \code{"sig"} is selected, the calibration is only applied 
#' in those points where the correlation between the ensemble mean and the observations is statistically significant.
#' @param alpha Significance (0.1 by default) of the ensemble mean correlation (i.e. \code{"alpha = 0.05"} would correspond to a 95\% confidence level). Only works if \code{"apply.to = sig"}.  
#' @return climate4R grid. Calibrated forecasts.					
#' @importFrom transformeR getShape
#' @importFrom stats cor.test sd var
#' @export
#' @references
#' \itemize{
#' \item Eade R., D. Smith, A. Scaife, et al, 2014: Do seasonal-to-decadal climate predictions underestimate the predictability of the real world? Geophys. Res. Lett., 41(15):5620-5628, doi:10.1002/2014GL061146
#' }
#' @author R. Manzanas.
#' @family calibration
#' @examples{
#' ## loading seasonal forecasts (CFS) and observations (NCEP) of boreal winter temperature over Iberia
#' require(climate4R.datasets)
#' data("CFS_Iberia_tas"); fcst = CFS_Iberia_tas
#' data("NCEP_Iberia_tas"); obs = NCEP_Iberia_tas
#' ## passing from daily data to seasonal averages
#' fcst = aggregateGrid(fcst, aggr.y = list(FUN = "mean", na.rm = TRUE))
#' obs = aggregateGrid(obs, aggr.y = list(FUN = "mean", na.rm = TRUE))
#' ## interpolating forecasts to the observations' resolution
#' fcst = interpGrid(fcst, new.coordinates = getGrid(obs))
#' ## applying calibration
#' fcst.cal = calRPC(fcst, obs, crossval = TRUE, apply.to = "all")
#' ## plotting climatologies
#' library(visualizeR)
#' spatialPlot(makeMultiGrid(climatology(obs),
#'                          climatology(fcst, by.member = FALSE),
#'                          climatology(fcst.cal, by.member = FALSE)),
#'            backdrop.theme = "coastline",
#'            layout = c(3, 1),
#'            names.attr = c("NCEP", "CFS (raw)", "CFS (calibrated)"))
#' }

calRPC <- function(fcst.grid, obs.grid, crossval = TRUE, apply.to = c("all", "sig"), alpha = 0.1) {
  # Eade et al. 2014: http://onlinelibrary.wiley.com/doi/10.1002/2014GL061146/abstract
  
  apply.to = match.arg(apply.to, choices = c("all","sig"))
  
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
      tryCatch({
        if (crossval) {
          ## leave-one-out cross-validation
          aux = sapply(1:ntimes, function(x) {
            fcst.train = fcst[,-x,ilat,ilon]
            fcst.test = fcst[,x,ilat,ilon]
            obs.train = obs[-x,ilat,ilon]
            clim.obs = mean(obs.train, na.rm = T)
            obs.train = obs.train - clim.obs
            clim.fcst = mean(fcst.train, na.rm = T)
            fcst.train = fcst.train - clim.fcst
            fcst.test = fcst.test - clim.fcst
            ens.mean = colMeans(fcst.train, na.rm = T)
            sigma.em = sd(ens.mean, na.rm = T)
            sigma.ref = sd(obs.train, na.rm = T)
            
            rho = cor.test(obs.train, ens.mean, method = "pearson", alternative = "greater")
            
            if (apply.to == "sig") {
              if (rho$p.value < alpha) {  # statistically significant (alpha*100(%)) correlation
                adj.em.test = ((mean(fcst.test, na.rm = T) - mean(ens.mean, na.rm = T))*(sigma.ref*rho$estimate/sigma.em)) + mean(ens.mean, na.rm = T)
                sigma.noi = sqrt(var(fcst.test - mean(fcst.test, na.rm = T), na.rm = T))      
                ((fcst.test - mean(fcst.test, na.rm = T))*((sigma.ref*sqrt(1-(rho$estimate^2))) / (sigma.noi))) + adj.em.test + clim.obs
              } else {
                fcst.test + clim.fcst
              }
            } else if (apply.to == "all") {
              adj.em.test = ((mean(fcst.test, na.rm = T) - mean(ens.mean, na.rm = T))*(sigma.ref*rho$estimate/sigma.em)) + mean(ens.mean, na.rm = T)
              sigma.noi = sqrt(var(fcst.test - mean(fcst.test, na.rm = T), na.rm = T))      
              ((fcst.test - mean(fcst.test, na.rm = T))*((sigma.ref*sqrt(1-(rho$estimate^2))) / (sigma.noi))) + adj.em.test + clim.obs
            }
          })
          fcst.cal[,,ilat,ilon] = aux; rm(aux)
        } else {
          fcst.train = fcst[,,ilat,ilon]
          fcst.test = fcst.train
          obs.train = obs[,ilat,ilon]
          clim.obs = mean(obs.train, na.rm = T)
          obs.train = obs.train - clim.obs
          clim.fcst = mean(fcst.train, na.rm = T)
          fcst.train = fcst.train - clim.fcst
          fcst.test = fcst.test - clim.fcst
          
          ens.mean = colMeans(fcst.train, na.rm = T)
          sigma.em = sd(ens.mean, na.rm = T)
          sigma.ref = sd(obs.train, na.rm = T)
          
          rho = cor.test(obs.train, ens.mean, method = "pearson", alternative = "greater")
          
          if (apply.to == "sig") {
            if (rho$p.value < alpha) {  # statistically significant (alpha*100(%)) correlation
              adj.em.test = ((colMeans(fcst.test, na.rm = T) - mean(ens.mean, na.rm = T))*(sigma.ref*rho$estimate/sigma.em)) + mean(ens.mean, na.rm = T)
              sigma.noi = sqrt(apply(fcst.test - matrix(colMeans(fcst.test, na.rm = T), nrow = nmemb, ncol = ntimes, byrow = T), 2, "var", na.rm = T))   
              fcst.cal[,,ilat,ilon] = ((fcst.test - matrix(colMeans(fcst.test, na.rm = T), nrow = nmemb,  ncol = ntimes, byrow = T))*((sigma.ref*sqrt(1-(rho$estimate^2))) / (sigma.noi))) + adj.em.test + clim.obs
            } else {
              fcst.cal[,,ilat,ilon] = fcst.test + clim.fcst
            }
          } else if (apply.to == "all") {
            adj.em.test = ((colMeans(fcst.test, na.rm = T) - mean(ens.mean, na.rm = T))*(sigma.ref*rho$estimate/sigma.em)) + mean(ens.mean, na.rm = T)
            sigma.noi = sqrt(apply(fcst.test - matrix(colMeans(fcst.test, na.rm = T), nrow = nmemb, ncol = ntimes, byrow = T), 2, "var", na.rm = T))   
            fcst.cal[,,ilat,ilon] = ((fcst.test - matrix(colMeans(fcst.test, na.rm = T), nrow = nmemb,  ncol = ntimes, byrow = T))*((sigma.ref*sqrt(1-(rho$estimate^2))) / (sigma.noi))) + adj.em.test + clim.obs
          }
        }
      }, error = function(x) {
        # NA data
      })
    }
  }
  fcst.out = fcst.grid
  fcst.out$Data = fcst.cal
  attributes(fcst.out$Data)$dimensions = attributes(fcst.grid$Data)$dimensions
  return(fcst.out)
}
