##     calLM.R Calibration of monthly/seasonal forecasts
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
#' @description This function performs an EMOS-like linear regression between the ensemble mean and the corresponding observations. 
#' To correct the forecast variance, the standardized anomalies are rescaled by the standard deviation of the predictive distribution from the linear fitting.
#' @note Ensemble Model Output Statistics (EMOS) methods use the correspondence between the ensemble mean and the observations in the calibration process.
#' @param fcst.grid climate4R grid. Forecasts to be calibrated (typically on a monthly/seasonal basis). At the moment, only gridded data are supported.
#' @param obs.grid climate4R grid. Reference observations the forecasts are calibrated towards (typically on a monthly/seasonal basis).
#' @param crossval Logical. TRUE (default) for leave-one-out-out cross-validation. FALSE for not cross-validation.
#' @param apply.to Character. If \code{"all"} is selected, all forecasts are calibrated. Alternatively, if \code{"sig"} is selected, the calibration is only applied 
#' in those points where the correlation between the ensemble mean and the observations is statistically significant.
#' @param alpha Significance (0.1 by default) of the ensemble mean correlation (i.e. \code{"alpha = 0.05"} would correspond to a 95\% confidence level). Only works if \code{"apply.to = sig"}.  
#' @return climate4R grid. Calibrated forecasts.
#' @importFrom transformeR getShape
#' @importFrom stats cor.test lm predict sd	pnorm	
#' @export
#' @author R. Manzanas and J. Bhend.
#' @seealso \code{\link{calDebias}}, \code{\link{calInflation}}, \code{\link{calRPC}}, \code{\link{calNGR}}
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
#' fcst.cal = calLM(fcst, obs, crossval = TRUE, apply.to = "all")
#' ## plotting climatologies
#' library(visualizeR)
#' spatialPlot(makeMultiGrid(climatology(obs),
#'                          climatology(fcst, by.member = FALSE),
#'                          climatology(fcst.cal, by.member = FALSE)),
#'            backdrop.theme = "coastline",
#'            layout = c(3, 1),
#'            names.attr = c("NCEP", "CFS (raw)", "CFS (calibrated)"))
#' }

calLM <- function(fcst.grid, obs.grid, crossval = TRUE, apply.to = c("all", "sig"), alpha = 0.1) {
  
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
      print(sprintf("... lat %d of %d ...", ilat, nlat))
    }
    for (ilon in 1:nlon) {
      tryCatch({
        if (crossval) {  
          ## leave-one-out cross-validation
          aux = sapply(1:ntimes, function(x) {
            obs.train = obs[-x,ilat,ilon]
            fcst.train = fcst[,-x,ilat,ilon]
            ens.mean = colMeans(fcst.train)
            clim.fcst = mean(fcst.train, na.rm = T)
            
            fcst.test = fcst[,x,ilat,ilon]
            fcst.test.anom = (fcst.test - mean(fcst.test, na.rm = T)) / sd(fcst.test, na.rm = 2)
            
            rho = cor.test(obs.train, ens.mean, method = "pearson", alternative = "greater")
            
            if (apply.to == "sig") {
              if (rho$p.value < alpha) {  # statistically significant (alpha*100(%)) correlation
                flm <- lm(obs ~ fcst, data.frame(obs = obs.train, fcst = ens.mean))
                
                plm = predict(flm, data.frame(fcst = mean(fcst.test)), 
                              interval = "prediction", level =  pnorm(1) - pnorm(-1))
                plm[,1] + -apply(plm[,1:2, drop=F], 1, diff) * fcst.test.anom
              } else {
                fcst.test
              }
            } else if (apply.to == "all") {
              flm <- lm(obs ~ fcst, data.frame(obs = obs.train, fcst = ens.mean))
              
              plm = predict(flm, data.frame(fcst = mean(fcst.test)), 
                            interval = "prediction", level =  pnorm(1) - pnorm(-1))
              plm[,1] + -apply(plm[,1:2, drop=F], 1, diff) * fcst.test.anom
            }
          })
          fcst.cal[,,ilat,ilon] = aux; rm(aux)
        } else {
          obs.train = obs[,ilat,ilon]
          fcst.train = fcst[,,ilat,ilon]
          ens.mean = colMeans(fcst.train)
          sigma.e.train = apply(fcst.train, 2, sd, na.rm = T)
          clim.fcst = mean(fcst.train, na.rm = T)
          
          fcst.test = fcst.train
          fcst.test.anom = (fcst.train - matrix(ens.mean, nrow = nmemb, ncol = ntimes, byrow = T)) / matrix(sigma.e.train, nrow = nmemb, ncol = ntimes, byrow = T)
          
          rho = cor.test(obs.train, ens.mean, method = "pearson", alternative = "greater")
          if (apply.to == "sig") {
            if (rho$p.value < alpha) {  # statistically significant (alpha*100(%)) correlation
              flm <- lm(obs ~ fcst, data.frame(obs = obs.train, fcst = ens.mean))
              plm = predict(flm, data.frame(fcst = ens.mean), 
                            interval = "prediction", level =  pnorm(1) - pnorm(-1))
              
              aux = sapply(1:ntimes, function(x){
                plm[x,1] + -apply(plm[x,1:2, drop=F], 1, diff) * fcst.test.anom[,x]
              })
            } else {
              aux = fcst.test
            }
          } else if (apply.to == "all") {
            flm <- lm(obs ~ fcst, data.frame(obs = obs.train, fcst = ens.mean))
            plm = predict(flm, data.frame(fcst = ens.mean), 
                          interval = "prediction", level =  pnorm(1) - pnorm(-1))
            
            aux = sapply(1:ntimes, function(x){
              plm[x,1] + -apply(plm[x,1:2, drop=F], 1, diff) * fcst.test.anom[,x]
            })
          }
          fcst.cal[,,ilat,ilon] = aux; rm(aux)
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
