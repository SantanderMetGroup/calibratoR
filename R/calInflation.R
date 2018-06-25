##     calInlation.R Calibration of monthly/seasonal forecasts
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

#' @title Inflation calibration method for seasonal forecasts
#' @description This function implements the EMOS method introduced in Doblas-Reyes et al. 2005
#' and recently applied in Torralba et al. 2017 (method 2) to produce reliable operational seasonal forecasts of wind speed.
#' After Weigel et al. 2009, this method is sometimes referred to as climate conserving recalibration.
#' @note Ensemble Model Output Statistics (EMOS) methods use the correspondence between the ensemble mean and the observations in the calibration process.
#' @param fcst.grid climate4R grid. Forecasts to be calibrated (typically on a monthly/seasonal basis). At the moment, only gridded data are supported.
#' @param obs.grid climate4R grid. Reference observations the forecasts are calibrated towards (typically on a monthly/seasonal basis).
#' @param crossval Logical. TRUE (default) for leave-one-out cross-validation. FALSE for not cross-validation.
#' @param apply.to Character. If \code{"all"} is selected, all forecasts are calibrated. Alternatively, if \code{"sig"} is selected, the calibration is only applied 
#' in those points where the correlation between the ensemble mean and the observations is statistically significant.
#' @param alpha Significance (0.1 by default) of the ensemble mean correlation (i.e. \code{"alpha = 0.05"} would correspond to a 95\% confidence level). Only works if \code{"apply.to = sig"}.
#' @return climate4R grid. Calibrated forecasts.
#' @importFrom transformeR getShape
#' @importFrom stats cor.test sd
#' @export
#' @references
#' \itemize{
#' \item Doblas-Reyes, F.J., R. Hagedorn, and T.N. Palmer, 2005: The rationale behind the success of multi-model ensembles in seasonal forecasting-II. Calibration and combination. Tellus, 57A, 234–252, doi:10.3402/tellusa.v57i3.14658
#' \item Torralba, V., F.J. Doblas-Reyes, D. MacLeod, I. Christel, and M. Davis, 2017: Seasonal Climate Prediction: A New Source of Information for the Management of Wind Energy Resources. J. Appl. Meteor. Climatol., 56, 1231–1247, https://doi.org/10.1175/JAMC-D-16-0204.1 
#' \item Weigel, A.P., M.A. Liniger, and C. Appenzeller, 2009: Seasonal ensemble forecasts: Are recalibrated single models better than multimodels? Mon. Wea. Rev., 137, 1460-1479, https://doi.org/10.1175/2008MWR2773.1
#' }
#' @seealso \code{\link{calMVA}}, \code{\link{calLM}}, \code{\link{calRPC}}, \code{\link{calNGR}}
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
#' fcst.cal = calInflation(fcst, obs, crossval = TRUE, apply.to = "all")
#' ## plotting climatologies
#' library(visualizeR)
#' spatialPlot(makeMultiGrid(climatology(obs),
#'                          climatology(fcst, by.member = FALSE),
#'                          climatology(fcst.cal, by.member = FALSE)),
#'            backdrop.theme = "coastline",
#'            layout = c(3, 1),
#'            names.attr = c("NCEP", "CFS (raw)", "CFS (calibrated)"))
#' }

calInflation <- function(fcst.grid, obs.grid, crossval = TRUE, apply.to = c("all", "sig"), alpha = 0.1) {
  ## Method 2 in Torralba et al. 2017: http://www.bsc.es/ESS/sites/default/files/imce/amspaper_final.pdf
  
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
            # sigma.e = mean(apply(fcst.train, 1, sd))
            sigma.e = sd(fcst.train - matrix(ens.mean, nrow = nmemb, ncol = ntimes-1, byrow = T), na.rm = T)
            sigma.ref = sd(obs.train, na.rm = T)
            
            rho = cor.test(obs.train, ens.mean, method = "pearson", alternative = "greater")
            
            if (apply.to == "sig") {
              if (rho$p.value < alpha) {  # statistically significant (alpha*100(%)) correlation
                a = rho$estimate*(sigma.ref/sigma.em)
                b = sqrt(1-(rho$estimate^2))*(sigma.ref/sigma.e)
                
                zeta = fcst.test - mean(fcst.test, na.rm = T)
                # fcst.cal[,x,ilat,ilon] = (a*(mean(fcst.test, na.rm = T))) + (b*zeta) + clim.obs
                (a*(mean(fcst.test, na.rm = T))) + (b*zeta) + clim.obs
              } else {
                # fcst.cal[,x,ilat,ilon] = fcst.test + clim.fcst
                fcst.test + clim.fcst
              }
            } else if (apply.to == "all") {
              a = rho$estimate*(sigma.ref/sigma.em)
              b = sqrt(1-(rho$estimate^2))*(sigma.ref/sigma.e)
              
              zeta = fcst.test - mean(fcst.test, na.rm = T)
              # fcst.cal[,x,ilat,ilon] = (a*(mean(fcst.test, na.rm = T))) + (b*zeta) + clim.obs
              (a*(mean(fcst.test, na.rm = T))) + (b*zeta) + clim.obs
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
          # sigma.e = mean(apply(fcst.train, 1, sd))
          sigma.e = sd(fcst.train - matrix(ens.mean, nrow = nmemb, ncol = ntimes, byrow = T), na.rm = T)
          sigma.ref = sd(obs.train, na.rm = T)
          
          rho = cor.test(obs.train, ens.mean, method = "pearson", alternative = "greater")
          
          if (apply.to == "sig") {
            if (rho$p.value < alpha) {  # statistically significant (alpha*100(%)) correlation
              a = rho$estimate*(sigma.ref/sigma.em)
              b = sqrt(1-(rho$estimate^2))*(sigma.ref/sigma.e)
              
              zeta = fcst.test - mean(fcst.test, na.rm = T)
              fcst.cal[,,ilat,ilon] = (a*(mean(fcst.test, na.rm = T))) + (b*zeta) + clim.obs
            } else {
              fcst.cal[,,ilat,ilon] = fcst.test + clim.fcst
            }
          } else if (apply.to == "all") {
            a = rho$estimate*(sigma.ref/sigma.em)
            b = sqrt(1-(rho$estimate^2))*(sigma.ref/sigma.e)
            
            zeta = fcst.test - mean(fcst.test, na.rm = T)
            fcst.cal[,,ilat,ilon] = (a*(mean(fcst.test, na.rm = T))) + (b*zeta) + clim.obs
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


