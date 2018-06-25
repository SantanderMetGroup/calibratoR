##     calNGR.R Calibration of monthly/seasonal forecasts
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

#' @title Non-homogeneous Gaussian Regression calibration for seasonal forecasts
#' @description This function implements an EMOS method that adopts the Non-homogeneous Gaussian Regression (NGR; Gneiting et al. 2005) for use with ensemble forecasts. 
#' In particular, it uses a constant term and the ensemble mean signal as predictors for the calibrated forecast mean and a constant term and the ensemble spread for the inflation (shrinkage) 
#' of the ensemble spread. NGR approaches have been applied in many previous works, but mostly in the context of short-term forecasts. 
#' Only Tippet and Barnston 2008 have used NGRs in the context of seasonal forecasting.
#' @note Ensemble Model Output Statistics (EMOS) methods use the correspondence between the ensemble mean and the observations in the calibration process.
#' @param fcst.grid climate4R grid. Forecasts to be calibrated (typically on a monthly/seasonal basis). At the moment, only gridded data are supported.
#' @param obs.grid climate4R grid. Reference observations the forecasts are calibrated towards (typically on a monthly/seasonal basis).
#' @param crossval Logical. TRUE (default) for leave-one-out-out cross-validation. FALSE for not cross-validation.
#' @param type Logical. NGR or ensNGR. In the former (latter) case, the parameters are optimized by minimizing the CRPS (fair CRPS).
#' @param apply.to Character. If \code{"all"} is selected, all forecasts are calibrated. Alternatively, if \code{"sig"} is selected, the calibration is only applied 
#' in those points where the correlation between the ensemble mean and the observations is statistically significant.
#' @param alpha Significance (0.1 by default) of the ensemble mean correlation (i.e. \code{"alpha = 0.05"} would correspond to a 95\% confidence level). Only works if \code{"apply.to = sig"}.  
#' @return climate4R grid. Calibrated forecasts.		
#' @importFrom transformeR getShape
#' @importFrom SpecsVerification FairCrps GaussCrps
#' @importFrom stats cor.test optim sd coef lm
#' @export
#' @references
#' \itemize{
#' \item Gneiting, T., A.E. Raftery, A.H. Westveld, and T.Goldman, 2005: Calibrated Probabilistic Forecasting Using Ensemble Model Output Statistics and Minimum CRPS Estimation. Mon. Wea. Rev., 133 (5): 1098–1118. doi:10.1175/MWR2904.1
#' \item Tippett, M.K. and A.G. Barnston, 2008: Skill of Multimodel ENSO Probability Forecasts. Mon. Wea. Rev., 136, 3933–3946, https://doi.org/10.1175/2008MWR2431.1 
#' }
#' @author R. Manzanas and J. Bhend.
#' @seealso \code{\link{calMVA}}, \code{\link{calInflation}}, \code{\link{calLM}}, \code{\link{calRPC}}
#' @examples{
#' ## loading seasonal forecasts (CFS) and observations (NCEP) of boreal winter temperature over Iberia
#' require(transformeR)
#' data("CFS_Iberia_tas"); fcst = CFS_Iberia_tas
#' data("NCEP_Iberia_tas"); obs = NCEP_Iberia_tas
#' ## passing from daily data to seasonal averages
#' fcst = aggregateGrid(fcst, aggr.y = list(FUN = "mean", na.rm = TRUE))
#' obs = aggregateGrid(obs, aggr.y = list(FUN = "mean", na.rm = TRUE))
#' ## interpolating forecasts to the observations' resolution
#' fcst = interpGrid(fcst, new.coordinates = getGrid(obs))
#' ## applying calibration
#' fcst.cal = calNGR(fcst, obs, crossval = TRUE, type = "NGR", apply.to = "all")
#' ## plotting climatologies
#' library(visualizeR)
#' spatialPlot(makeMultiGrid(climatology(obs),
#'                          climatology(fcst, by.member = FALSE),
#'                          climatology(fcst.cal, by.member = FALSE)),
#'            backdrop.theme = "coastline",
#'            layout = c(3, 1),
#'            names.attr = c("NCEP", "CFS (raw)", "CFS (calibrated)"))
#' }

calNGR <- function(fcst.grid, obs.grid, crossval = TRUE, type = c("NGR", "ensNGR"), apply.to = c("all", "sig"), alpha = 0.1) {
  
  apply.to = match.arg(apply.to, choices = c("all","sig"))
  type = match.arg(type, choices = c("NGR","ensNGR"))
  
  fcst = fcst.grid$Data
  obs = obs.grid$Data
  
  stopifnot(identical(dim(fcst)[-1], dim(obs)))  # check for equality of dimensions between fcst and obs
  
  nmemb = getShape(fcst.grid, "member")
  ntimes = getShape(fcst.grid, "time")
  nlat = getShape(fcst.grid, "lat")
  nlon = getShape(fcst.grid, "lon")
  
  # internal optimization function
  if (type == "ensNGR") {
    optfun <- function(pars, fmn, fsd, fanom, obs) {
      mean(FairCrps(pars[1] + pars[2]*fmn + sqrt(pars[3]**2 + pars[4]**2 * fsd**2) * fanom, obs))
    }
  } else if (type == "NGR") {
    optfun <- function(pars, fmn, fsd, fanom, obs) {
      mean(GaussCrps(pars[1] + pars[2]*fmn, sqrt(pars[3]**2 + pars[4]**2*fsd**2), obs))
    }
  }
  
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
            obs.train = obs[-x,ilat,ilon]
            fcst.train = fcst[,-x,ilat,ilon]
            ens.mean.train = colMeans(fcst.train)
            sigma.e.train = apply(fcst.train, 2, sd, na.rm = T)
            fcst.train.anom = (fcst.train - matrix(ens.mean.train, nrow = nmemb, ncol = ntimes-1, byrow = T)) / matrix(sigma.e.train, nrow = nmemb, ncol = ntimes-1, byrow = T)
            clim.fcst = mean(fcst.train, na.rm = T)
            
            fcst.test = fcst[,x,ilat,ilon]
            ens.mean.test = mean(fcst.test, na.rm = T)
            sigma.e.test = sd(fcst.test, na.rm = T)
            fcst.test.anom = (fcst.test - ens.mean.test) / sigma.e.test
            
            rho = cor.test(obs.train, ens.mean.train, method = "pearson", alternative = "greater")
            if (apply.to == "sig") {
              if (rho$p.value < alpha) {  # statistically significant (alpha*100(%)) correlation
                opars = optim(c(coef(lm(obs.train ~ ens.mean.train)),0,1), optfun,
                              fmn = ens.mean.train, fsd = sigma.e.train, fanom = t(fcst.train.anom), obs = obs.train,
                              control = list(maxit = 1e5), method = "BFGS")
                opars$par[1] + opars$par[2]*ens.mean.test +
                  sqrt(opars$par[3]**2 + opars$par[4]**2 * sigma.e.test**2) * fcst.test.anom
              } else {
                fcst.test
              }
            } else if (apply.to == "all") {
              opars = optim(c(coef(lm(obs.train ~ ens.mean.train)),0,1), optfun,
                            fmn = ens.mean.train, fsd = sigma.e.train, fanom = t(fcst.train.anom), obs = obs.train,
                            control = list(maxit = 1e5), method = "BFGS")
              opars$par[1] + opars$par[2]*ens.mean.test +
                sqrt(opars$par[3]**2 + opars$par[4]**2 * sigma.e.test**2) * fcst.test.anom
            }
            
          })
          fcst.cal[,,ilat,ilon] = aux; rm(aux)
        } else {
          obs.train = obs[,ilat,ilon]
          fcst.train = fcst[,,ilat,ilon]
          ens.mean.train = colMeans(fcst.train)
          sigma.e.train = apply(fcst.train, 2, sd, na.rm = T)
          fcst.train.anom = (fcst.train - matrix(ens.mean.train, nrow = nmemb, ncol = ntimes, byrow = T)) / matrix(sigma.e.train, nrow = nmemb, ncol = ntimes, byrow = T)
          clim.fcst = mean(fcst.train, na.rm = T)
          
          fcst.test = fcst.train
          sigma.e.test = sigma.e.train
          fcst.test.anom = fcst.train.anom
          
          rho = cor.test(obs.train, ens.mean.train, method = "pearson", alternative = "greater")
          
          if (apply.to == "sig") {
            if (rho$p.value < alpha) {  # statistically significant (alpha*100(%)) correlation
              opars = optim(c(coef(lm(obs.train ~ ens.mean.train)),0,1), optfun,
                            fmn = ens.mean.train, fsd = sigma.e.train, fanom = t(fcst.train.anom), obs = obs.train,
                            control = list(maxit = 1e5), method = "BFGS")
              aux = sapply(1:ntimes, function(x) {
                opars$par[1] + opars$par[2]*ens.mean.train[x] + 
                  sqrt(opars$par[3]**2 + opars$par[4]**2 * sigma.e.test[x]**2) * fcst.test.anom[,x]
              })
            } else {
              aux = fcst.test
            }
          } else if (apply.to == "all") {
            opars = optim(c(coef(lm(obs.train ~ ens.mean.train)),0,1), optfun,
                          fmn = ens.mean.train, fsd = sigma.e.train, fanom = t(fcst.train.anom), obs = obs.train,
                          control = list(maxit = 1e5), method = "BFGS")
            aux = sapply(1:ntimes, function(x) {
              opars$par[1] + opars$par[2]*ens.mean.train[x] + 
                sqrt(opars$par[3]**2 + opars$par[4]**2 * sigma.e.test[x]**2) * fcst.test.anom[,x]
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

