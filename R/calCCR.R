##     calCCR.R Climate-conserving recalibration of seasonal forecasts
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

#' @title Inflation calibration method for seasonal forecasts (a.k.a. climate conserving recalibration)
#' @description This function implements the EMOS method introduced in Doblas-Reyes et al. 2005
#' and recently applied in Torralba et al. 2017 (method 2) to produce reliable operational seasonal forecasts of wind speed.
#' After Weigel et al. 2009, this method is sometimes referred to as climate conserving recalibration (CCR).
#' @note Ensemble Model Output Statistics (EMOS) methods use the correspondence between the ensemble mean and the observations in the calibration process.
#' @param fcst.grid Training data
#' @param obs.grid Reference observations
#' @param newfcst.grid Optional. Forecast to correct (usually the operative) using the parameters from \code{fcst.grid} (typically the hindcast)
#' @param crossval Logical, applied only when \code{newfcst.grid = NULL}. Performs a leave-one-out cross validation.
#' @param apply.to Character. If \code{"all"} is selected, all forecasts are calibrated. Alternatively, if \code{"sig"} is selected, the calibration is only applied 
#' in those points where the correlation between the ensemble mean and the observations is statistically significant. This correlation
#' is computed upong the hindcast ensemble mean and the reference observations.
#' @param alpha Significance (0.1 by default) of the ensemble mean correlation (i.e. \code{"alpha = 0.05"} would correspond to a 95\% confidence level). Only works if \code{"apply.to = sig"}.
#' Note that the global correlation analysis will be skipped when \code{apply.to = "all"}.
# #' @param detrend In case the correlation significance test is undertaken, should the data be (linearly) detrended first?
# #' Detrending is internally achieved by \code{\link{detrendGrid}}.
# #' @template templateParallelParams
# #' @template templateParallel
#' @return A climate4R grid of calibrated forecasts
#' @importFrom transformeR getShape
#' @importFrom stats cor.test sd
#' @importFrom magrittr %>% extract2 extract multiply_by
#' @importFrom easyVerification veriApply
#' @export
#' @references
#' \itemize{
#' \item Doblas-Reyes, F.J., R. Hagedorn, and T.N. Palmer, 2005: The rationale behind the success of multi-model ensembles in seasonal forecasting-II. Calibration and combination. Tellus, 57A, 234–252, doi:10.3402/tellusa.v57i3.14658
#' \item Torralba, V., F.J. Doblas-Reyes, D. MacLeod, I. Christel, and M. Davis, 2017: Seasonal Climate Prediction: A New Source of Information for the Management of Wind Energy Resources. J. Appl. Meteor. Climatol., 56, 1231–1247, https://doi.org/10.1175/JAMC-D-16-0204.1 
#' \item Weigel, A.P., M.A. Liniger, and C. Appenzeller, 2009: Seasonal ensemble forecasts: Are recalibrated single models better than multimodels? Mon. Wea. Rev., 137, 1460-1479, https://doi.org/10.1175/2008MWR2773.1
#' }
#' @family calibration
#' @author R. Manzanas and V. Torralba.
#' 

# load("data/cfs.forecast.rda", verbose = TRUE)
# load("data/cfs.hindcast.rda", verbose = TRUE)
# load("data/ncep.ref.rda", verbose = TRUE)
# 
# y <- ncep.ref
# x <- cfs.hindcast
# newdata <- cfs.forecast
# # newdata <- NULL
# crossval = TRUE
# apply.to = "all"
# alpha = 0.1
# parallel = TRUE
# max.ncores = 16
# ncores = 2
# detrend = FALSE

calCCR <- function(fcst.grid, obs.grid, newfcst.grid = NULL, crossval = TRUE, apply.to = c("all", "sig"), alpha = 0.1) {
    ## Torralba et al. 2017: http://www.bsc.es/ESS/sites/default/files/imce/amspaper_final.pdf
    
    apply.to = match.arg(apply.to, choices = c("all","sig"))
    
    fcst = fcst.grid$Data
    obs = obs.grid$Data
    
    stopifnot(identical(dim(fcst)[-1], dim(obs)))  # check for equality of dimensions between fcst and obs
    
    nmemb = dim(fcst)[1]
    ntimes = dim(fcst)[2]
    nlat = dim(fcst)[3]
    nlon = dim(fcst)[4]
    
    if (is.null(newfcst.grid)) {
        fcst.cal = NA*fcst
        for (ilat in 1:nlat) {
            if (!(ilat/10) - trunc(ilat/10)) {
                #print(sprintf("... lat %d of %d ...", ilat, nlat))
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
                                if (rho$p.value < alpha) {  #rho significativa a una confianza del alpha*100(%)
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
                        fcst.out = fcst.grid
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
                            if (rho$p.value < alpha) {  #rho significativa a una confianza del alpha*100(%)
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
    } else {
        # calibrating new forecast
        newfcst = newfcst.grid$Data
        stopifnot(dim(fcst)[c(1,3,4)] == dim(newfcst))  # check for equality of dimensions between hind and fcst
        
        fcst.cal = NA*newfcst
        for (ilat in 1:nlat) {
            if (!(ilat/10) - trunc(ilat/10)) {
                print(sprintf("... lat %d of %d ...", ilat, nlat))
            }
            for (ilon in 1:nlon) {
                tryCatch({
                    fcst.train = fcst[,,ilat,ilon]
                    fcst.test = newfcst[,ilat,ilon]
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
                        if (rho$p.value < alpha) {  #rho significativa a una confianza del alpha*100(%)
                            a = rho$estimate*(sigma.ref/sigma.em)
                            b = sqrt(1-(rho$estimate^2))*(sigma.ref/sigma.e)
                            
                            zeta = fcst.test - mean(fcst.test, na.rm = T)
                            fcst.cal[,ilat,ilon] = (a*(mean(fcst.test, na.rm = T))) + (b*zeta) + clim.obs
                        } else {
                            fcst.cal[,ilat,ilon] = fcst.test + clim.fcst
                        }
                    } else if (apply.to == "all") {
                        a = rho$estimate*(sigma.ref/sigma.em)
                        b = sqrt(1-(rho$estimate^2))*(sigma.ref/sigma.e)
                        
                        zeta = fcst.test - mean(fcst.test, na.rm = T)
                        fcst.cal[,ilat,ilon] = (a*(mean(fcst.test, na.rm = T))) + (b*zeta) + clim.obs
                    }
                }, error = function(x) {
                    # NA data
                })
            }
        }
        fcst.out = newfcst.grid
    }
    fcst.out$Data = fcst.cal
    return(fcst.out)
}



# calCCR <- function(y, x, newdata = NULL, crossval = TRUE,
#                    detrend = TRUE, apply.to = "all", alpha = 0.1,
#                    parallel = FALSE, max.ncores = 16, ncores = NULL) {
#     stopifnot(is.logical(detrend))
#     apply.to <- match.arg(apply.to, choices = c("all", "sig"))
#     inp <- prepareCalInputs(x, y, newdata, crossval) 
#     x <- inp$x
#     y <- inp$y
#     newdata <- inp$newdata
#     crossval <- inp$crossval
#     inp <- NULL
#     n.mem <- getShape(newdata, dimension = "member")
#     coords <- getCoordinates(x)
#     x.coords <- coords$x
#     y.coords <- coords$y
#     n.points <- ifelse(typeofGrid(y) == "station", nrow(coords), length(x.coords) * length(y.coords))
#     if (apply.to == "sig") {
#         if (detrend) {
#             message("[", Sys.time(), "] Detrending ...")
#             xrho <- suppressMessages(detrendGrid(grid = x, parallel = parallel, ncores = ncores))
#             yrho <- suppressMessages(detrendGrid(grid = y, parallel = parallel, ncores = ncores)) 
#         } else {
#             xrho <- x
#             yrho <- y %>% redim(member = FALSE)
#         }
#         # Correlation is applied to (optionally detrended) x and y.
#         message("[", Sys.time(), "] Calculating ensemble mean vs. obs correlation ...")
#         rho.pval <- suppressMessages(veriApply(verifun = "corfun",
#                                                method = "pearson",
#                                                alternative = "greater",
#                                                fcst = xrho$Data,
#                                                obs = yrho$Data,
#                                                tdim = 2,
#                                                ensdim = 1,
#                                                parallel = parallel,
#                                                ncpus = ncores)) %>% extract2("pval") %>% easyVeri2grid(obs.grid = yrho,
#                                                                                                        verifun = "calibratoR:::corfun") %>% extract2("Data") %>% array3Dto2Dmat() %>% extract(1,)
#         xrho <- yrho <- NULL
#     }
#     ind <- if (apply.to == "all") {
#         1:n.points
#     } else {
#         which(rho.pval < alpha)
#     }
#     if (length(ind) == 0) {
#         warning("No significant correlations found for any of the points in the domain: no correction was applied")  
#         return(y)
#     }
#     rho.pval <- NULL
#     if (crossval) {
#         message("[", Sys.time(), "] Calculating parameters (leave-one-year out cross-validation mode) ...")
#         x <- NULL
#         ppars <- parallelCheck(parallel = parallel, max.ncores = max.ncores, ncores = ncores)
#         lapplyfun <- selectPar.pplyFun(parallel.pars = ppars, .pplyFUN = "lapply")
#         if (ppars$hasparallel) on.exit(parallel::stopCluster(ppars$cl))
#         yrs <- getYearsAsINDEX(newdata) %>% unique()
#         l1 <- lapplyfun(yrs, function(j) {
#             message("[", Sys.time(), "] Taking out ", j, " (", yrs[length(yrs)] - j, " years remaining ...)")
#             yr.ind <- setdiff(yrs, j)
#             obs.train <- subsetGrid(y, years = yr.ind) 
#             clim.obs <- suppressMessages(climatology(obs.train) %>% extract2("Data") %>% array3Dto2Dmat() %>% extract(1,))
#             suppressMessages(obs.train %<>% scaleGrid(time.frame = "monthly") %>% redim(drop = TRUE) %>% redim(member = FALSE))
#             sigma.ref <- getPooledMemberStat(obs.train, fun = "sd", na.rm = TRUE)
#             f.train <- subsetGrid(newdata, years = yr.ind)
#             f.test <- suppressMessages(subsetGrid(newdata, years = j) %>% scaleGrid(base = f.train, by.member = FALSE, time.frame = "monthly") %>% redim(drop = TRUE) %>% redim(member = TRUE))
#             suppressMessages(f.train %<>% scaleGrid(time.frame = "monthly", by.member = FALSE) %>% redim(drop = TRUE) %>% redim(member = TRUE)) 
#             ens.mean <- suppressMessages(aggregateGrid(f.train, aggr.mem = list(FUN = "mean", na.rm = TRUE)))
#             sigma.em <- getPooledMemberStat(ens.mean, fun = "sd", na.rm = TRUE)
#             # Sospechoso: (anomalias de anomalias...)
#             sigma.e <- suppressMessages(scaleGrid(f.train, base = ens.mean) %>% getPooledMemberStat(fun = "sd", na.rm = TRUE))
#             # sigma.e <- getPooledMemberStat(f.train, fun = "sd")
#             rho <- suppressMessages(veriApply(verifun = "corfun",
#                                               method = "pearson", alternative = "greater",
#                                               fcst = f.train$Data,
#                                               obs = obs.train$Data,
#                                               tdim = 2, ensdim = 1,
#                                               parallel = FALSE,
#                                               ncpus = ncores,
#                                               maxncpus = max.ncores)) %>% extract2("coef") %>% easyVeri2grid(obs.grid = obs.train,
#                                                                                                              verifun = "corfun") %>% extract2("Data") %>% array3Dto2Dmat() %>% extract(1,)
#             a <- rho * (sigma.ref/sigma.e)
#             b <- sqrt(1 - rho^2) * (sigma.ref/sigma.em)
#             clim.test <- suppressMessages(climatology(f.test, by.member = FALSE) %>% extract2("Data") %>% array3Dto2Dmat() %>% extract(1, ))
#             aux1 <- f.train
#             l <- lapply(1:n.mem, function(i) {
#                 z <- subsetGrid(f.train, members = i) %>% extract2("Data") %>% array3Dto2Dmat() %>% t()
#                 z[ind, ] <- (a * clim.test + (b * z) + clim.obs) %>% extract(ind, )
#                 # aux1$Data <- (a * clim.test + (b * z) + clim.obs) %>% t() %>% mat2Dto3Darray(x = x.coords, y = y.coords)  
#                 aux1$Data <- mat2Dto3Darray(t(z), x = x.coords, y = y.coords)  
#                 return(aux1)
#             })
#             return(do.call("bindGrid", c(l, dimension = "member")))
#         })
#         y <- newdata <- NULL
#         cal <- do.call("bindGrid", c(l1, dimension = "time")) %>% redim(drop = TRUE)
#     } else {
#         message("[", Sys.time(), "] Calculating parameters ...")
#         clim.obs <- suppressMessages(climatology(y) %>% extract2("Data") %>% array3Dto2Dmat() %>% extract(1,))
#         obs.train <- suppressMessages(scaleGrid(y, time.frame = "monthly") %>% redim(drop = TRUE) %>% redim(member = FALSE))
#         sigma.ref <- getPooledMemberStat(obs.train, fun = "sd", na.rm = TRUE)
#         f.train <- suppressMessages(scaleGrid(grid = x, by.member = FALSE, time.frame = "monthly") %>% redim(drop = TRUE) %>% redim(member = TRUE))
#         f.test <- suppressMessages(scaleGrid(grid = newdata, base = x, by.member = FALSE, time.frame = "monthly") %>% redim(drop = TRUE) %>% redim(member = TRUE))
#         ens.mean <- suppressMessages(aggregateGrid(f.train, aggr.mem = list(FUN = "mean", na.rm = TRUE)))
#         sigma.em <- getPooledMemberStat(ens.mean, fun = "sd", na.rm = TRUE)
#         # Sospechoso: (anomalias de anomalias...)
#         # sigma.e <- suppressMessages(scaleGrid(f.train, base = ens.mean) %>% getPooledMemberStat(fun = "sd", na.rm = TRUE))
#         sigma.e <- getPooledMemberStat(grid = f.train, fun = "sd", na.rm = TRUE)
#         rho <- suppressMessages(veriApply(verifun = "corfun",
#                                           method = "pearson", alternative = "greater",
#                                           fcst = f.train$Data,
#                                           obs = obs.train$Data,
#                                           tdim = 2, ensdim = 1,
#                                           parallel = parallel,
#                                           ncpus = ncores,
#                                           maxncpus = max.ncores)) %>% extract2("coef") %>% easyVeri2grid(obs.grid = obs.train,
#                                                                                                          verifun = "calibratoR::corfun") %>% extract2("Data") %>% array3Dto2Dmat() %>% extract(1,)
#         a <- rho * (sigma.ref/sigma.e)
#         b <- sqrt(1 - rho^2) * (sigma.ref/sigma.em)
#         clim.test <- suppressMessages(climatology(f.test, by.member = FALSE) %>% extract2("Data") %>% array3Dto2Dmat() %>% extract(1, ))
#         aux1 <- f.test
#         message("[", Sys.time(), "] Applying calibration ...")
#         ppars <- parallelCheck(parallel, max.ncores, ncores)
#         lapplyfun <- selectPar.pplyFun(parallel.pars = ppars, .pplyFUN = "lapply")
#         if (ppars$hasparallel) on.exit(parallel::stopCluster(ppars$cl))
#         l <- lapplyfun(1:n.mem, function(i) {
#             z <- subsetGrid(f.test, members = i) %>% extract2("Data") %>% array3Dto2Dmat() %>% t()
#             z[ind, ] <- (a * clim.test + (b * z) + clim.obs) %>% extract(ind, )
#             aux1$Data <- mat2Dto3Darray(t(z), x = x.coords, y = y.coords)  
#             return(aux1)
#         })
#         cal <- do.call("bindGrid", c(l, dimension = "member"))
#     }
#     message("[", Sys.time(), "] Done.")
#     if (length(unique(cal$Members)) != n.mem) cal$Members <- paste0("Member_", 1:n.mem)
#     return(cal)
# }
#     
    
#' @keywords internal
    
calCCR_old <- function(fcst.grid, obs.grid, crossval = TRUE, apply.to = c("all", "sig"), alpha = 0.1) {
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







calInflation <- function(fcst.grid, obs.grid, crossval = TRUE, apply.to = c("all", "sig"), alpha = 0.1) {
  .Deprecated(new = "calCCR", old = "calInflation")
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






