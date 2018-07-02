#' @param x climate4R grid. Forecasts used as predictors (typically a hindcast on a monthly/seasonal basis).
#' @param y climate4R grid. Reference observations the forecasts are calibrated towards (typically on a monthly/seasonal basis).
#' @param newdata Data to be corrected using the parameters from \code{x}. Typically the operational forecast.
#' By default \code{newdata} is \code{NULL}, so the hindcast (\code{x}) is corrected. 
#' @param crossval Logical. Default to \code{TRUE}, performs a leave-one-year-out cross-validation, so
#' the parameters used for the calibration are computed using all years but the current one.
#' Unused when an operational forecast is passed to \code{newdata}.
