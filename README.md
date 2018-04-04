# What is calibratoR?

[calibratoR](https://github.com/SantanderMetGroup/calibratoR) implements several methods for statistical calibration of climate forecasts (typically on a monthly/seasonal basis), which range from simple mean/variance adjustment to more sophisticated Ensemble Model Output Statistics (EMOS) options which take into account the existing correspondence between the ensemble mean and the observations in the calibration process. This R package works with the *grid* object developed for the [climate4R](http://www.meteo.unican.es/climate4r) bundle (see [Iturbide et al. 2018](***) for details). Currently, only gridded data are supported (point-wise stations will be also supported in a future release).

The recommended installation procedure is to use the `install_github` command from the devtools R package:

```r
devtools::install_github(c("SantanderMetGroup/transformeR", "SantanderMetGroup/calibratoR"))
```
**IMPORTANT:** Note that `transformeR` and `SpecsVerification` must be previously installed on your system. The latter can be retrieved from CRAN:

```r
install.packages("SpecsVerification")
```
---

**Reference and further information on the methods implemented:** 
Manzanas et al. (2018): 
