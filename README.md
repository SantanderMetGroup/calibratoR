# What is calibratoR?

[calibratoR](https://github.com/SantanderMetGroup/calibratoR) implements several methods for statistical calibration of climate forecasts (e.g. seasonal forecasts) on a monthly/seasonal basis. The implemented methods range from simple mean/variance adjustment to more sophisticated Ensemble Model Output Statistics (EMOS) options which take into account the temporal correspondence between the ensemble mean and the observations in the calibration process. This package is complementary with [downscaleR](https://github.com/SantanderMetGroup/downscaleR), which implements standard bias correction and (perfect prog) statistical downscaling on a daily basis. 

This package is fully compatible with the [climate4R](http://www.meteo.unican.es/climate4r) framework (see [Iturbide et al. 2018](https://doi.org/10.1016/j.envsoft.2018.09.009) for details). Currently, only gridded data are supported (point-wise stations will be also supported in a future release).

The recommended installation procedure is to use the `install_github` command from the devtools R package:

```r
devtools::install_github(c("SantanderMetGroup/transformeR", "SantanderMetGroup/calibratoR"))
```
**IMPORTANT:** Note that `transformeR` and `SpecsVerification` must be previously installed on your system. The latter can be retrieved from CRAN:

```r
install.packages("SpecsVerification")
```
---

---
**Reference and further information:**

**[Application to the C3S dataset]** Manzanas et al. (2018) Bias adjustment and ensemble recalibration methods for seasonal forecasting: A comprehensive intercomparison using the C3S dataset. Submitted to Climate Dynamics.

**[General description of the climate4R framework]** Iturbide et al. (2018) climate4R: An R-based Open Framework for Reproducible Climate Data Access and Post-processing. **Environmental Modelling and Software**. https://doi.org/10.1016/j.envsoft.2018.09.009
Check out the companion notebooks for the two examples [GitHub](https://github.com/SantanderMetGroup/notebooks).

