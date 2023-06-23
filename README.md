# What is calibratoR?

[calibratoR](https://github.com/SantanderMetGroup/calibratoR) implements several methods for the statistical calibration of climate forecasts (e.g. seasonal forecasts) on a monthly/seasonal basis; from simple mean/variance adjustment to more sophisticated Ensemble Model Output Statistics (EMOS) options which take into account the temporal correspondence between the ensemble mean and the observations in the calibration process. This package is complementary with [downscaleR](https://github.com/SantanderMetGroup/downscaleR), which implements standard bias correction and (perfect prog) statistical downscaling on a daily basis. It is fully compatible with the [climate4R](http://www.meteo.unican.es/climate4r) framework (see [Iturbide et al. 2018](https://doi.org/10.1016/j.envsoft.2018.09.009) for details). Currently, only gridded data are supported (point-wise stations will be also supported in a future release).

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

**[Comparison of bias adjustment, ensemble recalibration, MOS and PP techniques]** Manzanas et al. (2020) Statistical adjustment, calibration and downscaling of seasonal forecasts: a case-study for Southeast Asia. **Climate Dynamics**. https://doi.org/10.1007/s00382-020-05145-1

**[Application to the C3S dataset]** Manzanas et al. (2019) Bias adjustment and ensemble recalibration methods for seasonal forecasting: A comprehensive intercomparison using the C3S dataset. **Climate Dynamics**. https://doi.org/10.1007/s00382-019-04640-4

**[General description of the climate4R framework]** Iturbide et al. (2018) climate4R: An R-based Open Framework for Reproducible Climate Data Access and Post-processing. **Environmental Modelling and Software**. https://doi.org/10.1016/j.envsoft.2018.09.009
Check out the companion notebooks for the two examples [GitHub](https://github.com/SantanderMetGroup/notebooks).

