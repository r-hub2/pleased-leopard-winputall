
<!-- README.md is generated from README.Rmd. Please edit that file -->

# winputall

<!-- badges: start -->
<!-- badges: end -->

The goal of winputall is to provide cost input per crop. Using a
time-varying random parameters model developed in Koutchade et al.,
(2024) <https://hal.science/hal-04318163>, this package allows
allocating variable input costs among crops produced by farmers based on
panel data including information on input expenditure aggregated at the
farm level and acreage shares. It also considers in fairly way the
weighting data and can allow integrating time-varying and time-constant
control variables.

## Installation

You can install the development version of winputall like so:

``` r
devtools::install_gitlab("obafemi-philippe.koutchade/winputalloc", host = "https://forgemia.inra.fr")
```

## Example

This is a basic example which shows you how to estimate the population’s
parameters and calibrate the input uses per crop:

``` r
library("winputall")
data(my_winputall_data)
fit <- rpinpallEst(data = my_winputall_data,
                  id_time = c("id","year"),
                  total_input = "tx",
                  crop_acreage = c("s_crop1","s_crop2","s_crop3"),
                  distrib_method = "lognormal",
                  sim_method = "map_imh",
                  calib_method = "cmode",
                  saem_control = list(nb_SA = 10, nb_smooth = 10, estim_rdraw = 100))
#>   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |=====                                                                 |   7%  |                                                                              |=======                                                               |  10%  |                                                                              |=========                                                             |  13%  |                                                                              |============                                                          |  17%  |                                                                              |==============                                                        |  20%  |                                                                              |================                                                      |  23%  |                                                                              |===================                                                   |  27%  |                                                                              |=====================                                                 |  30%  |                                                                              |=======================                                               |  33%  |                                                                              |==========================                                            |  37%  |                                                                              |============================                                          |  40%  |                                                                              |==============================                                        |  43%  |                                                                              |=================================                                     |  47%  |                                                                              |===================================                                   |  50%  |                                                                              |=====================================                                 |  53%  |                                                                              |========================================                              |  57%  |                                                                              |==========================================                            |  60%  |                                                                              |============================================                          |  63%  |                                                                              |===============================================                       |  67%  |                                                                              |=================================================                     |  70%  |                                                                              |===================================================                   |  73%  |                                                                              |======================================================                |  77%  |                                                                              |========================================================              |  80%  |                                                                              |==========================================================            |  83%  |                                                                              |=============================================================         |  87%  |                                                                              |===============================================================       |  90%  |                                                                              |=================================================================     |  93%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%
```

### Summary of The following command chows

You can get the summary of estimations (estimated population’s
parameters and corresponding standard errors and p-value ) using this
command:

``` r
summary(fit)
#> $beta_b
#>        Estimate     StdErr   t_value
#> Crop_1 0.769346 0.05763893 13.347682
#> Crop_2 0.558724 0.06724847  8.308353
#> Crop_3 1.214241 0.08331344 14.574373
#> 
#> $diag_omega_b
#>          Estimate     StdErr   t_value
#> Crop_1 0.03926041 0.02491322 1.5758870
#> Crop_2 0.02356368 0.02915231 0.8082954
#> Crop_3 0.01298999 0.02738246 0.4743908
#> 
#> $diag_omega_e
#>          Estimate      StdErr  t_value
#> Crop_1 0.01719986 0.006280618 2.738561
#> Crop_2 0.03308994 0.018630540 1.776113
#> Crop_3 0.02193888 0.017672999 1.241378
#> 
#> $omega_u
#>         Estimate     StdErr  t_value
#> omega 0.02830819 0.00684317 4.136707
#> 
#> $beta_d
#>                  Estimate     StdErr     t_value
#> fid_time21994 -0.15648189 0.08738422 -1.79073401
#> fid_time21995  0.09288702 0.11742020  0.79106512
#> fid_time21996  0.04716402 0.17208261  0.27407781
#> fid_time21997  0.12670658 0.15256684  0.83049879
#> fid_time21998  0.19082133 0.25175487  0.75796480
#> fid_time21999  0.05932563 0.15944949  0.37206533
#> fid_time21994  0.01077295 0.22963650  0.04691308
#> fid_time21995 -0.08312569 0.29451137 -0.28224951
#> fid_time21996  0.11944555 0.25740906  0.46403011
#> fid_time21997  0.05480202 0.24569143  0.22305221
#> fid_time21998  0.15184941 0.60086514  0.25271795
#> fid_time21999  0.01147808 0.28331506  0.04051349
#> fid_time21994 -0.05895420 0.24913695 -0.23663370
#> fid_time21995 -0.01731513 0.38444578 -0.04503919
#> fid_time21996  0.07550077 0.36221331  0.20844285
#> fid_time21997  0.07543108 0.18242641  0.41348769
#> fid_time21998  0.17971862 0.22887088  0.78524023
#> fid_time21999  0.25175808 0.16375859  1.53737334
#> 
#> $sim_r_squared
#>         Estimate
#> sim_r2 0.9769992
```

### Check the convergence of estimation

To check the convergence, run for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

### Distribution of estimated crop input uses

The following code shows you the distribution of crop input uses:

``` r
print(fit)
#> Call:
#> rpinpallEst(data = my_winputall_data, id_time = c("id", "year"), 
#>     total_input = "tx", crop_acreage = c("s_crop1", "s_crop2", 
#>         "s_crop3"), distrib_method = "lognormal", sim_method = "map_imh", 
#>     calib_method = "cmode", saem_control = list(nb_SA = 10, nb_smooth = 10, 
#>         estim_rdraw = 100))
#> 
#> Distribution of crop input uses:
#>          s_crop1_sh_X s_crop2_sh_X s_crop3_sh_X
#> Mean        2.3653712    1.8483589    3.6464656
#> Std.dev.    0.5269317    0.2598688    0.4714726
#> Min.        1.3377069    1.2837059    2.8221989
#> 1st Qu      1.9746842    1.6656265    3.3064640
#> Median      2.2825371    1.8156327    3.5617256
#> 3rd Qu      2.6707134    1.9910962    3.9608675
#> Max.        4.2429049    2.5695979    5.9153713
#> 
#> Mean by year of crop input uses:
#>   year s_crop1_sh_X s_crop2_sh_X s_crop3_sh_X
#> 1 1993     2.196506     1.743965     3.269636
#> 2 1994     1.920565     1.814173     3.183703
#> 3 1995     2.462400     1.618560     3.317188
#> 4 1996     2.342072     2.009337     3.650486
#> 5 1997     2.553509     1.877619     3.648641
#> 6 1998     2.701774     2.076605     4.050690
#> 7 1999     2.380772     1.798252     4.404915
```

### Estimated crop input uses

The following code shows the first parts of the estimated crops’ input
uses:

``` r
library("Matrix")
head(fit$xit_pred)
#>   id year s_crop1_sh_X s_crop2_sh_X s_crop3_sh_X
#> 1  1 1993     2.325827     1.758998     3.171023
#> 2  1 1994     1.880519     1.702855     2.837764
#> 3  1 1995     2.468293     1.591250     3.065120
#> 4  1 1996     2.611693     2.049696     3.536285
#> 5  1 1997     2.834446     1.892086     3.561735
#> 6  1 1998     2.949463     2.076206     3.915467
```
