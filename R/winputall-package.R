#' The 'winputall' package.
#'
#' @description Using a time-varying random parameters model, this package allows
#'   allocating variable input costs among crops produced by farmers based on panel
#'   data including information on input expenditure aggregated at the farm level
#'   and acreage shares. It also considers in fairly way the weighting data and
#'   can allow integrating time-varying and time-constant control variables.
#'
#' @docType package
#' @name winputall-package
#' @aliases winputall
#' @useDynLib winputall, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom rstan vb
#' @importFrom rstan optimizing
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#'
#' @references
#' Koutchade Obafèmi Philippe, Fabienne Femenia, Alain Carpentier (2024).Variable Input Allocation Among Crops: A Time-Varying Random Parameters Approach https://hal.science/hal-04318163
#'
"_PACKAGE"
