#' @title Random Parameters Transformation
#' @description rp.model.fun is used to tranform the random parameters.
#' @param beta vector of random parameters. beta follows multivariate normal distribution.
#' @param opt a list of control parameters. See rpinpall.
#' @return the transformed beta
#' @export
rp_model_fun <- function(beta, opt){
  if(opt$distrib=="normal") {
    beta <- as.matrix(beta)
  }
  if(opt$distrib=="lognormal") {
    beta <- as.matrix(exp(beta))
  }
  if(opt$distrib == "censored-normal") {
    beta <- pmax(0,beta)
  }
  return(beta)
}


#' @title Convergence function
#' @description fun_converge is used to compute the relative error.
#' @param x vector of parameters at iteration n-1.
#' @param y vector of parameters at iteration n.
#' @param tol tolerance value.
#' @return the relative error of x.
#' @noRd
fun_converge <- function(x, y, tol){
  val <- max(abs((x-y )/(abs(y) + tol)))
  return(val)
}

## Allassonni?re and  Chevallier 2021
#' @title tempering function
#' @description tempering function
#' @param iter_max integer.
#' @param param vector of parameters at iteration n.
#' @return a vector.
#' @noRd
tempering <- function(iter_max,param=c(0.5,-4,2,4)){
  a <- param[1]
  b <- param[2]
  c <- param[3]
  r <- param[4]
  iter <- 1:iter_max
  kap <- (iter + c*r)/r
  temp <- ifelse(1 + a^kap + b*sin(kap)/kap <0, 1, 1 + a^kap + b*sin(kap)/kap)
  return(temp)
}

# Selection Matrix: observed
#' @noRd
mat_select_obs <- function(regime){
  mat_tp     <- diag(1,length(regime))
  val        <- mat_tp[regime > 0,]
  if(sum(regime)==1) val <- matrix(val,nrow = 1)
  return(val)
}
