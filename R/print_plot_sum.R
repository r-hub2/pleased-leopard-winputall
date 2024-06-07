#' @importFrom stats formula
#' @importFrom stats as.formula
#' @importFrom stats coef
#' @importFrom stats lm
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats terms
#' @importFrom stats printCoefmat
#' @importFrom stats model.response
#' @importFrom stats dnorm
#' @importFrom stats runif
#' @importFrom stats pt
#' @importFrom stats na.omit
#' @importFrom stats aggregate
#' @importFrom stats pnorm
#' @importFrom stats plogis
#' @importFrom stats var
#' @importFrom dplyr arrange
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr n
#' @importFrom dplyr mutate
#' @importFrom dplyr across
#' @importFrom dplyr all_of
#' @importFrom dplyr %>%
#' @importFrom Matrix bdiag
#' @importFrom MASS mvrnorm
#' @importFrom MASS mvrnorm
#' @importFrom stats optim
#' @importFrom matrixcalc vec
#' @importFrom matrixcalc duplication.matrix
#' @importFrom matrixcalc is.positive.definite
#' @importFrom matrixStats colWeightedMeans
#' @importFrom matrixStats logSumExp
#' @importFrom matrixStats colVars
#' @importFrom matrixStats colQuantiles
#' @importFrom ks invvech
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' @importFrom LearnBayes dmnorm
#' @importFrom future plan
#' @importFrom future availableCores
#' @importFrom future.apply future_lapply
#' @importFrom future multisession
#' @importFrom graphics abline
#' @importFrom plm make.pbalanced


rpinpall <- function(x, ...) UseMethod("rpinpall")
#summary.rpinpall <- function(x, ...) UseMethod("summary.rpinpall")

#' @describeIn rpinpallEst Displays the distribution of estimated crop input uses accounting for error by default
#' @param  x an object produced by the function \code{rpinpallEst}, to be displayed
#' @param  error logical. If TRUE, residuals are considered in variable input prediction
#' @param  \dots other arguments
#' @return Distribution of estimated crop input uses.
#' @export
print.rpinpall <- function(x, error = FALSE, ...) {
  message("Call:\n")
  print(x$call)
  if(error){
    xit       <- as.data.frame(x$xit_pred_with_err)
  }else{
    xit       <- as.data.frame(x$xit_pred)
  }

  m_mean   <- colMeans(xit[,-c(1,2)])
  m_std    <- sqrt(matrixStats::colVars(as.matrix(xit[,-c(1,2)])))
  m_dist   <- t(matrixStats::colQuantiles(as.matrix(xit[,-c(1,2)])))
  m_dist   <- rbind(m_mean,m_std,m_dist)
  rownames(m_dist) <- c("Mean","Std.dev.", "Min.", "1st Qu", "Median", "3rd Qu", "Max.")
  m_mean_y <- aggregate(xit[,-c(1,2)],list(xit[,2]), mean, na.rm=TRUE)
  colnames(m_mean_y)[1] <- c(colnames(xit)[2])
  message("\nDistribution of crop input uses: \n")
  print(m_dist)
  message("\nMean by year of crop input uses: \n")
  print(m_mean_y)
}

#' @describeIn rpinpallEst Displays a summary of estimated parameters
#' @param  object An object produced by the function \code{rpinpallEst}, to be displayed
#' @param  \dots Other arguments
#' @export
summary.rpinpall <- function(object, ...){

  coef_beta_b <- c(object$est_pop$par$beta_b)
  coef_beta_d <- c(object$est_pop$par$beta_d)
  coef_beta_m <- c(object$est_pop$par$beta_m)
  coef_domega_b <- c(diag(object$est_pop$par$omega_b))
  coef_domega_e <- c(diag(object$est_pop$par$omega_e))
  coef_omega_u <- c(object$est_pop$par$omega_u)
  se_beta_b <- c(object$est_stde$ecart_b)
  se_beta_d <- c(object$est_stde$ecart_d)
  se_beta_m <- c(object$est_stde$ecart_m)
  se_domega_b <- diag(object$est_stde$ecart_bb)
  se_domega_e <- diag(object$est_stde$ecart_e)
  se_omega_u <- diag(object$est_stde$ecart_u)
  tval_beta_b <- coef_beta_b / se_beta_b
  tval_beta_d <- coef_beta_d / se_beta_d
  tval_beta_m <- coef_beta_m / se_beta_m
  tval_domega_b <- coef_domega_b / se_domega_b
  tval_domega_e <- coef_domega_e / se_domega_e
  tval_omega_u <- coef_omega_u / se_omega_u

  TAB_beta_b <- cbind(Estimate = coef_beta_b,
                      StdErr = se_beta_b,
                      t_value = tval_beta_b)
  row.names(TAB_beta_b) <- paste0("Crop", sep="_", 1:nrow(TAB_beta_b))

  TAB_omega_b <- cbind(Estimate = coef_domega_b,
                       StdErr = se_domega_b,
                       t_value = tval_domega_b)
  row.names(TAB_omega_b) <- paste0("Crop", sep="_", 1:nrow(TAB_omega_b))

  TAB_omega_e <- cbind(Estimate = coef_domega_e,
                       StdErr = se_domega_e,
                       t_value = tval_domega_e)
  row.names(TAB_omega_e) <- paste0("Crop", sep="_", 1:nrow(TAB_omega_e))

  TAB_omega_u <- cbind(Estimate = coef_omega_u,
                       StdErr = se_omega_u,
                       t_value = tval_omega_u)
  row.names(TAB_omega_u) <- "omega"

  TAB_beta_d <- cbind(Estimate = coef_beta_d,
                      StdErr = se_beta_d,
                      t_value = tval_beta_d)
  row.names(TAB_beta_d) <- row.names(object$est_pop$par$beta_d)

  TAB_beta_m <- cbind(Estimate = coef_beta_m,
                      StdErr = se_beta_m,
                      t_value = tval_beta_m)
  row.names(TAB_beta_m) <- row.names(object$est_pop$par$beta_m)

  # R squared
  yit <- merge(object$yit_obs, object$yit_predict, by = colnames(object$yit_predict[,1:2]))
  ols <- lm(yit[,3] ~ yit[,4] )
  #ols <- lm(object$yit_obs[,3] ~ object$yit_pred[,3] )
  sim_r_squared <- matrix(summary(ols)$r.squared,1,1)
  rownames(sim_r_squared) <- "sim_r2"
  colnames(sim_r_squared) <- "Estimate"
  if(is.null(object$est_pop$par$beta_m)){
    res <- list(#call=object$call,
                beta_b=TAB_beta_b,
                diag_omega_b=TAB_omega_b,
                diag_omega_e=TAB_omega_e,
                omega_u=TAB_omega_u,
                beta_d=TAB_beta_d,
                sim_r_squared=sim_r_squared)
  }else{
    res <- list(#call=object$call,
                beta_b=TAB_beta_b,
                diag_omega_b=TAB_omega_b,
                diag_omega_e=TAB_omega_e,
                omega_u=TAB_omega_u,
                beta_d=TAB_beta_d,
                beta_m=TAB_beta_m,
                sim_r_squared=sim_r_squared)
  }
  res
}

#' @describeIn rpinpallEst Plot the "global" convergence indicator
#' @param  x An object produced by the function \code{rpinpallEst}, to be displayed
#' @param  \dots Other arguments
#' @export
plot.rpinpall <- function(x, ... ) {
  plot(x$conv_ind_cll, type="l")
  graphics::abline(v=c(x$opt$nb_SA, x$opt$nb_RS), col=c("red" , "blue"), lty=c(1,2), lwd=c(1, 3))
}

