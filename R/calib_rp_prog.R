#' @title Random Parameters Calibration
#' @description Calibration of individual parameters for each i given conditional
#'  distribution and estimated population parameters
#' @param data_list List of data
#' @param beta_list a list of simulated data obtained at previous iteration_
#' @param par List of population parameters
#' @param opt Options
#' @noRd
calib_rp <- function(data_list, beta_list, par, opt){
  nb_core <- future::availableCores() - 1
  nb_core <- ifelse(nb_core > 2, nb_core - 1, nb_core)
  nb_K    <- nrow(par$omega_b)
  nb_id   <- length(data_list)
  nb_obs  <- sum(sapply(1:nb_id, function(i) nrow(data_list[[i]]$y)))
  seq_i   <- seq_len(nb_id)
  boot_rp <- boot_rp_fact_calib(data_list, beta_list, par, opt)
  if(opt$doParallels==FALSE){
    results <- lapply(seq_i, function(i) boot_rp(i))
  }else{
    future::plan(future::multisession,workers=nb_core)
    results <- future.apply::future_lapply(seq_i, boot_rp,future.seed=TRUE)
  }

  ## Treatment of results
  beta_b_MAP_list   <- lapply(seq_i, function(x) results[[x]][["beta_b_MAP"]])
  beta_b_CM_list    <- lapply(seq_i, function(x) results[[x]][["beta_b_CM"]])
  beta_b_CSE_list   <- lapply(seq_i, function(x) results[[x]][["beta_b_CSE"]])
  beta_b_RSCD_list  <- lapply(seq_i, function(x) results[[x]][["beta_b_RSCD"]])
  ##############################################################################
  beta_b_calib_list <- NULL
  if(opt$calib_method == "cmean")
    beta_b_calib_list <- beta_b_CM_list
  if(opt$calib_method == "cmode")
    beta_b_calib_list <- beta_b_MAP_list
  if(opt$calib_method == "rscd")
    beta_b_calib_list <- beta_b_RSCD_list

  ## Treatment of results
  beta_i     <- beta_b_calib_list
  yit        <- lapply(1:nb_id, function(i){
    mydata_i      <- data_list[[i]]
    nb_T          <- nrow(mydata_i$y)
    trans_lambda  <- rp_model_fun(beta_i[[i]],opt)
    mat_SD        <- as.matrix(Matrix::bdiag(lapply(1:nb_T, function(t) matrix(mydata_i$x[t,], nrow=1))))
    val_temp      <- mat_SD %*% trans_lambda
    val           <- cbind(mydata_i$id,as.data.frame(val_temp))
    return(val)
  })
  ##############
  xit_without_err <- lapply(1:nb_id, function(i){
    mydata_i      <- data_list[[i]]
    trans_lambda  <- rp_model_fun(beta_i[[i]],opt)
    trans_lambda  <- matrix(trans_lambda, byrow = FALSE, nrow = nb_K)
    trans_lambda  <- t(trans_lambda)
    nb_T <- nrow(mydata_i$y)
    xi            <- cbind(mydata_i$id,as.data.frame(trans_lambda ))
    return(xi)
  })

  xit_with_err <- lapply(1:nb_id, function(i){
    mydata_i      <- data_list[[i]]
    trans_lambda  <- rp_model_fun(beta_i[[i]],opt)
    trans_lambda  <- matrix(trans_lambda, byrow = FALSE, nrow = nb_K)
    trans_lambda  <- t(trans_lambda)

    nb_T <- nrow(mydata_i$y)

    error <- mydata_i$y-yit[[i]][,3]
    blu_eps <- t(sapply(1:nb_T, function(t) error[t]*rep(1,ncol(mydata_i$x))))
    trans_lambda <- trans_lambda  + blu_eps
    xi           <- cbind(mydata_i$id,as.data.frame(trans_lambda))
    return(xi)
  })

  xit_approx_without_err <- Reduce("rbind", lapply(seq_i, function(i) xit_without_err[[i]]))
  xit_approx_with_err    <- Reduce("rbind", lapply(seq_i, function(i) xit_with_err[[i]]))
  yit_approx  <- Reduce("rbind", lapply(seq_i, function(i) yit[[i]]))

  if(opt$calib_method == "cmode" || opt$calib_method == "rscd"){
    result <- list(xit_pred_without_err = xit_approx_without_err, xit_pred_with_err=xit_approx_with_err,
                   yit_predict=yit_approx,
                   beta_b_calib_list = beta_b_calib_list)
  }else{
    result <- list(xit_pred_without_err = xit_approx_without_err, xit_pred_with_err=xit_approx_with_err,
                   yit_predict=yit_approx,
                   beta_b_calib_list = beta_b_calib_list, beta_b_CSE_list = beta_b_CSE_list)
  }
  return(result)
}
