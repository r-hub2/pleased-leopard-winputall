#' @title boot_rp_factory
#' @description boot_rp_factory is used to simulate the conditional distribution of random parameters
#' @param data_list List of individual data
#' @param beta_list a list of simulated data obtained at previous iteration_
#' @param par a list of values of parameters of population obtained at current iteration
#' @param opt a list of control parameters_ See rpinpall
#' @param temp temperature value at current iteration
#' @return boot_rp_factory returns a list of matrix
#' @noRd
boot_rp_factory <- function(data_list, beta_list, par , opt, temp) function(i){
  mydata_i    <- data_list[[i]]
  start_value <- beta_list[[i]]
  nbatch      <- opt$estim_rdraw
  w_i         <- c(mydata_i$w)
  nb_K        <- nrow(par$omega_b)
  nb_T        <- nrow(mydata_i$y)
  cvec_I      <- matrix(1,ncol=1, nrow = nb_T)
  R_i_tp      <- kronecker(diag(nb_T), par$omega_e*temp) + kronecker(tcrossprod(cvec_I), par$omega_b*temp)
  if(det(R_i_tp)==0){
    R_i   <-  R_i_tp + diag(.000001, nrow(R_i_tp) )
  }else{
    R_i   <-  R_i_tp
  }
  mat_ZD    <- Reduce("rbind", mydata_i$mat_ZD)
  if(is.null(par$beta_m)){
    m_beta_b  <- par$beta_b
    m_beta    <- kronecker(cvec_I,m_beta_b) + mat_ZD %*% par$beta_d
  }else{
    m_beta_b  <- par$beta_b + mydata_i$mat_MD%*%par$beta_m
    m_beta    <- kronecker(cvec_I,m_beta_b) + mat_ZD %*% par$beta_d
  }
  omega_u   <- c(par$omega_u * temp / w_i)
  ##########################################
  if(opt$distrib_method == "lognormal"){
    distrib_method <- 1
  }else if(opt$distrib_method == "normal"){
    distrib_method <- 2
  }else if(opt$distrib_method == "censored-normal"){
    distrib_method <- 3
  }

  standata <- list(T = nb_T, K = nb_K, y = c(mydata_i$y), alpha = c(m_beta),
                   omega = R_i, sigma = diag(omega_u) , big_s = mydata_i$mat_SD,
                   distrib_method = distrib_method)

  if(opt$sim_method == "nuts"){#nuts_stan
    sim_out  <- rstan::sampling(stanmodels$winputall_mod_ind, data = standata, chains = opt$estim_chains,
                                warmup = opt$nb_burn_sim, iter = nbatch + opt$nb_burn_sim , init = start_value,
                                refresh = 0, show_messages = FALSE, verbose = FALSE)
    list_of_draws <- rstan::extract(sim_out)
    sim_i     <- list_of_draws$mu
    loglike_c <- mean(list_of_draws$lp__)
  }else if(opt$sim_method == "variat"){#variat
    sim_out  <- rstan::vb(stanmodels$winputall_mod_ind, data = standata, output_samples = nbatch, refresh = 0 , importance_resampling = TRUE)
    list_of_draws <- rstan::extract(sim_out)
    sim_i     <- list_of_draws$mu
    loglike_c <- mean(Matrix::tail(list_of_draws$lupost, nbatch))
  }else if(opt$sim_method == "marg_imh"){

    standata_sim_imh <- c(standata, list(nb_sim = nbatch + opt$nb_burn_sim,
                                         m_beta_importance = c(m_beta),
                                         omega_importance = R_i,
                                         start_value = c(start_value)))

    out_sim <- rstan::sampling(stanmodels$imh_allx_ind,
                               standata_sim_imh,
                               algorithm = "Fixed_param" ,
                               warmup  = 0, chains = 1, iter = 1,refresh = 0,
                               show_messages = FALSE, verbose = FALSE)

    out_sim_list <- rstan::extract(out_sim)
    sim_i        <- Matrix::tail(apply(out_sim_list$beta_accept, 3, "rbind"), c(nbatch, opt$nb_rp))
    loglike_c    <- mean(Matrix::tail(out_sim_list$log_p_mean, nbatch))

  }else if(opt$sim_method == "lapl_approx"){
    out       <- rstan::optimizing(stanmodels$winputall_mod_ind, data = standata, init = "start_value",  draws = nbatch, importance_resampling = TRUE, as_vector = FALSE)
    sim_i     <- out$theta_tilde[,startsWith(colnames(out$theta_tilde), "mu")]
    loglike_c <- mean(out$log_p)
  }else if(opt$sim_method == "map_imh"){
    out_map_tp  <- try(rstan::optimizing(stanmodels$winputall_mod_ind, data = standata, init = "start_value",  hessian = TRUE, as_vector = FALSE), silent = TRUE)
    if(inherits(out_map_tp,"try-error")){
      m_beta_importance <- c(m_beta)
      omega_importance  <- R_i
    }else if(out_map_tp$return_code == 0){
      m_beta_importance <- out_map_tp$par$mu
      omega_importance_tp <- try(- solve(out_map_tp$hessian), silent = TRUE)
      if(inherits(omega_importance_tp,"try-error")){
        omega_importance <- R_i
      }else{
        omega_importance <- omega_importance_tp
      }

    }else if(out_map_tp$return_code != 0){
      m_beta_importance <- c(m_beta)
      omega_importance  <- R_i
    }

    # check if omega_importance is positive definite
    omega_importance <- (omega_importance + t(omega_importance))/2
    if(matrixcalc::is.positive.definite(omega_importance) == FALSE)
      omega_importance <- diag(diag(omega_importance))

    standata_sim_imh <- c(standata, list(nb_sim = nbatch + opt$nb_burn_sim,
                                         m_beta_importance = m_beta_importance,
                                         omega_importance = omega_importance,
                                         start_value = c(start_value)))

    out_sim <- rstan::sampling(stanmodels$imh_allx_ind,
                               standata_sim_imh,
                               algorithm = "Fixed_param" ,
                               warmup  = 0, chains = 1, iter = 1,refresh = 0,
                               show_messages = FALSE, verbose = FALSE)

    out_sim_list <- rstan::extract(out_sim)
    sim_i        <- Matrix::tail(apply(out_sim_list$beta_accept, 3, "rbind"), c(nbatch, opt$nb_rp))
    loglike_c    <- mean(Matrix::tail(out_sim_list$log_p_mean, nbatch))

  }else if(opt$sim_method == "mhrw"){

    standata_sim_mhrw <- c(standata, list(nb_sim = nbatch + opt$nb_burn_sim,
                                         omega_importance = R_i,
                                         start_value = c(start_value)))

    out_sim <- rstan::sampling(stanmodels$mhrw_allx_ind,
                               standata_sim_mhrw,
                               algorithm = "Fixed_param" ,
                               warmup  = 0, chains = 1, iter = 1,refresh = 0,
                               show_messages = FALSE, verbose = FALSE)

    out_sim_list <- rstan::extract(out_sim)
    sim_i        <- Matrix::tail(apply(out_sim_list$beta_accept, 3, "rbind"), c(nbatch, opt$nb_rp))
    loglike_c    <- mean(Matrix::tail(out_sim_list$log_p_mean, nbatch))
  }else if(opt$sim_method == "mhrw_imh"){
    low <- 1
    upper <- nbatch
    nbatch_rw <- low + (upper - low) * stats::plogis(stats::rnorm(1))
    nbatch_rw <- round(nbatch_rw)
    nbatch_ind <- nbatch - nbatch_rw

    standata_sim_mhrw <- c(standata, list(nb_sim = nbatch_rw + opt$nb_burn_sim,
                                          omega_importance = R_i,
                                          start_value = c(start_value)))
    out_sim <- rstan::sampling(stanmodels$mhrw_allx_ind,
                               standata_sim_mhrw,
                               algorithm = "Fixed_param" ,
                               warmup  = 0, chains = 1, iter = 1,refresh = 0,
                               show_messages = FALSE, verbose = FALSE)

    out_sim_list <- rstan::extract(out_sim)
    sim_i_rw        <- Matrix::tail(apply(out_sim_list$beta_accept, 3, "rbind"), c(nbatch_rw, opt$nb.rp))
    loglike_c_rw    <- Matrix::tail(out_sim_list$log_p_mean, nbatch_rw)

    standata_sim_imh <- c(standata, list(nb_sim = nbatch_ind + opt$nb_burn_sim,
                                         m_beta_importance = c(m_beta),
                                         omega_importance = R_i,
                                         start_value = c(start_value)))

    out_sim <- rstan::sampling(stanmodels$imh_allx_ind,
                               standata_sim_imh,
                               algorithm = "Fixed_param" ,
                               warmup  = 0, chains = 1, iter = 1,refresh = 0,
                               show_messages = FALSE, verbose = FALSE)


    out_sim_list   <- rstan::extract(out_sim)
    sim_i_ind      <- Matrix::tail(apply(out_sim_list$beta_accept, 3, "rbind"), c(nbatch_ind, opt$nb.rp))
    sim_i          <- rbind(sim_i_rw,sim_i_ind)
    loglike_c_ind  <- Matrix::tail(out_sim_list$log_p_mean, nbatch_ind)
    loglike_c      <- mean(c(loglike_c_ind, loglike_c_rw))
  }

  mm        <- matrix(1:ncol(sim_i),byrow = T, ncol = nb_K)
  sim_l_i   <- lapply(1:nb_T, function(t) sim_i[,mm[t,1]:mm[t,nb_K]])

  # Sufficient statistic
  V_i <- solve(solve(par$omega_b) + nb_T * solve(par$omega_e))
  m_i <- sapply(1:nrow(sim_i), function(r) {
    m_1 <- Reduce("+", lapply(1:nb_T, function(t) {
      m_1 <- sim_l_i[[t]][r,] - mydata_i$mat_ZD[[t]] %*% par$beta_d
    }))
    val <- V_i %*% (solve(par$omega_e) %*% m_1 + solve(par$omega_b) %*% m_beta_b )
  })
  m_i <- t(m_i)

  ss_rp_1   <- matrixStats::colMeans2(m_i)
  ss_rp_2   <- Matrix::crossprod(m_i)/nrow(m_i)
  ss_rp_3   <- Reduce("rbind",lapply(1:nb_T, function(t)  matrixStats::colMeans2(sim_l_i[[t]]- m_i)))
  ss_rp_4   <- Reduce("+",lapply(1:nb_T, function(t) Matrix::crossprod(sim_l_i[[t]]- m_i)/nrow(m_i)))

  ss_rp_5   <- 0
  for(r in 1:nrow(sim_i)){
    trans_lambda    <- rp_model_fun(sim_i[r,],opt)
    err_e         <- (mydata_i$y - mydata_i$mat_SD %*% trans_lambda)
    ss_rp_5       <- ss_rp_5 + t(err_e) %*% diag(w_i) %*% err_e/nrow(m_i)
  }
  ss_rp_6   <- matrixStats::colMeans2(sim_i)
  return(list( ss_rp_1             = ss_rp_1,
               ss_rp_2             = ss_rp_2,
               ss_rp_3             = ss_rp_3,
               ss_rp_4             = ss_rp_4,
               ss_rp_5             = ss_rp_5,
               ss_rp_6             = ss_rp_6,
               loglike_c           = loglike_c))
}

#' @title boot_rp_fact_calib
#' @description boot_rp_fact_calib is used to estimate the standard eror of estimated parameters_
#' @param data_list List of individual data_
#' @param beta_list a list of simulated data obtained at previous iteration_
#' @param par a list of values of parameters of population obtained at current iteration
#' @param opt a list of control parameters_ See rpinpall_
#' @return boot_rp_fact_calib returns a list of matrix
#' @noRd
boot_rp_fact_calib <- function(data_list, beta_list, par , opt) function(i) {
  temp           <- 1
  mydata_i       <- data_list[[i]]
  start_value    <- beta_list[[i]]
  nbatch         <- opt$calib_rdraw
  nb_K           <- nrow(par$omega_b)
  nb_T           <- nrow(mydata_i$y)
  w_i            <- c(mydata_i$w)
  cvec_I         <- matrix(1,ncol=1, nrow = nb_T)
  R_i_tp    <- kronecker(diag(nb_T), par$omega_e) + kronecker(tcrossprod(cvec_I),par$omega_b)
  if(det(R_i_tp)==0){
    R_i   <-  R_i_tp + diag(.000001, nrow(R_i_tp) )
  }else{
    R_i   <-  R_i_tp
  }
  mat_ZD    <- Reduce("rbind", mydata_i$mat_ZD)
  if(is.null(par$beta_m)){
    m_beta_b <- par$beta_b
    m_beta    <- kronecker(cvec_I,m_beta_b) + mat_ZD %*% par$beta_d
  }else{
    m_beta_b  <- par$beta_b + mydata_i$mat_MD%*%par$beta_m
    m_beta    <- kronecker(cvec_I,m_beta_b) + mat_ZD %*% par$beta_d
  }
  invw_i    <- c(1/w_i)
  omega_u   <- c(par$omega_u * invw_i)
  ##########################################
  if(opt$distrib_method == "lognormal"){
    distrib_method <- 1
  }else if(opt$distrib_method == "normal"){
    distrib_method <- 2
  }else if(opt$distrib_method == "censored-normal"){
    distrib_method <- 3
  }

  if(opt$doParallels == TRUE){
    calib_chains <- 1
  }else if(opt$calib_chains <= future::availableCores() - 1){
    calib_chains <- opt$calib_chains
  }else{
    calib_chains <- 1
    warning("calib_chains is upper than available cores. 1 chain is used")
  }

  if(opt$calib_method == "cmean" |  opt$calib_method == "rscd"){

    standata <- list(T = nb_T, K = nb_K, y = c(mydata_i$y), alpha = c(m_beta),
                     omega = R_i, sigma = diag(omega_u) , big_s = mydata_i$mat_SD, distrib_method = distrib_method)

    sim_out  <- rstan::sampling(stanmodels$winputall_mod_ind, data = standata, chains = calib_chains, cores = calib_chains,
                                iter = nbatch + opt$warmup, init = "start_value", warmup = opt$warmup,
                                refresh = 0, show_messages = FALSE, verbose = FALSE)

    list_of_draws <- rstan::extract(sim_out)
    sim_i     <- list_of_draws$mu
    selec    <- (1:calib_chains)*nbatch
    # Conditional mean and variance
    beta_b_RSCD <- colMeans(matrix(sim_i[selec,], nrow = calib_chains ))
    beta_b_CM   <- colMeans(sim_i)
    beta_b_CSE  <- sqrt(diag(crossprod(sim_i))/(calib_chains*nbatch))
    beta_b_MAP  <- NULL
  }else if(opt$calib_method=="cmode"){
    ## Conditional mode using optimizing of rstan
    standata <- list(T = nb_T, K = nb_K, y = c(mydata_i$y), alpha = c(m_beta),
                     omega = R_i, sigma = diag(omega_u) , big_s = mydata_i$mat_SD, distrib_method = distrib_method)
    out      <- rstan::optimizing(stanmodels$winputall_mod_ind, data = standata, init = "start_value",  hessian = TRUE, as_vector = FALSE)
    mu_map   <- out$par$mu
    var_map  <- solve(-out$hessian)
    beta_b_MAP      <- mu_map
    beta_b_CM       <- NULL
    beta_b_CSE      <- NULL
    beta_b_RSCD     <- NULL
  }
  return(list(beta_b_MAP = beta_b_MAP,
              beta_b_CM  = beta_b_CM,
              beta_b_CSE = beta_b_CSE ,
              beta_b_RSCD=beta_b_RSCD))
}


#' @title boot_rp_factory_stde
#' @description boot_rp_factory_stde is used to estimate the standard error of estimated parameters_
#' @param data_list List of individual data_
#' @param beta_list a list of simulated data obtained at previous iteration_
#' @param par a list of values of parameters of population obtained at current iteration
#' @param opt a list of control parameters_ See rpinpall_
#' @param dup_omega a duplication matrix_
#' @return boot_rp_factory_stde returns a list of matrix
#' @noRd
boot_rp_factory_stde <- function(data_list, beta_list, par, opt,dup_omega) function(i){
  mydata_i     <- data_list[[i]]
  start_value  <- beta_list[[i]]
  nbatch       <- opt$stde_rdraw
  nb_K         <- nrow(par$omega_b)
  nb_T         <- nrow(mydata_i$y)
  temp         <- 1
  cvec_I       <- matrix(1,ncol=1, nrow = nb_T)
  R_i_tp       <- kronecker(diag(nb_T), par$omega_e*temp) + kronecker(tcrossprod(cvec_I),par$omega_b*temp)
  if(det(R_i_tp)==0){
    R_i   <-  R_i_tp + diag(.000001, nrow(R_i_tp) )
  }else{
    R_i   <-  R_i_tp
  }
  mat_ZD    <- Reduce("rbind", mydata_i$mat_ZD)
  if(is.null(par$beta_m)){
    m_beta_b <- par$beta_b
    m_beta    <- kronecker(cvec_I,m_beta_b) + mat_ZD %*% par$beta_d
  }else{
    m_beta_b  <- par$beta_b + mydata_i$mat_MD%*%par$beta_m
    m_beta    <- kronecker(cvec_I,m_beta_b) + mat_ZD %*% par$beta_d
  }
  ###
  invw_i    <- c(1/mydata_i$w)
  omega_u   <- c(par$omega_u * temp * invw_i)
  ##########################################
  if(opt$distrib_method == "lognormal"){
    distrib_method <- 1
  }else if(opt$distrib_method == "normal"){
    distrib_method <- 2
  }else if(opt$distrib_method == "censored-normal"){
    distrib_method <- 3
  }

  standata <- list(T = nb_T, K = nb_K, y = c(mydata_i$y), alpha = c(m_beta),
                   omega = R_i, sigma = diag(omega_u) , big_s = mydata_i$mat_SD,
                   distrib_method = distrib_method)

  if(opt$sim_method == "nuts"){#nuts_stan
    sim_out  <- rstan::sampling(stanmodels$winputall_mod_ind, data = standata, chains = opt$estim_chains,
                                warmup = opt$nb_burn_sim, iter = nbatch + opt$nb_burn_sim , init = start_value,
                                refresh = 0, show_messages = FALSE, verbose = FALSE)
    list_of_draws <- rstan::extract(sim_out)
    sim_i     <- list_of_draws$mu
  }else if(opt$sim_method == "variat"){#variat
    sim_out  <- rstan::vb(stanmodels$winputall_mod_ind, data = standata, output_samples = nbatch, refresh = 0 , importance_resampling = TRUE)
    list_of_draws <- rstan::extract(sim_out)
    sim_i     <- list_of_draws$mu
  }else if(opt$sim_method == "marg_imh"){

    standata_sim_imh <- c(standata, list(nb_sim = nbatch + opt$nb_burn_sim,
                                         m_beta_importance = c(m_beta),
                                         omega_importance = R_i,
                                         start_value = c(start_value)))

    out_sim <- rstan::sampling(stanmodels$imh_allx_ind,
                               standata_sim_imh,
                               algorithm = "Fixed_param" ,
                               warmup  = 0, chains = 1, iter = 1,refresh = 0,
                               show_messages = FALSE, verbose = FALSE)

    out_sim_list <- rstan::extract(out_sim)
    sim_i        <- Matrix::tail(apply(out_sim_list$beta_accept, 3, "rbind"), c(nbatch, opt$nb_rp))
  }else if(opt$sim_method == "lapl_approx"){
    out       <- rstan::optimizing(stanmodels$winputall_mod_ind, data = standata, init = "start_value",  draws = nbatch, importance_resampling = TRUE, as_vector = FALSE)
    sim_i     <- out$theta_tilde[,startsWith(colnames(out$theta_tilde), "mu")]
  }else if(opt$sim_method == "map_imh"){
    out_map_tp  <- try(rstan::optimizing(stanmodels$winputall_mod_ind, data = standata, init = "start_value",  hessian = TRUE, as_vector = FALSE), silent = TRUE)
    if(inherits(out_map_tp,"try-error")){
      m_beta_importance <- c(m_beta)
      omega_importance  <- R_i
    }else if(out_map_tp$return_code == 0){
      m_beta_importance <- out_map_tp$par$mu
      omega_importance_tp <- try(- solve(out_map_tp$hessian), silent = TRUE)
      if(inherits(omega_importance_tp,"try-error")){
        omega_importance <- R_i
      }else{
        omega_importance <- omega_importance_tp
      }
    }else if(out_map_tp$return_code != 0){
      m_beta_importance <- c(m_beta)
      omega_importance  <- R_i
    }

    standata_sim_imh <- c(standata, list(nb_sim = nbatch + opt$nb_burn_sim,
                                         m_beta_importance = m_beta_importance,
                                         omega_importance = omega_importance,
                                         start_value = c(start_value)))

    out_sim <- rstan::sampling(stanmodels$imh_allx_ind,
                               standata_sim_imh,
                               algorithm = "Fixed_param" ,
                               warmup  = 0, chains = 1, iter = 1,refresh = 0,
                               show_messages = FALSE, verbose = FALSE)

    out_sim_list <- rstan::extract(out_sim)
    sim_i        <- Matrix::tail(apply(out_sim_list$beta_accept, 3, "rbind"), c(nbatch, opt$nb_rp))

  }else if(opt$sim_method == "mhrw"){

    standata_sim_mhrw <- c(standata, list(nb_sim = nbatch + opt$nb_burn_sim,
                                          omega_importance = R_i,
                                          start_value = c(start_value)))

    out_sim <- rstan::sampling(stanmodels$mhrw_allx_ind,
                               standata_sim_mhrw,
                               algorithm = "Fixed_param" ,
                               warmup  = 0, chains = 1, iter = 1,refresh = 0,
                               show_messages = FALSE, verbose = FALSE)

    out_sim_list <- rstan::extract(out_sim)
    sim_i        <- Matrix::tail(apply(out_sim_list$beta_accept, 3, "rbind"), c(nbatch, opt$nb_rp))
    loglike_c    <- mean(Matrix::tail(out_sim_list$log_p_mean, nbatch))
  }else if(opt$sim_method == "mhrw_imh"){
    low <- 1
    upper <- nbatch
    nbatch_rw <- low + (upper - low) * stats::plogis(stats::rnorm(1))
    nbatch_rw <- round(nbatch_rw)
    nbatch_ind <- nbatch - nbatch_rw

    standata_sim_mhrw <- c(standata, list(nb_sim = nbatch_rw + opt$nb_burn_sim,
                                          omega_importance = R_i,
                                          start_value = c(start_value)))
    out_sim <- rstan::sampling(stanmodels$mhrw_allx_ind,
                               standata_sim_mhrw,
                               algorithm = "Fixed_param" ,
                               warmup  = 0, chains = 1, iter = 1,refresh = 0,
                               show_messages = FALSE, verbose = FALSE)

    out_sim_list <- rstan::extract(out_sim)
    sim_i_rw        <- Matrix::tail(apply(out_sim_list$beta_accept, 3, "rbind"), c(nbatch_rw, opt$nb.rp))

    standata_sim_imh <- c(standata, list(nb_sim = nbatch_ind + opt$nb_burn_sim,
                                         m_beta_importance = c(m_beta),
                                         omega_importance = R_i,
                                         start_value = c(start_value)))

    out_sim <- rstan::sampling(stanmodels$imh_allx_ind,
                               standata_sim_imh,
                               algorithm = "Fixed_param" ,
                               warmup  = 0, chains = 1, iter = 1,refresh = 0,
                               show_messages = FALSE, verbose = FALSE)


    out_sim_list   <- rstan::extract(out_sim)
    sim_i_ind      <- Matrix::tail(apply(out_sim_list$beta_accept, 3, "rbind"), c(nbatch_ind, opt$nb.rp))
    sim_i          <- rbind(sim_i_rw,sim_i_ind)
  }

  mm        <- matrix(1:ncol(sim_i),byrow = T, ncol = nb_K)
  sim_l_i   <- lapply(1:nb_T, function(t) sim_i[,mm[t,1]:mm[t,nb_K]])

  V_i <- solve(solve(par$omega_b) + nb_T * solve(par$omega_e))
  m_i <- sapply(1:nrow(sim_i), function(r) {
    m_1 <- Reduce("+", lapply(1:nb_T, function(t) sim_l_i[[t]][r,] - mydata_i$mat_ZD[[t]] %*% par$beta_d ))
    val <- V_i %*% (solve(par$omega_e) %*% m_1 + solve(par$omega_b) %*% m_beta_b )
  })
  m_i <- t(m_i)

  ######
  score_beta_b  <- lapply(1:nrow(m_i), function(r)  solve(par$omega_b) %*% as.matrix(m_i[r,] - m_beta_b) )
  score_beta_b  <- Reduce("+",score_beta_b)/nrow(m_i)

  if(is.null(par$beta_m)){
    score_beta_m <- NULL
  }else{
    score_beta_m  <- lapply(1:nrow(m_i), function(r)  t(mydata_i$mat_MD ) %*% solve(par$omega_b) %*% as.matrix(m_i[r,]-par$beta_b) )
    score_beta_m  <- Reduce("+",score_beta_m)/nrow(m_i)
  }

  score_omega_b <- lapply(1:nrow(m_i), function(r){
    m1 <- solve(par$omega_b) - solve(par$omega_b) %*% tcrossprod(m_i[r,] - m_beta_b) %*% solve(par$omega_b) - solve(par$omega_b) %*% V_i %*% solve(par$omega_b)
    m1 <- -0.5 * t(dup_omega) %*% vec(m1)
  })
  score_omega_b <- Reduce("+",score_omega_b)/nrow(m_i)
  ######
  score_beta_d  <- lapply(1:nrow(m_i), function(r){
    m0 <- Reduce("+",lapply(1:nb_T, function(t) t(mydata_i$mat_ZD[[t]]) %*% solve(par$omega_e) %*% (sim_l_i[[t]][r,]- m_i[r,] - mydata_i$mat_ZD[[t]] %*% par$beta_d) ))
    m1 <- m0
  })
  score_beta_d <- Reduce("+",score_beta_d)/nrow(m_i)

  score_omega_e <- lapply(1:nrow(m_i), function(r){
    m0 <- Reduce("+",lapply(1:nb_T, function(t) (sim_l_i[[t]][r,]- m_i[r,] - mydata_i$mat_ZD[[t]] %*% par$beta_d) %*% t(sim_l_i[[t]][r,] - m_i[r,] - mydata_i$mat_ZD[[t]] %*% par$beta_d)))
    m1 <- nb_T * solve(par$omega_e) - solve(par$omega_e) %*% m0 %*% solve(par$omega_e) - nb_T * solve(par$omega_e) %*% V_i %*% solve(par$omega_e)
    m1 <- -0.5 * t(dup_omega) %*% vec(m1)
  })
  score_omega_e <- Reduce("+",score_omega_e)/nrow(m_i)

  # invw_i    <- c(1/mydata_i$w)# inverse of data weight
  # omega_u <- par$omega_u * invw_i
  score_omega_u <- lapply(1:nrow(m_i), function(r){
    trans_lambda  <- rp_model_fun(sim_i[r,],opt)
    err_r         <- mydata_i$y - mydata_i$mat_SD %*% trans_lambda
    m0            <- sum(err_r^2)
    m1            <- -0.5 * (nb_T /c(par$omega_u) - m0 / c(par$omega_u^2))
  })
  score_omega_u <- Reduce("+",score_omega_u)/nrow(m_i)


  return(list(  score_beta_b    =  score_beta_b,
                score_beta_m    =  score_beta_m,
                score_beta_d    =  score_beta_d,
                score_omega_b   =  score_omega_b,
                score_omega_u   =  score_omega_u,
                score_omega_e   =  score_omega_e))
}
