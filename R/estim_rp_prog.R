
#' @title Estimation of Variable Input Use per Crop
#'
#' @description This function allows estimating the parameters of the random
#'    parameters distribution: the population's parameters.
#' @param data_list list of individual data created in `rpinpallEst`.
#' @param par a list of values of parameters of population.
#' @param opt a list of control parameters. See 'rpinpall'
#' @return estim_rp returns a list of matrix
#' @noRd
estim_rp <- function(data_list, par, opt){

  start_func <- Sys.time()
  #
  nb_id      <- length(data_list)
  nb_obs     <- sum(sapply(1:nb_id, function(i) nrow(data_list[[i]]$y)))
  nb_K       <- nrow(par$omega_b)


  # initialization
  ss_rp_1_current   <- 0
  ss_rp_2_current   <- 0
  ss_rp_3_current   <- lapply(1:nb_id, function(i) 0)
  ss_rp_4_current   <- 0
  ss_rp_5_current   <- 0
  ss_rp_6_current   <- lapply(1:nb_id, function(i) 0)
  loglike_c_current <- 0
  iter              <- 0L
  nb_conv           <- 0L

  ## SA step length
  alpha_SA <- rep(0,opt$maxit)
  nb_SA <- opt$nb_SA + opt$nb_burn_saem
  alpha_SA[1:(nb_SA-1)] <- 1
  alpha_SA[nb_SA:opt$maxit] <- (1:(opt$maxiter- nb_SA + 1))^(-opt$p_SA)

  #
  seq_i <- seq_len(nb_id)
  value_current_par   <- 10

  #
  beta_b_all      <- NULL
  beta_d_all      <- NULL
  omega_b_all     <- NULL
  omega_e_all     <- NULL
  omega_u_all     <- NULL
  loglike_c_all   <- NULL

  beta_b  <- par$beta_b
  beta_m  <- par$beta_m
  beta_d  <- par$beta_d
  omega_b <- par$omega_b
  omega_e <- par$omega_e
  omega_u <- par$omega_u

  if (opt$showProgress) {
    pb <- txtProgressBar(min = 0, max = opt$maxit , initial = 0, style = 3)
  }
  #initialization
  beta_list <- lapply(1:nb_id, function(i){
    mydata_i  <- data_list[[i]]
    nb_T      <- nrow(mydata_i$y)
    cvec_I    <- matrix(1,ncol=1, nrow = nb_T)
    mat_ZD    <- Reduce("rbind", mydata_i$mat_ZD)

    if(is.null(par$beta_m)){
      m_beta_b  <- beta_b
      m_beta    <- kronecker(cvec_I,m_beta_b) +  mat_ZD %*% beta_d
    }else{
      m_beta_b  <- beta_b +  mydata_i$mat_MD%*%beta_m
      m_beta    <- kronecker(cvec_I,m_beta_b) +  mat_ZD %*% beta_d
    }
    sim_l <- m_beta
  })

  temp_all     <- tempering(opt$maxit,opt$doTemp_param) # Allassonni?re and  Chevallier 2021
  temp         <- 1
  max_nb_core  <- future::availableCores()

  while(nb_conv < 3L && iter < opt$maxiter){
    iter       <- iter + 1L
    alpha_iter <- alpha_SA[iter]
    if(opt$doTempering == TRUE && iter < opt$nb_RS){
      temp <- temp_all[iter]
    }else{
      temp <- 1
    }

    ##
    boot_rp   <- boot_rp_factory(data_list, beta_list, par, opt, temp)
    if(opt$doParallels==FALSE){
      results <- lapply(seq_i, function(i) boot_rp(i) )
    }else{
      nb_core   <- ifelse(opt$nb_core > max_nb_core, max_nb_core, opt$nb_core)
      future::plan(future::multisession, workers=nb_core)
      results <- future.apply::future_lapply(seq_i, boot_rp, future.seed = TRUE)
    }

    ## Treatment of results
    list_ss_rp_1     <- lapply(seq_i, function(x) results[[x]][["ss_rp_1"]])
    list_ss_rp_2     <- lapply(seq_i, function(x) results[[x]][["ss_rp_2"]])
    list_ss_rp_3     <- lapply(seq_i, function(x) results[[x]][["ss_rp_3"]])
    list_ss_rp_4     <- lapply(seq_i, function(x) results[[x]][["ss_rp_4"]])
    list_ss_rp_5     <- lapply(seq_i, function(x) results[[x]][["ss_rp_5"]])
    list_ss_rp_6     <- lapply(seq_i, function(x) results[[x]][["ss_rp_6"]])

    #
    list_loglike     <- lapply(seq_i, function(x) results[[x]][["loglike_c"]])
    loglike_c        <- loglike_c_current + alpha_iter * (Reduce("+",  list_loglike) - loglike_c_current)

    ss_rp_1          <- ss_rp_1_current  +  alpha_iter * (Reduce("rbind", list_ss_rp_1) - ss_rp_1_current)
    ss_rp_2          <- ss_rp_2_current  +  alpha_iter * (Reduce("+",list_ss_rp_2) - ss_rp_2_current)
    ss_rp_3          <- lapply(1:nb_id, function(i) ss_rp_3_current[[i]]  +  alpha_iter * (list_ss_rp_3[[i]]  - ss_rp_3_current[[i]]) )
    ss_rp_4          <- ss_rp_4_current  +  alpha_iter * (Reduce("+",list_ss_rp_4)  - ss_rp_4_current)
    ss_rp_5          <- ss_rp_5_current  +  alpha_iter * (Reduce("+",list_ss_rp_5)  - ss_rp_5_current)
    ss_rp_6          <- lapply(1:nb_id, function(i) ss_rp_6_current[[i]]  +  alpha_iter * (list_ss_rp_6[[i]]  - ss_rp_6_current[[i]]) )

    ###
    ss_rp_1_current  <- ss_rp_1
    ss_rp_2_current  <- ss_rp_2
    ss_rp_3_current  <- ss_rp_3
    ss_rp_4_current  <- ss_rp_4
    ss_rp_5_current  <- ss_rp_5
    ss_rp_6_current  <- ss_rp_6
    loglike_c_current <- loglike_c
    #Update of beta_list

    beta_list        <- ss_rp_6
    # b_i <- ss_rp_1

    # Updating of parameters_ Cconvergence control is applied after K_A

    if(iter > opt$nb_burn_saem){

      V_i_b   <- Reduce("+", lapply(1:nb_id, function(i) {
        nb_T  <- nrow(data_list[[i]]$y)
        val   <- solve(solve(omega_b) + nb_T * solve(omega_e))
      }))

      V_i_e   <- Reduce("+", lapply(1:nb_id, function(i) {
        nb_T  <- nrow(data_list[[i]]$y)
        val   <- nb_T*solve(solve(omega_b) + nb_T * solve(omega_e))
      }))

      for(j in 1:2){

        # Update of beta_b, beta_m and omega_b

        if(is.null(par$beta_m)){
          beta_b  <- as.matrix(colMeans(ss_rp_1))
          beta_m <- NULL
          omega_b <- ss_rp_2/nb_id - tcrossprod(beta_b) + V_i_b/nb_id
        }else{
          beta_b  <- as.matrix(colMeans(ss_rp_1)) - Reduce("+", lapply(1:nb_id, function(i) data_list[[i]]$mat_MD%*%beta_m))/nb_id
          sum1 <- Reduce("+", lapply(1:nb_id, function(i) t(data_list[[i]]$mat_MD)%*% solve(omega_b) %*% data_list[[i]]$mat_MD))
          sum2 <- Reduce("+", lapply(1:nb_id, function(i) t(data_list[[i]]$mat_MD)%*% solve(omega_b) %*% (as.matrix(ss_rp_1[i,]) - beta_b)))
          beta_m <- solve(sum1)%*%sum2
          sum1 <- Reduce("+", lapply(1:nb_id, function(i) (beta_b + data_list[[i]]$mat_MD%*%beta_m) %*% ss_rp_1[i,]))
          sum2 <- Reduce("+", lapply(1:nb_id, function(i) tcrossprod(beta_b + data_list[[i]]$mat_MD%*%beta_m)))
          omega_b <- (ss_rp_2 - sum1 -t(sum1) + sum2 )/nb_id +  V_i_b/nb_id
        }

        # Update of beta_d
        sum1 <- 0
        sum2 <- 0
        for (i in 1:nb_id) {
          nb_T     <- nrow(data_list[[i]]$y)
          for(t in 1:nb_T){
            sum1 <- sum1 + t(data_list[[i]]$mat_ZD[[t]])%*%solve(omega_e)%*%data_list[[i]]$mat_ZD[[t]]
            sum2 <- sum2 + t(data_list[[i]]$mat_ZD[[t]])%*%solve(omega_e)%*%as.matrix(ss_rp_3[[i]][t,])
          }
        }
        beta_d      <- solve(sum1)%*%sum2

        # Update of omega_e
        sum1 <- 0
        sum2 <- 0
        # sum3 <- 0
        for (i in 1:nb_id) {
          nb_T   <- nrow(data_list[[i]]$y)
          for(t in 1:nb_T){
            sum1   <- sum1   + data_list[[i]]$mat_ZD[[t]] %*% beta_d %*% t(as.matrix(ss_rp_3[[i]][t,]))
            sum2   <- sum2   + data_list[[i]]$mat_ZD[[t]] %*% beta_d %*% t(data_list[[i]]$mat_ZD[[t]] %*% beta_d)
          }
        }
        omega_e <- (ss_rp_4 - t(sum1) - sum1 + sum2)/nb_obs  + V_i_e/nb_obs

        if(opt$cov_rp_ISO)
          omega_e    <- mean(diag(omega_e))*diag(nb_K)
      }

      # check if omega_b is positive definite
      omega_b <- (omega_b + t(omega_b))/2
      if(matrixcalc::is.positive.definite(omega_b) == FALSE)
        omega_b <- diag(diag(omega_b))

      # check if omega_e is positive definite
      omega_e <- (omega_e + t(omega_e))/2
      if(matrixcalc::is.positive.definite(omega_e) == FALSE)
        omega_e <- diag(diag(omega_e))

      # update of omega_u
      omega_u  <- c(ss_rp_5)/nb_obs

      ##
      if(opt$doDiagEps=="mod0"){
        omega_b <- omega_b
        omega_e <- omega_e
      }else if(opt$doDiagEps=="mod1"){
        omega_b <- diag(diag(omega_b))
        omega_e <- omega_e
      }else if(opt$doDiagEps=="mod2"){
        omega_b <- omega_b
        omega_e <- diag(diag(omega_e))
      }else if(opt$doDiagEps=="mod3"){
        omega_b <- diag(diag(omega_b))
        omega_e <- diag(diag(omega_e))
      }

      #
      beta_b_all  <- rbind(beta_b_all,c(beta_b))
      beta_d_all  <- rbind(beta_d_all, c(beta_d))
      omega_b_all <- rbind(omega_b_all,diag(omega_b))
      omega_e_all <- rbind(omega_e_all,diag(omega_e))
      omega_u_all <- rbind(omega_u_all,omega_u)
      loglike_c_all <- rbind(loglike_c_all,loglike_c)

      par$beta_b  <- beta_b
      par$beta_m  <- beta_m
      par$beta_d  <- beta_d
      par$omega_u <- omega_u
      par$omega_b <- omega_b
      par$omega_e <- omega_e

      if(iter > opt$nb_SA){
        #value_new_par       <- loglike_c
        value_new_par       <- c(par$beta_b, diag(par$omega_b),diag(par$omega_e), par$beta_d, par$beta_m)

        conv_value_par      <- fun_converge(value_new_par, value_current_par, opt$tol )
        value_current_par   <- value_new_par
        if (conv_value_par  < opt$tol ){
          nb_conv   <- nb_conv + 1L
        }else {
          nb_conv   <- 0L
        }
      }else{
        conv_value_par <- NA
      }
      if(opt$showIterConvLL)
        message("iteration: ", iter,
                ", conv_value_par:  ", format(round(conv_value_par,3) , digits=2, nsmall=3),
                #", omega_u:  ", format(round(omega_u,3) , digits=2, nsmall=3),
                #", omega_b:  ", format(round(diag(omega_b)[1:2],3) , digits=2, nsmall=3),
                #", omega_e:  ", format(round(diag(omega_e),3) , digits=2, nsmall=3),
                ", loglike_compl: ", format(round(loglike_c,3) , digits=2, nsmall=3))
    }

    if (opt$showProgress)
      setTxtProgressBar(pb, iter)
  }
  if (opt$showProgress)
    close(pb)

  beta_i  <- beta_list
  yit     <- lapply(1:nb_id, function(i){
    mydata_i      <- data_list[[i]]
    nb_T          <- nrow(mydata_i$y)
    trans_lambda  <- rp_model_fun(beta_i[[i]],opt)
    mat_SD        <- as.matrix(Matrix::bdiag(lapply(1:nb_T, function(t) matrix(mydata_i$x[t,], nrow=1))))
    val_temp      <- mat_SD %*% trans_lambda
    val           <- cbind(mydata_i$id, as.data.frame(val_temp))
    return(val)
  })

  # prediction of crop input use

  xit_without_err <- lapply(1:nb_id, function(i){
    mydata_i      <- data_list[[i]]
    trans_lambda  <- rp_model_fun(beta_i[[i]],opt)
    trans_lambda  <- matrix(trans_lambda, byrow = FALSE, nrow = nb_K)
    trans_lambda  <- t(trans_lambda)
    xi            <- cbind(mydata_i$id,as.data.frame(trans_lambda ))

    return(xi)
  })

  xit_with_err <- lapply(1:nb_id, function(i){
    mydata_i      <- data_list[[i]]
    trans_lambda  <- rp_model_fun(beta_i[[i]],opt)
    trans_lambda  <- matrix(trans_lambda, byrow = FALSE, nrow = nb_K)
    trans_lambda  <- t(trans_lambda)

    nb_T <- nrow(mydata_i$y)

    error   <- mydata_i$y-yit[[i]][,3]
    blu_eps <- t(sapply(1:nb_T, function(t) error[t]*rep(1,ncol(mydata_i$x))))
    trans_lambda <- trans_lambda + blu_eps
    xi           <- cbind(mydata_i$id,as.data.frame(trans_lambda))
    return(xi)
  })

  xit_approx_without_err <- Reduce("rbind", lapply(seq_i, function(i) xit_without_err[[i]]))
  xit_approx_with_err    <- Reduce("rbind", lapply(seq_i, function(i) xit_with_err[[i]]))
  yit_approx             <- Reduce("rbind", lapply(seq_i, function(i) yit[[i]]))

  convergence_indicator <- list(beta_b_all=beta_b_all, beta_d_all=beta_d_all, omega_u_all=omega_u_all,
                                omega_e_all=omega_e_all, omega_b_all=omega_b_all,loglike_c_all=loglike_c_all)

  times <- Sys.time() - start_func

  results              <- NULL
  results$par          <- par
  results$beta_list    <- beta_list
  results$loglike_comp <- loglike_c
  results$conv_value_par <- conv_value_par
  results$nb_iter        <- iter
  results$xit_pred_without_err  <- xit_approx_without_err
  results$xit_pred_with_err     <- xit_approx_with_err
  results$yit_pred              <- yit_approx
  results$conv_indicators       <- convergence_indicator
  results$times                 <- times

  return(results)

}


#' @title stde_rp
#' @description stde_rp is used to estimate the standard eror of estimated parameters_
#' @param data_list List of individual data_
#' @param beta_list a list of simulated data obtained at previous iteration_
#' @param par a list of values of parameters of population obtained at current iteration
#' @param opt a list of control parameters_ See rpinpall_
#' @return stde_rp returns a list of matrix
#' @noRd
stde_rp <- function(data_list, beta_list, par, opt){
  nb_core     <- future::availableCores() - 1
  nb_core     <- ifelse(nb_core > 2, nb_core - 1, nb_core)
  nb_K        <- nrow(par$omega_b)
  nb_id       <- length(data_list)    #
  nb_obs      <- sum(sapply(1:nb_id, function(i) nrow(data_list[[i]]$y)))
  seq_i       <- seq_len(nb_id)
  dup_omega   <- matrixcalc::duplication.matrix(nb_K)
  ## Calibration de beta_b for each i given conditional distribution
  boot_rp     <- boot_rp_factory_stde(data_list,beta_list, par, opt,dup_omega)
  if(opt$doParallels==FALSE){
    results <- lapply(seq_i, function(i) boot_rp(i))
  }else{
    future::plan(future::multisession,workers=nb_core)
    results <- future.apply::future_lapply(seq_i, boot_rp,future.seed=TRUE)
  }
  ## Treatment of results
  list_score_beta_b     <- lapply(seq_i, function(x) results[[x]][["score_beta_b"]])
  list_score_beta_d     <- lapply(seq_i, function(x) results[[x]][["score_beta_d"]])
  list_score_omega_b    <- lapply(seq_i, function(x) results[[x]][["score_omega_b"]])
  list_score_omega_e    <- lapply(seq_i, function(x) results[[x]][["score_omega_e"]])
  list_score_omega_u    <- lapply(seq_i, function(x) results[[x]][["score_omega_u"]])
  #
  score_beta_b          <- t(do.call(cbind,list_score_beta_b))
  score_beta_d          <- t(do.call(cbind,list_score_beta_d))
  score_omega_b         <- t(do.call(cbind,list_score_omega_b))
  score_omega_e         <- t(do.call(cbind,list_score_omega_e))
  score_omega_u         <- t(do.call(cbind,list_score_omega_u))

  mean_score_b     <- colMeans(score_beta_b)
  variance_score_b <- crossprod(score_beta_b)/nb_id - tcrossprod(mean_score_b)
  variance_b       <- try( solve(variance_score_b)/nb_id, silent = TRUE)
  if(inherits(variance_b, "try-error")){
    ecart_b <- NULL
  }

  else
    ecart_b          <- sqrt(diag(variance_b))

  mean_score_d     <- colMeans(score_beta_d)
  variance_score_d <- crossprod(score_beta_d)/nb_id - tcrossprod(mean_score_d)
  variance_d       <- try(solve(variance_score_d)/nb_id, silent = TRUE)
  if(inherits(variance_d, "try-error"))
    ecart_d <- NULL
  else
    ecart_d          <- sqrt(diag(variance_d))

  mean_score_bb     <- colMeans(score_omega_b)
  variance_score_bb <- crossprod(score_omega_b)/nb_id - tcrossprod(mean_score_bb)
  variance_bb       <- try(solve(variance_score_bb)/nb_id, silent = TRUE)
  if(inherits(variance_bb, "try-error"))
    ecart_bb <- NULL
  else
    ecart_bb          <- sqrt(diag(variance_bb))

  ecart_bb          <- ks::invvech(ecart_bb)

  mean_score_e     <- colMeans(score_omega_e)
  variance_score_e <- crossprod(score_omega_e)/nb_id - tcrossprod(mean_score_e)
  variance_e       <- try(solve(variance_score_e)/nb_id, silent = TRUE)
  if(inherits(variance_e, "try-error"))
    ecart_e <- NULL
  else
    ecart_e          <- sqrt(diag(variance_e))
  ecart_e          <- ks::invvech(ecart_e)


  mean_score_u     <- mean(score_omega_u)
  variance_score_u <- crossprod(score_omega_u)/nb_id - tcrossprod(mean_score_u)
  variance_u       <- solve(variance_score_u)/nb_id
  ecart_u          <- sqrt(variance_u)

  ecart_m <- NULL
  if(!is.null(par$beta_m)){
    list_score_beta_m   <- lapply(seq_i, function(x) results[[x]][["score_beta_m"]])
    score_beta_m     <- t(do.call(cbind,list_score_beta_m))
    mean_score_m     <- colMeans(score_beta_m)
    variance_score_m <- crossprod(score_beta_m)/nb_id - tcrossprod(mean_score_m)
    variance_m       <- try(solve(variance_score_m)/nb_id, silent = TRUE)
    if(inherits(variance_m, "try-error"))
      ecart_m <- NULL
    else
      ecart_m          <- sqrt(diag(variance_m))
  }


  return(list(ecart_b=ecart_b, ecart_m=ecart_m, ecart_d=ecart_d,ecart_bb=ecart_bb, ecart_e=ecart_e, ecart_u=ecart_u))
}


