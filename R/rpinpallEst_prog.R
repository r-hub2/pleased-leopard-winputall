
#' @title Fit Input Allocation Random Parameters Model
#'
#' @description Designed to fit a random parameters input allocation
#'   model proposed in Koutchade et al., (2024) <https://hal.science/hal-04318163>.
#'   It provides crops input cost for each individual at each time and can account
#'   for Weighted Panel Data.
#'
#' @param data name of the data frame or matrix containing all the variables
#' included in the model.
#' @param id_time  first (individual) and second (time) level variables allowing
#' characterizing panel data.
#' @param total_input variable (name) containing the total input used  at farm
#' level per ha to be allocated to the different crops.
#' @param crop_acreage  list of variables containing the acreage of the different crops.
#' @param crop_indvar optional list of vector of (time-varying) variables
#' specific to each crop used to control for observed (individual and/or temporal)
#' characteristics in the estimation process. Default=NULL.
#' @param crop_rp_indvar optional list of vector of (time-constant) variables
#' specific to each crop used to control for observed time-constant
#' characteristics in the estimation process. Default=NULL.
#' @param weight optional variable containing weights of individual sample farms.
#' Default=NULL (equal weight is given to each farm). Default=NULL.
#' @param distrib_method assumption on the distribution of input use
#' per crop (x_kit): "normal", "lognormal" or "censored-normal". Default="lognormal".
#' @param sim_method method used to draw the random parameters in the simulation
#' step of the SAEM algorithm in the estimation process: "map_imh" (independant
#' Metropolis Hasting with Laplace approximation as proposal distribution),
#' "marg_imh" (independant Metropolis Hasting with marginal distribution of
#' random pararameters as proposal distribution), "mhrw" (Metropolis Hasting Random Walk),
#' "mhrw_imh" (combined "imh" and "mhrw"), "nuts","variat" and"lapl_approx". Default= "map_imh".
#' @param calib_method method used: "cmode" (conditional mode),
#' "cmean" (conditional mean). Default="cmode".
#' @param saem_control list of options for the SAEM algorithm. See 'Details
#' @param par_init list of some parameters' initialization.
#'
#' @details
#' An SAEM algorithm is used to perform the estimation of input uses per crop.
#' Different options can be specified by the user for this algorithm in the saem_control argument.
#' The saem_control argument is list that can supply any of the following component.
#'
#' * `nb_burn_saem`: Number of iterations of the burn-in phase where individual
#'    parameters are sampled from their conditional distribution using
#'    sim_method and the initial values for model parameters without update
#'    these parameters. Default=20.
#' * `nb_SA`: Number of iterations in the exploration phase where algorithm
#'    explore parameters space without memory. The parameter that controls the
#'    convergence of the algorithm is set to 1. Default=200.
#' * `nb_smooth`: Number of iterations in the smoothing phase.Default=200 and
#'    the parameter that controls the convergence of the algorithm is set to 0.85 by default.
#' *  `nb_RS`: Number of iterations where tempering approach is used
#' * `tol`: Tolerance value for the convergence. Default 1.10-3.
#' * `estim_rdraw`: Number of random draws using in the estimation process. Default=100
#' * `calib_rdraw`: Number of random draws using in the calibration process. Default=100
#' * `stde_rdraw`:  Number of random draws using for computation of estimation
#'    standard errors. Default=100
#' * `p_SA`: Parameter determining  step sizes in the Stochastic
#'    Approximation (SA) step. Must be comprise between 0 and 1. Default=0.85
#' * `doParallels`:  Logical.If TRUE a parallel processing is used when more
#'    than 2 cores are available. Default=FALSE
#' * `doTempering`: Logical. If TRUE the tempering approach proposed by
#'    (Allassonnière and Chevallier, 2021) is used to avoid convergence to local
#'    maxima. Default=TRUE
#' * `doDiagEps`: = "2",
#' * `showProgress`: Logical. If TRUE the evolution of the estimation process is
#'    displayed graphically at the bottom of the screen. Default=TRUE
#' * `showIterConvLL`: Logical. If TRUE iteration number and convergence value
#'    are displayed during the estimation process. Default=FALSE
#' @md
#'
#'
#'
#' @return This function returns a list with the following components:
#' \itemize{
#' \item {`xit_pred`: matrix of predicted crop input used per ha.}
#' \item {`xit_pred_with_error`: matrix of predicted crop input used per ha.}
#' \item {`yit_predict`: vector of predicted totat input used.}
#' \item {`est_pop` list of results of estimation: estimated parameters.}
#' \item {`est_stde`list of  parameters standard errors.}
#' \item {`call`: a copy of the function call.}
#' \item {`opt`: a list of saem algorithm control parameters.}
#' \item {`conv_ind_cll`: vector of convergence indicator.}
#' \item {`data_list`: a list of individual data used for estimation.}
#' }
#'
#' @references Koutchade, O. P., Carpentier A. and Femenia F. (2024).
#'
#' @examples
#' \donttest{
#' data(my_winputall_data)
#' mydata <- my_winputall_data
#' fit <- rpinpallEst(data = my_winputall_data,
#'                    id_time = c("id","year"),
#'                    total_input = "tx",
#'                    crop_acreage = c("s_crop1","s_crop2","s_crop3"),
#'                    distrib_method = "lognormal",
#'                    sim_method = "map_imh",
#'                    calib_method = "cmode",
#'                    saem_control = list(nb_SA = 10, nb_smooth = 10, estim_rdraw = 10))
#' print(fit)
#' plot(fit)
#' summary(fit)
#' head(fit$xit_pred)
#' }
#'
#' @export
rpinpallEst <-   function(data,
                          id_time,
                          total_input,
                          crop_acreage,
                          crop_indvar = NULL,
                          crop_rp_indvar = NULL,
                          weight = NULL,
                          distrib_method = c("lognormal", "normal", "censored-normal"),
                          sim_method = c("map_imh","mhrw","marg_imh","mhrw_imh","nuts","variat","lapl_approx"),
                          calib_method = c("cmode","cmean","rscd","estim-sim"),
                          saem_control = list(),
                          par_init = list()){

  # check if 'data' argument has been specified

  if (missing(data))
    data <- NULL

  no.data <- is.null(data)

  if (no.data) {
    data <- sys.frame(sys.parent())
  } else {
    if (!is.data.frame(data))
      data <- data.frame(data)
  }
  ID_1 <- ID_2 <- NULL
  data[,c("ID_1","ID_2")] <- data[,id_time]
  data$ID_1 <- as.factor(data$ID_1)
  data$ID_2 <- as.factor(data$ID_2)
  noDms <- names(data) # get variable name in data


  # create crop acreage share
  data$total_acreage <- rowSums(data[, crop_acreage])
  crop_acreage_sh <- paste0(crop_acreage, sep="_", "sh")
  if(any(crop_acreage_sh %in% names(data))){
    data = data[,!(names(data) %in% crop_acreage_sh)]
    data[, crop_acreage_sh] <- data[, crop_acreage] / rowSums(data[, crop_acreage])
  }else{
    data[, crop_acreage_sh] <- data[, crop_acreage] / rowSums(data[, crop_acreage])
  }

  # check that 'id_time', 'total_input', 'crop_acreage' arguments have been specified

  if(is.null(id_time))
    stop(" argument 'id_time' must be specified")
  if(length(noNms <- id_time[!id_time %in% noDms]))
    stop("unknown names in argument 'id_time': ", paste(noNms, collapse = ", "))

  if(is.null(total_input))
    stop(" argument 'total_input' must be specified")
  if(length(noNms <- total_input[!total_input %in% noDms]))
    stop("unknown names in argument 'total_input': ", paste(noNms, collapse = ", "))

  if(is.null(crop_acreage))
    stop(" argument 'crop_acreage' must be specified")
  if(length(noNms <- crop_acreage[!crop_acreage %in% noDms]))
    stop("unknown names in argument 'crop_acreage': ", paste(noNms, collapse = ", "))

  nb.K   <- length(crop_acreage)

  # check that 'weight' argument has been specified

  if(!is.null(weight) && length(noNms <- weight[!weight %in% noDms]))
    stop(" unknown names of argument 'weight': ", paste(noNms, collapse = ", "))


  # check that crop_indvar and crop_rp_indvar arguments have been specified

  unlistexogenous_i <- unlist(crop_rp_indvar)
  if(is.list(crop_rp_indvar) && length(noNms <- unlistexogenous_i[!unlistexogenous_i %in% noDms]))
    stop(" unknown names of 'crop_rp_indvar:  ",  paste(noNms, collapse = ", "))

  if(is.list(crop_rp_indvar) && length(crop_rp_indvar)!=nb.K)
    stop(" argument crop_rp_indvar is a list and its length may be equal to the number of crops", paste(nb.K))

  unlistexogenous <- unlist(crop_indvar)
  if(is.list(crop_indvar) && length(noNms <- unlistexogenous[!unlistexogenous %in% noDms]))
    stop(" unknown names of 'crop_indvar:  ",  paste(noNms, collapse = ", "))

  if(is.list(crop_indvar) && length(crop_indvar)!=nb.K)
    stop(" argument crop_indvar is a list and its length may be equal to the number of crops", paste(nb.K))

  # check that dim of total_input argument

  if(length(total_input)!=1)
    stop(" dim of 'total_input' should be 1  ")

  # match distrib_method and  sim_method

  distrib_method <- match.arg(distrib_method)
  sim_method     <- match.arg(sim_method)
  calib_method   <- match.arg(calib_method)

  # other arguments definition

  opt    <- list(nb_burn_saem = 10,
                 nb_SA = 200,
                 nb_smooth = 200,
                 nb_burn_sim = 50,
                 tol = 1.e-3,
                 nb_core = 4,
                 estim_rdraw = 100,
                 estim_thinning  = 1,
                 estim_chains = 1, warmup = 50,
                 calib_rdraw = 200,
                 calib_chains = 1,
                 stde_rdraw = 100,
                 p_SA = 0.85,
                 showIterConvLL = FALSE,
                 doTempering = FALSE,
                 doTemp_param = c(0.5,-4,2,4),
                 showProgress = TRUE,
                 doParallels = FALSE,
                 doDiagEps = "mod2",
                 cov_rp_ISO = FALSE)


  doDiagEps <- opt$doDiagEps
  if(length(noNms <- doDiagEps[!doDiagEps %in% c("mod0","mod1","mod2","mod3")]))
    stop(" unknown names of argument 'doDiagEps' in saem_control: ", paste(noNms, collapse = ", "))

  # check for saem_control argument and update opt

  nmsO <- names(opt)
  opt[(namo <- names(saem_control))] <- saem_control
  if (length(noNms <- namo[!namo %in% nmsO]))
    warning("unknown names in 'saem_control' arguments: ", paste(noNms, collapse = ", "))

  opt$maxiter <- opt$nb_burn_saem + opt$nb_SA + opt$nb_smooth
  opt$nb_RS <- opt$nb_burn_saem + round(opt$nb_SA/2)
  # check for parallel processing, doParallels and opt$nb_core arguments
  nb_core_c <- future::availableCores()
  if(is.numeric(opt$nb_core) && opt$nb_core > nb_core_c && opt$doParallels == TRUE)
    stop("nb_core must be lower than: ", paste(nb_core_c, collapse = ", "))

  if(opt$doParallels==TRUE && nb_core_c<=2) {
    opt$doParallels <- FALSE
    warning("serial processing is used: nb_core <=2")
  }

  opt$distrib_method  <- distrib_method
  opt$sim_method     <- sim_method
  opt$calib_method   <- calib_method


  # arrange data given id_time,

  data   <- dplyr::arrange(data, ID_1, ID_2) # unquote variable name

  # check if each individual is present at least 3 times in data

  data <- dplyr::mutate(dplyr::group_by(data, dplyr::across(all_of(id_time[1]))), rowCount = dplyr::n())

  if(any(data$rowCount < 3)){
    data    <- subset(data, data$rowCount >= 3 )
    warning("\nonly individual presents at least 3 times are considered \n")
  }

  # check that 'id_time', 'total_input', 'crop_acreage' arguments have no missings

  if (any(is.na(data[,id_time])))
    stop("argument 'id_time' should not contain any NAs.")

  # check that 'total_input' arguments has no missings

  if (any(is.na(data[,total_input])))
    stop("argument total_input should not contain any NAs.")

  # check that  'crop_acreage' arguments has no missings

  if (any(is.na(data[,crop_acreage])))
    stop("argument 'crop_acreage'  should not contain any NAs.")

  # check that exogenous argument has no missings

  if (any(is.na(data[,unlist(crop_indvar)])))
    stop("argument crop_indvar should not contain any NAs.")

  # check that crop_indvar argument has no missings

  if (any(is.na(data[,unlist(crop_indvar)])))
    stop("argument 'crop_indvar' should not contain any NAs.")

  if (any(is.na(data[,unlist(crop_rp_indvar)])))
    stop("argument 'crop_rp_indvar' should not contain any NAs.")

  # check that exogenous argument has no missings


  # # make sure id is a character variable
  #
  # id <- as.character(id)

  # check that


  # create dummy year

  Ndid_time2 <- NULL
  fid_time2  <- factor(as.matrix(data[,"ID_2"]))
  did_time2  <-  model.matrix(~fid_time2) # Dummy id_time2
  did_time2  <-  did_time2[,-1]
  data       <- data.frame(data,did_time2)
  Ndid_time2 <- colnames(did_time2)

  # check if 'x' is time-invariant within subjects

  is.constant <- function(x, na.rm=TRUE) {
    if (all(is.na(x))) {
      TRUE
    } else {
      res <- all(x == na.omit(x)[1], na.rm=na.rm)
      res[is.na(res)] <- FALSE
      res
    }
  }

  if(!is.null(crop_indvar)){

    vconst <- sapply(1:length(unlistexogenous), function(t){
      const <- tapply(data[,unlistexogenous[t]], data[,id_time[1]], is.constant, na.rm=TRUE)

      # check if 'x' is time-invariant for all subjects
      all.const <- all(const)
    })

    # get indvar time constant exotcont
    if(any(vconst))
      stop("arguments 'crop_indvar' should be time-varying or NULL.")
  }

  if(!is.null(crop_rp_indvar)){

    vconst <- sapply(1:length(unlistexogenous_i), function(t){
      const <- tapply(data[,unlistexogenous_i[t]], data[,id_time[1]], is.constant, na.rm=TRUE)

      # check if 'x' is time-invariant for all subjects
      all.const <- all(const)
    })

    # get indvar time constant exotcont
    if(any(!vconst))
      stop("arguments 'crop_rp_indvar' should be time-constant or NULL.")

  }

  # get names used to create individual data

  y.name  <- total_input
  x.name  <- crop_acreage_sh
  m_name  <- crop_rp_indvar
  z_name  <- lapply(1:nb.K, function(k)  {
    if(is.null(crop_indvar)){
      c(Ndid_time2)
    }else{
      c(crop_indvar[[k]])
    }
  })


  # Creation of data for each individual

  # message("\nData transformation step  \n")
  data      <- dplyr::arrange(data,ID_1, ID_2)
  hh        <- levels(factor(as.matrix(data[,id_time[1]])))
  nb.id     <- length(hh)
  data_list <- NULL
  X.mat     <- NULL
  Y.mat     <- NULL
  data_balanced <- data[,id_time]
  data_balanced <- dplyr::mutate(data_balanced, IND=1)
  data_balanced <- plm::make.pbalanced(data_balanced)
  data_balanced[is.na(data_balanced[,"IND"]),"IND"] <- 0
  for (i in 1:nb.id) {
    id    <- as.matrix(data[data[,id_time[1]]==hh[i], id_time])
    ireg  <- as.matrix(data_balanced[data_balanced[,id_time[1]]==hh[i], "IND"])
    y     <- as.matrix(data[data[,id_time[1]]==hh[i], y.name])
    x     <- as.matrix(data[data[,id_time[1]]==hh[i], x.name])
    z     <- lapply(1:nb.K, function(k) as.matrix(data[data[,id_time[1]]==hh[i], z_name[[k]]]))

    S_D   <- as.matrix(Matrix::bdiag(lapply(1:nrow(x), function(t) matrix(x[t,], nrow = 1))))
    Z_D   <- lapply(1:nrow(y), function (t) t(as.matrix(Matrix::bdiag(lapply(1:nb.K, function(k) z[[k]][t,]))) ))

    m     <- NULL
    M_D   <- NULL
    if(!is.null(crop_rp_indvar)){
      m   <- lapply(1:nb.K, function(k) as.matrix(data[data[,id_time[1]]==hh[i], m_name[[k]]])[1,])
      M_D <- t(as.matrix(Matrix::bdiag(lapply(1:nb.K, function(k) m[[k]]))) )
    }

    if(!is.null(weight)){
      w   <- as.matrix(data[data[,id_time[1]]==hh[i], weight])
    }else{
      nb.T <- nrow(y)
      w    <- rep(1,nb.T)
    }

    data_list[[i]] <-  list(y      = y,
                            x      = x,
                            w      = w,
                            z      = z,
                            m      = m,
                            mat_SD = S_D,
                            mat_ZD = Z_D,
                            mat_MD = M_D,
                            id     = id,
                            ireg   = ireg)
  }

  # Initialization
  Y.mat   <- as.matrix(data[,y.name])
  X.mat   <- as.matrix(data[,x.name])
  est_ols <- lm(Y.mat ~  -1 +  X.mat)
  #summary(est_ols)
  beta.b.0  <- est_ols$coef
  mean_all <- mean(Y.mat)
  var_all <- stats::var(Y.mat)

  beta.b.0.min <- min(beta.b.0[which(beta.b.0>0)])
  for(l in 1:length(x.name)){
    beta.b.0[l] <- ifelse(beta.b.0[l]<=0,beta.b.0.min/2,beta.b.0[l])
  }

  if(opt$distrib_method %in% c("normal", "censored-normal")){
    beta.b  <- as.matrix(c(beta.b.0))
    omega.e <- diag(nb.K)*(beta.b.0/2)^2
    omega.b <- diag(nb.K)*(beta.b.0/2)^2
  }


  if(opt$distrib_method == "lognormal"){
    var_norm <- (beta.b.0/2)^2
    mean_norm <- beta.b.0
    var_log     <- log(1+var_norm/(mean_norm^2))
    mean_log    <- log(beta.b.0) - 0.5*var_log
    beta.b  <- mean_log#as.matrix(log(beta.b.0))
    omega.e <- 0.05*diag(nb.K)#diag(var_log)#0.1*diag(nb.K) #
    omega.b <- 0.05*diag(nb.K)#diag(var_log)#0.1*diag(nb.K)
  }

  beta.d  <- matrix(0,nrow = sum(lengths(z_name)))
  rownames(beta.d) <- unlist(z_name)
  beta_m <- NULL
  if(!is.null(crop_rp_indvar)){
    beta_m <- matrix(0,nrow = sum(lengths(m_name)))
    rownames(beta_m) <- unlist(m_name)
  }

  omega.u <- sum(est_ols$residuals^2)/nrow(Y.mat) #0.1


  par <- NULL
  par$beta_b    <- beta.b
  par$beta_m    <- beta_m
  par$beta_d    <- beta.d
  par$omega_b   <- omega.b
  par$omega_e   <- omega.e
  par$omega_u   <- omega.u

  ## check for par_init argument and update par
  nmsO <- names(par)
  par[(namo <- names(par_init))] <- par_init
  if (length(noNms <- namo[!namo %in% nmsO]))
    stop("unknown names in 'par_init' arguments: ", paste(noNms, collapse = ", "))

  ###########################
  est     <- list()
  est$call <- match.call() # get the original call
  #
  # message("\nEstimation Population Parameters step  \n")
  # est_pop_tp  <- try(estim_rp(data_list, par, opt))
  # if(inherits(est_pop_tp, "try-error")){
  #   est_pop <- NULL
  #   stop( "Problem in estimation   \n" )
  est_pop_tp  <- suppressWarnings(tryCatch(estim_rp(data_list, par, opt)))
  if(inherits(est_pop_tp, "error")){
    est_pop <- NULL
    stop( "Problem in estimation   \n" )
  }else{
    # message("Standard Error Calculation step \n")
    est_pop      <- est_pop_tp
    rownames(est_pop$par$beta_b) <- paste0("Beta_crop",sep="_", 1:nb.K)
    rownames(est_pop$par$beta_m) <- unlist(m_name)
    rownames(est_pop$par$beta_d) <- unlist(z_name)
    par          <- est_pop$par
    beta_list    <- est_pop$beta_list
    est_stde_tp  <- suppressWarnings(try(stde_rp(data_list, beta_list, par, opt), silent = TRUE))
    if(inherits(est_stde_tp, "try-error")) {
      est_stde <- NULL
      warning( "Problem in standard error computation  \n" )
    }else{
      est_stde <- est_stde_tp
    }

    #
    # message("Calibration step  \n")
    if(calib_method != "estim-sim"){
      est_calib_rp_tp <- suppressWarnings(try(calib_rp(data_list, beta_list, par, opt), silent = TRUE))
      if(inherits(est_calib_rp_tp, "try-error")) {
        xit_pred_without_err <- est_pop$xit_pred_without_err
        xit_pred_with_err    <- est_pop$xit_pred_with_err
        colnames(xit_pred_with_err) <- colnames(xit_pred_without_err) <- c(id_time,paste0(x.name, sep="_", "X"))
        yit_predict <- est_pop$yit_pred
        colnames(yit_predict) <- c(id_time,y.name)
        warning( "Problem in calibration step: last simulation in the estimation is reported \n" )
      }else{
        xit_pred_without_err <- est_calib_rp_tp$xit_pred_without_err
        xit_pred_with_err    <- est_calib_rp_tp$xit_pred_with_err
        colnames(xit_pred_with_err) <- colnames(xit_pred_without_err) <- c(id_time,paste0(x.name, sep="_", "X"))
        yit_predict <- est_calib_rp_tp$yit_pred
        colnames(yit_predict) <- c(id_time,y.name)
      }
    }else{
      xit_pred_without_err        <- est_pop$xit_pred_without_err
      xit_pred_with_err           <- est_pop$xit_pred_with_err
      colnames(xit_pred_with_err) <- colnames(xit_pred_without_err) <- c(id_time,paste0(x.name, sep="_", "X"))
      yit_predict                 <- est_pop$yit_pred
      colnames(yit_predict)       <- c(id_time,y.name)
    }

  }

  xit_pred_without_err   <- merge(data[,id_time],xit_pred_without_err, by = id_time)
  xit_pred_with_err      <- merge(data[,id_time],xit_pred_with_err, by = id_time)
  yit_predict            <- merge(data[,id_time],yit_predict, by=id_time)

  est$xit_pred             <- xit_pred_without_err
  est$xit_pred_without_err <- xit_pred_without_err
  est$xit_pred_with_err    <- xit_pred_with_err
  est$yit_predict          <- yit_predict
  est$yit_obs              <- data[,c(id_time,total_input)]
  est$est_pop              <- est_pop
  est$conv_ind_cll         <- est$est_pop$conv_indicators$loglike_c_all
  est$est_stde             <- est_stde
  est$data_list            <- data_list
  est$data                 <- data
  est$opt                  <- opt
  # Name
  class(est) <- "rpinpall"
  est
}


#' @title Predict random parameter
#' @description  rpinpallPred is used to predict
#' @param x an object of class 'rpinpall' that contains data_list, beta_list, est_pop, xit_predict, yit_predict and opt.
#' @param calib_method  method used to calibrate individual input uses per crop: conditional mean ("cmean") or conditional mode ("cmode")
#' @param saem_control a list of saem_control parameters for prediction: rdraw.calib
#' @return a list of matrix
#' \itemize{
#' \item {`xit_predict`: matrix of predicted crop input used per ha.}
#' \item {`xit_predict_with_error`: matrix of predicted crop input used per ha.}
#' \item {`yit_predict`: vector of predicted totat input used.}
#' }
#' @noRd
rpinpallPred <- function(x, calib_method=c("cmode","cmean","rscd"),
                         saem_control=list()){

  # check if x belong class rpinpall

  if(!inherits(x, "rpinpall"))
    stop(" argument must belong to class 'rpinpall': ", paste(x, collapse = ", "))

  # get from the rpinpall object
  data_list <- x$data_list
  beta_list <- x$est_pop$beta_list
  par       <- x$est_pop$par

  #
  opt <- x$opt

  nmsO <- names(opt)
  opt[(namo <- names(saem_control))] <- saem_control
  if (length(noNms <- namo[!namo %in% nmsO]))
    stop("unknown names in 'saem_control' arguments: ", paste(noNms, collapse = ", "))

  calib_method  <- match.arg(calib_method)
  opt$calib_method <- calib_method

  # message("Calibration Step \n")
  est_calib_tp    <- suppressWarnings(try(calib_rp(data_list, beta_list, par, opt) ))
  if(inherits(est_calib_tp, "try-error")) {
    est_calib <- NULL
    warning( "Problem in calibration \n" )
  }else{
    est_calib <- est_calib_tp
  }

  est_calib$xit_pred                       <- est_calib$xit_pred_without_err
  colnames(est_calib$xit_pred)             <- colnames(x$xit_pred_without_err)
  colnames(est_calib$xit_pred_without_err) <- colnames(x$xit_pred_without_err)
  colnames(est_calib$xit_pred_with_err)    <- colnames(x$xit_pred_with_err)
  colnames(est_calib$yit_predict)          <- colnames(x$yit_predict)
  class(est_calib) <- "rpinpall"
  est_calib
}
#' ````

