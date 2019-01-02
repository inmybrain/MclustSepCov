##--------------------------------------------------------------------------
## Miscellaneous functions for tidy coding
##--------------------------------------------------------------------------

## (190102) source from https://stackoverflow.com/a/47045368/4008219
if(require("rstudioapi", quietly = TRUE)){
  getScriptTitle <- function(){
    path <- rstudioapi::getSourceEditorContext()$path
    filename <- substr(path, rev(gregexpr("/", path)[[1]])[1] + 1, nchar(path))
    # res <- substr(filename, 1, rev(gregexpr("_v", filename)[[1]])[1] - 1)
    return(filename)
  }
} else{
  getScriptTitle <- function(){
    path <- commandArgs() %>% 
      tibble::as.tibble() %>%
      tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
      dplyr::filter(key == "--file") %>%
      dplyr::pull(value)
    
    filename <- substr(path, rev(gregexpr("/", path)[[1]])[1] + 1, nchar(path))
    return(filename)
  }
}


CheckFile <-
  function(filename,
           fileext,
           filepath,
           cnt = 1,
           today = format(Sys.Date(), "%y%m%d")) {
    if (nchar(fileext) > 0) {
      fullname <- paste0(today, "_", filename, "_", cnt, ".", fileext)
      while (any(fullname %in% list.files(path = filepath))) {
        cnt <- cnt + 1
        fullname <-
          paste0(today, "_", filename, "_", cnt, ".", fileext)
      }
    } else{
      # file should be a folder
      fullname <- paste0(today, "_", filename, "_", cnt)
      while (any(fullname %in% list.files(path = filepath))) {
        cnt <- cnt + 1
        fullname <- paste0(today, "_", filename, "_", cnt)
      }
    }
    return(paste0(filepath, fullname))
  }


cat_parameter <- function(df_para){
  df_para <- data.frame(df_para)
  cat(sprintf("---------------------------------\n"))
  cat(sprintf("Parameters are varying as follows\n"))
  for(i in 1:ncol(df_para)){
    if(is.numeric(df_para[,i])){
      value <- paste(sprintf("%.3f", unique(df_para[,i])), collapse = ", ")
    } else{
      value <- paste(unique(df_para[,i]), collapse = ", ")
    }
    cat(sprintf("%s = %s\n", 
                colnames(df_para)[i],
                value
    ))
  }
  cat(sprintf("---------------------------------\n"))
  invisible()
}

##--------------------------------------------------------------------------
## Data generation
##--------------------------------------------------------------------------
getSynData <- function(r,
                       p,
                       q,
                       K,
                       wt_cluster1,
                       type_mean,
                       mean_len,
                       type_vari,
                       type_time,
                       rho_vari = NULL,
                       rho_time = NULL)
{
  ## Check type_vari, type_time and generate the covariance matrix
  if(type_vari == type_time & type_vari == "INDEP"){
    Sigma <- lapply(1:K, function(noarg) diag(rep(1, p * q)))
  } else if(type_vari == type_time & type_vari == "O.AR"){
    Sigma <- lapply(1:K, function(noarg) getCovariance(p * q, rho = 0.3, "AR"))
  } else if(type_vari == type_time & type_vari == "O.CS"){
    Sigma <- lapply(1:K, function(noarg) getCovariance(p * q, rho = 0.3, "CS"))
  } else{ # Separable covariance case
    
    ## Check 0 < rho_vari, rho_time < 1
    if(any(rho_vari > 1 | rho_vari < 0 | rho_time > 1 | rho_time < 0)){
      stop("rho_vari or rho_time has a parameter not in (0,1).\n")
    }
    
    U <- lapply(rho_vari, function(rho) getCovariance(p, rho, "AR"))
    if(type_time == "VAR"){# V = "VAR"
      V <- lapply(rho_time, function(rho) getCovariance(q, rho, "AR"))
    } else if(type_time == "VCS"){# V = "VCS"
      V <- lapply(rho_time, function(rho) getCovariance(q, rho, "CS"))
    } else if(type_time == "EAR"){# V = "EAR"
      V <- lapply(rep(rho_time, K), function(rho) getCovariance(q, rho, "AR"))
    } else if(type_time == "ECS"){# V = "ECS"
      V <- lapply(rep(rho_time, K), function(rho) getCovariance(q, rho, "CS"))
    } else{
      stop(sprintf("type_time=\"%s\" is not available in data generation.\n", type_time))
    }
    Sigma <- Map(kronecker, U, V)
  }
  
  if(K == 2){
    wt_cluster <- c(wt_cluster1, 1 - wt_cluster1)
    if(type_mean == "one-vector"){
      mu <- list(rep(0, p*q), mean_len / sqrt(p*q) * rep(1, p*q))
    } else if (type_mean == "step-vector"){
      stepvec <- rep(c(1,0), times = c(round(p*q/4), p*q - round(p*q/4)))
      mu <- list(rep(0, p*q), mean_len / sqrt(sum(stepvec)) * stepvec)
    } else if (type_mean == "maxvar-vector"){
      mu <- list(rep(0, p*q), mean_len * eigen(Sigma)$vectors[,1])
    }
  } else{
    stop("K > 2 is not supported.\n")
  }
  if(sum(round(r * wt_cluster)) != r){
    stop("Please set wt_cluster * r be integers.\n")
  }
  
  id_cluster <- rep(1:K, times = r * wt_cluster)
  id_cluster_tab <- table(id_cluster)
  
  Y <- vector(mode = "list", length = K)
  for(i in 1:K){
    Y[[i]] <- mvtnorm::rmvnorm(n = id_cluster_tab[names(id_cluster_tab) == i], 
                      mean = mu[[i]], 
                      sigma = Sigma[[i]])
  }
  return(list(
    list_Y = Y,
    Y = Reduce(rbind, Y),
    info = list(mu = mu,
                Sigma = Sigma,
                id_cluster = id_cluster,
                r = r,
                p = p,
                q = q,
                K = K,
                wt_cluster1 = wt_cluster1,
                type_mean = type_mean,
                mean_len = mean_len,
                type_vari = type_vari,
                type_time = type_time,
                rho_vari = rho_vari,
                rho_time = rho_time
    )
  )
  )
}

##--------------------------------------------------------------------------
## Model-based clustering 
##--------------------------------------------------------------------------

#' @name Mclust_SEP_cpp
#' @title The model-based clustering for longitudinal data
#' @description This is a wrapper function of \code{\link{Mclust_SEP_each_cpp}}. All arguments except \code{save_fit} will be passed to \code{\link{Mclust_SEP_each_cpp}}.
#' @param Y A r-by-(\code{p}*\code{q}) matrix where r is the sample size, and ordering of columns should be carefully set (see `Details').
#' @param p,q An integer value for the number of multi-variables and the number of time points, respectively.
#' @param Ks A sequence of positive integers indicating the number of mixture components, each of which will be used in \code{K} of \code{\link{Mclust_SEP_each_cpp}}.
#' @param type_cov A sequence of character strings indicating covariance structures, each of which will be used in \code{\link{Mclust_SEP_each_cpp}}. Default is \code{'all'}, which runs all available models. See `Details'.
#' @param tol Tolerance constant for convergence. Default is \code{1e-3}.
#' @param maxit Maximum number of iterations. Default is \code{500}.
#' @param save_fit A logical value indicating whether to save all fitted mixture models. If \code{FALSE}, the best model is only available by \code{best_model}. Default is \code{TRUE}.
#' @return A list with components:
#' \item{best_model}{A list of the mixture models with the largest BIC. If there is a tie, they are all returned.}
#' \item{bic_table}{Table filled with BIC values. Type of covariance models are given by rows and values in \code{Ks} by columns.}
#' \item{res_Mclust_SEP}{If \code{save_fit} is TRUE, all fitted models are stored. It is a nested list with the first layer corresponding to covariance models specified in \code{type_cov} and the second to values of \code{Ks}.}
#' @details 
#' The first \code{q} components from each row of \code{Y} denote \code{q} variables at time point 1, the second \code{q} are those at time point 2, and so on until time point \code{p}. Under separability, the covariance matrix of row vectors of \code{Y} is represented by \eqn{U_{p\times p} \otimes V_{q \times q}} for some covariance factors \eqn{U_{p\times p}, V_{q \times q}}.\cr\cr
#' \code{type_cov} should be in ``XXX-YYY'' format. ``XXX'' is for the multivariable covariance \eqn{U_{p\times p}}, and ``YYY'' for the temporal covariance \eqn{V_{q \times q}}. They will be passed respectively to \code{type_vari} and \code{type_time} in \code{\link{Mclust_SEP_each_cpp}}. Available options are as follows;
#' \itemize{
#'   \item Heteroscadatsic covariance structure : \code{VVV-VUN} (unstructured), \code{VVV-VAR} (AR), \code{VVV-VCS} (CS).
#'   \item Homoscadatsic covariance structure : \code{EEE-EUN} (unstructured), \code{EEE-EAR} (AR), \code{EEE-ECS} (CS).
#' }
#' For initialization of cluster membership, see `Details' in \code{\link{Mclust_SEP_each_cpp}}.
#' @seealso \code{\link{Mclust_SEP_each_cpp}}
#' @examples 
#' # Gaussian mixture model with two components
#' K <- 2
#' p <- 2
#' q <- 3
#' U <- lapply(1:K, function(noarg) getCovariance(p, 0.3, "AR"))
#' V <- lapply(1:K, function(noarg) getCovariance(q, 0.2, "CS"))
#' Sigma <- Map(kronecker, U, V) # separable covariance matrix
#' mu <- list(rep(0, p * q), 5 / sqrt(p*q) * rep(1, p * q)) # distinct mean vectors
#' Y <- vector(mode = "list", length = K)
#' for(i in 1:K){
#'   Y[[i]] <- mvtnorm::rmvnorm(n = 20, mean = mu[[i]], sigma = Sigma[[i]])
#' }
#' fit <- Mclust_SEP_cpp(Y = Reduce(rbind, Y), p = p, q = q, Ks = 2, type_cov = "EEE-ECS")
#' @export
Mclust_SEP_cpp <- function(Y, 
                           p, 
                           q,
                           Ks, 
                           type_cov,
                           tol = 1e-3, 
                           maxit = 500,
                           save_fit = TRUE){
  if(type_cov == "all"){
    type_cov <- expand.grid(type_vari = c("VVV"),
                            type_time = c("VUN", "VAR", "VCS")
    )
    type_cov <- rbind(type_cov, 
                      expand.grid(type_vari = c("EEE"),
                                  type_time = c("EUN", "EAR", "ECS")
                      )
    )
    type_cov <- paste0(type_cov$type_vari, "-", type_cov$type_time)
  }
  list_type_cov <- strsplit(type_cov, split = "-")
  type_vari <- unlist(lapply(list_type_cov, function(x) x[1]))
  type_time <- unlist(lapply(list_type_cov, function(x) x[2]))
  
  ## Check out available inputs
  if(any(!type_vari %in% c("VVV", "EEE"))){
    stop("Please check available options for type_time\n");
  }
  if(any(!type_time %in% c("VUN", "VAR", "VCS", "EUN", "EAR", "ECS"))){
    stop("Please check available options for type_time\n");
  }
  if(any(type_vari == "VVV" & !type_time %in% c("VUN", "VAR", "VCS"))){
    stop("type_vari=\"VVV\" can be matched with one of type_time=\"VUN\", \"VAR\", \"VCS\"\n");
  }
  if(any(type_vari == "EEE" & !type_time %in% c("EUN", "EAR", "ECS"))){
    stop("type_vari=\"EEE\" can be matched with one of type_time=\"EUN\", \"EAR\", \"ECS\"\n");
  }
  
  ## Run Mclust_SEP_cpp for each combination of separable factors
  res_Mclust_SEP <- list()
  bic_table <- matrix(0, nrow = length(type_vari), ncol = length(Ks))
  rownames(bic_table) <- type_cov
  colnames(bic_table) <- as.character(Ks)
  
  for(i in seq_along(type_vari)){
    res_Mclust_SEP[[i]] <- list()
    names(res_Mclust_SEP)[i] <- type_cov[i]
    for(k in seq_along(Ks)){
      # err_flag <- TRUE # until err_flag = FALSE, perform the model-based clustering
      # while(err_flag){
      res_Mclust_SEP[[i]][[k]] <- 
        # try(
        Mclust_SEP_each_cpp(Y = Y, 
                            p = p, 
                            q = q, 
                            K = Ks[k], 
                            type_vari = type_vari[i],
                            type_time = type_time[i],
                            tol = tol, 
                            maxit = maxit)
      # , silent = T)
      # if(inherits(res_Mclust_SEP[[i]][[k]], "try-error")){
      #   print(res_Mclust_SEP[[i]][[k]]) # show error
      #   cat(sprintf("Refit Mclust_SEP (K=%d,type_cov=%s...)\n", Ks[k], type_cov[i]))
      # } else{
      #   err_flag <- FALSE # model is fitted without error
      # }
      # }
      names(res_Mclust_SEP[[i]])[k] <- Ks[k]
      bic_table[i,k] <- res_Mclust_SEP[[i]][[k]]$BIC
    }
  }
  cat("Mclust_SEP_cpp is done!\n")
  
  ## Extract best models
  idx_best_model <- which(abs(bic_table - max(bic_table)) < .Machine$double.eps, arr.ind = T)
  best_model <- vector("list", length = nrow(idx_best_model))
  for(i in seq_along(best_model)){
    best_model[[i]] <- res_Mclust_SEP[[idx_best_model[i,1]]][[idx_best_model[i,2]]]
  }
  
  ## Object to return
  res <- list(best_model = best_model,
              bic_table = bic_table
  )
  if(save_fit){
    res$res_Mclust_SEP <- res_Mclust_SEP
  }
  return(res)
}