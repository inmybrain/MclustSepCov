##------- Source from Fun_MclustSepCov_v5.0.R: do not edit by hand



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
#' U <- lapply(1:K, function(noarg) getCovariance(p, 0.3, "AR")) # getCovariance in Rcpp
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
  # Y = data_mix_gaussian$Y
  # p = data_mix_gaussian$info$p
  # q = data_mix_gaussian$info$q
  # Ks = Ks
  # type_cov = "all"
  # save_fit = T
  # tol = 1e-3 
  # maxit = 500
  # 
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
      err_flag <- TRUE # until err_flag = FALSE, perform the model-based clustering
      while(err_flag){
        res_Mclust_SEP[[i]][[k]] <- 
          try(
            Mclust_SEP_each_cpp(Y = Y, 
                                p = p, 
                                q = q, 
                                K = Ks[k], 
                                type_vari = type_vari[i],
                                type_time = type_time[i],
                                tol = tol, 
                                maxit = maxit), silent = T)
        if(inherits(res_Mclust_SEP[[i]][[k]], "try-error")){
          print(res_Mclust_SEP[[i]][[k]]) # show error
          cat(sprintf("Refit Mclust_SEP (K=%d,type_cov=%s...)\n", Ks[k], type_cov[i]))
        } else{
          err_flag <- FALSE # model is fitted without error
        }
      }
      names(res_Mclust_SEP[[i]])[k] <- Ks[k]
      bic_table[i,k] <- res_Mclust_SEP[[i]][[k]]$BIC
    }
  }
  cat("Mclust_SEP_cpp is done!\n")
  
  ## Extract best models
  bic_table[bic_table %in% c(Inf, -Inf)] <- NA
  idx_best_model <- which(abs(bic_table - max(bic_table, na.rm = TRUE)) < .Machine$double.eps, arr.ind = T)
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
