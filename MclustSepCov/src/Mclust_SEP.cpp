//------- Source from Fun_MclustSepCov_v1.0.cpp: do not edit by hand
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

#include "misc.h"
#include "em_step.h"
#include "optim_temporal.h"
//' @name Mclust_SEP_each_cpp
//' @title The maximum likelihood estimation of the mixture distribution
//' @description Perform the EM algorithm for fitting the finite Gaussian mixture distribution with covariance separability.
//' @param Y Same as in \code{\link{Mclust_SEP_cpp}}.
//' @param p,q Same as in \code{\link{Mclust_SEP_cpp}}.
//' @param K A positive integer indicating the number of mixture components.
//' @param type_vari,type_time A character string indicating a structure of covariance factors, passed from \code{\link{Mclust_SEP_cpp}}. See `Details' for available options.
//' @param tol Same as in \code{\link{Mclust_SEP_cpp}}.
//' @param maxit Same as in \code{\link{Mclust_SEP_cpp}}.
//' @return A list with components:
//' \item{loglik}{The log-likelihood function.}
//' \item{df}{The degrees of freedom or the number of parameters of the mixture model.}
//' \item{BIC}{The Bayesian information crieteria.}
//' \item{K}{\code{K} from the input.}
//' \item{id_cluster}{Cluster membership of samples.}
//' \item{wt_cluster}{A matrix of dimension r-by-\code{K} whose row represents the maximum a posteriori of a sample.}
//' \item{EM_iter}{The number of iterations in the EM algorithm.}
//' \item{mu}{A matrix of the estimated mean whose row is the mean vector of a mixture component.}
//' \item{U,V}{A cube containing \code{K} slices of the estimated covariance matrix.}
//' \item{type_vari, type_time}{\code{type_vari}, \code{type_time} from the input.}
//' @details 
//' Cluster membership from 1 to \code{K} is randomly assigned to each sample.\cr\cr
//' \code{type_vari} specifies a type of the multivariable covariance \eqn{U_{p \times p}};
//' \itemize{
//'   \item Heteroscadatsic : \code{VVV} (unstructured)
//'   \item Homoscadatsic : \code{EEE} (unstructured)
//' }
//' and \code{type_time} the temporal covariance \eqn{V_{q\times q}};
//' \itemize{
//'   \item Heteroscadatsic :  \code{VUN} (unstructured), \code{VAR} (AR), \code{VCS} (CS) 
//'   \item Homoscadatsic :  \code{EUN} (unstructured), \code{EAR} (AR), \code{ECS} (CS)
//' }
//' @seealso \code{\link{Mclust_SEP_cpp}}
//' @examples
//' # Gaussian mixture model with two components
//' K <- 2
//' p <- 2
//' q <- 3
//' U <- lapply(1:K, function(noarg) getCovariance(p, 0.3, "AR"))
//' V <- lapply(1:K, function(noarg) getCovariance(q, 0.2, "CS"))
//' Sigma <- Map(kronecker, U, V) # separable covariance matrix
//' mu <- list(rep(0, p * q), 5 / sqrt(p*q) * rep(1, p * q)) # distinct mean vectors
//' Y <- vector(mode = "list", length = K)
//' for(i in 1:K){
//'   Y[[i]] <- mvtnorm::rmvnorm(n = 20, mean = mu[[i]], sigma = Sigma[[i]])
//' }
//' fit <- Mclust_SEP_each_cpp(Y = Reduce(rbind, Y), p = p, q = q, K = 2, type_vari = "EEE", type_time = "ECS")
//' @export
// [[Rcpp::export]]
Rcpp::List Mclust_SEP_each_cpp(
    arma::mat Y,
    int p,
    int q,
    int K,
    std::string type_vari,
    std::string type_time,
    double tol = 1e-4,
    int maxit = 500){
  // Input :
  // - Y : r-by-(pq) matrix
  // Initial
  //// - id_clust : division of elements into Ks clusters
  //// - wt_clust  : mixture weights or probabilities
  //// - mu, U, V : parameters used in mixture Gaussian models
  
  int r = Y.n_rows;
  int iter = 0;
  
  IntegerVector id_cluster0 = Rcpp::sample(K, r, true);
  vec wt_cluster0 = zeros(K);
  for(int i = 0; i < r; i++){
    wt_cluster0(id_cluster0(i) - 1) += 1;
  }
  wt_cluster0 /= (double) r;
  
  cube U0(p, p, K), V0(q, q, K);
  U0.each_slice() = getSampleCov(Y.cols(regspace<uvec>(0,  q,  Y.n_cols - 1)));
  V0.each_slice() = getSampleCov(Y.cols(0, q - 1));
  // U0.each_slice() = eye(p, p);
  // V0.each_slice() = eye(q, q);
  
  // mu0 : K-by-pq matrix, each row corresponds to a mean vector of length pq in each group
  mat mu0 = zeros(K, p * q);
  vec n_cluster = zeros(K);
  
  for (int i = 0; i < r; i++){
    mu0.row(id_cluster0(i) - 1) += Y.row(i);
    n_cluster(id_cluster0(i) - 1) += 1;
  }
  for (int k = 0; k < K; k++){
    mu0.row(k) /= (double) n_cluster(k);
  }
  
  List es = Estep_cpp(Y,
                      mu0,
                      U0,
                      V0,
                      wt_cluster0
  );
  mat W0 = es["W"]; // r-by-K
  double L0 = es["likelihood"];
  
  vec wt_cluster;
  mat mu;
  cube U, V;
  // EM algorithm
  while(iter < maxit){
    iter += 1;
    
    // M-step
    List ms = Mstep_cpp(Y,
                        W0,
                        p,
                        q,
                        type_vari,
                        type_time,
                        tol,
                        200
    );
    mu = as<mat>(ms["M"]);
    U = as<cube>(ms["U"]);
    V = as<cube>(ms["V"]);
    wt_cluster = as<vec>(ms["wt_cluster"]);
    // Rcout <<  iter << "th, wt_cluster0=\n" <<  wt_cluster0 << "\n";
    if(any(wt_cluster < 1.0 / r)){
      Rcout << "On the boundary of parameter space... adjusting it...\n";
      iter -= 1;
      W0 += 1e-3;
      continue; 
    }
    // E-step
    es = Estep_cpp(Y,
                   mu,
                   U,
                   V,
                   wt_cluster
    );
    
    if(any(iter == regspace(0, 100, maxit))){
      Rcout <<  iter << "th, L=" <<  as<double>(es["likelihood"]) << "\n";
    }
    // stopping criteria
    if((abs(as<double>(es["likelihood"]) - L0) < tol) & (iter > maxit * 0.01)){
      break;
    }
    W0 = as<mat>(es["W"]);
    mu0 = mu;
    U0 = U;
    V0 = V;
    wt_cluster0 = wt_cluster;
    L0 = as<double>(es["likelihood"]);
  } // end of the EM algorithm
  
  uvec id_cluster(r);
  for(int ii = 0; ii < r; ii++){
    id_cluster(ii) = as<mat>(es["W"]).row(ii).index_max() + 1;
  }
  
  // BIC
  double LL = accu(log(sum(exp(as<mat>(es["logden"])), 1))); // log-likelihood of Y_1,obs, ..., Y_n,obs given final parameter
  
  double df = 0;
  if(type_vari == "VVV"){
    df = df + K * p * (p + 1.0) / 2.0;
    if(type_time == "VUN"){// V = "VUN"
      df = df + K * (q * (q + 1) / 2.0 - 1);
    } else if(type_time == "VAR"){// V = "VAR"
      df = df + K * 1.0;
    } else if(type_time == "VCS"){// V = "VCS"
      df = df + K * 1.0;
    } else{
      stop("Please check available options for type_time\n");
    }
  } else if(type_vari == "EEE"){
    df = df + p * (p + 1) / 2.0;
    if(type_time == "EUN"){// V = "EUN"
      df = df + q * (q + 1) / 2.0 - 1;
    } else if(type_time == "EAR"){// V = "EAR"
      df = df + 1.0;
    } else if(type_time == "ECS"){// V = "ECS"
      df = df + 1.0;
    } else{
      stop("Please check available options for type_time\n");
    }
  } else{
    stop("Please check available options for type_vari.\n");
  }
  df = df + K - 1.0 + K * p * q ; // parameters used to model mean and mixing probability
  double BIC = 2.0 * LL - log(r) * df;
  
  return List::create(
    Named("loglik") = LL,
    Named("df") = df,
    Named("BIC") = BIC,
    Named("K") = K,
    Named("id_cluster") = id_cluster,
    Named("wt_cluster") = wt_cluster,
    Named("EM_iter") = iter,
    Named("mu") = mu,
    Named("U") = U,
    Named("V") = V,
    Named("type_vari") = type_vari,
    Named("type_time") = type_time
  );
}
