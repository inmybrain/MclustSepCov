//------- Source from Fun_MclustSepCov_v1.0.cpp: do not edit by hand
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

#include "misc.h"
#include "optim_temporal.h"
Rcpp::List Estep_cpp(arma::mat Y, 
                     arma::mat mu, 
                     arma::cube U, 
                     arma::cube V, 
                     arma::vec wt_cluster){
  // U : p x p x K, K slices have a p-by-p covariance factor, respectively.
  // V : q x q x K, K slices have a q-by-q covariance factor, respectively.
  
  int r = Y.n_rows, K = wt_cluster.size();
  // density in log scale
  mat logden(r, K);
  mat W(r, K);
  
  for (int k = 0; k < K; k++){
    logden.col(k) = dmvnrm_arma(Y, mu.row(k), kron(U.slice(k), V.slice(k)), true) + log(wt_cluster(k));
  }
  for(int i = 0; i < r; i++){
    W.row(i) = exp(logden.row(i) - log_sum_exp(conv_to<vec>::from(logden.row(i))));
  }
  double lk = accu(logden % W);
  
  return List::create(Named("likelihood") = lk,
                      Named("W") = W,
                      Named("logden") = logden
  );
}

Rcpp::List Mstep_cpp(arma::mat Y,
                     arma::mat W,
                     int p,
                     int q,
                     std::string type_vari,
                     std::string type_time,
                     double eps,
                     int maxit){
  // Y : r-by-pq matrix
  // W : r-by-K matrix
  int r = Y.n_rows;
  int K = W.n_cols;
  
  // Mean
  cube M_cube(p, q, K, fill::zeros);
  mat M_mat(K, p * q);
  for(int kk = 0; kk < K; kk++){
    for(int ii = 0; ii < r; ii++){
      M_cube.slice(kk) += reshape(Y.row(ii) * W(ii, kk) / accu(W.col(kk)), q, p).t();
    }
  }
  for(int kk = 0; kk < K; kk++){
    M_mat.row(kk) = vectorise(M_cube.slice(kk), 1);
  }
  
  // Group probability
  rowvec wt_cluster = mean(W, 0);
  
  // Covariance factors, U (covariance for variables) and V (covariance for time)
  cube U(p, p, K), U_new(p, p, K), V(q, q, K), V_new(q, q, K);
  
  // Initialization
  for(int kk = 0; kk < K; kk++){
    U.slice(kk) = eye(p, p);
    V.slice(kk) = eye(q, q);
  }
  U_new = U; 
  V_new = V;
  
  if(type_vari == "VVV"){ // U = "VVV"
    for(int kk = 0; kk < K; kk++){
      mat Z = zeros(p * p, q * q);
      for(int ii = 0; ii < r; ii++){
        mat Yr = reshape(Y.row(ii), q, p).t();
        Z += W(ii, kk) * kron(Yr - M_cube.slice(kk), Yr - M_cube.slice(kk));
      }
      if(type_time == "VUN"){// V = "VUN"
        int iter = 0;
        while(iter < maxit){
          iter += 1;
          U_new.slice(kk) = reshape(Z * vectorise(inv(V.slice(kk))) / (q * accu(W.col(kk))), p, p);
          V_new.slice(kk) = reshape(Z.t() * vectorise(inv(U_new.slice(kk))) / (p * accu(W.col(kk))), q, q);
          
          if((accu(square(U_new.slice(kk) - U.slice(kk))) < pow(eps, 2.0)) * (accu(square(V_new.slice(kk) - V.slice(kk))) < pow(eps, 2.0))){
            U.slice(kk) = U_new.slice(kk);
            V.slice(kk) = V_new.slice(kk);
            break;
          }
          U.slice(kk) = U_new.slice(kk);
          V.slice(kk) = V_new.slice(kk);
        }
      } else if(type_time == "VAR"){// V = "VAR"
        int iter = 0;
        mat Y_tilde = zeros(r * p, q);
        List res_LB;
        while(iter < maxit){
          iter += 1;
          U_new.slice(kk) = reshape(Z * vectorise(inv(V.slice(kk))) / (q * accu(W.col(kk))), p, p);
          mat U_new_inv_sqrt = sqrtmat_sympd(inv(U_new.slice(kk)));
          int idx = 0;
          for(int ii = 0; ii < r; ii++){
            Y_tilde.rows(idx * p, (idx + 1) * p - 1) = U_new_inv_sqrt * pow(W(ii, kk), 0.5) * (reshape(Y.row(ii), q, p).t() - M_cube.slice(kk));
            idx += 1;
          }
          res_LB = LB_algorithm_cpp(accu(W.col(kk)) * p, Y_tilde, "AR", 0.001, 1, maxit);
          V_new.slice(kk) = getCovariance(q, as<double>(res_LB["rho"]), "AR");
          
          
          if((accu(square(U_new.slice(kk) - U.slice(kk))) < pow(eps, 2.0)) & (accu(square(V_new.slice(kk) - V.slice(kk))) < pow(eps, 2.0))){
            U.slice(kk) = U_new.slice(kk);
            V.slice(kk) = V_new.slice(kk);
            break;
          }
          U.slice(kk) = U_new.slice(kk);
          V.slice(kk) = V_new.slice(kk);
        }
      } else if(type_time == "VCS"){// V = "VCS"
        int iter = 0;
        mat Y_tilde = zeros(r * p, q);
        List res_LB;
        while(iter < maxit){
          iter += 1;
          U_new.slice(kk) = reshape(Z * vectorise(inv(V.slice(kk))) / (q * accu(W.col(kk))), p, p);
          mat U_new_inv_sqrt = sqrtmat_sympd(inv(U_new.slice(kk)));
          int idx = 0;
          for(int ii = 0; ii < r; ii++){
            Y_tilde.rows(idx * p, (idx + 1) * p - 1) = U_new_inv_sqrt * pow(W(ii, kk), 0.5) * (reshape(Y.row(ii), q, p).t() - M_cube.slice(kk));
            idx += 1;
          }
          res_LB = LB_algorithm_cpp(accu(W.col(kk)) * p, Y_tilde, "CS", 0.5, 1, maxit);
          V_new.slice(kk) = getCovariance(q, as<double>(res_LB["rho"]), "CS");
          
          if((accu(square(U_new.slice(kk) - U.slice(kk))) < pow(eps, 2.0)) * (accu(square(V_new.slice(kk) - V.slice(kk))) < pow(eps, 2.0))){
            U.slice(kk) = U_new.slice(kk);
            V.slice(kk) = V_new.slice(kk);
            break;
          }
          U.slice(kk) = U_new.slice(kk);
          V.slice(kk) = V_new.slice(kk);
        }
      } else{
        stop("Please check available options for type_time.\n");
      } // type_time
    }
  } else if(type_vari == "EEE"){ // U = "EEE"
    mat Z = zeros(p * p, q * q);
    for(int ii = 0; ii < r; ii++){
      for(int kk = 0; kk < K; kk++){
        mat Yr = reshape(Y.row(ii), q, p).t();
        Z += W(ii, kk) * kron(Yr - M_cube.slice(kk), Yr - M_cube.slice(kk));
      }
    }
    mat Y_tilde = zeros<mat>(K * r * p, q);
    if(type_time == "EUN"){// V = "EUN"
      int iter = 0;
      while(iter < maxit){
        iter += 1;
        U_new.slice(0) = reshape(Z * vectorise(inv(V.slice(0))) / (q * r), p, p);
        V_new.slice(0) = reshape(Z.t() * vectorise(inv(U_new.slice(0))) / (p * r), q, q);
        
        if((accu(square(U_new.slice(0) - U.slice(0))) < pow(eps, 2.0)) * (accu(square(V_new.slice(0) - V.slice(0))) < pow(eps, 2.0))){
          U.slice(0) = U_new.slice(0);
          V.slice(0) = V_new.slice(0);
          break;
        }
        U.slice(0) = U_new.slice(0);
        V.slice(0) = V_new.slice(0);
      }
    } else if(type_time == "EAR"){// V = "EAR"
      int iter = 0;
      List res_LB;
      while(iter < maxit){
        iter += 1;
        U_new.slice(0) = reshape(Z * vectorise(inv(V.slice(0))) / (q * r), p, p);
        mat U_new_inv_sqrt = sqrtmat_sympd(inv(U_new.slice(0)));
        int idx = 0;
        for(int ii = 0; ii < r; ii++){
          for(int kk = 0; kk < K; kk++){
            Y_tilde.rows(idx * p, (idx + 1) * p - 1) = U_new_inv_sqrt * pow(W(ii, kk), 0.5) * (reshape(Y.row(ii), q, p).t() - M_cube.slice(kk));
            idx += 1;
          }
        }
        res_LB = LB_algorithm_cpp(r * p, Y_tilde, "AR", 0.001, 1, maxit);
        V_new.slice(0) = getCovariance(q, as<double>(res_LB["rho"]), "AR");
        
        if((accu(square(U_new.slice(0) - U.slice(0))) < pow(eps, 2.0)) * (accu(square(V_new.slice(0) - V.slice(0))) < pow(eps, 2.0))){
          U.slice(0) = U_new.slice(0);
          V.slice(0) = V_new.slice(0);
          break;
        }
        U.slice(0) = U_new.slice(0);
        V.slice(0) = V_new.slice(0);
      }
    } else if(type_time == "ECS"){// V = "ECS"
      int iter = 0;
      List res_LB;
      while(iter < maxit){
        iter += 1;
        U_new.slice(0) = reshape(Z * vectorise(inv(V.slice(0))) / (q * r), p, p);
        mat U_new_inv_sqrt = sqrtmat_sympd(inv(U_new.slice(0)));
        int idx = 0;
        for(int ii = 0; ii < r; ii++){
          for(int kk = 0; kk < K; kk++){
            Y_tilde.rows(idx * p, (idx + 1) * p - 1) = U_new_inv_sqrt * pow(W(ii, kk), 0.5) * (reshape(Y.row(ii), q, p).t() - M_cube.slice(kk));
            idx += 1;
          }
        }
        res_LB = LB_algorithm_cpp(r * p, Y_tilde, "CS", 0.5, 1, maxit);
        V_new.slice(0) = getCovariance(q, as<double>(res_LB["rho"]), "CS");
        
        if((accu(square(U_new.slice(0) - U.slice(0))) < pow(eps, 2.0)) * (accu(square(V_new.slice(0) - V.slice(0))) < pow(eps, 2.0))){
          U.slice(0) = U_new.slice(0);
          V.slice(0) = V_new.slice(0);
          break;
        }
        U.slice(0) = U_new.slice(0);
        V.slice(0) = V_new.slice(0);
      }
    }else{
      stop("Please check available options for type_time.\n");
    } // type_time
    for(int kk = 0; kk < K; kk++){
      U.slice(kk) = U.slice(0);
      V.slice(kk) = V.slice(0);
    }
  } else{
    stop("Please check available options for type_vari.\n");
  } // type_vari
  
  return List::create(
    Named("wt_cluster") = wt_cluster,
    Named("M") = M_mat,
    Named("U") = U,
    Named("V") = V
  );
}

