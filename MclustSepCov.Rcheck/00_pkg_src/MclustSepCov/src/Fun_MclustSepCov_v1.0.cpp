// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//--------------------------------------------------------------------------
// Miscellaneous functions
//--------------------------------------------------------------------------
// (190102) source from https://github.com/helske/seqHMM/blob/master/src/logSumExp.cpp
double log_sum_exp(const arma::vec& x) {
  unsigned int maxi = x.index_max();
  double maxv = x(maxi);
  if (!(maxv > -arma::datum::inf)) {
    return -arma::datum::inf;
  }
  double cumsum = 0.0;
  for (unsigned int i = 0; i < x.n_elem; i++) {
    if ((i != maxi) & (x(i) > -arma::datum::inf)) {
      cumsum += exp(x(i) - maxv);
    }
  }
  return maxv + log1p(cumsum);
}
// (190102) source from https://github.com/RcppCore/rcpp-gallery/blob/gh-pages/src/2013-07-13-dmvnorm_arma.Rmd
arma::vec dmvnrm_arma(arma::mat x,  
                      arma::rowvec mean,  
                      arma::mat sigma, 
                      bool logd = false) { 
  int n = x.n_rows;
  int xdim = x.n_cols;
  arma::vec out(n);
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * std::log(2.0 * M_PI);
  
  for (int i=0; i < n; i++) {
    arma::vec z = rooti * arma::trans( x.row(i) - mean) ;    
    out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;     
  }  
  
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

//--------------------------------------------------------------------------
// Data generation
//--------------------------------------------------------------------------
//' @name getCovariance
//' @title Generate temporal covariance matrices
//' @description Return a covariance matrix with temporal structure.
//' @param q Dimension of a covariance matrix.
//' @param rho Temporal correlation in \eqn{(-1, 1)}. For \code{type='CS'}, \eqn{-1/\sqrt{q-1}} is a lower bound of \code{rho} for the returned matrix to be positive definite.
//' @param type A character string indicating one of types of temporal structure; autoregressive model if \code{type='AR'} and compound symmetry model if \code{type='CS'}. See `Details' for their structures.
//' @return A \code{q}-by-\code{q} temporal covariance matrix with diagonals 1.
//' @details Following temporal structures are available. 
//' \enumerate{
//'   \item Autogressive structure: \eqn{V = \big(\rho^{|i-j|}; 1\le i,j \le q\big)}.
//'   \item Compound symmetry structure: \eqn{V = \big({\rm I}(i=j) + \rho {\rm I}(i \neq j); 1\le i,j \le q\big)}.
//' }
//' @examples
//' getCovariance(3, 0.3, "AR") # AR
//' getCovariance(3, 0.3, "CS") # CS
//' 
//' # AR structure with heterogeneous variances
//' hvar <- c(1, 2, 3) # variances
//' diag(sqrt(hvar)) %*% getCovariance(3, 0.3, "AR") %*% diag(sqrt(hvar))
//' @export
// [[Rcpp::export]]
arma::mat getCovariance(int q, 
                        double rho, 
                        std::string type){
  mat Sig(q, q);
  if(type == "AR"){
    for(int i = 0; i < q; i++){
      for(int j = 0; j < q; j++){
        Sig(i,j) = pow(rho, abs(i - j));
      }
    }
  } else if(type == "CS"){
    Sig.fill(rho);
    for(int i = 0; i < q; i++){
      Sig(i,i) = 1;
    }
  }
  return Sig;
}

//--------------------------------------------------------------------------
// Optimization for temporal correlation
//--------------------------------------------------------------------------
Rcpp::List obj_AR_cpp(double rho,
                      double a,
                      arma::mat Z
){
  int d = Z.n_cols;
  
  double ssZ1 = accu(pow(Z.col(0), 2.0));
  double ssZcross = accu(Z.cols(1, d - 1) % Z.cols(0, d - 2));
  double ssZ_1st = accu(Z.cols(0, d - 2) % Z.cols(0, d - 2));
  double ssZ_last = accu(Z.cols(1, d - 1) % Z.cols(1, d - 1));
  
  double val = a * (d - 1) * log(1.0 - pow(rho, 2.0)) + ssZ1 + 
    (ssZ_last - 2.0 * rho * ssZcross + pow(rho, 2.0) * ssZ_1st) / (1.0 - pow(rho, 2.0));
  
  double der = -(2.0 * a * (d - 1) * rho) / (1.0 - pow(rho, 2.0)) - 
    2.0 * (1.0 + pow(rho, 2.0)) / pow((1.0 - pow(rho, 2.0)), 2.0) * ssZcross +
    2.0 * rho / pow((1.0 - pow(rho, 2.0)), 2.0) * (ssZ_1st + ssZ_last);
  
  double hes = - 2.0 * a * (d - 1) * (1.0 + pow(rho, 2.0)) / pow(1.0 - pow(rho, 2.0), 2.0) - 
    4.0 * rho * (3.0 + pow(rho, 2.0)) / pow((1 - pow(rho, 2.0)), 3.0) * ssZcross + 
    2.0 * (1.0 + 3.0 * pow(rho, 2.0)) / pow((1 - pow(rho, 2.0)), 3.0) * (ssZ_1st + ssZ_last);
  
  return List::create(
    Named("objective") = val,
    Named("derivative") = der,
    Named("hessian") = hes
  );
}
Rcpp::List obj_CS_cpp(double rho,
                      double a,
                      arma::mat Z
){
  int d = Z.n_cols;
  
  double ssZ1 = accu(pow(Z * ones<vec>(d), 2.0));
  double ssZ = accu(Z % Z);
  
  double val = a * (d * log(1 - rho) + log(1 + rho * d / (1 - rho))) +
    (ssZ - rho * ssZ1 / (1 + (d - 1) * rho))/ (1 - rho);
  
  double der1 = a * (d - 1.0) * (1.0 / (1.0 + (d - 1.0) * rho) - 1.0 / (1.0 - rho));
  double der2 = ssZ / pow(1.0 - rho, 2.0);
  double der3 = -pow(1.0 - rho, -2.0) * ssZ1 * (1.0 + (d - 1.0) * pow(rho, 2.0)) / pow(1.0 + (d - 1.0) * rho, 2.0);
  
  double hes1 = - a * (d - 1.0) * ( (d - 1.0) / pow(1.0 + (d - 1.0) * rho, 2.0) +
                       1.0 / pow(1.0 - rho, 2.0)
  );
  double hes2 = 2.0 * ssZ / pow(1.0 - rho, 3.0);
  double hes31 = (d - 1.0) * rho * (1.0 - rho) * (1.0 + (d - 1.0) * rho);
  double hes32 = (1.0 + (d - 1.0) * pow(rho, 2.0)) * (1.0 + (d - 1.0) * rho);
  double hes33 = -(d - 1.0) * (1.0 + (d - 1.0) * pow(rho, 2.0)) * (1.0 - rho);
  double hes3 = -2.0 * ssZ1 / (pow(1.0 - rho, 3.0) * pow(1.0 + (d - 1.0) * rho, 3.0) ) * (hes31 + hes32 + hes33);
  return List::create(
    Named("objective") = val,
    Named("derivative") = der1 + der2 + der3,
    Named("hessian") = hes1 + hes2 + hes3
  );
}
Rcpp::List obj_cpp(double rho,
                   double a,
                   arma::mat Z,
                   std::string type){
  if(type == "AR"){
    return obj_AR_cpp(rho, a, Z);
  } else if (type == "CS"){
    return obj_CS_cpp(rho, a, Z);
  } else {
    stop("Please check available options for type.\n");
  }
}
Rcpp::List log_barrier_cpp(double rho, 
                           double lambda, 
                           int p, 
                           std::string type){
  if (type == "AR"){
    return List::create(Named("objective") = lambda * (log(rho + 1.0) + log(1.0 - rho)),
                        Named("derivative") = lambda * (1.0 / (rho + 1.0) - 1.0 / (1.0 - rho)),
                        Named("hessian") = lambda * (1.0 / pow(rho + 1.0, 2.0) - 1.0/pow(1.0 - rho, 2.0))
    );
  } else if (type == "CS"){
    return List::create(Named("objective") = lambda * (log(rho + 1.0 / (p - 1.0)) + log(1.0 - rho)),
                        Named("derivative") = lambda * (1.0 / (rho + 1.0 / (p - 1.0)) - 1.0 / (1.0 - rho)),
                        Named("hessian") = lambda * (1.0 / pow((rho + 1.0 / (p - 1.0)), 2.0) - 1.0 / pow(1.0 - rho, 2.0))
    );
  } else {
    stop("Please check available options for type.\n");
  }
}

//' @name Optimization
//' @title Newton-Raphson's algorithm to find the optimal temporal correlation
//' @description Solve the constrained minimization problem using the log-barrier method to find the maximum likelihood estimator (MLE) of temporal correlation. The objective function is described in `Details'.
//' @param a A positive constant. See `Details'.
//' @param Z A matrix of sample vectors at row.
//' @param type A character string indicating a type of temporal covariance matrix. Available options are \code{'AR'} and \code{'CS'}.
//' @param rho0 An initial value for the temporal correlation coefficient. We empirically found that 0.001 works well for \code{type='AR'} and 0.5 for \code{type='CS'}.
//' @param lambda A positive constant multiplied to the log-barrier term. Default is 1.
//' @param maxit The maximum number for iterations. Default is 500.
//' @details The objective function is divided into two parts; the Gaussian log-likelihood function (up to constant multiplication) with mean 0 and covariance matrix \eqn{\Sigma = \Sigma(\rho)} and the log-barrier function. The former is written by 
//' \deqn{h(\rho; a, Z) = a \log|\Sigma| + {\rm tr}(\Sigma^{-1} S),}
//' where \eqn{a>0} and \eqn{S = Z^{\rm T} Z}, and the latter is
//' \deqn{b(\rho; u, l) = \log(u - \rho) + \log(\rho - l),}
//' where \eqn{u,l} is an upper and a lower bound of \eqn{\rho}, respectively. These quantities depend on \code{type} as follows;
//' \itemize{
//'   \item If \code{type='AR'}, \eqn{\Sigma= \big(\rho^{|i-j|}; 1\le i,j \le q\big)} and  \eqn{l=-1, u=1},
//'   \item if \code{type='CS'}, \eqn{\Sigma = \big({\rm I}(i=j) + \rho {\rm I}(i \neq j); 1\le i,j \le q\big)} and \eqn{l=-1/\sqrt{q-1}, u=1},
//' }
//' where \eqn{q=}\code{ncol(Z)}. The objective function is, hence,
//' \deqn{h(\rho; a, Z) - \lambda ~ b(\rho; u, l).}
//' @examples
//' q <- 10
//' # AR model
//' set.seed(6)
//' Y <- mvtnorm::rmvnorm(100, rep(0, q), getCovariance(q, 0.3, "AR"))
//' LB_algorithm_cpp(a = nrow(Y), Z = Y, rho0 = 1e-3, type = "AR")
//' 
//' # CS model
//' set.seed(6)
//' Y <- mvtnorm::rmvnorm(100, rep(0, q), getCovariance(q, 0.3, "CS"))
//' LB_algorithm_cpp(a = nrow(Y), Z = Y, rho0 = 1e-3, type = "CS")
//' @exports
// [[Rcpp::export]]
Rcpp::List LB_algorithm_cpp(double a,
                            arma::mat Z,
                            std::string type,
                            double rho0,
                            double lambda = 1.0,
                            int maxit = 500){
  double lb_rho, ub_rho;
  if(type == "AR"){
    lb_rho = -1.0;
    ub_rho =  1.0;
  } else if(type == "CS"){
    lb_rho = -1.0 / (Z.n_cols - 1.0);
    ub_rho =  1.0;
  } else{
    stop("No other types are suppored; AR and CS only.\n");
  }
  
  int iter = 0;
  int conv = 0;
  double rho1, val0, val1, der0, hes0;
  // double der1, hes1;
  
  List obj = obj_cpp(rho0, a, Z, type);
  List barr = log_barrier_cpp(rho0, lambda, Z.n_cols, type);
  val0 = as<double>(obj["objective"]) - as<double>(barr["objective"]);
  der0 = as<double>(obj["derivative"]) - as<double>(barr["derivative"]);
  hes0 = as<double>(obj["hessian"]) - as<double>(barr["hessian"]);
  
  while(iter < maxit){
    rho1 = rho0 - der0 / hes0;
    
    // If an estimate is out of boundary, then force it to be around the boundary.
    if(lb_rho > rho1){
      rho1 = lb_rho + 1e-7;
    }
    if(ub_rho < rho1){
      rho1 = ub_rho - 1e-7;
    }
    lambda *= R::runif(0, 1);
    
    obj = obj_cpp(rho1, a, Z, type);
    barr = log_barrier_cpp(rho1, lambda, Z.n_cols, type);
    val1 = as<double>(obj["objective"]) - as<double>(barr["objective"]);
    if(abs(val1 - val0 / val0) < 1e-3){
      conv = 1;
      break;
    } else{
      val0 = val1;
      der0 = as<double>(obj["derivative"]) - as<double>(barr["derivative"]);
      hes0 = as<double>(obj["hessian"]) - as<double>(barr["hessian"]);
      rho0 = rho1;
      iter += 1;
    }
  }
  val1 = as<double>(obj["objective"]) - as<double>(barr["objective"]);
  // der1 = as<double>(obj["derivative"]) - as<double>(barr["derivative"]);
  // hes1 = as<double>(obj["hessian"]) - as<double>(barr["hessian"]);
  
  return List::create(Named("value") = val1,
                      Named("rho") = rho1,
                      Named("iter") = iter,
                      Named("conv") = conv
  );
}

//--------------------------------------------------------------------------
// EM algorithm
//--------------------------------------------------------------------------
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

//--------------------------------------------------------------------------
// Model-based clustering
//--------------------------------------------------------------------------
// mat getSampleCov(mat X){
//   int n = X.n_rows;
//   
//   return X.t() * (diagmat(ones(n)) - ones(n, n) / (double) n) * X / (n - 1.0);
// }
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
    double tol = 1e-3,
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
  // U0.each_slice() = getSampleCov(Y.cols(regspace<uvec>(0,  q,  Y.n_cols - 1)));
  // V0.each_slice() = getSampleCov(Y.cols(0, q - 1));
  U0.each_slice() = eye(p, p);
  V0.each_slice() = eye(q, q);
  
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
                        maxit
    );
    mu = as<mat>(ms["M"]);
    U = as<cube>(ms["U"]);
    V = as<cube>(ms["V"]);
    wt_cluster = as<vec>(ms["wt_cluster"]);
    
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
