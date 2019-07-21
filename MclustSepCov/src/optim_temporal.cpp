//------- Source from Fun_MclustSepCov_v1.0.cpp: do not edit by hand
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


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

