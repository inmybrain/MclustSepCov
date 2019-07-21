//------- Source from Fun_MclustSepCov_v1.0.cpp: do not edit by hand
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


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
  return out;
}

arma::mat getSampleCov(arma::mat X){
  int n = X.n_rows;

  return X.t() * (diagmat(ones(n)) - ones(n, n) / (double) n) * X / (n - 1.0);
}

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


