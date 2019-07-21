//------- Source from Fun_MclustSepCov_v1.0.cpp: do not edit by hand
#ifndef misc
#define misc

double log_sum_exp(const arma::vec& x) 
;
arma::vec dmvnrm_arma(arma::mat x,
                      arma::rowvec mean,
                      arma::mat sigma,
                      bool logd = false) 
;
arma::mat getSampleCov(arma::mat X)
;
arma::mat getCovariance(int q,
                        double rho,
                        std::string type)
;
#endif
