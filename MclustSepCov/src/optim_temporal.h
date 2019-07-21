//------- Source from Fun_MclustSepCov_v1.0.cpp: do not edit by hand
#ifndef optim_temporal
#define optim_temporal

Rcpp::List obj_AR_cpp(double rho,
                      double a,
                      arma::mat Z
)
;
Rcpp::List obj_CS_cpp(double rho,
                      double a,
                      arma::mat Z
)
;
Rcpp::List obj_cpp(double rho,
                   double a,
                   arma::mat Z,
                   std::string type)
;
Rcpp::List log_barrier_cpp(double rho, 
                           double lambda, 
                           int p, 
                           std::string type)
;
Rcpp::List LB_algorithm_cpp(double a,
                            arma::mat Z,
                            std::string type,
                            double rho0,
                            double lambda = 1.0,
                            int maxit = 500)
;
#endif
