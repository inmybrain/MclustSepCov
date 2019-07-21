//------- Source from Fun_MclustSepCov_v1.0.cpp: do not edit by hand
#ifndef em_step
#define em_step

Rcpp::List Estep_cpp(arma::mat Y, 
                     arma::mat mu, 
                     arma::cube U, 
                     arma::cube V, 
                     arma::vec wt_cluster)
;
Rcpp::List Mstep_cpp(arma::mat Y,
                     arma::mat W,
                     int p,
                     int q,
                     std::string type_vari,
                     std::string type_time,
                     double eps,
                     int maxit)
;
#endif
