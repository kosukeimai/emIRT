#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "estimate_poisIRT.h"

RcppExport SEXP poisIRT_estimate(SEXP alpha_startSEXP, 
                                 SEXP psi_startSEXP, 
                                 SEXP beta_startSEXP,
                                 SEXP x_startSEXP,
                                 SEXP ySEXP,
                                 SEXP iSEXP,
                                 SEXP niSEXP,
                                 SEXP psi_muSEXP, 
                                 SEXP psi_sigmaSEXP,
                                 SEXP alpha_muSEXP, 
                                 SEXP alpha_sigmaSEXP, 
                                 SEXP beta_muSEXP, 
                                 SEXP beta_sigmaSEXP, 
                                 SEXP x_muSEXP, 
                                 SEXP x_sigmaSEXP, 
                                 SEXP threadsSEXP, 
                                 SEXP verboseSEXP,
                                 SEXP maxitSEXP,
                                 SEXP threshSEXP,
                                 SEXP checkfreqSEXP
                                 ) {
  BEGIN_RCPP
    SEXP resultSEXP ;
  {
    Rcpp::RNGScope __rngScope ;
    Rcpp::traits::input_parameter<arma::mat>::type alpha_start(alpha_startSEXP) ;
    Rcpp::traits::input_parameter<arma::mat>::type psi_start(psi_startSEXP) ;
    Rcpp::traits::input_parameter<arma::mat>::type beta_start(beta_startSEXP) ;
    Rcpp::traits::input_parameter<arma::mat>::type x_start(x_startSEXP) ;
    Rcpp::traits::input_parameter<arma::mat>::type y(ySEXP) ;
    Rcpp::traits::input_parameter<arma::mat>::type i(iSEXP) ;
    Rcpp::traits::input_parameter<int>::type ni(niSEXP) ;
    Rcpp::traits::input_parameter<double>::type psi_mu(psi_muSEXP) ;
    Rcpp::traits::input_parameter<double>::type psi_sigma(psi_sigmaSEXP) ;
    Rcpp::traits::input_parameter<double>::type alpha_mu(alpha_muSEXP) ;
    Rcpp::traits::input_parameter<double>::type alpha_sigma(alpha_sigmaSEXP) ;
    Rcpp::traits::input_parameter<double>::type beta_mu(alpha_muSEXP) ;
    Rcpp::traits::input_parameter<double>::type beta_sigma(beta_sigmaSEXP) ;
    Rcpp::traits::input_parameter<double>::type x_mu(x_muSEXP) ;
    Rcpp::traits::input_parameter<double>::type x_sigma(x_sigmaSEXP) ;
    Rcpp::traits::input_parameter<int>::type threads(threadsSEXP) ;
    Rcpp::traits::input_parameter<bool>::type verbose(verboseSEXP) ;
    Rcpp::traits::input_parameter<int>::type maxit(maxitSEXP) ;
    Rcpp::traits::input_parameter<double>::type thresh(threshSEXP) ;
    Rcpp::traits::input_parameter<int>::type checkfreq(checkfreqSEXP) ;
    
    Rcpp::List result = estimate_poisIRT(alpha_start,
                                 psi_start,
                                 beta_start,
                                 x_start,
                                 y, 
                                 i, 
                                 ni, 
                                 psi_mu,
                                 psi_sigma,
                                 alpha_mu,
                                 alpha_sigma, 
                                 beta_mu,
                                 beta_sigma, 
                                 x_mu,
                                 x_sigma, 
                                 threads,
                                 verbose,
                                 maxit,
                                 thresh,
                                 checkfreq
                                 ) ;
    PROTECT(resultSEXP = Rcpp::wrap(result)) ;
  }
  UNPROTECT(1);
  return(resultSEXP) ;
  END_RCPP
    }
