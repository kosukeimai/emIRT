// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "estimate_hierIRT.h"

RcppExport SEXP hierIRT_estimate(SEXP alpha_startSEXP, 
                                 SEXP beta_startSEXP, 
                                 SEXP gamma_startSEXP, 
                                 SEXP sigma_startSEXP, 
                                 SEXP eta_startSEXP, 
                                 SEXP ySEXP,
                                 SEXP zSEXP,
                                 SEXP gSEXP,
                                 SEXP iSEXP,
                                 SEXP jSEXP,
                                 SEXP gammamuSEXP, 
                                 SEXP gammasigmaSEXP,
                                 SEXP betamuSEXP, 
                                 SEXP betasigmaSEXP, 
                                 SEXP sigmavSEXP, 
                                 SEXP sigmasSEXP, 
                                 SEXP NDSEXP, 
                                 SEXP NGSEXP, 
                                 SEXP NISEXP, 
                                 SEXP NJSEXP, 
                                 SEXP NLSEXP, 
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
    Rcpp::traits::input_parameter<arma::mat>::type beta_start(beta_startSEXP) ;
    Rcpp::traits::input_parameter<arma::mat>::type gamma_start(gamma_startSEXP) ;
    Rcpp::traits::input_parameter<arma::mat>::type sigma_start(sigma_startSEXP) ;
    Rcpp::traits::input_parameter<arma::mat>::type eta_start(eta_startSEXP) ;
    Rcpp::traits::input_parameter<arma::mat>::type y(ySEXP) ;
    Rcpp::traits::input_parameter<arma::mat>::type z(zSEXP) ;
    Rcpp::traits::input_parameter<arma::mat>::type g(gSEXP) ;
    Rcpp::traits::input_parameter<arma::mat>::type i(iSEXP) ;
    Rcpp::traits::input_parameter<arma::mat>::type j(jSEXP) ;
    Rcpp::traits::input_parameter<arma::mat>::type gammamu(gammamuSEXP) ;
    Rcpp::traits::input_parameter<arma::mat>::type gammasigma(gammasigmaSEXP) ;
    Rcpp::traits::input_parameter<arma::mat>::type betamu(betamuSEXP) ;
    Rcpp::traits::input_parameter<arma::mat>::type betasigma(betasigmaSEXP) ;
    Rcpp::traits::input_parameter<arma::mat>::type sigmav(sigmavSEXP) ;
    Rcpp::traits::input_parameter<arma::mat>::type sigmas(sigmasSEXP) ;
    Rcpp::traits::input_parameter<int>::type ND(NDSEXP) ;
    Rcpp::traits::input_parameter<int>::type NG(NGSEXP) ;
    Rcpp::traits::input_parameter<int>::type NI(NISEXP) ;
    Rcpp::traits::input_parameter<int>::type NJ(NJSEXP) ;
    Rcpp::traits::input_parameter<int>::type NL(NLSEXP) ;
    Rcpp::traits::input_parameter<int>::type threads(threadsSEXP) ;
    Rcpp::traits::input_parameter<bool>::type verbose(verboseSEXP) ;
    Rcpp::traits::input_parameter<int>::type maxit(maxitSEXP) ;
    Rcpp::traits::input_parameter<double>::type thresh(threshSEXP) ;
    Rcpp::traits::input_parameter<int>::type checkfreq(checkfreqSEXP) ;
    
    Rcpp::List result = estimate_hierIRT(alpha_start,
                                 beta_start,
                                 gamma_start,
                                 sigma_start,
                                 eta_start,
                                 y, 
                                 z, 
                                 g, 
                                 i, 
                                 j, 
                                 gammamu,
                                 gammasigma,
                                 betamu,
                                 betasigma, 
                                 sigmav, 
                                 sigmas, 
                                 ND,
                                 NG,
                                 NI,
                                 NJ,
                                 NL,
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
