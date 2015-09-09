#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "estimate_ordIRT.h"

RcppExport SEXP ordIRT_estimate(SEXP tau_startSEXP, 
				 				SEXP DD_startSEXP, 
                                 SEXP beta_startSEXP, 
                                 SEXP x_startSEXP,
                                 SEXP ySEXP,
                                 SEXP xmuSEXP, 
                                 SEXP xsigmaSEXP,
                                 SEXP betamuSEXP, 
                                 SEXP betasigmaSEXP, 
                                 SEXP DSEXP, 
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
    Rcpp::traits::input_parameter<arma::mat>::type tau_start(tau_startSEXP);
    Rcpp::traits::input_parameter<arma::mat>::type DD_start(DD_startSEXP);
    Rcpp::traits::input_parameter<arma::mat>::type beta_start(beta_startSEXP) ;
    Rcpp::traits::input_parameter<arma::mat>::type x_start(x_startSEXP) ;
    Rcpp::traits::input_parameter<arma::mat>::type y(ySEXP) ;
    Rcpp::traits::input_parameter<arma::mat>::type xmu(xmuSEXP) ;
    Rcpp::traits::input_parameter<arma::mat>::type xsigma(xsigmaSEXP) ;
    Rcpp::traits::input_parameter<arma::mat>::type betamu(betamuSEXP) ;
    Rcpp::traits::input_parameter<arma::mat>::type betasigma(betasigmaSEXP) ;
    Rcpp::traits::input_parameter<int>::type D(DSEXP) ;
    Rcpp::traits::input_parameter<int>::type threads(threadsSEXP) ;
    Rcpp::traits::input_parameter<bool>::type verbose(verboseSEXP) ;
    Rcpp::traits::input_parameter<int>::type maxit(maxitSEXP) ;
    Rcpp::traits::input_parameter<double>::type thresh(threshSEXP) ;
    Rcpp::traits::input_parameter<int>::type checkfreq(checkfreqSEXP) ;
    
    Rcpp::List result = estimate_ordIRT(tau_start,
                                 DD_start,
                                 beta_start,
                                 x_start,
                                 y, 
                                 xmu,
                                 xsigma,
                                 betamu,
                                 betasigma, 
                                 D,
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
