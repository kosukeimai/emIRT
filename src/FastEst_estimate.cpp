// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "estimate.h"

RcppExport SEXP FastEst_estimate(SEXP alpha_startSEXP,
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
                                 SEXP convtypeSEXP,
                                 SEXP threshSEXP,
                                 SEXP checkfreqSEXP,
                                 SEXP withLBSEXP,
                                 SEXP withProbsSEXP,
                                 SEXP asEMSEXP
                                 ) {
  BEGIN_RCPP
    SEXP resultSEXP ;
  {
    Rcpp::RNGScope __rngScope ;
    Rcpp::traits::input_parameter<arma::mat>::type alpha_start(alpha_startSEXP) ;
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
    Rcpp::traits::input_parameter<int>::type convtype(convtypeSEXP) ;
    Rcpp::traits::input_parameter<int>::type maxit(maxitSEXP) ;
    Rcpp::traits::input_parameter<double>::type thresh(threshSEXP) ;
    Rcpp::traits::input_parameter<int>::type checkfreq(checkfreqSEXP) ;
    Rcpp::traits::input_parameter<bool>::type withLB(withLBSEXP) ;
    Rcpp::traits::input_parameter<bool>::type withProbs(withProbsSEXP) ;
    Rcpp::traits::input_parameter<bool>::type asEM(asEMSEXP) ;

    Rcpp::List result = estimate(alpha_start,
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
                                 convtype,
                                 maxit,
                                 thresh,
                                 checkfreq,
                                 withLB,
                                 withProbs,
                                 asEM
                                 ) ;
    PROTECT(resultSEXP = Rcpp::wrap(result)) ;
  }
  UNPROTECT(1);
  return(resultSEXP) ;
  END_RCPP
    }
