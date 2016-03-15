// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "estimate_endorseIRT.h"

RcppExport SEXP endorseIRT_estimate(SEXP alpha_startSEXP,
                                    SEXP beta_startSEXP,
                                    SEXP w_startSEXP,
                                    SEXP theta_startSEXP,
                                    SEXP gamma_startSEXP,
                                    SEXP ySEXP,
                                    SEXP alphamuSEXP,
                                    SEXP alphasigmaSEXP,
                                    SEXP betamuSEXP,
                                    SEXP betasigmaSEXP,
                                    SEXP wmuSEXP,
                                    SEXP wsigmaSEXP,
                                    SEXP thetamuSEXP,
                                    SEXP thetasigmaSEXP,
                                    SEXP gammamuSEXP,
                                    SEXP gammasigmaSEXP,
                                    SEXP threadsSEXP,
                                    SEXP verboseSEXP,
                                    SEXP maxitSEXP,
                                    SEXP threshSEXP,
                                    SEXP checkfreqSEXP,
                                    SEXP convtypeSEXP
                                    ) {
    BEGIN_RCPP
        SEXP resultSEXP ;
    {
        Rcpp::RNGScope __rngScope ;
        Rcpp::traits::input_parameter<arma::mat>::type alpha_start(alpha_startSEXP) ;
        Rcpp::traits::input_parameter<arma::mat>::type beta_start(beta_startSEXP) ;
        Rcpp::traits::input_parameter<arma::mat>::type w_start(w_startSEXP) ;
        Rcpp::traits::input_parameter<arma::mat>::type theta_start(theta_startSEXP) ;
        Rcpp::traits::input_parameter<arma::mat>::type gamma_start(gamma_startSEXP) ;
        Rcpp::traits::input_parameter<arma::mat>::type y(ySEXP) ;
        Rcpp::traits::input_parameter<arma::mat>::type alphamu(alphamuSEXP) ;
        Rcpp::traits::input_parameter<arma::mat>::type alphasigma(alphasigmaSEXP) ;
        Rcpp::traits::input_parameter<arma::mat>::type betamu(betamuSEXP) ;
        Rcpp::traits::input_parameter<arma::mat>::type betasigma(betasigmaSEXP) ;
        Rcpp::traits::input_parameter<arma::mat>::type wmu(wmuSEXP) ;
        Rcpp::traits::input_parameter<arma::mat>::type wsigma(wsigmaSEXP) ;
        Rcpp::traits::input_parameter<arma::mat>::type thetamu(thetamuSEXP) ;
        Rcpp::traits::input_parameter<arma::mat>::type thetasigma(thetasigmaSEXP) ;
        Rcpp::traits::input_parameter<arma::mat>::type gammamu(gammamuSEXP) ;
        Rcpp::traits::input_parameter<arma::mat>::type gammasigma(gammasigmaSEXP) ;
        Rcpp::traits::input_parameter<int>::type threads(threadsSEXP) ;
        Rcpp::traits::input_parameter<bool>::type verbose(verboseSEXP) ;
        Rcpp::traits::input_parameter<int>::type maxit(maxitSEXP) ;
        Rcpp::traits::input_parameter<double>::type thresh(threshSEXP) ;
        Rcpp::traits::input_parameter<int>::type checkfreq(checkfreqSEXP) ;
        Rcpp::traits::input_parameter<int>::type convtype(convtypeSEXP) ;

        // Rcpp::List result ;

        Rcpp::List result = estimate_endorseIRT (alpha_start,
                                                 beta_start,
                                                 w_start,
                                                 theta_start,
                                                 gamma_start,
                                                 y,
                                                 alphamu,
                                                 alphasigma,
                                                 betamu,
                                                 betasigma,
                                                 wmu,
                                                 wsigma,
                                                 thetamu,
                                                 thetasigma,
                                                 gammamu,
                                                 gammasigma,
                                                 threads,
                                                 verbose,
                                                 maxit,
                                                 thresh,
                                                 checkfreq,
                                                 convtype
                                                 ) ;
        PROTECT(resultSEXP = Rcpp::wrap(result)) ;
    }
    UNPROTECT(1);
    return(resultSEXP) ;
    END_RCPP
        }
