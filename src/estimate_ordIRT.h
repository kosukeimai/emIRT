// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef ESTIMATE_ORDIRT_H
#define ESTIMATE_ORDIRT_H

#include <RcppArmadillo.h>

Rcpp::List estimate_ordIRT(arma::mat tau_start,
                    arma::mat DD_start,
                    arma::mat beta_start,
                    arma::mat x_start,
                    arma::mat y,
                    arma::mat xmu,
                    arma::mat xsigma,
                    arma::mat betamu,
                    arma::mat betasigma,
                    unsigned int D = 1,
                    unsigned int threads = 1,
                    bool verbose = true,
                    unsigned int maxit = 2500,
                    double thresh = 1e-9,
                    unsigned int checkfreq = 50
                    ) ;

#endif
