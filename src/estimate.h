// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef ESTIMATE_H
#define ESTIMATE_H

#include <RcppArmadillo.h>

Rcpp::List estimate(arma::mat alpha_start,
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
                    int convtype = 1,
                    unsigned int maxit = 2500,
                    double thresh = 1e-9,
                    unsigned int checkfreq = 50,
                    bool withLB = false,
                    bool withProbs = false,
                    bool asEM = false
                    ) ;

#endif
