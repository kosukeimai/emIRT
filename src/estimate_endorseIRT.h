// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef ESTIMATE_ENDORSEIRT_H
#define ESTIMATE_ENDORSEIRT_H

#include <RcppArmadillo.h>

Rcpp::List estimate_endorseIRT (arma::mat alpha_start,
                                arma::mat beta_start,
                                arma::mat w_start,
                                arma::mat theta_start,
                                arma::mat gamma_start,
                                arma::mat y,
                                arma::mat alphamu,
                                arma::mat alphasigma,
                                arma::mat betamu,
                                arma::mat betasigma,
                                arma::mat wmu,
                                arma::mat wsigma,
                                arma::mat thetamu,
                                arma::mat thetasigma,
                                arma::mat gammamu,
                                arma::mat gammasigma,
                                unsigned int threads = 0,
                                bool verbose = true,
                                unsigned int maxit = 2500,
                                double thresh = 1e-9,
                                unsigned int checkfreq = 50,
                                int convtype = 1
                                ) ;
#endif
