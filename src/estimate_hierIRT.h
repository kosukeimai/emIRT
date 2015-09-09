// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef ESTIMATE_HIERIRT_H
#define ESTIMATE_HIERIRT_H

#include <RcppArmadillo.h>

Rcpp::List estimate_hierIRT(arma::mat alpha_start,
               arma::mat beta_start,
               arma::mat gamma_start,
               arma::mat sigma_start,
               arma::mat eta_start,
               arma::mat y,
               arma::mat z,
               arma::mat g,
               arma::mat i,
               arma::mat j,
               arma::mat gammamu,
               arma::mat gammasigma,
               arma::mat betamu,
               arma::mat betasigma,
               arma::mat sigmav,
               arma::mat sigmas,
               unsigned int ND,
               unsigned int NG,
               unsigned int NI,
               unsigned int NJ,
               unsigned int NL,
               unsigned int threads = 0,
               bool verbose = true,
               unsigned int maxit = 2500,
               double thresh = 1e-9,
               unsigned int checkfreq = 50
               );

#endif
