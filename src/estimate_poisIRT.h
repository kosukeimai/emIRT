// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef ESTIMATE_POISIRT_HPP
#define ESTIMATE_POISIRT_HPP

#include <RcppArmadillo.h>

Rcpp::List estimate_poisIRT (arma::mat alpha_start,
               arma::mat psi_start,
               arma::mat beta_start,
               arma::mat x_start,
               arma::mat Y,
               arma::mat i,
               int ni,
               double psi_mu,
               double psi_sigma,
               double alpha_mu,
               double alpha_sigma,
               double beta_mu,
               double beta_sigma,
               double x_mu,
               double x_sigma,
               unsigned int threads = 0,
               bool verbose = true,
               unsigned int maxit = 2500,
               double thresh = 1e-9,
               unsigned int checkfreq = 50
               ) ;

#endif
