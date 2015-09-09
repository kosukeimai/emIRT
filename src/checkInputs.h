// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

void checkInputs (const arma::mat &alpha_start,
                  const arma::mat &beta_start,
                  const arma::mat &x_start,
                  const arma::mat &y,
                  const arma::mat &xmu,
                  const arma::mat &xsigma,
                  const arma::mat &betamu,
                  const arma::mat &betasigma,
                  bool verbose,
                  unsigned int maxit,
                  double thresh,
                  unsigned int checkfreq,
                  unsigned int D,
                  unsigned int threads,
                  unsigned int N,
                  unsigned int J
                  ) ;
