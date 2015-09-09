// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef CALCLB_H
#define CALCLB_H

#include <RcppArmadillo.h>

double calcLB (const arma::mat &y,
               const arma::mat &Eystar,
               const arma::mat &Ex,
               const arma::mat &Vx,
               const arma::mat &xmu,
               const arma::mat &xsigma,
               const arma::mat &Eb2,
               const arma::mat &Vb2,
               const arma::mat &betamu,
               const arma::mat &betasigma
               ) ;

# endif
