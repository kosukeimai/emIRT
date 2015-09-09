// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef CALCFITSTATS_H
#define CALCFITSTATS_H

#include <RcppArmadillo.h>

arma::mat calcProb1 (const arma::mat &alpha,
                     const arma::mat &beta,
                     const arma::mat &x,
                     const int N,
                     const int J
                     ) ;

arma::mat calcProbObs (const arma::mat &probs1,
                       const arma::mat &y,
                       const int N,
                       const int J
                       ) ;

arma::mat calcCS (const arma::mat &probs1,
                  const arma::mat &y,
                  const double thresh,
                  const int N,
                  const int J
                  ) ;

double calcCSR (const arma::mat &cs,
                const int N,
                const int J,
                const int nYY,
                const int nYN
                ) ;

double calcGMP (const arma::mat &probsobs,
                const unsigned int nYnil,
                const unsigned int nYna
                ) ;

#endif
