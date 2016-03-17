// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef CHECKCONV_H
#define CHECKCONV_H

#include <RcppArmadillo.h>

bool checkConv(const arma::mat &oldx,
               const arma::mat &curx,
               const arma::mat &olda,
               const arma::mat &cura,
               const arma::mat &oldb,
               const arma::mat &curb,
               const int D,
               const int counter,
               const double thresh,
               arma::mat &convtrace,
               const int convtype
               ) ;

bool checkConv_endorseIRT (const arma::mat &oldalpha,
                           const arma::mat &curalpha,
                           const arma::mat &oldbeta,
                           const arma::mat &curbeta,
                           const arma::mat &oldtheta,
                           const arma::mat &curtheta,
                           const arma::mat &oldw,
                           const arma::mat &curw,
                           const arma::mat &oldgamma,
                           const arma::mat &curgamma,
                           const double thresh,
                           const int convtype
                           ) ;


#endif
