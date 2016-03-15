// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef UPDATEHIGHERMOMENTS_ENDORSEIRT_H
#define UPDATEHIGHERMOMENTS_ENDORSEIRT_H

#include <RcppArmadillo.h>

void  updateHigherMoments_endorseIRT (const arma::mat &Etheta,
                                      const arma::mat &Vtheta,
                                      const arma::mat &Ew,
                                      const arma::mat &Vw,
                                      const arma::mat &Egamma,
                                      const arma::mat &Vgamma,
                                      arma::mat &Etheta2,
                                      arma::mat &Etheta3,
                                      arma::mat &Etheta4,
                                      arma::mat &Ew2,
                                      arma::mat &Ew3,
                                      arma::mat &Ew4,
                                      arma::mat &Egamma2
                                      ) ;

#endif
