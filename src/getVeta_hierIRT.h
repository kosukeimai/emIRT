// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETVETA_HIERIRT_HPP
#define GETVETA_HIERIRT_HPP

#include <RcppArmadillo.h>

arma::mat getVeta_hierIRT(const arma::mat &i,
                    const arma::mat &j,
                    const arma::mat &g,
                    const arma::mat &curEsigma,
                    const arma::mat &curEbb,
                    const int NI,
                    const int NL
                    );

#endif
