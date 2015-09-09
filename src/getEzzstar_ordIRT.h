// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETEZZSTAR_ORDIRT_H
#define GETEZZSTAR_ORDIRT_H

#include <RcppArmadillo.h>

arma::mat getEzzstar_ordIRT(const arma::mat &Ezstar,
                    const arma::mat &Vzstar,
                    const int N,
                    const int J
                    );

#endif
