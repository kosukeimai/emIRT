// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETEBB_HIERIRT_H
#define GETEBB_HIERIRT_H

#include <RcppArmadillo.h>

void getEbb_hierIRT(arma::mat &Ebb,
					const arma::mat &Eb,
					const arma::cube &Vb2,
                    const int NJ
                    ) ;

#endif
