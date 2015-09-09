// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETEBA_HIERIRT_H
#define GETEBA_HIERIRT_H

#include <RcppArmadillo.h>

void getEba_hierIRT(arma::mat &Eba,
					const arma::mat &Eb,
					const arma::mat &Ea,
					const arma::cube &Vb2,
                    const int NJ
                    );

#endif
