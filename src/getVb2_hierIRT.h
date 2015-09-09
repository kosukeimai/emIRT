// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETVB2_HIERIRT_HPP
#define GETVB2_HIERIRT_HPP

#include <RcppArmadillo.h>

void getVb2_hierIRT(arma::cube &Vb2,
					arma::mat &betasigma,
					arma::cube &Ex2x2,
					const arma::mat &i,
                    const arma::mat &j,
                    const int NL,
                    const int NJ
                    );

#endif
