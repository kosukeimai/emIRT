// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifndef GETEGG_HIERIRT_H
#define GETEGG_HIERIRT_H

#include <RcppArmadillo.h>

void getEgg_hierIRT(arma::cube &Egg,
					const arma::mat &Egamma,
					const arma::cube &Vgamma,
                    const int NG
                    );

#endif
