// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

double entN (const arma::mat &sigma
             ) {
    //
    int D = sigma.n_rows ;

    double q1 = (D / 2) * (1 + std::log(2 * arma::datum::pi)) ;
    double q2 = (1/2) * std::log(arma::det(sigma)) ;

    return(q1 + q2) ;
}
