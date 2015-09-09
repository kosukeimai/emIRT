// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
arma::mat getVx_ordIRT(const arma::mat &Ebb,
                const arma::mat &Edd,
                const arma::mat &sigma2_x) {

	arma::mat Bdd(1,1);
	arma::mat Vx(1,1);

	Bdd = trans(Ebb)*Edd;
	Vx(0,0) = 1/( (1/sigma2_x(0,0)) + Bdd(0,0));

    return(Vx) ;
}
