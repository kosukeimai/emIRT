// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
arma::mat getEzzstar_ordIRT(const arma::mat &Ezstar,
                    const arma::mat &Vzstar,
                    const int N,
                    const int J
                    ) {

  int i,j;
  arma::mat Ezzstar(N,J);

#pragma omp parallel for private(i,j)
  for(i=0; i < N; i++){
	for(j=0; j<J; j++){
	  Ezzstar(i,j) = Ezstar(i,j)*Ezstar(i,j) + Vzstar(i,j);
	}
  }

  return(Ezzstar) ; 
} 



