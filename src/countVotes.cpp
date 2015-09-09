// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

int countVotes(const arma::mat &y,
               const int code
               ) {
    int cnt = 0 ;
    for (unsigned int n = 0 ; n < y.n_rows ; n++) {
        for (unsigned int j = 0 ; j < y.n_cols ; j++) {
            if (y(n, j) == code) cnt++ ;
        }
    }
    return(cnt) ;
}
