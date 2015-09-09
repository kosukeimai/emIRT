// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

arma::mat calcProb1 (const arma::mat &alpha,
                     const arma::mat &beta,
                     const arma::mat &x,
                     const int N,
                     const int J
                     ) {
    arma::mat probs(N, J, arma::fill::zeros) ;
    arma::mat xbeta = x * trans(beta) ;

#pragma omp parallel for
    for (int n = 0 ; n < N ; n++) {
        for (int j = 0 ; j < J ; j++) {
            probs(n, j) = R::pnorm(alpha(j, 0) + xbeta(n, j),
                                   0.0, 1.0,
                                   1, 0
                                   ) ;
        }
    }
    return(probs) ;
}

arma::mat calcProbObs (const arma::mat &probs1,
                       const arma::mat &y,
                       const int N,
                       const int J
                       ) {
    arma::mat probobs = probs1 ;
    probobs.fill(0.0) ;

#pragma omp parallel for
    for (int n = 0 ; n < N ; n++) {
        for (int j = 0 ; j < J ; j++) {
            if (y(n, j) == 1) {
                // yea
                probobs(n, j) = probs1(n, j) ;
            } else if (y(n, j) == -1) {
                //nay
                probobs(n, j) = 1 - probs1(n, j) ;
            } else {
                // missing, NIL
                probobs(n, j) = 1 ;
            }
        }
    }
    return(probobs) ;
}

arma::mat calcCS (const arma::mat &probs1,
                  const arma::mat &y,
                  const double thresh,
                  const int N,
                  const int J
                  ) {
    // Init Qtys
    arma::mat cs = probs1;
    cs.zeros() ;

    // Check classification for each vote
#pragma omp parallel for
    for (int n = 0 ; n < N ; n++) {
        for (int j = 0 ; j < J ; j++) {
            int ypred = 11 ;

            // Prediction
            if (probs1(n, j) > thresh) {
                ypred = 1 ;
            } else if (probs1(n, j) <= thresh) {
                ypred = -1 ;
            }

            // Assume a class error
            cs(n, j) = -1 ;
            if (y(n, j) == 0) {
                // if missing, missing class
                cs(n, j) = 0 ;
            } else if (y(n, j) == 9) {
                // if notInLeg, nil class
                cs(n, j) = 9 ;
            } else if (y(n, j) == 1) {
                if (ypred == 1) {
                    // obs 1, pred 1, = success
                    cs(n, j) = 1 ;
                }
            } else if (y(n, j) == -1) {
                if (ypred == -1) {
                    cs(n, j) = 1 ;
                }
            }
        }
    }

    // Returns
    return(cs) ;
}

double calcCSR (const arma::mat &cs,
                const int N,
                const int J,
                const int nYY,
                const int nYN
                ) {
    double cssum = 0 ;
    for (int n = 0 ; n < N ; n++) {
        for (int j = 0 ; j < J ; j++) {
            if (cs(n, j) == 1) {
                cssum++ ;
            }
        }
    }
    return((cssum) / (nYY + nYN)) ;
}

double calcGMP (const arma::mat &probsobs,
                const unsigned int nYnil,
                const unsigned int nYna
                ) {

    Rcpp::NumericVector allp(probsobs.begin(), probsobs.end()) ;
    Rcpp::NumericVector alllogp = -log(allp) ;
    std::sort(alllogp.begin(), alllogp.end()) ;

    double tmp = (-1 * std::accumulate(alllogp.begin(), alllogp.end(), 0.0)) ;
    double gmp = exp(tmp / (alllogp.size() - nYnil - nYna)) ;

    return(gmp) ;
}
