// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

bool checkConv(const arma::mat &oldx,
               const arma::mat &curx,
               const arma::mat &olda,
               const arma::mat &cura,
               const arma::mat &oldb,
               const arma::mat &curb,
               const int D,
               const int counter,
               const double thresh,
               arma::mat &convtrace,
               const int convtype
               ) {
    // Init Qtys
    bool isconv = false ;
    double devX = 100.0;
    double devA = 100.0;
    double devB = 100.0 ;

    if (convtype == 1) {
        // correlation
        devX = 1 - (cor(oldx, curx)).min() ;
        devA = 1 - (cor(olda, cura)).min() ;
        devB = 1 - (cor(oldb, curb)).min() ;
    }
    if (convtype == 2) {
        // maximum absolute deviation
        devX = (abs(curx - oldx)).max() ;
        devA = (abs(cura - olda)).max() ;
        devB = (abs(curb - oldb)).max() ;
    }

    // Rcout << devX << " " << devA << " " << devB << std::endl ;

    convtrace(counter - 2, 0) = devX ;
    convtrace(counter - 2, 1) = devA ;
    convtrace(counter - 2, 2) = devB ;

    bool check = (devX < thresh) & (devA < thresh) & (devB < thresh) ;

    if (check) {
        return(true) ;
    } else {
        return(false) ;
    }

    return(isconv) ;
}




bool checkConv_endorseIRT (const arma::mat &oldalpha,
                           const arma::mat &curalpha,
                           const arma::mat &oldbeta,
                           const arma::mat &curbeta,
                           const arma::mat &oldtheta,
                           const arma::mat &curtheta,
                           const arma::mat &oldw,
                           const arma::mat &curw,
                           const arma::mat &oldgamma,
                           const arma::mat &curgamma,
                           const double thresh,
                           const int convtype
                           ) {
    // Init Qtys
    bool isconv = false ;
    double devalpha = 100.0;
    double devbeta = 100.0;
    double devtheta = 100.0;
    double devw = 100.0;
    double devgamma = 100.0;

    if (convtype == 1) {
        // correlation
        devalpha = 1 - (cor(oldalpha, curalpha)).min() ;
        devbeta = 1 - (cor(oldbeta, curbeta)).min() ;
        devtheta = 1 - (cor(oldtheta, curtheta)).min() ;
        devw = 1 - (cor(oldw, curw)).min() ;
        devgamma = 1 - (cor(oldgamma, curgamma)).min() ;
    }
    if (convtype == 2) {
        // maximum absolute deviation
        devalpha = (abs(curalpha - oldalpha)).max() ;
        devbeta = (abs(curbeta - oldbeta)).max() ;
        devtheta = (abs(curtheta - oldtheta)).max() ;
        devw = (abs(curw - oldw)).max() ;
        devgamma = (abs(curgamma - oldgamma)).max() ;
    }

    // Rcout << devX << " " << devA << " " << devB << std::endl ;

    bool check = (devalpha < thresh) &
        (devbeta < thresh) &
        (devtheta < thresh) &
        (devw < thresh) &
        (devgamma < thresh) ;

    if (check) {
        return(true) ;
    } else {
        return(false) ;
    }

    return(isconv) ;
}
