// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

// Copyright (C) 2012-2015 Jonathan Olmsted
// RcppTN: <https://github.com/olmjo/RcppTN/>
//
// This source code comes from the RcppTN project and is subject to the terms of
// the GNU Public License (>= 2.0). A copy of this license can be obtained at
// <http://www.gnu.org/licenses/gpl-2.0.html>.

const double PI = 3.141592653589793238463 ;

#include <Rcpp.h>

double enttn1(const double mean,
              const double sd,
              const double low,
              const double high
              ) {
    // Normalize bounds
    double alpha = (low - mean) / sd ;
    double beta = (high - mean) / sd ;

    // Calculate Intermediate Qtys
    double q1 = R::pnorm(beta, 0.0, 1.0, true, false) ;
    double q2 = R::pnorm(alpha, 0.0, 1.0, true, false) ;
    double q3 = q1 - q2 ;

    double q4 = R::dnorm(alpha, 0.0, 1.0, false) ;
    double q5 = R::dnorm(beta, 0.0, 1.0, false) ;

    double num1 = alpha * q4 ;
    double num2 = beta * q5 ;
    if (alpha == R_NegInf) {
        num1 = 0.0 ;
    }
    if (beta == R_PosInf) {
        num2 = 0.0 ;
    }

    double num = num1 - num2 ;
    double denom = 2 * q3 ;
    double temp = (.5 * log(2 * PI * exp(1.0)) +
                   log(sd * q3)
                   ) ;

    double res = (temp +
                  num / denom
                  ) ;

    return(res) ;
}
