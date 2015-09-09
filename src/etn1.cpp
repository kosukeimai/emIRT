// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

// Copyright (C) 2012-2015 Jonathan Olmsted
// RcppTN: <https://github.com/olmjo/RcppTN/>
//
// This source code comes from the RcppTN project and is subject to the terms of
// the GNU Public License (>= 2.0). A copy of this license can be obtained at
// <http://www.gnu.org/licenses/gpl-2.0.html>.

# include <Rcpp.h>

using namespace Rcpp ;

/// No Truncation
double e0 (const double mean,
           const double sd,
           const double low,
           const double high
           ) {
    return(mean) ;
}

/// Truncated Below and Above
double e1 (const double mean,
           const double sd,
           const double low,
           const double high
           ) {
    double s_low = (low - mean) / sd ;
    double s_high = (high - mean) / sd ;

    double q1 = R::dnorm(s_low, 0.0, 1.0, false) ;
    double q2 = R::dnorm(s_high, 0.0, 1.0, false) ;
    double q3 = R::pnorm(s_low, 0.0, 1.0, true, false) ;
    double q4 = R::pnorm(s_high, 0.0, 1.0, true, false) ;

    return(mean + sd * ((q1 - q2) / (q4 - q3))) ;
}

/// Truncated Below
double e2 (const double mean,
           const double sd,
           const double low
           ) {
    double s_low = (low - mean) / sd ;

    double q1 = R::dnorm(s_low, 0.0, 1.0, false) ;
    double q3 = R::pnorm(s_low, 0.0, 1.0, true, false) ;

    return(mean + sd * ((q1) / (1 - q3))) ;
}

/// Truncated Above
double e3 (const double mean,
           const double sd,
           const double high
           ) {
    double s_high = (high - mean) / sd ;

    double q2 = R::dnorm(s_high, 0.0, 1.0, false) ;
    double q4 = R::pnorm(s_high, 0.0, 1.0, true, false) ;

    return(mean - sd * ((q2) / (q4))) ;
}


/// Main Function
double etn1(const double mean,
            const double sd,
            const double low,
            const double high
            ) {
    // Init Useful Values
    double out = NA_REAL ;
    //

    if (low == R_NegInf &&
        high == R_PosInf
        ) {
      // No Truncation
        out = e0(mean, sd, low, high) ;
    } else if (low == R_NegInf) {
        // Truncated Above
        out = e3(mean, sd, high) ;
    } else if (high == R_PosInf) {
        // Truncation Below
        out = e2(mean, sd, low) ;
    } else {
        // Truncation Above and Below
        out = e1(mean, sd, low, high) ;
    }
    // Check if Mirror Problem is Numerically Better
    if (!std::isfinite(out)) {
        return(-e1(-mean, sd, -high, -low)) ;
    }
    //
    return(out) ;
}
