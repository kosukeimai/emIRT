// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

void checkInputs (const arma::mat &alpha_start,
                  const arma::mat &beta_start,
                  const arma::mat &x_start,
                  const arma::mat &y,
                  const arma::mat &xmu,
                  const arma::mat &xsigma,
                  const arma::mat &betamu,
                  const arma::mat &betasigma,
                  bool verbose,
                  unsigned int maxit,
                  double thresh,
                  unsigned int checkfreq,
                  unsigned int D,
                  unsigned int threads,
                  unsigned int N,
                  unsigned int J
                  ) {
    if (verbose) {
        Rcpp::Rcout << "Checking for Valid Inputs:" << std::endl ;
    }
    // Simples Checks
    if (verbose) {
        Rcpp::Rcout << "- Control Parameters" << std::endl ;
    }

    if (thresh <= 0) {
        throw std::runtime_error("Threshold not positive.") ;
    }

    if (maxit < 2) {
        throw std::runtime_error("Max. iterations not > 1.") ;
    }

    if (checkfreq < 1) {
        throw std::runtime_error("Check frequency not positve.") ;
    }

    if (threads < 0) {
        throw std::runtime_error("Number of threads not non-negative.") ;
    }
    if (D <= 0) {
        throw std::runtime_error("Number of dimensions not positive.") ;
    }

    // Check Dimensions
    if (verbose) {
        Rcpp::Rcout << "-" << D << " Dimensional Inputs" << std::endl ;
    }

    //// Priors X
    if ((xmu.n_rows != D) |
        (xmu.n_cols != 1)
        ) {
        throw std::runtime_error("X prior mean not D x 1.") ;
    }
    if ((xsigma.n_rows != D) |
        (xsigma.n_cols != D)
        ) {
        throw std::runtime_error("X prior covariance not D x D.") ;
    }
    //// Priors Alpha, Beta
    if ((betamu.n_rows != (D + 1)) |
        (betamu.n_cols != 1)
        ) {
        throw std::runtime_error("Beta prior mean not (D + 1) x 1.") ;
    }
    if ((betasigma.n_rows != (D + 1)) |
        (betasigma.n_cols != (D + 1))
        ) {
        throw std::runtime_error("Beta prior covariance not (D + 1) x (D  + 1)") ;
    }

    //// Starts X
    if ((x_start.n_rows != N) |
        (x_start.n_cols != D)
        ) {
        throw std::runtime_error("X starts not N x D.") ;
    }
    //// Starts Alpha, Beta
    if ((beta_start.n_rows != J) |
        (beta_start.n_cols != D)
        ) {
        throw std::runtime_error("Beta starts not J X D.") ;
    }
    if ((alpha_start.n_rows != J) |
        (alpha_start.n_cols != 1)
        ) {
        throw std::runtime_error("Alpha starts not J X 1.") ;
    }

    // Check Positive-Definiteness of Prior Variances
    arma::mat R ;
    bool test = arma::chol(R, xsigma) ;
    if (!test) {
        throw std::runtime_error("X prior covariance not positive-definite.") ;
    }
    test = arma::chol(R, betasigma) ;
    if (!test) {
        throw std::runtime_error("Beta prior covariance not positive-definite.") ;
    }
}
