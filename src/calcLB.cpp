// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>
//#include <RcppTN.h>

#include "entN.h"
#include "enttn1.h"
#include "vtn1.h"

double calcLB (const arma::mat &y,
               const arma::mat &Eystar,
               const arma::mat &Ex,
               const arma::mat &Vx,
               const arma::mat &xmu,
               const arma::mat &xsigma,
               const arma::mat &Eb2,
               const arma::mat &Vb2,
               const arma::mat &betamu,
               const arma::mat &betasigma
               ) {
    /// /////////////////
    /// Common Quantities
    /// /////////////////

    const int nJ = y.n_cols ;
    const int nN = y.n_rows ;
    const int nD = Ex.n_cols ;
    arma::mat ones(nN, 1) ;
    ones.ones() ;
    arma::mat Ex2 = Ex ;
    Ex2.insert_cols(0, ones) ;
    const arma::mat  mu = Ex2 * Eb2.t() ;

    ///
    ///
    ///


    /// /////////////////////////
    /// Variational Distributions
    /// /////////////////////////

    // Latent Propensities
    double elqystar = (0
                       // + -std::log(nN * nJ) // WRONG
                       ) ;
    {
#pragma omp parallel for reduction(+:elqystar)
        for (int n = 0 ; n < nN ; n++) {
            for (int j = 0 ; j < nJ ; j++) {
                double low = y(n,j) == 1 ? 0.0 : R_NegInf ;
                double high = y(n,j) == -1 ? 0.0 : R_PosInf ;
//                double tmpB = RcppTN::enttn1(mu(n, j), 1.0, low, high) ;
                double tmpB = enttn1(mu(n, j), 1.0, low, high) ;
                elqystar += -tmpB ;
            }
        }
    }
    //

    // Ideal Points
    double elqx = (0
                   // + -std::log(nN) // WRONG
                   - nN * entN(Vx)
                   ) ;
    //

    // Bill Parameters
    double elqb2 = (0
                    // + -std::log(nJ) // WRONG
                    - nJ * entN(Vb2)
                    ) ;
    //

    ///
    ///
    ///

    /// //////
    /// Priors
    /// //////

    // Ideal Points
    double elpx = (0
                   // -std::log(nN) // WRONG
                   -(nN * nD) * std::log(2 * M_PI) / 2
                   -(nN / 2) * std::log(arma::det(xsigma))
                   ) ;
    {
#pragma omp parallel for reduction(+:elpx)
        for (int n = 0 ; n < nN ; n++) {
            double tmp = arma::as_scalar(Ex.row(n) * inv(xsigma) * Ex.row(n).t() + // optim
                                         trace(xsigma * Vx) + //optim
                                         xmu.t() * inv(xsigma) * xmu - //optim
                                         2 * xmu.t() * inv(xsigma) * Ex.row(n).t() // optim
                                         ) ;
            elpx += -(1/2) * tmp ;
        }
    }
    //

    // Bill Parameters
    double elpb2 = (0
                    // -std::log(nJ) // WRONG
                    -(nJ * (nD + 1)) * std::log(2 * M_PI) / 2
                    -(nJ / 2) * std::log(arma::det(betasigma))
                    ) ;
    {
#pragma omp parallel for reduction(+:elpb2)
        for (int j = 0 ; j < nJ ; j++) {
            double tmp = arma::as_scalar(Eb2.row(j) * inv(betasigma) * Eb2.row(j).t() + // optim
                                         trace(betasigma * Vb2) + //optim
                                         betamu.t() * inv(betasigma) * betamu - //optim
                                         2 * betamu.t() * inv(betasigma) * Eb2.row(j).t() //optim
                                         ) ;
            elpb2 += -(1/2) * tmp ;
        }
    }
    //

    /// //////////
    /// Likelihood
    /// //////////

    // Indicator
    double elpy = 0.0 ;
    // Summation of all 0's
    //

    // Augmented Likelihood
    double elpystar = (0
                       // -std::log(nN * nJ) // WRONG
                       -((nN * nJ) / 2) * std::log(2 * M_PI)
                       ) ;
    {
#pragma omp parallel for reduction(+:elpystar)
        for (int n = 0 ; n < nN ; n++) {
            for (int j = 0 ; j < nJ ; j++) {
                double low = y(n,j) == 1 ? 0.0 : R_NegInf ;
                double high = y(n,j) == -1 ? 0.0 : R_PosInf ;

                // E[ystar_ij^2]
//                double q1 = (RcppTN::vtn1(mu(n, j), 1.0, low, high) +
                double q1 = (vtn1(mu(n, j), 1.0, low, high) +
                             pow(Eystar(n, j), 2) // optim
                             ) ;

                // E[x_i] E[beta_j alpha_j]
                double q2_ = arma::as_scalar(Ex.row(n) *
                                             (Vb2.submat(1, 0, nD + 1, 0) +
                                              Eb2.submat(j, 1, j, nD).t() * Eb2(j, 0)
                                              )
                                             ) ; //optim

                // E[x2_i beta2_j beta2_j x2_i]
                double q2 = (trace((Vx + Ex.row(n).t() * Ex.row(n)) *
                                   (Vb2.submat(1, 1, nJ, nJ) +
                                    Eb2.submat(j, 1, j, nD + 1).t() * Eb2.submat(j, 1, j, nD + 1))
                                   ) +
                             Vb2(0, 0) + pow(Eb2(j, 0), 2) +
                             2 * q2_
                             ) ;

                // 2 E[ystar_ij] E[x2_i] E[b2_j]
                double q3 = 2 * Eystar(n, j) * mu(n, j) ; // optim

                // .5 * E[istar_ij ^ 2] + E[x2_i beta2_j beta2_j x2_i] - 2 E[ystar_ij] E[x2_i] E[b2_j]
                elpystar += -(1/2) * (q1 + q2 - q3) ;
            }
        }
    }
    //

    ///
    ///
    ///

    /// ///////
    /// Returns
    /// ///////

    const double lb = (elpy
                       + elpystar
                       + elpx
                       + elpb2
                       - elqystar
                       - elqx
                       - elqb2
                       ) ;
    return(lb) ;

    ///
    ///
    ///
}
