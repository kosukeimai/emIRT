// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#define DEBUG false

#ifdef _OPENMP
#include <omp.h>
#endif

#include <RcppArmadillo.h>

#include "getEystar_endorseIRT.h"

#include "getEalpha_endorseIRT.h"
#include "getEbeta_endorseIRT.h"
#include "getEgamma_endorseIRT.h"

#include "getEtheta_endorseIRT.h"
#include "getEw_endorseIRT.h"

#include "updateHigherMoments_endorseIRT.h"

#include "checkConv.h"

using namespace Rcpp ;


void checkContainer (arma::mat container,
                     char* note
                     ) {
    if (!container.is_finite()) {
        Rcout << note << std::endl ;
    }

}

List estimate_endorseIRT (arma::mat alpha_start,
                          arma::mat beta_start,
                          arma::mat w_start,
                          arma::mat theta_start,
                          arma::mat gamma_start,
                          arma::mat y,
                          arma::mat alphamu,
                          arma::mat alphasigma,
                          arma::mat betamu,
                          arma::mat betasigma,
                          arma::mat wmu,
                          arma::mat wsigma,
                          arma::mat thetamu,
                          arma::mat thetasigma,
                          arma::mat gammamu,
                          arma::mat gammasigma,
                          unsigned int threads = 0,
                          bool verbose = true,
                          unsigned int maxit = 2500,
                          double thresh = 1e-9,
                          unsigned int checkfreq = 50,
                          int convtype = 1
                          ) {

    // Init Qtys

    //// Data Parameters
    unsigned int nJ = y.n_cols ;
    unsigned int nN = y.n_rows ;

    //// Admin
    unsigned int threadsused = 0 ;
    unsigned int counter = 0 ;
    bool isconv = false ;

    //// Init "Current" Containers
    arma::mat curEystar = y ;

    arma::mat curEalpha = alpha_start ;
    arma::mat curValpha = alpha_start * 0 + alphasigma(0, 0) ;

    arma::mat curEbeta = beta_start ;
    arma::mat curVbeta = beta_start * 0 + betasigma(0, 0) ;

    arma::mat curEw = w_start ;
    arma::mat curVw = w_start * 0 + wsigma(0, 0) ;

    arma::mat curEw2 = pow(curEw, 2) ;
    arma::mat curEw3 = pow(curEw, 3) ;
    arma::mat curEw4 = pow(curEw, 4) ;

    arma::mat curEtheta = theta_start ;
    arma::mat curVtheta = theta_start * 0 + thetasigma(0, 0) ;

    arma::mat curEtheta2 = pow(curEtheta, 2) ;
    arma::mat curEtheta3 = pow(curEtheta, 3) ;
    arma::mat curEtheta4 = pow(curEtheta, 4) ;

    arma::mat curEgamma = gamma_start ;
    arma::mat curVgamma = gammasigma ;

    arma::mat curEgamma2 = pow(curEgamma, 2) ;

    //// Init "Old" Containers
    arma::mat oldEalpha = alpha_start ;
    arma::mat oldEbeta = beta_start ;

    arma::mat oldEw = w_start ;
    arma::mat oldEtheta = theta_start ;

    arma::mat oldEgamma = gamma_start ;

    // OpenMP Support
#ifdef _OPENMP
    omp_set_num_threads(1) ;
    if (threads > 0) {
        omp_set_num_threads(threads) ;
        threadsused = omp_get_max_threads() ;
    }
#endif

    // Rcout << "par: " << threads << " " << threadsused << std::endl ;

    // Update Higher Moments 1x for consistency

    // updateHigherMoments_endorseIRT(curEtheta,
    //                                curVtheta,
    //                                curEw,
    //                                curVw,
    //                                curEgamma,
    //                                curVgamma,
    //                                curEtheta2,
    //                                curEtheta3,
    //                                curEtheta4,
    //                                curEw2,
    //                                curEw3,
    //                                curEw4,
    //                                curEgamma2
    //                                ) ;

    // Main Loop Until Convergence
    while (counter < maxit) {
        counter++ ;

        // Update ystar

        getEystar_endorseIRT(curEalpha,
                             curEbeta,
                             curEw,
                             curEtheta,
                             curEgamma,
                             y,
                             nN,
                             nJ,
                             curEystar,
                             curEtheta2,
                             curEw2
                             ) ;


        // Update alpha

        // checkContainer(curEystar, "ystar") ;

        getEalpha_endorseIRT(curEystar,
                             curEbeta,
                             curEtheta,
                             curEw,
                             curEgamma,
                             alphamu,
                             alphasigma,
                             nN,
                             nJ,
                             curEalpha,
                             curValpha,
                             curEtheta2,
                             curEw2
                             ) ;


        // checkContainer(curEalpha, "alpha") ;

        // Update beta

        getEbeta_endorseIRT(curEystar,
                            curEalpha,
                            curEtheta,
                            curEw,
                            curEgamma,
                            betamu,
                            betasigma,
                            nN,
                            nJ,
                            curEbeta,
                            curVbeta,
                            curEtheta2,
                            curEw2
                            ) ;

        // checkContainer(curEalpha, "beta") ;


        // update w

        getEw_endorseIRT(curEystar,
                         curEalpha,
                         curEbeta,
                         curEtheta,
                         curEgamma,
                         wmu,
                         wsigma,
                         nN,
                         nJ,
                         oldEw,
                         curEw,
                         curVw,
                         curEgamma2,
                         curEtheta2,
                         curEtheta3
                         ) ;

        // checkContainer(curEalpha, "w") ;

        // updateHigherMoments_endorseIRT(curEtheta,
        //                                curVtheta,
        //                                curEw,
        //                                curVw,
        //                                curEgamma,
        //                                curVgamma,
        //                                curEtheta2,
        //                                curEtheta3,
        //                                curEtheta4,
        //                                curEw2,
        //                                curEw3,
        //                                curEw4,
        //                                curEgamma2
        //                                ) ;

        // Update theta

        getEtheta_endorseIRT(curEystar,
                             curEalpha,
                             curEbeta,
                             curEw,
                             curEgamma,
                             thetamu,
                             thetasigma,
                             nN,
                             nJ,
                             oldEtheta,
                             curEtheta,
                             curVtheta,
                             curEw2,
                             curEw3,
                             curEgamma2
                             ) ;

        // checkContainer(curEalpha, "theta") ;

        // updateHigherMoments_endorseIRT(curEtheta,
        //                                curVtheta,
        //                                curEw,
        //                                curVw,
        //                                curEgamma,
        //                                curVgamma,
        //                                curEtheta2,
        //                                curEtheta3,
        //                                curEtheta4,
        //                                curEw2,
        //                                curEw3,
        //                                curEw4,
        //                                curEgamma2
        //                                ) ;

        // Update gamma

        // curEgamma = getEgamma_endorseIRT(curEystar,
        //                                  curEalpha,
        //                                  curEbeta,
        //                                  curEtheta,
        //                                  curEw,
        //                                  gammamu,
        //                                  gammasigma,
        //                                  nN,
        //                                  nJ,
        //                                  curVgamma,
        //                                  curEtheta2,
        //                                  curEtheta3,
        //                                  curEtheta4,
        //                                  curEw2,
        //                                  curEw3,
        //                                  curEw4
        //                                  ) ;

        // checkContainer(curEalpha, "gamma") ;

        // updateHigherMoments_endorseIRT(curEtheta,
        //                                curVtheta,
        //                                curEw,
        //                                curVw,
        //                                curEgamma,
        //                                curVgamma,
        //                                curEtheta2,
        //                                curEtheta3,
        //                                curEtheta4,
        //                                curEw2,
        //                                curEw3,
        //                                curEw4,
        //                                curEgamma2
        //                                ) ;



        // Check for Interupt & Update Progress
        if (counter % checkfreq == 0) {
            R_CheckUserInterrupt() ;

            if (verbose) {
                Rcout << "Iteration: " << counter << std::endl ;
            }
        }

        // Check Convergence on Its >= 2
        if (counter > 1) {
            isconv = checkConv_endorseIRT(oldEalpha, curEalpha,
                                          oldEbeta, curEbeta,
                                          oldEtheta, curEtheta,
                                          oldEw, curEw,
                                          oldEgamma, curEgamma,
                                          thresh,
                                          convtype
                                          ) ;
        }


        // Stop if Converged
        if (isconv) {
            // Update Misc Parameters Conditional on Converged Params

            break ;
        }

        // Update Old Values If Not Converged
        oldEalpha = curEalpha ;
        oldEbeta = curEbeta ;
        oldEw = curEw ;
        oldEtheta = curEtheta ;
        oldEgamma = curEgamma ;
    }

    // Prepare Output

    // Rescale to RM Over Parametrization

    // curEw = curEw * pow(curEgamma(0, 0), 1/2) ;
    // curEtheta = curEtheta * pow(curEgamma(0, 0), 1/2) ;

    // curVw = curVw * curEgamma(0, 0) ;
    // curVtheta = curVtheta * curEgamma(0, 0) ;


    // curEgamma = curEgamma / curEgamma(0, 0) ;
    // curVgamma = curVgamma / pow(curEgamma(0, 0), 2) ;



    // Returns

    List ret ;
    List means ;
    List vars ;
    List runtime ;
    List fit ;

    means["alpha"] = curEalpha ;
    means["beta"] = curEbeta ;
    means["w"] = curEw ;
    means["theta"] = curEtheta ;

    // means["gamma"] = curEgamma ;

    vars["w"] = curVw ;
    vars["theta"] = curVtheta ;

    // vars["gamma"] = curVgamma ;

    vars["alpha"] = curValpha ;
    vars["beta"] = curVbeta ;

    runtime["iters"] = counter ;
    runtime["conv"] = isconv ;
    runtime["threads"] = threadsused ;
    runtime["tolerance"] = thresh ;

    ret["means"] = means ;
    ret["vars"] = vars ;
    ret["runtime"] = runtime ;

    ret["N"] = nN ;
    ret["J"] = nJ ;

    return(ret) ;
}
