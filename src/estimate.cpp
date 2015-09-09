// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#define DEBUG false

#ifdef _OPENMP
#include <omp.h>
#endif

#include <RcppArmadillo.h>

#include "getEystar.h"

#include "getEb2.h"
#include "getVb2.h"
#include "getEbb.h"
#include "getEba.h"

#include "getEx.h"
#include "getEx2x2.h"
#include "getVx.h"

#include "checkInputs.h"
#include "checkConv.h"

#include "calcFitStats.h"
#include "countVotes.h"

#include "calcLB.h"

using namespace Rcpp ;

List estimate (arma::mat alpha_start,
               arma::mat beta_start,
               arma::mat x_start,
               arma::mat y,
               arma::mat xmu,
               arma::mat xsigma,
               arma::mat betamu,
               arma::mat betasigma,
               unsigned int D = 1,
               unsigned int threads = 0,
               bool verbose = true,
               int convtype = 1,
               unsigned int maxit = 2500,
               double thresh = 1e-9,
               unsigned int checkfreq = 50,
               bool withLB = false,
               bool withProbs = false,
               bool asEM = false
               ) {
    // Init Qtys
    //// Data Parameters
    unsigned int nD = D ;
    unsigned int nJ = y.n_cols ;
    unsigned int nN = y.n_rows ;
    unsigned int nYY = countVotes(y, 1) ;
    unsigned int nYN = countVotes(y, -1) ;
    unsigned int nYnil = countVotes(y, 9) ;
    unsigned int nYna = countVotes(y, 0) ;

    //// Admin
    unsigned int threadsused = 0 ;
    bool withstats = false ;
    unsigned int counter = 0 ;
    bool isconv = false ;
    arma::mat convtrace(maxit - 1,
                        3,
                        arma::fill::zeros
                        ) ;
    arma::mat statstrace(maxit - 1,
                         2,
                         arma::fill::zeros
                         ) ; // GMP, Cor

    arma::mat lbtrace(maxit - 1,
                      1,
                      arma::fill::zeros
                      ) ;

    arma::mat probs1 ;
    arma::mat probsobs ;

    if (asEM) {
        withLB = false ;
    }

    //// Check Inputs
    checkInputs(alpha_start,
                beta_start,
                x_start,
                y,
                xmu,
                xsigma,
                betamu,
                betasigma,
                verbose,
                maxit,
                thresh,
                checkfreq,
                nD,
                threads,
                nN,
                nJ
                ) ;

    //// Init "Current" Containers
    arma::mat curEystar ;

    arma::mat curEa = alpha_start ;
    arma::mat curEb = beta_start ;
    arma::mat curEba = beta_start.t() * alpha_start ;
    arma::mat curEb2 ;

    arma::mat curEbb = arma::eye(D, D) ;
    arma::mat curVb2 = betasigma ;
    arma::mat curVb = curVb2.submat(1, 1, D, D) ;

    arma::mat curEx = x_start ;

    arma::mat curEx2x2 = arma::eye(D + 1, D + 1);
    arma::mat curVx = xsigma ;

    double curgmp = 0.0 ;
    double curlb = 0.0 ;

    //// Init "Old" Containers
    arma::mat oldEa = alpha_start ;
    arma::mat oldEb = beta_start ;
    arma::mat oldEx = x_start ;

    // OpenMP Support
#ifdef _OPENMP
    omp_set_num_threads(1) ;
    if (threads > 0) {
        omp_set_num_threads(threads) ;
        threadsused = omp_get_max_threads() ;
    }
#endif

    // Main Loop Until Convergence
    while (counter < maxit) {

        counter++ ;

        // Update Ystar
        if (DEBUG) {
            Rcout << counter << std::endl ;
            Rcout << "y" << std::endl ;
        }

        curEystar = getEystar(curEa,
                              curEb,
                              curEx,
                              y,
                              nD,
                              nN,
                              nJ
                              ) ;


        // Update [Alpha, Beta]
        if (DEBUG) {
            Rcout << "b" << std::endl ;
        }

        curEx2x2 = getEx2x2(curEx,
                            curVx,
                            nN,
                            nD
                            ) ;
        curVb2 = getVb2(curEx,
                        curEx2x2,
                        betasigma
                        ) ;
        curVb = curVb2.submat(1, 1, D, D) ;

        curEb2 = getEb2(curEystar,
                        curEx,
                        curVb2,
                        betamu,
                        betasigma,
                        nJ,
                        nD,
                        asEM
                        ) ;

        // curEb2.print("curEb2 vi") ;
        curEb2 = getEb2(curEystar,
                        curEx,
                        curVb2,
                        betamu,
                        betasigma,
                        nJ,
                        nD,
                        asEM
                        ) ;

        curEa = curEb2.col(0) ;
        curEb = curEb2.cols(1, D) ;
        curEba = getEba(curEb,
                        curEa,
                        curVb2,
                        nJ,
                        nD,
                        asEM
                        ) ;


        // Update X
        if (DEBUG) {
            Rcout << "x" << std::endl ;
        }

        curEbb = getEbb(curEb,
                        curVb,
                        nJ,
                        nD
                        ) ;
        curVx = getVx(curEb,
                      trans(curEb),
                      curEbb,
                      xsigma
                      ) ;

        // EM and VI
        curEx = getEx(curEystar,
                      curEb,
                      curVx,
                      curEba,
                      xmu,
                      xsigma,
                      nN,
                      nJ,
                      nD,
                      asEM
                      ) ;

        if (DEBUG) {
            curVb.print("Vb") ;
            curEbb.print("Ebb") ;
            curVx.print("Vx") ;
            curEx2x2.print("Ex2x2") ;
        }


        // Check for Interupt & Update Progress
        if (counter % checkfreq == 0) {
            R_CheckUserInterrupt() ;

            if (verbose) {
                Rcout << "Iteration: " << counter << std::endl ;
            }
        }

        // Check Convergence on Its >= 2
        if (counter > 1) {
            isconv = checkConv(oldEx, curEx,
                               oldEa, curEa,
                               oldEb, curEb,
                               nD,
                               counter,
                               thresh,
                               convtrace,
                               convtype
                               ) ;

            // Calculate LB only if Requested
            if (withLB) {
                curlb = calcLB(y,
                               curEystar,
                               curEx,
                               curVx,
                               xmu,
                               xsigma,
                               curEb2,
                               curVb2,
                               betamu,
                               betasigma
                               ) ;
            }
            lbtrace(counter - 2, 0) = curlb ;

            // Store other stats
            if (withstats) {
                // probs1 = calcProb1(curEa, curEb, curEx, nN, nJ) ;
                // probsobs = calcProbObs(probs1, y, nN, nJ) ;
                // curgmp = calcGMP(probsobs, nYnil, nYna) ;
                // statstrace(counter - 2, 0) = curgmp ;
            }
        }


        // Stop if Converged
        if (isconv) {
            // Update Misc Parameters Conditional on Converged X,Beta

            curEx2x2 = getEx2x2(curEx,
                                curVx,
                                nN,
                                nD
                                ) ;
            curVb2 = getVb2(curEx,
                            curEx2x2,
                            betasigma
                            ) ;
            curVb = curVb2.submat(1, 1, D, D) ;

            curEb2 = getEb2(curEystar,
                            curEx,
                            curVb2,
                            betamu,
                            betasigma,
                            nJ,
                            nD,
                            asEM
                            ) ;
            curEa = curEb2.col(0) ;
            curEb = curEb2.cols(1, D) ;
            curEba = getEba(curEb,
                            curEa,
                            curVb2,
                            nJ,
                            nD,
                            asEM
                            ) ;

            curEbb = getEbb(curEb,
                            curVb,
                            nJ,
                            nD
                            ) ;
            curVx = getVx(curEb,
                          trans(curEb),
                          curEbb,
                          xsigma
                          ) ;

            curEx = getEx(curEystar,
                          curEb,
                          curVx,
                          curEba,
                          xmu,
                          xsigma,
                          nN,
                          nJ,
                          nD,
                          asEM
                          ) ;
            break ;
        }

        // Update Old Values If Not Converged
        oldEx = curEx ;
        oldEa = curEa ;
        oldEb = curEb ;
    }

    // Prepare Output
    probs1 = calcProb1(curEa, curEb, curEx, nN, nJ) ;
    probsobs = calcProbObs(probs1, y, nN, nJ) ;
    curgmp = calcGMP(probsobs, nYnil, nYna) ;

    // Returns
    arma::mat convtrace2 = convtrace.rows(0, counter - 2) ;
    arma::mat lbtrace2 = lbtrace.rows(0, counter - 2) ;

    // if (withstats) {
    //     arma::mat statstrace2 = statstrace.rows(0, counter - 2) ;
    // }

    List ret ;
    List means ;
    List vars ;
    List runtime ;
    List fit ;

    arma::mat cs = calcCS(probs1, y, 0.5, nN, nJ) ;

    fit["csr"] = calcCSR(cs, nN, nJ, nYY, nYN) ;
    fit["gmp"] = curgmp ;

    means["x"] = curEx ;
    means["a"] = curEa ;
    means["b"] = curEb ;

    vars["x"] = curVx ;
    vars["beta"] = curVb2 ;

    runtime["iters"] = counter ;
    runtime["conv"] = isconv ;
    runtime["threads"] = threadsused ;
    runtime["tolerance"] = thresh ;

    ret["convstats"] = convtrace2 ;

    if (withLB) {
        ret["lbtrace"] = lbtrace2 ;
    } else {
        ret["lbtrace"] = R_NilValue ;
    }

    if (withProbs) {
        ret["probs"] = probsobs ;
    } else {
        ret["probs"] = R_NilValue ;
    }

    ret["means"] = means ;
    ret["vars"] = vars ;
    ret["runtime"] = runtime ;
    ret["fit"] = fit ;

    ret["n"] = nN ;
    ret["j"] = nJ ;
    ret["d"] = nD ;

    ret["type"] = "vi" ;
    if (asEM) {
        ret["type"] = "em" ;
    }

    return(ret) ;
}
