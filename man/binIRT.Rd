\name{binIRT}
\alias{binIRT}
\title{ Two-parameter Binary IRT estimation via EM }
\description{
  \code{binaryIRT} estimates a binary IRT model with two response categories. Estimation
  is conducted using the EM algorithm described in the reference paper below. The algorithm will
  produce point estimates that are comparable to those of \code{\link[pscl]{ideal}},
  but will do so much more rapidly and also scale better with larger data sets.
}

\usage{
  binIRT(.rc, .starts = NULL, .priors = NULL, .D = 1L, .control = NULL,
  .anchor_subject = NULL, .anchor_outcomes = FALSE)
}

\arguments{
  \item{.rc}{  a list object, in which .rc$votes is a matrix of numeric values containing the data
    to be scaled.  Respondents are assumed to be on rows, and items assumed to be on columns, so the
    matrix is assumed to be of dimension (N x J).  For each item, `1', and `-1' represent different
    responses (i.e. yes or no votes) with `0' as a missing data record.}
  \item{.starts}{ a list containing several matrices of starting values for the parameters.
      The list should contain the following matrices:
    \itemize{
    \item{\code{alpha}}{  A (J x 1) matrix of starting values for the item difficulty parameter \eqn{alpha}.}
    \item{\code{beta}}{  A (J x D) matrix of starting values for the item discrimination parameter \eqn{\beta}.}
    \item{\code{x}}{  An (N x D) matrix of starting values for the respondent ideal points \eqn{x_i}.}
    }
  }
  \item{.priors}{  list, containing several matrices of starting values for the parameters.
     The list should contain the following matrices:
    \itemize{
    \item{\code{x$mu}}{  A (D x D) prior means matrix for respondent ideal points \eqn{x_i}.}
    \item{\code{x$sigma}}{  A (D x D) prior covariance matrix for respondent ideal points \eqn{x_i}.}
    \item{\code{beta$mu}}{  A (D+1 x 1) prior means matrix for \eqn{\alpha_j} and \eqn{\beta_j}.}
    \item{\code{beta$sigma}}{  A (D+1 x D+1) prior covariance matrix for \eqn{\alpha_j} and \eqn{\beta_j}.}
     }
  }
  \item{.D}{ integer, indicates number of dimensions to estimate. Only a
    1 dimension is currently supported. If a higher dimensional model is
    requested, \code{binIRT} exits with an error.}
  \item{.control}{list, specifying some control functions for estimation.  Options include the following:
    \itemize{
    \item{\code{threads}}{  integer, indicating number of cores to use. Default is to use a single core, but more can be
    	supported if more speed is desired.}
    \item{\code{verbose}}{  boolean, indicating whether output during estimation should be verbose or not.  Set FALSE by default.}
    \item{\code{thresh}}{  numeric. Algorithm will run until all
      parameters have a correlation greater than (1 - threshold) across
      consecutive iterations. Set at 1e-6 by default.}
    \item{\code{maxit}}{  integer. Sets the maximum number of iterations the algorithm can run. Set at 500 by default.}
    \item{\code{checkfreq}}{  integer. Sets frequency of verbose output by number of iterations. Set at 50 by default.}
    }
  }
  \item{.anchor_subject}{ integer, the index of the subect to be used in
  anchoring the orientation/polarity of the underlying latent
  dimensions. Defaults to \code{NULL} and no anchoring is done.
}
\item{.anchor_outcomes}{ logical, should an outcomes-based metric be
  used to anchor the orientation of the underlying space. The
  outcomes-based anchoring uses a model-free/non-parametric
  approximation to quantify each item's difficulty and each subject's
  ability. The post-processing then rotates the model-dependent results
  to match the model-free polarity. Defaults to \code{FALSE} and no
  anchoring is done.}
}

\value{
    An object of class \code{binIRT}.

  \item{means}{ list, containing several matrices of point estimates for the parameters corresponding
    to the inputs for the priors.  The list should contain the following matrices.
    \itemize{
    \item{\code{x}}{  A (N x 1) matrix of point estimates for the respondent ideal points \eqn{x_i}.}
    \item{\code{beta}}{   A (J x D+1 ) matrix of point estimates for the item parameters \eqn{\alpha} and \eqn{\beta}.}
    }
   }

  \item{vars}{ list, containing several matrices of variance estimates for parameters corresponding
    to the inputs for the priors.  Note that these variances are those recovered via variational approximation,
    and in most cases they are known to be far too small and generally unusable.  Better estimates of variances
    can be obtained manually via the parametric bootstrap.  The list should contain the following matrices:
    \itemize{
    \item{\code{x}}{  A (N x 1) matrix of variances for the respondent ideal points \eqn{x_i}.}
    \item{\code{beta}}{  A (J x D+1 ) matrix of variances for the item parameters \eqn{\alpha} and \eqn{\beta}.}
    }
   }

  \item{runtime}{ A list of fit results, with elements listed as follows:}
    \itemize{
    \item{\code{iters}}{  integer, number of iterations run.}
    \item{\code{conv}}{  integer, convergence flag. Will return 1 if threshold reached, and 0 if maximum number of iterations reached.}
    \item{\code{threads}}{  integer, number of threads used to estimated model.}
    \item{\code{tolerance}}{  numeric, tolerance threshold for convergence.  Identical to thresh argument in input to .control list.}
    }
  \item{n}{ Number of respondents in estimation, should correspond to number of rows in roll call matrix.}
  \item{j}{ Number of items in estimation, should correspond to number of columns in roll call matrix.}
  \item{d}{ Number of dimensions in estimation.}
  \item{call}{ Function call used to generate output.}

}


\author{
Kosuke Imai \email{kimai@princeton.edu}

James Lo \email{jameslo@princeton.edu}

Jonathan Olmsted \email{jpolmsted@gmail.com}
}

\references{

Kosuke Imai, James Lo, and Jonathan Olmsted ``Fast Estimation of Ideal Points with Massive Data.''
Working Paper. Available at \url{http://imai.princeton.edu/research/fastideal.html}.

}

\seealso{
  '\link{convertRC}', '\link{makePriors}', '\link{getStarts}'.
}

\examples{

## Data from 109th US Senate
data(s109)

## Convert data and make starts/priors for estimation
rc <- convertRC(s109)
p <- makePriors(rc$n, rc$m, 1)
s <- getStarts(rc$n, rc$m, 1)

## Conduct estimates
lout <- binIRT(.rc = rc,
                .starts = s,
                .priors = p,
                .control = {
                    list(threads = 1,
                         verbose = FALSE,
                         thresh = 1e-6
                         )
                }
                )

## Look at first 10 ideal point estimates
lout$means$x[1:10]


lout2 <- binIRT(.rc = rc,
                .starts = s,
                .priors = p,
                .control = {
                    list(threads = 1,
                         verbose = FALSE,
                         thresh = 1e-6
                         )
                },
                .anchor_subject = 2
                )
                                        # Rotates so that Sen. Sessions (R AL)
                                        # has more of the estimated trait

lout3 <- binIRT(.rc = rc,
                .starts = s,
                .priors = p,
                .control = {
                    list(threads = 1,
                         verbose = FALSE,
                         thresh = 1e-6
                         )
                },
                .anchor_subject = 10
                )
                                        # Rotates so that Sen. Boxer (D CA)
                                        # has more of the estimated trait

cor(lout2$means$x[, 1],
    lout3$means$x[, 1]
    )
                                        # = -1 --> same numbers, flipped
                                        # orientation

}
\keyword{ multivariate }
