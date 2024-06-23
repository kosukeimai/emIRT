\name{makePriors}
\alias{makePriors}
\title{ Generate Priors for \code{binIRT} }
\description{ 
  \code{makePriors} generates diffuse priors for \code{binIRT}.
}

\usage{
  makePriors(.N = 20, .J = 100, .D = 1)
}

\arguments{
  \item{.N}{ integer, number of subjects/legislators to generate priors for.}
  \item{.J}{ integer, number of items/bills to generate priors for.}
  \item{.D}{ integer, number of dimensions.}

}

\value{
    \describe{
    \item{\code{x$mu}}{  A (D x D) prior means matrix for respondent ideal points \eqn{x_i}.}
    \item{\code{x$sigma}}{  A (D x D) prior covariance matrix for respondent ideal points \eqn{x_i}.}
    \item{\code{beta$mu}}{  A ( D+1 x 1) prior means matrix for \eqn{\alpha_j} and \eqn{\beta_j}.}
    \item{\code{beta$sigma}}{  A ( D+1 x D+1 ) prior covariance matrix for \eqn{\alpha_j} and \eqn{\beta_j}.}
     }
}

\author{ 
Kosuke Imai \email{imai@harvard.edu}

James Lo \email{jameslo@princeton.edu}

Jonathan Olmsted \email{jpolmsted@gmail.com}
}

\references{

Kosuke Imai, James Lo, and Jonathan Olmsted. (2016).  ``Fast Estimation of Ideal Points with Massive Data.''
Working Paper. American Political Science Review, Vol. 110, No. 4 (December), pp. 631-656.

}

\seealso{
  '\link{binIRT}', '\link{getStarts}', '\link{convertRC}'.
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

}
\keyword{ multivariate }
