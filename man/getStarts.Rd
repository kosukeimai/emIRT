\name{getStarts}
\alias{getStarts}
\title{ Generate Starts for \code{binIRT} }
\description{ 
  \code{getStarts} generates starting values for \code{binIRT}.
}

\usage{
  getStarts(.N, .J, .D, .type = "zeros")
}

\arguments{
  \item{.N}{ integer, number of subjects/legislators to generate starts for.}
  \item{.J}{ integer, number of items/bills to generate starts for.}
  \item{.D}{ integer, number of dimensions.}
  \item{.type}{ ``zeros'' and ``random'' are the only valid types, will generate starts accordingly.}
}

\value{
  \item{\code{alpha}}{  A (J x 1) matrix of starting values for the item difficulty parameter \eqn{alpha}.}
  \item{\code{beta}}{  A (J x D) matrix of starting values for the item discrimination parameter \eqn{\beta}.}
  \item{\code{x}}{  An (N x D) matrix of starting values for the respondent ideal points \eqn{x_i}.}
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
  '\link{binIRT}', '\link{makePriors}', '\link{convertRC}'.
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
