\name{convertRC}
\alias{convertRC}
\title{ Convert Roll Call Matrix Format }
\description{ 
  \code{convertRC} takes a \code{\link[pscl]{rollcall}} object and converts it
  into a format suitable for estimation with '\link{binIRT}'.
}

\usage{
  convertRC(.rc, type = "binIRT") 
}

\arguments{
  \item{.rc}{  a \code{\link[pscl]{rollcall}} object containing votes to be scaled using '\link{binIRT}'.}
  \item{type}{ string, only ``binIRT'' is supported for now, and argument is ignored. }
}

\value{
    An object of class \code{rollcall}, with votes recoded such that yea=1, nay=-1, missing data = 0.
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
  '\link{binIRT}', '\link{makePriors}', '\link{getStarts}'.
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
