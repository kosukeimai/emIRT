\name{manifesto}
\alias{manifesto}

\title{
  German Manifesto Data
}

\description{
  A word frequency matrix containing word frequencies from 25 German party manifestos between 1990-2005.
  Obtained from Slapin and Proksch AJPS paper, also used in Lo, Slapin and Proksch.
}

\usage{
data(manifesto)
}

\format{list, containing the following elements:
    \describe{
    \item{data.manif}{ Term-document matrix, formatted for input to \code{ordIRT()}.}
    \item{starts.manif}{ Start values, formatted for input to \code{ordIRT()}. }
    \item{priors.manif}{ Priors, formatted for input to \code{ordIRT()}.}
  }
}


\references{

Jonathan Slapin and Sven-Oliver Proksch. 2009. ``A Scaling Model for Estimating Time-Series Party
Positions from Texts.'' American Journal of Political Science 52(3), 705-722

James Lo, Jonathan Slapin, and Sven-Oliver Proksch. 2016. ``Ideological Clarify in Multiparty
Competition: A New Measure and Test Using Election Manifestos.'' British Journal of Political
Science, 1-20

Kosuke Imai, James Lo, and Jonathan Olmsted. (2016). ``Fast Estimation of Ideal
Points with Massive Data.'' American Political Science Review, Vol. 110, No. 4 (December), pp. 631-656.

}

\seealso{
  '\link{poisIRT}'.
}

\examples{

\dontrun{
## Load German Manifesto data
data(manifesto)

## Estimate variational Wordfish model
lout <- poisIRT(.rc = manifesto$data.manif,
		    i = 0:(ncol(manifesto$data.manif)-1),
		    NI=ncol(manifesto$data.manif),
                    .starts = manifesto$starts.manif,
                    .priors = manifesto$priors.manif,
                    .control = {list(
                        threads = 1,
                        verbose = TRUE,
                        thresh = 1e-6,
						maxit=1000
                        )})

## Positional Estimates for Parties
lout$means$x  
}

}

\keyword{datasets}
