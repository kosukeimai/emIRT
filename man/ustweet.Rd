\name{ustweet}
\alias{ustweet}

\title{
  U.S. Twitter Following Data
}

\description{
  Data from U.S. Twitter follower data, obtained from Barbera's (2015) replication archive.
}

\usage{data(ustweet)}

\format{list, containing the following elements:
    \describe{
    \item{data}{ Term-document matrix, formatted for input to \code{networkIRT()}.}
    \item{starts}{ Start values, formatted for input to \code{networkIRT()}. }
    \item{priors}{ Priors, formatted for input to \code{networkIRT()}.}
  }
}

\references{

Pablo Barbera. 2015. ``Birds of the Same Feather Tweet Together: Bayesian Ideal Point Estimation
Using Twitter Data.'' Political Analysis 23(1), 76-91

Kosuke Imai, James Lo, and Jonathan Olmsted. (2016). ``Fast Estimation of Ideal Points with Massive Data.''
American Political Science Review, Vol. 110, No. 4 (December), pp. 631-656.

}

\seealso{
  '\link{networkIRT}'.
}

\examples{

\dontrun{
data(ustweet)

## A ridiculously short run to pass CRAN
## For a real test, set maxit to a more reasonable number to reach convergence
lout <- networkIRT(.y = ustweet$data,
                   .starts = ustweet$starts,
                   .priors = ustweet$priors,
                   .control = {list(verbose = TRUE,
                                    maxit = 3,
                                    convtype = 2,
                                    thresh = 1e-6,
                                    threads = 1
                                    )
                           },
                   .anchor_item = 43
                   )

}

}

\keyword{datasets}
