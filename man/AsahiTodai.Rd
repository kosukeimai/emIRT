\name{AsahiTodai}
\alias{AsahiTodai}

\title{
  Asahi-Todai Elite Survey
}

\description{
  The Asahi-Todai Elite survey was conducted by the University of Tokyo in collaboration with
  a major national newspaper, the Asahi Shimbun, covering all candidates (both incumbents and
  challengers) for the eight Japanese Upper and Lower House elections that occurred between 2003
  and 2013. In six out of eight waves, the survey was also administered to a nationally representative
  sample of voters with the sample size ranging from approximately 1,100 to about 2,000. The novel
  feature of the data is that there are a set of common policy questions, which can be used to scale
  both politicians and voters over time on the same dimension.

  All together, the data set contains a total of N = 19,443 respondents, including 7,734 politicians
  and 11,709 voters. There are J = 98 unique questions in the survey, most of which consisted of
  questions asking for responses on a 5-point Likert scale. However, these scales were collapsed
  into a 3-point Likert scale for estimation with \code{ordIRT()}.  In the data set, we include
  estimates obtained via MCMC using both the 3 and 5-point scale data. See Hirano et al. 2011 for
  more details.

}

\usage{
data(AsahiTodai)
}

\value{

  \item{AsahiTodai}{ list, containing the following elements:
    \itemize{
    \item{\code{dat.all}}{ Survey data, formatted for input to \code{ordIRT()}.}
    \item{\code{start.values}}{ Start values, formatted for input to \code{ordIRT()}. }
    \item{\code{priors}}{ Priors, formatted for input to \code{ordIRT()}.}
    \item{\code{ideal3}}{ Ideal point estimates with data via MCMC, using collapsed 3-category data.}
    \item{\code{ideal5}}{ Ideal point estimates with data via MCMC, using original 5-category data.}
    \item{\code{obs.attri}}{ Attribute data of the respondents.}
    }
  }
}

\references{

Shigeo Hirano, Kosuke Imai, Yuki Shiraito and Masaaki Taniguchi. 2011. ``Policy Positions in
Mixed Member Electoral Systems:Evidence from Japan.'' Working Paper available at \url{http://imai.princeton.edu/research/japan.html}.

Kosuke Imai, James Lo, and Jonathan Olmsted. ``Fast Estimation of Ideal Points with Massive Data.''
Working Paper available at \url{http://imai.princeton.edu/research/fastideal.html}.

}

\seealso{
  '\link{ordIRT}'.
}

\examples{

\dontrun{
### Real data example: Asahi-Todai survey (not run)
## Collapses 5-category ordinal survey items into 3 categories for estimation
data(AsahiTodai)
out.varinf <- ordIRT(.rc = AsahiTodai$dat.all, .starts = AsahiTodai$start.values,
					.priors = AsahiTodai$priors, .D = 1,
					.control = {list(verbose = TRUE,
                     thresh = 1e-6, maxit = 500)})

## Compare against MCMC estimates using 3 and 5 categories
cor(ideal3, out.varinf$means$x)
cor(ideal5, out.varinf$means$x)  
}

}

\keyword{datasets}
