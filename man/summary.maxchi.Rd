\name{summary.maxchi}
\alias{summary.maxchi}
\title{ Summary of a maxchi object  }
\description{
Print a summary of a maxchi object
}
\usage{ \method{summary}{maxchi}(object, \dots) }
\arguments{
  \item{object}{ A list object of class \code{maxchi} output by the \code{\link{maxchi}} function }
  \item{\dots}{ Additional arguments to the summary function (currently unused)}
}
\details{
The output of \code{\link{maxchi}} includes information on all site-
and pair-specific chi-squares that exceed the 90th percentile of
the permutation distribution. The \code{summary.maxchi} function computes
maximum site-specific measures (i.e. maxima over all sequence-pair
chi-squares from the same site) and reports all sequence pairs
that tied for these maximum values.
}
\value{
  \item{siteWinlocs}{ Window locations with significant recombination signal }
  \item{siteChisqs}{ The maximum chi-square at these window locations }
  \item{pairs}{ The pair(s) that give the maximum chi-square statistic at each window location }
}
\author{ Brad McNeney <mcneney@stat.sfu.ca>, Jinko Graham <jgraham@stat.sfu.ca>,
Sigal Blay <sblay@sfu.ca> }
\seealso{ \code{\link{maxchi}} }
\examples{ demo(maxchi) }
\keyword{ print }
