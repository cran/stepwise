\name{summary.phylpro}
\alias{summary.phylpro}
\title{ Summary of a phylpro object }
\description{
Print a summary of a phylpro object
}
\usage{ \method{summary}{phylpro}(object, \dots) }
\arguments{
  \item{object}{ A list object of class \code{phylpro} output by the \code{\link{phylpro}} function }
  \item{\dots}{ Additional arguments to the summary function (currently unused)}
}
\details{
The output of \code{\link{phylpro}} includes information on all 
site- and target sequence-specific correlations
that exceed the 90th percentile of
the permutation distribution. The \code{summary.phylpro} function computes
minimum site-specific measures (i.e. minima over all target-sequence
correlations from the same site) and reports all target sequences
that tied for these minimum values.
}
\value{
  \item{siteWinlocs}{Window locations with significant recombination signal}
  \item{siteCorrs}{The minimum correlation at these window locations}
  \item{target.seqs}{The target sequence(s) that give the minimum correlation at each window location}
}
\author{ Brad McNeney <mcneney@stat.sfu.ca>, Jinko Graham <jgraham@stat.sfu.ca>,
Sigal Blay <sblay@sfu.ca> }
\seealso{ \code{\link{phylpro}} }
\examples{ demo(phylpro) }
\keyword{ print }
