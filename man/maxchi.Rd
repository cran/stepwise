\name{maxchi}
\alias{maxchi}
\title{ Stepwise detection of recombination breakpoints using MaxChi }
\description{
Stepwise detection of recombination breakpoints using the maximum chi-square (MaxChi) method at each step 
}
\usage{
maxchi(input_file, breaks, winHalfWidth, permReps)
}
\arguments{
  \item{input_file}{character string indicating the name of a Phylip format data input file}
  \item{breaks}{an integer vector of ordered site(s) just before the previously declared breakpoints}
  \item{winHalfWidth}{the window half width to use}
  \item{permReps}{the number of Monte Carlo replicates to use for the permutation distribution}
}
\details{
The maxchi function implements the maximum chi-square (MaxChi) method
for detecting recombination breakpoints
(Maynard Smith 1992) using a moving window of fixed width. 
Breakpoints detected in previous steps of a
stepwise search may be conditioned upon. 

For a given position of the moving window on the sequence alignment,
and for a given pair of sequences, a chi-square statistic is
computed to compare two proportions: the proportion
of sites at which the sequences agree
in the left half-window and the proportion of sites in
at which the sequences agree in the right half-window.
Discordance between the two proportions may
reflect a recombination event, located at the window centre,
in the history of the two sequences.
The maximum chi-square over all sequence
pairs is regarded as a summary of the evidence for recombination at
the window centre. The individual chi-square statistics may
also be of interest for suggesting pairs of sequences segments that descend from
historical recombination events.
Significance of observed chi-square statistics
is assessed by a Monte Carlo permutation test.
When conditioning on breakpoints proposed at previous steps of
a stepwise search, permutation is restricted to sites within blocks
defined by the previously proposed breakpoints, as described by
Graham et al. (2004).
}
\value{
  \item{polyposn}{The site numbers of all ungapped polymorphic sites in the alignment}
  \item{chisqs}{ Observed chi-square statistics that exceed the 90th percentile of the permutation null distribution  }
  \item{winlocs}{ Window centres corresponding to the chi-square statistics in \code{chisqs} }
  \item{pairmem1}{ First member of each pair that lead to a significant chi-square statistic in \code{chisqs} }
  \item{pairmem2}{ Second member of each pair }  
  \item{quants}{ 90th and 95th percentiles of the permutation distribution }
}
\references{ 

Graham J, McNeney B and Seillier-Moiseiwitsch F (2004). Stepwise detection of 
recombination breakpoints in sequence alignments. \emph{Bioinformatics} Sep 23; 
[Epub ahead of print]

Maynard Smith J (1992). Analyzing the mosaic structure
of genes. \emph{J Mol Evol}, \bold{34}:126-129.

\url{http://stat-db.stat.sfu.ca/stepwise}
}
\author{ Brad McNeney <mcneney@stat.sfu.ca>, Jinko Graham <jgraham@stat.sfu.ca>,
Sigal Blay <sblay@sfu.ca> }
\seealso{ \code{\link{summary.maxchi}}, \code{\link{phylpro}}  }
\examples{ demo(maxchi) }
\keyword{ htest }
