\name{stepwise.seqs}
\alias{stepwise.seqs}
\alias{stepwise}
\docType{data}
\title{Simulated sequence alignment from the stepwise package}
\description{
    Simulated sequence alignment comprised of 30 sequences, each of length 1000 bases.  
}
\usage{data(stepwise.seqs)}
\format{
  A 30-by-1000 matrix of mode character.
  Each row represents a genetic sequence.
  Each column is a nucleotide in that genetic sequence. 
  }
\details{
A simulated data set analysed by Graham et al. (2004) to illustrate the stepwise recombination detection methods. The alignment was simulated using Treevolve version 1.3 (Grassly and Holmes 1997). For details on the parameter values used in the simulation, see
Graham et al. (2004).
}
\source{
Graham J, McNeney B and Seillier-Moiseiwitsch F (2004). Stepwise detection of 
recombination breakpoints in sequence alignments. 
\emph{Bioinformatics} Sep 23; [Epub ahead of print]
}
\references{
Grassley NC and Holmes EC (1997). A likelihood method for the detection
of selection and recombination using nucleotide sequences. 
\emph{Mol Biol Evol}, \bold{14(3)}: 239-247
}
\examples{ 
data(stepwise.seqs) 
dim(stepwise.seqs) # 30 x 1000
\dontrun{
# write the 30 sequences to a Phylip-format input file 
cat(paste(nrow(stepwise.seqs),ncol(stepwise.seqs),"\n"),
file="stepwise.phy")
write.table(stepwise.seqs,file="stepwise.phy",append=TRUE,quote=FALSE,
sep="", row.names=TRUE,col.names=FALSE)
}
}
\keyword{datasets}
