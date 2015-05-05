\name{Ecoli}
\alias{Ecoli}

\title{Ecoli gene expression and connectivity data from Kao et al. (2003)}
\usage{
data(Ecoli)
}

\description{
  Gene expression data obtained during Escherichia coli carbon source transition
  and connectivity data from the RegulonDB data base (Salgado et al., 2001). 
  The experiments  and data sets are described in Kao et al. (2003).
}



\value{
  A list with the following components:
  \item{CONNECdata}{a (100 x 16) matrix containing the connectivity data for 100 genes and 16 
  regulators. The data are coded as 1 (positive interaction), 0 (no interaction) and -1 (negative
  interaction).}
  \item{GEdata}{a (100 x 23) matrix containing gene expression data for 100 genes and 23 
  samples corresponding to differents times during carbon source transition.}
  \item{timepoint}{a numeric vector of length 23 containing the time points (in
  hours) for the 23 samples.}
  
}

\references{
 K. Kao, Y.-L. Yang, R. Boscolo, C. Sabatti, V. Roychowdhury and J. C. Liao
 (2004). Transcriptome-based determination of multiple transcription regulator
 activities in Escherichia coli by using network component analysis, PNAS
 \bold{101}, 641--646. 
 
 H. Salgado, A. Santos-Zavaleta, S. Gama-Castro, D. Millan-Zarate, E.
 Diaz-Peredo, F. Sanchez-Solano, E. Perez-Rueda, C. Bonavides-Martinez and J.
 Collado-Vides (2001). RegulonDB (version 3.2): transcriptional regulation and
 operon organization in Escherichia coli K-12, Nucleic Acids Research \bold{29},
 72--74. 
}

\source{The data are described in Kao et al. (2004) and can be freely downloaded from 
\url{http://www.seas.ucla.edu/~liaoj/downloads.html}.}

\examples{
# load plsgenomics library
library(plsgenomics)

# load data set
data(Ecoli)

# how many genes and how many transcription factors ?
dim(Ecoli$CONNECdata)


}
\keyword{datasets}
