\name{Colon}
\alias{Colon}

\title{Gene expression data from Alon et al. (1999)}
\usage{
data(Colon)
}

\description{
  Gene expression data (2000 genes for 62 samples) from the
 microarray experiments of Colon tissue samples of Alon et al. (1999). 
  
}

\details{ This data set contains 62 samples 
 with 2000 genes: 40 tumor tissues, coded 2 and 22 normal tissues, coded 1.  
}

\value{
A list with the following elements:
  \item{X}{a (62 x 2000) matrix giving the expression levels of 2000 
  genes for the 62 Colon tissue samples. Each row corresponds to a patient, each column to a
  gene.}
  \item{Y}{a numeric vector of length 62 giving the type of tissue sample (tumor or normal).} 
  \item{gene.names}{a vector containing the names of the 2000 genes for the gene
          expression matrix \code{X}.}

}

\source{The data are described in Alon et al. (1999) and can be freely
downloaded from 
\url{http://microarray.princeton.edu/oncology/affydata/index.html}.}

\references{
Alon, U. and Barkai, N. and Notterman, D.A. and Gish, K. and Ybarra, S. and Mack, D. and Levine, A.J. (1999).
Broad patterns of gene expression revealed by clustering analysis of tumor and normal colon tissues probed by oligonucleotide arrays,
Proc. Natl. Acad. Sci. USA,\bold{96}(12), 6745--6750.
}


\examples{
# load plsgenomics library
library(plsgenomics)

# load data set
data(Colon)

# how many samples and how many genes ?
dim(Colon$X)

# how many samples of class 1 and 2 respectively ?
sum(Colon$Y==1)
sum(Colon$Y==2)
}


\keyword{datasets}
