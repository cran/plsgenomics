\name{leukemia}
\alias{leukemia}

\title{Gene expression data from Golub et al. (1999)}
\usage{
data(leukemia)
}

\description{
  Gene expression data (3051 genes and 38 tumor mRNA samples) from the
leukemia microarray study of Golub et al. (1999). 
  
}

\details{
  
}

\value{
A list with the following elements:
  \item{X}{a (38 x 3051) matrix giving the expression levels of 3051 
  genes for 38 leukemia patients. Each row corresponds to a patient, each column to a
  gene.}
  \item{Y}{a numeric vector of length 38 giving the cancer class of each
  patient.} 
  \item{gene.names}{a matrix containing the names of the 3051 genes for the gene
          expression matrix \code{X}. The three columns correspond to
          the gene 'index', 'ID', and 'Name', respectively.
}
  
}

\source{The dataset was taken from
the R package multtest. The data are described in Golub et al. (1999) and can be freely
downloaded from \url{http://www-genome.wi.mit.edu/MPR/}.}

\references{
 S. Dudoit, J. Fridlyand and T. P. Speed (2002). Comparison of discrimination
 methods for the classification of tumors using gene expression data, Journal of
 the American Statistical Association \bold{97}, 77--87.  
 
 Golub et al. (1999). Molecular classification of cancer: class discovery
 and class prediction by gene expression monitoring, Science \bold{286},
 531--537.
}

\examples{
# load plsgenomics library
library(plsgenomics)

# load data set
data(leukemia)

# how many samples and how many genes ?
dim(leukemia$X)

# how many samples of class 1 and 2, respectively ?
sum(leukemia$Y==1)
sum(leukemia$Y==2)}


\keyword{datasets}
