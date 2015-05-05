\name{SRBCT}
\alias{SRBCT}

\title{Gene expression data from Khan et al. (2001)}
\usage{
data(SRBCT)
}

\description{
  Gene expression data (2308 genes for 83 samples) from the
 microarray experiments of Small Round Blue Cell Tumors (SRBCT) 
 of childhood cancer study of Khan et al. (2001). 
  
}

\details{ This data set contains 83 samples 
 with 2308 genes: 29 cases of Ewing sarcoma (EWS), coded 1, 
 11 cases of Burkitt lymphoma (BL), coded 2, 18 cases of 
 neuroblastoma (NB), coded 3, 25 cases of rhabdomyosarcoma (RMS), coded 4.
 A total of 63 training samples and 25 test samples are 
 provided in Khan et al. (2001). Five of the test set are
 non-SRBCT and are not considered here. The training sample indexes correspond
 to 1:65 and the test sample indexes (without non-SRBCT sample) correspond to 66:83.
  
}

\value{
A list with the following elements:
  \item{X}{a (88 x 2308) matrix giving the expression levels of 2308 
  genes for 88 SRBCT patients. Each row corresponds to a patient, each column to a
  gene.}
  \item{Y}{a numeric vector of length 88 giving the cancer class of each
  patient.} 
  \item{gene.names}{a matrix containing the names of the 2308 genes for the gene
          expression matrix \code{X}. The two columns correspond to
          the gene 'Image.Id.' and 'Gene.Description', respectively.}

}

\source{The data are described in Khan et al. (2001) and can be freely
downloaded from 
\url{http://research.nhgri.nih.gov/microarray/Supplement/}.}

\references{
 Khan, J. and Wei, J. S. and Ringner, M. and Saal, L. H. and Ladanyi,
 M. and Westermann, F. and Berthold, F. and Schwab,
 M. and Antonescu, C. R. and Peterson, C. and Meltzer, P. S. (2001).
 Classification and diagnostic prediction of cancers using gene expression profiling and artificial
 neural networks, Nature Medecine, 7, 673--679.
}



\examples{
# load plsgenomics library
library(plsgenomics)

# load data set
data(SRBCT)

# how many samples and how many genes ?
dim(SRBCT$X)

# how many samples of class 1, 2, 3 and 4, respectively ?
sum(SRBCT$Y==1)
sum(SRBCT$Y==2)
sum(SRBCT$Y==3)
sum(SRBCT$Y==4)}


\keyword{datasets}
