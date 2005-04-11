\name{TFA.estimate}
\alias{TFA.estimate}

\title{Prediction of Transcription Factor Activities using PLS}
\usage{
TFA.estimate(CONNECdata, GEdata, ncomp=NULL, nruncv=0, alpha=2/3, unit.weights=TRUE)
}
\arguments{
 \item{CONNECdata}{a (n x p) matrix containing the ChIP data for the n genes and the
  p predictors. The n genes must be the same as the n genes of \code{GEdata} and the ordering
  of the genes must also be the same. Each row of \code{ChIPdata} corresponds to a gene, each
  column to a transcription factor. \code{CONNECdata} might have either binary (e.g. 0-1) or numeric
  entries.}
  \item{GEdata}{a (n x m) matrix containing the gene expression levels of the n
 considered genes for m samples. Each row of \code{GEdata} corresponds to a gene, each
  column to a sample.}
    \item{ncomp}{if \code{nruncv=0}, \code{ncomp} is the number of latent 
  components to be  constructed. If \code{nruncv>0}, the  number of 
  latent components to be used for PLS
  regression is chosen from 1,...,\code{ncomp} using the cross-validation procedure described
  in Boulesteix and Strimmer (2005). If \code{ncomp=NULL}, \code{ncomp}
  is set to min(n,p).}
  \item{nruncv}{the number of cross-validation iterations to be performed for the choice of
  the number of latent components. If \code{nruncv=0}, cross-validation is not performed and
  \code{ncomp} latent components are used.}
  \item{alpha}{the proportion of genes to be included in the training set for
  the cross-validation procedure.}
  \item{unit.weights}{If \code{TRUE} then the latent components
    will be constructed from weight vectors that are standardized to length 1, 
    otherwise the weight vectors do not have length 1 but the latent components have
    norm 1.}
}

\description{
 The function \code{TFA.estimate} estimates the transcription factor activities from gene expression data and
  ChIP data using the PLS multivariate regression approach described in Boulesteix and Strimmer
  (2005).
}

\details{
The gene expression data as well as the ChIP data are assumed to have been
properly normalized. However, they do not have to be centered or scaled, since
centering and scaling are performed by the function \code{TFA.estimate} as a
preliminary step. 

The matrix \code{ChIPdata} containing the ChIP data for the n genes and p transcription 
factors might be replaced by  any 'connectivity' matrix whose entries give the strength
of the interactions between the genes and transcription factors. For instance, a connectivity
matrix obtained by aggregating qualitative information from various genomic databases 
might be used as argument  in place of ChIP data.

}

\value{
  A list with the following components:
  \item{TFA}{a (p x m) matrix containing the estimated transcription factor
  activities for the p transcription factors and the m samples.}
  \item{metafactor}{a (m x \code{ncomp}) matrix containing the metafactors for the m
  samples. Each row corresponds to a sample, each column to a metafactor.}
  \item{ncomp}{the number of latent components used in the PLS regression.} 
  }

\references{

A. L. Boulesteix and K. Strimmer (2005). Predicting Transcription Factor Activities 
from Combined Analysis of Microarray and ChIP Data: A Partial Least Squares Approach.
 
 S. de Jong (1993). SIMPLS: an alternative approach to partial least squares 
 regression, Chemometrics Intell. Lab. Syst. \bold{18}, 251--263.

}

\author{
  Anne-Laure Boulesteix (\url{http://www.statistik.lmu.de/~socher/}) 
   and 
  Korbinian Strimmer (\url{http://www.statistik.lmu.de/~strimmer/}).


}

\keyword{regression}

\seealso{\code{\link{pls.regression}}, \code{\link{pls.regression.cv}}.}

\examples{
# load plsgenomics library
library(plsgenomics)

# load Ecoli data
data(Ecoli)

# estimate TFAs based on 3 latent components
TFA.estimate(Ecoli$CONNECdata,Ecoli$GEdata,ncomp=3,nruncv=0)

# estimate TFAs and determine the best number of latent components simultaneously
TFA.estimate(Ecoli$CONNECdata,Ecoli$GEdata,ncomp=1:5,nruncv=20)

}

