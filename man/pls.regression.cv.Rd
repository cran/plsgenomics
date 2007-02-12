\name{pls.regression.cv}
\alias{pls.regression.cv}

\title{Determination of the number of latent components to be used in PLS regression}
\usage{
pls.regression.cv(Xtrain, Ytrain, ncomp, nruncv=20, alpha=2/3)
}
\arguments{
  \item{Xtrain}{a (ntrain x p) data matrix containing the predictors for the training data set.
   \code{Xtrain} may be a matrix or a data frame. Each row is an observation and each column is a 
   predictor variable.}
  \item{Ytrain}{a (ntrain x m) data matrix of responses. \code{Ytrain} may be a vector (if m=1), a matrix or a data frame. 
  If \code{Ytrain} is a matrix or a data
  frame, each row is an observation and each column is a response variable. If
  \code{Ytrain} is a vector, it contains the unique response variable for each
  observation.}
  \item{ncomp}{the vector of integers from which the best number of latent 
  components has to be chosen by cross-validation. If \code{ncomp} is of length
  1, the best number of components is chosen from  1,...,\code{ncomp}.}
  \item{nruncv}{the number of cross-validation iterations to be performed for the choice of
  the number of latent components.}
  \item{alpha}{the proportion of observations to be included in the training set at
  each cross-validation iteration.}
}

\description{
 The function \code{pls.regression.cv} determines the best number of latent components to be used for PLS regression using the cross-validation
  approach described in Boulesteix and Strimmer (2005).
}

\details{

  The cross-validation procedure described in Boulesteix and Strimmer (2005) 
  is used to  determine the best number of latent components to be used for classification.
  At each cross-validation run, \code{Xtrain} is split into a pseudo training
  set and a pseudo test set and the squared error is determined for each
  number of latent components. Finally, the function \code{pls.regression.cv} returns
  the number of latent components for which the mean squared error over
  the \code{nrun} partitions is minimal. 
}

\value{
  The number of latent components to be used in PLS regression, as
  determined by cross-validation.
  }

\references{

 A. L. Boulesteix and K. Strimmer (2005). Predicting Transcription Factor Activities 
 from Combined Analysis of Microarray and ChIP Data: A Partial Least Squares Approach.
 
  A. L. Boulesteix, K. Strimmer (2007). Partial least squares: a versatile tool for the analysis of high-dimensional genomic data. Briefings in Bioinformatics 7:32-44.

 S. de Jong (1993). SIMPLS: an alternative approach to partial least squares 
 regression, Chemometrics Intell. Lab. Syst. \bold{18}, 251--263.

}

\author{
  Anne-Laure Boulesteix (\url{http://www.slcmsr.net/boulesteix}) 
   and 
  Korbinian Strimmer (\url{http://strimmerlab.org/}).


}

\keyword{multivariate}

\seealso{\code{\link{pls.regression}}, \code{\link{TFA.estimate}},
\code{\link{pls.lda.cv}}.}

\examples{
# load plsgenomics library
library(plsgenomics)

# load Ecoli data
data(Ecoli)

# determine the best number of components for PLS regression using the cross-validation approach
# choose the best number from 1,2,3,4
pls.regression.cv(Xtrain=Ecoli$CONNECdata,Ytrain=Ecoli$GEdata,ncomp=4,nruncv=20)
# choose the best number from 2,3
pls.regression.cv(Xtrain=Ecoli$CONNECdata,Ytrain=Ecoli$GEdata,ncomp=c(2,3),nruncv=20)

}

