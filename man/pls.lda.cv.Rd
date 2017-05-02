\name{pls.lda.cv}
\alias{pls.lda.cv}

\title{Determination of the number of latent components to be used for
classification with PLS and LDA}
\usage{
pls.lda.cv(Xtrain, Ytrain,  ncomp, nruncv=20, alpha=2/3, priors=NULL)
}
\arguments{
  \item{Xtrain}{a (ntrain x p) data matrix containing the predictors for the training data set.
   \code{Xtrain} may be a matrix or a data frame. Each row is an observation and each column is a 
   predictor variable.}
  \item{Ytrain}{a vector of length ntrain giving the classes of the ntrain
  observations. The classes must be coded as 1,...,K (K>=2).}
  \item{ncomp}{the vector of integers from which the best number of latent 
  components has to be chosen by cross-validation. If \code{ncomp} is of length
  1, the best number of components is chosen from  1,...,\code{ncomp}.}  \item{nruncv}{the number 
  of cross-validation iterations to be performed for the choice of the number of latent components.}
  \item{alpha}{the proportion of observations to be included in the training set at
  each cross-validation iteration.}
  \item{priors}{The class priors to be used for linear discriminant analysis. If
  unspecified, the class proportions in the training set are used. }
   
}
\description{
 The function \code{pls.lda.cv} determines the best number of latent components to be used for 
 classification with PLS dimension reduction and linear discriminant analysis as described in
 Boulesteix (2004). 
  
}

\details{
  The cross-validation procedure described in Boulesteix (2004) is used to
  determine the best number of latent components to be used for classification.
  At each cross-validation run, \code{Xtrain} is split into a pseudo training
  set and a pseudo test set and the classification error rate is determined for each
  number of latent components. Finally, the function \code{pls.lda.cv} returns
  the number of latent components for which the mean classification rate over
  the \code{nrun} partitions is minimal. 
  
 
}

\value{
  
  The number of latent components to be used for classification.
  
  
}

\references{
 
 A. L. Boulesteix (2004). PLS dimension reduction for classification with microarray data,
 Statistical Applications in Genetics and Molecular Biology \bold{3}, Issue 1, Article 33.

 A. L. Boulesteix, K. Strimmer (2007). Partial least squares: a versatile tool for the analysis 
 of high-dimensional genomic data. Briefings in Bioinformatics 7:32-44.

 S. de Jong (1993). SIMPLS: an alternative approach to partial least squares 
 regression, Chemometrics Intell. Lab. Syst. \bold{18}, 251--263.

}

\author{
Anne-Laure Boulesteix 
(\url{http://www.ibe.med.uni-muenchen.de/organisation/mitarbeiter/020_professuren/boulesteix/eng.html}) 
   
 

}
\seealso{\code{\link{pls.lda}}, \code{\link{pls.regression.cv}}.}

\examples{
\dontrun{
## between 5~15 seconds
# load plsgenomics library
library(plsgenomics)

# load leukemia data
data(leukemia)

# Determine the best number of components to be used for classification using the 
# cross-validation procedure
# choose the best number from 2,3,4
pls.lda.cv(Xtrain=leukemia$X,Ytrain=leukemia$Y,ncomp=2:4,nruncv=20)
# choose the best number from 1,2,3
pls.lda.cv(Xtrain=leukemia$X,Ytrain=leukemia$Y,ncomp=3,nruncv=20)

}
}
\keyword{multivariate}
