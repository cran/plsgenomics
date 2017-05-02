\name{pls.lda}
\alias{pls.lda}

\title{Classification with PLS Dimension Reduction and Linear Discriminant
Analysis}
\usage{
pls.lda(Xtrain, Ytrain, Xtest=NULL, ncomp, nruncv=0, alpha=2/3, priors=NULL)
}
\arguments{
  \item{Xtrain}{a (ntrain x p) data matrix containing the predictors for the training data set.
   Xtrain may be a matrix or a data frame. Each row is an observation and each column is a 
   predictor variable.}
  \item{Ytrain}{a vector of length ntrain giving the classes of the ntrain
  observations. The classes must be coded as 1,...,K (K>=2).}
  \item{Xtest}{a (ntest x p) data matrix containing the predictors for the test
  data set.  \code{Xtest} may also be a
  vector of length p (corresponding to only one test observation). If
  \code{Xtest=NULL}, the training data set is considered as test data set as
  well.}
  \item{ncomp}{if \code{nruncv=0}, \code{ncomp} is the number of latent components to 
  be used for PLS dimension reduction. If
  \code{nruncv>0}, the cross-validation procedure described in Boulesteix (2004) is used
  to choose the best number of components from the vector of integers \code{ncomp} or from 
  1,...,\code{ncomp} if \code{ncomp} is of length 1.}
  \item{nruncv}{the number of cross-validation iterations to be performed for the choice of
  the number of latent components. If \code{nruncv=0}, cross-validation is not performed 
  and \code{ncomp} latent components are used.}
  \item{alpha}{the proportion of observations to be included in the training set at
  each cross-validation iteration.}
  \item{priors}{The class priors to be used for linear discriminant analysis. If
  unspecified, the class proportions in the training set are used. }
   
}
\description{
The function \code{pls.lda}  performs binary or multicategorical classification using the method 
described in Boulesteix (2004) which consists in PLS dimension reduction and linear
discriminant analysis applied on the PLS components. 
  
}

\details{
  The function \code{pls.lda} proceeds as follows to predict the class of the
  observations from the test data set. 
  First, the SIMPLS algorithm is run on \code{Xtrain} and \code{Ytrain} to
  determine the new PLS components based on the training observations only. 
  The new PLS components are then computed for the test
  data set. Classification is performed by applying classical linear 
  discriminant analysis (LDA) to the new components. Of course, the LDA
  classifier is built using the training observations only. 
  
 
}

\value{
  A list with the following components:
  \item{predclass}{the vector containing the predicted classes of the ntest observations from
  \code{Xtest}.}
  \item{ncomp}{the number of latent components used for classification.}
  \item{pls.out}{an object containing the results from the call of the \code{pls.regression}
  function (from the \code{plsgenomics} package).}
  \item{lda.out}{an object containing the results from the call of the \code{lda} 
  function (from the \code{MASS} package).}
  \item{pred.lda.out}{an object containing the results from the call of the \code{predict.lda} 
  function (from the \code{MASS} package).}
  
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
\seealso{\code{\link{pls.regression}}, \code{\link{variable.selection}},
\code{\link{pls.lda.cv}}.}

\examples{
# load plsgenomics library
library(plsgenomics)

# load leukemia data
data(leukemia)

# Classify observations 1,2,3 (test set) using observations 4 to 38 (training set), 
# with 2 PLS components
pls.lda(Xtrain=leukemia$X[-(1:3),],Ytrain=leukemia$Y[-(1:3)],Xtest=leukemia$X[1:3,],
     	ncomp=2,nruncv=0)

# Classify observations 1,2,3 (test set) using observations 4 to 38 (training set), 
# with the best number of components as determined by cross-validation
pls.lda(Xtrain=leukemia$X[-(1:3),],Ytrain=leukemia$Y[-(1:3)],Xtest=leukemia$X[1:3,],
		ncomp=1:4,nruncv=20)

}
\keyword{multivariate}
