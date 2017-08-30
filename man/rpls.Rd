\name{rpls}
\alias{rpls}
\alias{logit.pls}

\title{Ridge Partial Least Square for binary data}
\usage{
rpls(Ytrain,Xtrain,Lambda,ncomp,Xtest=NULL,NbIterMax=50)
}
\arguments{
  \item{Xtrain}{a (ntrain x p) data matrix of predictors. \code{Xtrain} must be a matrix. 
  Each row corresponds to an observation and each column to a predictor variable.}
  \item{Ytrain}{a ntrain vector of responses. \code{Ytrain} must be a vector. 
  \code{Ytrain} is a \{0,1\}-valued vector and contains the response variable for each
  observation.}
  \item{Xtest}{a (ntest x p) matrix containing the predictors for the test data
  set. \code{Xtest} may also be a vector of length p (corresponding to only one 
  test observation).If \code{Xtest} is not equal to NULL, then the prediction 
  step is made for these new predictor variables.}
  \item{Lambda}{a positive real value. \code{Lambda} is the ridge regularization parameter.}
  \item{ncomp}{a positive integer. \code{ncomp} is the number of PLS components. 
  If \code{ncomp}=0,then the Ridge regression is performed without reduction 
  dimension. }
  \item{NbIterMax}{a positive integer. \code{NbIterMax} is the maximal number of iterations in the 
  Newton-Rapson parts.}
}
\description{
 The function \code{mrpls} performs prediction using Fort and Lambert-Lacroix (2005) RPLS 
 algorithm.
  
}

\details{
  The columns of the data matrices \code{Xtrain} and \code{Xtest} may not be standardized, 
  since standardizing is performed by the function \code{rpls} as a preliminary step
  before the algorithm is run. 

  The procedure described in Fort and Lambert-Lacroix (2005) is used to determine
  latent components to be used for classification and when \code{Xtest} 
  is not equal to NULL, the procedure predicts the labels for these new 
  predictor variables.  
 
}

\value{
  A list with the following components:
  \item{Coefficients}{the (p+1) vector containing the coefficients weighting the 
  design matrix.}
  \item{hatY}{the ntrain vector containing the estimated \{0,1\}-valued labels for the 
  observations from \code{Xtrain}.}
  \item{hatYtest}{the ntest vector containing the predicted \{0,1\}-valued labels for the 
  observations from \code{Xtest}.}
  \item{proba}{the ntrain vector containing the estimated probabilities for the 
  observations from \code{Xtrain}.}
  \item{proba.test}{the ntest vector containing the predicted probabilities for the 
  observations from \code{Xtest}.}
  \item{DeletedCol}{the vector containing the column number of \code{Xtrain} when the 
  variance of the corresponding predictor variable is null. Otherwise \code{DeletedCol}=NULL}
  \item{hatYtest_k}{If \code{ncomp} is greater than 1, \code{hatYtest_k} is a \{0,1\}-valued matrix of 
  size ntest x ncomp in such a way that the kth column corresponds to the predicted label 
  obtained with k PLS components.}

  }

\references{
 
 G. Fort and S. Lambert-Lacroix (2005). Classification using Partial Least Squares with 
 Penalized Logistic Regression, Bioinformatics, vol 21,  n 8, 1104-1111. 

}

\author{
  Sophie Lambert-Lacroix 
  (\url{http://membres-timc.imag.fr/Sophie.Lambert/}). 

}
\seealso{\code{\link{rpls.cv}}, \code{\link{mrpls}}, \code{\link{mrpls.cv}}.}

\examples{
# load plsgenomics library
library(plsgenomics)

# load Colon data
data(Colon)
IndexLearn <- c(sample(which(Colon$Y==2),12),sample(which(Colon$Y==1),8))

# preprocess data
res <- preprocess(Xtrain= Colon$X[IndexLearn,], Xtest=Colon$X[-IndexLearn,],
                    Threshold = c(100,16000),Filtering=c(5,500),
                    log10.scale=TRUE,row.stand=TRUE)
# the results are given in res$pXtrain and res$pXtest

# perform prediction by RPLS
resrpls <- rpls(Ytrain=Colon$Y[IndexLearn]-1,Xtrain=res$pXtrain,Lambda=0.6,ncomp=1,Xtest=res$pXtest)
resrpls$hatY
sum(resrpls$Ytest!=Colon$Y[-IndexLearn])

# prediction for another sample
Xnew <- res$pXtest[1,]
# Compute the linear predictor for each classes expect class 0
eta <- c(1,Xnew) \%*\% resrpls$Coefficients
Ypred <- which.max(c(0,eta))
Ypred+1



}

