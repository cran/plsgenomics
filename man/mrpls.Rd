\name{mrpls}
\alias{mrpls}

\title{Ridge Partial Least Square for categorical data}
\usage{
mrpls(Ytrain,Xtrain,Lambda,ncomp,Xtest=NULL,NbIterMax=50)
}
\arguments{
  \item{Xtrain}{a (ntrain x p) data matrix of predictors. \code{Xtrain} must be a matrix. 
  Each row corresponds to an observation and each column to a predictor variable.}
  \item{Ytrain}{a ntrain vector of responses. \code{Ytrain} must be a vector. 
  \code{Ytrain} is a \{1,...,c+1\}-valued vector and contains the response variable for each
  observation. c+1 is the number of classes.}
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
 The function \code{mrpls} performs prediction using Fort et al. (2005) MRPLS algorithm.
  
}

\details{
  The columns of the data matrices \code{Xtrain} and \code{Xtest} may not be standardized, 
  since standardizing is performed by the function \code{mrpls} as a preliminary step
  before the algorithm is run. 

  The procedure described in Fort et al. (2005) is used to determine
  latent components to be used for classification and when \code{Xtest} 
  is not equal to NULL, the procedure predicts the labels for these new 
  predictor variables.  
 
}

\value{
  A list with the following components:
  \item{Ytest}{the ntest vector containing the predicted labels for the observations from 
  \code{Xtest}.} 
  \item{Coefficients}{the (p+1) x c matrix containing the coefficients weighting the block 
  design matrix.}
  \item{DeletedCol}{the vector containing the column number of \code{Xtrain} when the 
  variance of the corresponding predictor variable is null. Otherwise \code{DeletedCol}=NULL}
  \item{hatY}{If \code{ncomp} is greater than 1, \code{hatY} is a matrix of size ntest x ncomp 
  in such a way that the kth column corresponds to the predicted label obtained with k PLS components.}

  }

\references{
 
 G. Fort, S. Lambert-Lacroix and Julie Peyre (2005). R�duction de dimension dans les mod�les 
 lin�aires g�n�ralis�s : application � la classification supervis�e de donn�es issues des biopuces.
 Journal de la SFDS, tome 146, n1-2, 117-152. 

}

\author{
  Sophie Lambert-Lacroix (\url{http://www-lmc.imag.fr/lmc-sms/Sophie.Lambert}). 

}
\seealso{\code{\link{mrpls.cv}}, \code{\link{rpls}}, \code{\link{rpls.cv}}.}

\examples{
# load plsgenomics library
library(plsgenomics)

# load SRBCT data
data(SRBCT)
IndexLearn <- c(sample(which(SRBCT$Y==1),10),sample(which(SRBCT$Y==2),4),sample(which(SRBCT$Y==3),7),sample(which(SRBCT$Y==4),9))

# perform prediction by MRPLS
res <- mrpls(Ytrain=SRBCT$Y[IndexLearn],Xtrain=SRBCT$X[IndexLearn,],Lambda=0.001,ncomp=2,Xtest=SRBCT$X[-IndexLearn,])
sum(res$Ytest!=SRBCT$Y[-IndexLearn])

# prediction for another sample
Xnew <- SRBCT$X[83,]
# Compute the linear predictor for each classes expect class 1
eta <- diag(t(cbind(c(1,Xnew),c(1,Xnew),c(1,Xnew))) \%*\% res$Coefficients)
Ypred <- which.max(c(0,eta))
Ypred
SRBCT$Y[83]

}
\keyword{}
