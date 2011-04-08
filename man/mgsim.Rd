\name{mgsim}
\alias{mgsim}

\title{GSIM for categorical data}
\usage{
mgsim(Ytrain,Xtrain,Lambda,h,Xtest=NULL,NbIterMax=50)
}
\arguments{
  \item{Xtrain}{a (ntrain x p) data matrix of predictors. \code{Xtrain} must be a matrix. 
  Each row corresponds to an observation and each column to a predictor variable.}
  \item{Ytrain}{a ntrain vector of responses. \code{Ytrain} must be a vector. 
  \code{Ytrain} is a \{1,...,c+1\}-valued vector and contains the response variable for each
  observation. c+1 is the number of classes.}
  \item{Xtest}{a (ntest x p) matrix containing the predictors for the test data
  set. \code{Xtest} may also be a vector of length p (corresponding to only one 
  test observation). If \code{Xtest} is not equal to NULL, then the prediction 
  step is made for these new predictor variables.}
  \item{Lambda}{a positive real value. \code{Lambda} is the ridge regularization parameter.}
  \item{h}{a strictly positive real value. \code{h} is the bandwidth for GSIM step A.}
  \item{NbIterMax}{a positive integer. \code{NbIterMax} is the maximal number of iterations in the 
  Newton-Rapson parts.}
}
\description{
 The function \code{mgsim} performs prediction using Lambert-Lacroix and Peyre's MGSIM algorithm.
  
}

\details{
  The columns of the data matrices \code{Xtrain} and \code{Xtest} may not be standardized, 
  since standardizing is performed by the function \code{mgsim} as a preliminary step
  before the algorithm is run. 

  The procedure described in Lambert-Lacroix and Peyre (2005) is used to estimate 
  the c projection directions and the coefficients of the parametric fit obtained 
  after projecting predictor variables onto the estimated directions. When \code{Xtest} 
  is not equal to NULL, the procedure predicts the labels for these new predictor variables.  
 
}

\value{
  A list with the following components:
  \item{Ytest}{the ntest vector containing the predicted labels for the observations from 
  \code{Xtest}.} 
  \item{beta}{the (p x c) matrix containing the c estimated projection directions.}
  \item{Coefficients}{the (2 x c) matrix containing the coefficients of the parametric fit obtained 
  after projecting predictor variables onto these estimated directions. }
  \item{DeletedCol}{the vector containing the column number of \code{Xtrain} when the 
  variance of the corresponding predictor variable is null. Otherwise \code{DeletedCol}=NULL}
  \item{Cvg}{the 0-1 value indicating convergence of the algorithm (1 for convergence, 0 otherwise).}
  }

\references{
 
S. Lambert-Lacroix,  J. Peyre . (2006) Local likelyhood regression in  generalized linear single-index models with applications to microarrays data. Computational Statistics and Data Analysis, vol 51, n 3, 2091-2113. 

}

\author{
  Sophie Lambert-Lacroix (\url{http://www-lmc.imag.fr/lmc-sms/Sophie.Lambert}) 
   and 
  Julie Peyre (\url{http://www-lmc.imag.fr/lmc-sms/Julie.Peyre/}).

}
\seealso{\code{\link{mgsim.cv}}, \code{\link{gsim}}, \code{\link{gsim.cv}}.}

\examples{
# load plsgenomics library
library(plsgenomics)

# load SRBCT data
data(SRBCT)
IndexLearn <- c(sample(which(SRBCT$Y==1),10),sample(which(SRBCT$Y==2),4),sample(which(SRBCT$Y==3),7),sample(which(SRBCT$Y==4),9))

# perform prediction by MGSIM
res <- mgsim(Ytrain=SRBCT$Y[IndexLearn],Xtrain=SRBCT$X[IndexLearn,],Lambda=0.001,h=19,Xtest=SRBCT$X[-IndexLearn,])
res$Cvg
sum(res$Ytest!=SRBCT$Y[-IndexLearn])

# prediction for another sample
Xnew <- SRBCT$X[83,]
# projection of Xnew onto the c estimated direction
Xproj <- Xnew \%*\% res$beta
# Compute the linear predictor for each classes expect class 1
eta <- diag(cbind(rep(1,3),t(Xproj)) \%*\% res$Coefficients)
Ypred <- which.max(c(0,eta))
Ypred
SRBCT$Y[83]

}


