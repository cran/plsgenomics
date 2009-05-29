\name{mgsim.cv}
\alias{mgsim.cv}

\title{Determination of the ridge regularization parameter and the bandwidth to be used for
classification with GSIM for categorical data}
\usage{
mgsim.cv(Ytrain,Xtrain,LambdaRange,hRange,NbIterMax=50)

}
\arguments{
  \item{Xtrain}{a (ntrain x p) data matrix of predictors. \code{Xtrain} must be a matrix. 
  Each row corresponds to an observation and each column to a predictor variable.}
  \item{Ytrain}{a ntrain vector of responses. \code{Ytrain} must be a vector. 
  \code{Ytrain} is a \{1,...,c+1\}-valued vector and contains the response variable for each
  observation. c+1 is the number of classes.}
  \item{LambdaRange}{the vector of positive real value from which the best ridge regularization parameter
  has to be chosen by cross-validation.}
  \item{hRange}{the vector of strictly positive real value from which the best bandwidth
  has to be chosen by cross-validation.}
  \item{NbIterMax}{a positive integer. \code{NbIterMax} is the maximal number of iterations in the 
  Newton-Rapson parts.}
}
\description{
 The function \code{mgsim.cv} determines the best ridge regularization parameter and bandwidth to be
 used for classification with MGSIM as described in Lambert-Lacroix and Peyre (2005).
  
}

\details{
  The cross-validation procedure described in Lambert-Lacroix and Peyre (2005)
  is used to determine the best ridge regularization parameter and bandwidth to be
  used for classification with GSIM for categorical data (for binary data see 
  \code{\link{gsim}} and \code{\link{gsim.cv}}).
  At each cross-validation run, \code{Xtrain} is split into a pseudo training
  set (ntrain-1 samples) and a pseudo test set (1 sample) and the 
  classification error rate is determined for each
  value of ridge regularization parameter and bandwidth. Finally, the function 
  \code{mgsim.cv} returns the values of the ridge regularization parameter and 
  bandwidth for which the mean classification error rate is minimal. 
  
 
}

\value{
  A list with the following components:
  \item{Lambda}{the optimal regularization parameter.} 
  \item{h}{the optimal bandwidth parameter.}
  }

\references{
 
S. Lambert-Lacroix,  J. Peyre . (2006) Local likelyhood regression in  generalized linear single-index models with applications to microarrays data. Computational Statistics and Data Analysis, vol 51, n 3, 2091-2113. 

}

\author{
  Sophie Lambert-Lacroix (\url{http://www-lmc.imag.fr/lmc-sms/Sophie.Lambert}) 
   and 
  Julie Peyre (\url{http://www-lmc.imag.fr/lmc-sms/Julie.Peyre/}).

}
\seealso{\code{\link{mgsim}}, \code{\link{gsim}}, \code{\link{gsim.cv}}.}

\examples{
# load plsgenomics library
library(plsgenomics)

# load SRBCT data
data(SRBCT)
IndexLearn <- c(sample(which(SRBCT$Y==1),10),sample(which(SRBCT$Y==2),4),sample(which(SRBCT$Y==3),7),sample(which(SRBCT$Y==4),9))

# Determine optimum h and lambda 
hl <- mgsim.cv(Ytrain=SRBCT$Y[IndexLearn],Xtrain=SRBCT$X[IndexLearn,],LambdaRange=c(0.1),hRange=c(7,20))

# perform prediction by MGSIM
res <- mgsim(Ytrain=SRBCT$Y[IndexLearn],Xtrain=SRBCT$X[IndexLearn,],Lambda=hl$Lambda,h=hl$h,Xtest=SRBCT$X[-IndexLearn,])
res$Cvg
sum(res$Ytest!=SRBCT$Y[-IndexLearn])

}
\keyword{}
