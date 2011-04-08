\name{rpls.cv}
\alias{rpls.cv}

\title{Determination of the ridge regularization parameter and the number of PLS 
components to be used for classification with RPLS for binary data}
\usage{
rpls.cv(Ytrain,Xtrain,LambdaRange,ncompMax,NbIterMax=50)
}
\arguments{
  \item{Xtrain}{a (ntrain x p) data matrix of predictors. \code{Xtrain} must be a matrix. 
  Each row corresponds to an observation and each column to a predictor variable.}
  \item{Ytrain}{a ntrain vector of responses. \code{Ytrain} must be a vector. 
  \code{Ytrain} is a \{1,2\}-valued vector and contains the response variable for each
  observation.}
  \item{LambdaRange}{the vector of positive real value from which the best ridge regularization parameter
  has to be chosen by cross-validation.}
  \item{ncompMax}{a positive integer. the best number of components is chosen from  1,...,\code{ncompMax}. 
  If \code{ncompMax}=0,then the Ridge regression is performed without reduction 
  dimension. }
  \item{NbIterMax}{a positive integer. \code{NbIterMax} is the maximal number of iterations in the 
  Newton-Rapson parts.}
}
\description{
 The function \code{rpls.cv} determines the best ridge regularization parameter and the best number of PLS
 components to be used for classification for Fort and Lambert-Lacroix (2005) RPLS algorithm.

}

\details{
  A cross-validation procedure is used to determine the best ridge regularization parameter and 
  number of PLS components to be used for classification with RPLS for binary data 
  (for categorical data see \code{\link{mrpls}} and \code{\link{mrpls.cv}}).
  At each cross-validation run, \code{Xtrain} is split into a pseudo training
  set (ntrain-1 samples) and a pseudo test set (1 sample) and the classification error rate is determined for each
  value of ridge regularization parameter and number of components. Finally, the function 
  \code{mrpls.cv} returns the values of the ridge regularization parameter and 
  bandwidth for which the mean classification error rate is minimal. 
 
}


\value{
  A list with the following components:
  \item{Lambda}{the optimal regularization parameter.} 
  \item{ncomp}{the optimal number of PLS components.}
  }

\references{
 
 G. Fort and S. Lambert-Lacroix (2005). Classification using Partial Least Squares with 
 Penalized Logistic Regression, Bioinformatics, vol 21,  n 8, 1104-1111. 

}

\author{
  Sophie Lambert-Lacroix (\url{http://www-lmc.imag.fr/lmc-sms/Sophie.Lambert}). 

}
\seealso{\code{\link{rpls}}, \code{\link{mrpls}}, \code{\link{mrpls.cv}}.}

\examples{
# load plsgenomics library
# load plsgenomics library
library(plsgenomics)

# load Colon data
data(Colon)
IndexLearn <- c(sample(which(Colon$Y==2),12),sample(which(Colon$Y==1),8))

# preprocess data
res <- preprocess(Xtrain= Colon$X[IndexLearn,], Xtest=Colon$X[-IndexLearn,],Threshold = c(100,16000),Filtering=c(5,500),log10.scale=TRUE,row.stand=TRUE)
# the results are given in res$pXtrain and res$pXtest

# Determine optimum ncomp and lambda
nl <- rpls.cv(Ytrain=Colon$Y[IndexLearn],Xtrain=res$pXtrain,LambdaRange=c(0.1,1),ncompMax=3)

# perform prediction by RPLS
resrpls <- rpls(Ytrain=Colon$Y[IndexLearn],Xtrain=res$pXtrain,Lambda=nl$Lambda,ncomp=nl$ncomp,Xtest=res$pXtest)
sum(resrpls$Ytest!=Colon$Y[-IndexLearn])

}

