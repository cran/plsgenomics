\name{gsim.cv}
\alias{gsim.cv}

\title{Determination of the ridge regularization parameter and the bandwidth to be used for
classification with GSIM for binary data}
\usage{
gsim.cv(Xtrain, Ytrain,LambdaRange,hARange,hB=NULL, NbIterMax=50)
}
\arguments{
  \item{Xtrain}{a (ntrain x p) data matrix of predictors. \code{Xtrain} must be a matrix. 
  Each row corresponds to an observation and each column to a predictor variable.}
  \item{Ytrain}{a ntrain vector of responses. \code{Ytrain} must be a vector. 
  \code{Ytrain} is a \{1,2\}-valued vector and contains the response variable for each
  observation.}
  \item{LambdaRange}{the vector of positive real value from which the best ridge regularization 
  parameter has to be chosen by cross-validation.}
  \item{hARange}{the vector of strictly positive real value from which the best bandwidth
  has to be chosen by cross-validation for GSIM step A.}
  \item{hB}{a strictly positive real value. \code{hB} is the bandwidth for 
  GSIM step B. if \code{hB} is equal to NULL, then hB value is chosen using a
  plug-in method.}
  \item{NbIterMax}{a positive integer. \code{NbIterMax} is the maximal number of
  iterations in the Newton-Rapson parts.}
}
\description{
 The function \code{gsim.cv} determines the best ridge regularization parameter and bandwidth to be
 used for classification with GSIM as described in Lambert-Lacroix and Peyre (2005).
  
}

\details{
  The cross-validation procedure described in Lambert-Lacroix and Peyre (2005)
  is used to determine the best ridge regularization parameter and bandwidth to be
  used for classification with GSIM for binary data (for categorical data see 
  \code{\link{mgsim}} and \code{\link{mgsim.cv}}).
  At each cross-validation run, \code{Xtrain} is split into a pseudo training
  set (ntrain - 1 samples) and a pseudo test set (1 sample) and the classification error rate is 
  determined for each value of ridge regularization parameter and bandwidth. Finally, the function 
  \code{gsim.cv} returns the values of the ridge regularization parameter and 
  bandwidth for which the mean classification error rate is minimal. 
  
 
}

\value{
  A list with the following components:
  \item{Lambda}{the optimal regularization parameter.} 
  \item{hA}{the optimal bandwidth parameter.}
  }

\references{
 
 S. Lambert-Lacroix,  J. Peyre . (2006) Local likelyhood regression in  generalized linear 
 single-index models with applications to microarrays data. Computational Statistics and 
 Data Analysis, vol 51, n 3, 2091-2113.

}

\author{
  Sophie Lambert-Lacroix 
  (\url{http://membres-timc.imag.fr/Sophie.Lambert/}) 
   and 
  Julie Peyre (\url{http://www-lmc.imag.fr/lmc-sms/Julie.Peyre/}).
}

\seealso{\code{\link{mgsim}}, \code{\link{gsim}}, \code{\link{gsim.cv}}.}

\examples{
\dontrun{
## between 5~15 seconds
# load plsgenomics library
library(plsgenomics)

# load Colon data
data(Colon)
IndexLearn <- c(sample(which(Colon$Y==2),12),sample(which(Colon$Y==1),8))

Xtrain <- Colon$X[IndexLearn,]
Ytrain <- Colon$Y[IndexLearn]
Xtest <- Colon$X[-IndexLearn,]

# preprocess data
resP <- preprocess(Xtrain= Xtrain, Xtest=Xtest,Threshold = c(100,16000),Filtering=c(5,500),
				log10.scale=TRUE,row.stand=TRUE)

# Determine optimum h and lambda
hl <- gsim.cv(Xtrain=resP$pXtrain,Ytrain=Ytrain,hARange=c(7,20),LambdaRange=c(0.1,1),hB=NULL)

# perform prediction by GSIM  
res <- gsim(Xtrain=resP$pXtrain,Ytrain=Ytrain,Xtest=resP$pXtest,Lambda=hl$Lambda,hA=hl$hA,hB=NULL)
res$Cvg
sum(res$Ytest!=Colon$Y[-IndexLearn])
}
}

