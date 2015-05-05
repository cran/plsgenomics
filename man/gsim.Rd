\name{gsim}
\alias{gsim}

\title{GSIM for binary data}
\usage{
gsim(Xtrain, Ytrain, Xtest=NULL, Lambda, hA, hB=NULL, NbIterMax=50)
}
\arguments{
  \item{Xtrain}{a (ntrain x p) data matrix of predictors. \code{Xtrain} must be a matrix. 
  Each row corresponds to an observation and each column to a predictor variable.}
  \item{Ytrain}{a ntrain vector of responses. \code{Ytrain} must be a vector. 
  \code{Ytrain} is a \{1,2\}-valued vector and contains the response variable for each
  observation.}
  \item{Xtest}{a (ntest x p) matrix containing the predictors for the test data
  set. \code{Xtest} may also be a vector of length p (corresponding to only one
  test observation). If \code{Xtest} is not equal to NULL, then the prediction 
  step is made for these new predictor variables.}
  \item{Lambda}{a positive real value. \code{Lambda} is the ridge regularization parameter.}
  \item{hA}{a strictly positive real value. \code{hA} is the bandwidth for GSIM step A.}
  \item{hB}{a strictly positive real value. \code{hB} is the bandwidth for 
  GSIM step B. if \code{hB} is equal to NULL, then hB value is chosen using a
  plug-in method.}
  \item{NbIterMax}{a positive integer. \code{NbIterMax} is the maximal number of
  iterations in the Newton-Rapson parts.}
}
\description{
 The function \code{gsim} performs prediction using Lambert-Lacroix and Peyre's GSIM algorithm.
  }

\details{
  The columns of the data matrices \code{Xtrain} and \code{Xtest} may not be standardized, 
  since standardizing is performed by the function \code{gsim} as a preliminary step
  before the algorithm is run. 

  The procedure described in Lambert-Lacroix and Peyre (2005) is used to estimate 
  the projection direction beta. When \code{Xtest} 
  is not equal to NULL, the procedure predicts the labels for these new predictor variables.  
 
}

\value{
  A list with the following components:
  \item{Ytest}{the ntest vector containing the predicted labels for the observations from 
  \code{Xtest}.} 
  \item{beta}{the p vector giving the projection direction estimated.}
  \item{hB}{the value of hB used in step B of GSIM (value given by the user or
  estimated by plug-in if the argument value was equal to NULL)}
  \item{DeletedCol}{the vector containing the column number of \code{Xtrain} when the 
  variance of the corresponding predictor variable is null. Otherwise \code{DeletedCol}=NULL}
  \item{Cvg}{the 0-1 value indicating convergence of the algorithm (1 for convergence, 
  0 otherwise).}
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

\seealso{\code{\link{gsim.cv}}, \code{\link{mgsim}}, \code{\link{mgsim.cv}}.}

\examples{
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

# perform prediction by GSIM
res <- gsim(Xtrain=resP$pXtrain,Ytrain= Ytrain,Xtest=resP$pXtest,Lambda=10,hA=50,hB=NULL)
   
res$Cvg
sum(res$Ytest!=Colon$Y[-IndexLearn])

}

