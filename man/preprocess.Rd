\name{preprocess}
\alias{preprocess}

\title{preprocess for microarray data}
\usage{
preprocess(Xtrain, Xtest=NULL,Threshold=c(100,16000),Filtering=c(5,500),
		log10.scale=TRUE,row.stand=TRUE)
}
\arguments{
\item{Xtrain}{a (ntrain x p) data matrix of predictors. \code{Xtrain} must be a matrix. 
  Each row corresponds to an observation and each column to a predictor variable.}
\item{Xtest}{a (ntest x p) matrix containing the predictors for the test data
  set. \code{Xtest} may also be a vector of length p (corresponding to only one 
  test observation).}
\item{Threshold}{ a vector of length 2 containing the values (threshmin,threshmax) for 
thresholding data in preprocess. Data is thresholded to value threshmin and ceiled to value 
threshmax. If \code{Threshold} is NULL then no thresholding is done. By default, if the value 
given for \code{Threshold} is not valid, no thresholding is done.}
\item{Filtering}{ a vector of length 2 containing the values (FiltMin,FiltMax) for filtering genes 
in preprocess. Genes with max/min$<= FiltMin$ and (max-min)$<= FiltMax$ are excluded. 
If \code{Filtering} is NULL then no thresholding is done. By default, if the value given for 
\code{Filtering} is not valid, no filtering is done.}
\item{log10.scale}{a logical value equal to TRUE if a log10-transformation has to be done.}
\item{row.stand}{a logical value equal to TRUE if a standardisation in row has to be done.}
}

\description{
 The function \code{preprocess} performs a preprocessing of microarray data.
  }
  
\details{
  The pre-processing steps recommended by Dudoit et al. (2002) are performed.
  The default values are those adapted for Colon data.
 
}  
\value{ A list with the following components:
  \item{pXtrain}{the (ntrain x p') matrix containing the preprocessed train data.}
  \item{pXtest}{the (ntest x p') matrix containing the preprocessed test data.}
}

\references{
Dudoit, S. and Fridlyand, J. and Speed, T. (2002). Comparison of discrimination methods 
for the classification of tumors using gene expression data, 
Journal of the American Statistical Association, 97, 77--87.
}

\author{
  Sophie Lambert-Lacroix 
  (\url{http://membres-timc.imag.fr/Sophie.Lambert/}) 
   and 
  Julie Peyre (\url{https://membres-ljk.imag.fr/Julie.Peyre/}).
}

\examples{
# load plsgenomics library
library(plsgenomics)

# load Colon data
data(Colon)
IndexLearn <- c(sample(which(Colon$Y==2),27),sample(which(Colon$Y==1),14))

Xtrain <- Colon$X[IndexLearn,]
Ytrain <- Colon$Y[IndexLearn]
Xtest <- Colon$X[-IndexLearn,]

# preprocess data
resP <- preprocess(Xtrain= Xtrain, Xtest=Xtest,Threshold = c(100,16000),Filtering=c(5,500),
				log10.scale=TRUE,row.stand=TRUE)

# how many genes after preprocess ?
dim(resP$pXtrain)[2]
}

