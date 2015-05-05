\name{rirls.spls}
\alias{rirls.spls}

\title{Classification by Ridge Iteratively Reweighted Least Squares followed by 
Adaptive Sparse PLS regression for binary response}
\usage{
rirls.spls(Xtrain, Ytrain, lambda.ridge, lambda.l1, ncomp, 
           Xtest=NULL, adapt=TRUE, maxIter=100, svd.decompose=TRUE)
}
\arguments{
	\item{Xtrain}{a (ntrain x p) data matrix of predictors. \code{Xtrain} must be a matrix. 
	Each row corresponds to an observation and each column to a predictor variable.}
	\item{Ytrain}{a ntrain vector of responses. \code{Ytrain} must be a vector or a one 
	column matrix. \code{Ytrain} is a \{0,1\}-valued vector and contains the response variable 
	for each observation.}
	\item{Xtest}{a (ntest x p) matrix containing the predictors for the test data set. 
	\code{Xtest} may also be a vector of length p (corresponding to only one test observation).
	If \code{Xtest} is not equal to NULL, then the prediction step is made for these new 
	predictor variables.}
	\item{lambda.ridge}{a positive real value. \code{lambda.ridge} is the ridge regularization 
	parameter for the RIRLS algorithm (see details).}
	\item{lambda.l1}{a positive real value, in [0,1]. \code{lambda.l1} is the sparse penalty 
	parameter for the dimension reduction step by sparse PLS (see details).}
	\item{ncomp}{a positive integer. \code{ncomp} is the number of PLS components. 
	If \code{ncomp}=0,then the Ridge regression is performed without dimension reduction.}
	\item{adapt}{a boolean value, indicating whether the sparse PLS selection step sould be 
	adaptive or nor.}
	\item{maxIter}{a positive integer. \code{maxIter} is the maximal number of iterations in the 
	Newton-Raphson parts in the RIRLS algorithm (see details).}
	\item{svd.decompose}{a boolean parameter. \code{svd.decompose} indicates wether or not should 
	the design matrix X be decomposed by SVD (singular values decomposition) for the RIRLS step 
	(see details).}
}
\description{
	The function \code{rirls.spls} performs compression, variable selection and classification (with possible prediction) 
	using Durif et al. (2015) RIRLS-SPLS algorithm.
  
}

\details{
	The columns of the data matrices \code{Xtrain} and \code{Xtest} may not be standardized, 
	since standardizing is performed by the function \code{rirls.spls} as a preliminary step
	before the algorithm is run. 
	
	The procedure described in Durif et al. (2015) is used to determine
	latent components to be used for classification and when \code{Xtest} 
	is not equal to NULL, the procedure predicts the labels for these new 
	predictor variables.  
 
}

\value{
	A list with the following components:
	\item{Coefficients}{the (p+1) vector containing the coefficients of the 
	design matrix and intercept in the logistic model explaining the response Y.}
	\item{hatY}{the (ntrain) vector containing the estimated reponse value on the train 
	set of predictors Xtrain.}
	\item{hatYtest}{the (ntest) vector containing the predicted labels for the observations from 
	\code{Xtest} if non null.}
	\item{DeletedCol}{the vector containing the column number of \code{Xtrain} when the 
	variance of the corresponding predictor variable is null. Otherwise \code{DeletedCol}=NULL}
	\item{A}{the active set of predictors selected by the procedures. \code{A} is a subset of 1:p}
	\item{converged}{a \{0,1\} value indicating whether the IRLS algorithm converged 
	in less than \code{maxIter} iterations or not.}
	\item{X.score}{a (n x ncomp) matrix being the observations coordinates or scores in the 
	new component basis produced by the compression step (sparse PLS). Each column t.k of 
	\code{X.score} is a new component.}
	\item{X.weight}{a (p x ncomp) matrix being the coefficients of predictors in each components 
	produced by sparse PLS. Each column w.k of \code{X.weight} verifies t.k = Xtrain x w.k (as a 
	matrix product).}
	\item{Xtrain}{the design matrix.}
	\item{sXtrain}{the scaled design matrix.}
	\item{Ytrain}{the response observations.}
	\item{sPseudoVar}{the scaled pseudo-response as produced by the RIRLS-algorithm and then 
	being scaled.}
	\item{lambda.ridge}{the ridge hyper-parameter used to fit the model.}
	\item{lambda.l1}{the sparse hyper-parameter used to fit the model.}
	\item{ncomp}{the number of components used to fit the model.}
	\item{V}{the (ntrain x ntrain) matrix used to weight the metric in the sparse PLS step. 
	\code{V} is the inverse of the covariance matrix of the pseudo-response produced by the 
	RIRLS step.}
	\item{proba.test}{the (ntest) vector of estimated probabilities for the observations in 
	\code{Xtest}, used to predict the \code{hatYtest} labels.}
}

\references{
 
	G. Durif, F. Picard, S. Lambert-Lacroix (2015). Adaptive sparse PLS for logistic regression, 
	(in prep), available on (\url{http://arxiv.org/abs/1502.05933}). 

}

\author{
	Ghislain Durif (\url{http://lbbe.univ-lyon1.fr/-Durif-Ghislain-.html}). 
     
     Adapted in part from rpls code by S. Lambert-Lacroix 
     (function available in this package).

}
\seealso{\code{\link{rirls.spls.tune}}.}

\examples{
### load plsgenomics library
library(plsgenomics)

### generating data
n <- 50
p <- 100
sample1 <- sample.bin(n=n, p=p, kstar=20, lstar=2, beta.min=0.25, beta.max=0.75, 
					mean.H=0.2, sigma.H=10, sigma.F=5)

X <- sample1$X
Y <- sample1$Y

### splitting between learning and testing set
index.train <- sort(sample(1:n, size=round(0.7*n)))
index.test <- (1:n)[-index.train]

Xtrain <- X[index.train,]
Ytrain <- Y[index.train,]

Xtest <- X[index.test,]
Ytest <- Y[index.test,]

### fitting the model, and predicting new observations
model1 <- rirls.spls(Xtrain=Xtrain, Ytrain=Ytrain, lambda.ridge=2, lambda.l1=0.5, ncomp=2, 
					Xtest=Xtest, adapt=TRUE, maxIter=100, svd.decompose=TRUE)
str(model1)

### prediction error rate
sum(model1$hatYtest!=Ytest) / length(index.test)

}

