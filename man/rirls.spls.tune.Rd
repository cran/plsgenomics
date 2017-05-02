\name{rirls.spls.tune}
\alias{rirls.spls.tune}

\title{Tuning parameters (ncomp, lambda.l1, lambda.ridge) for Ridge Iteratively Reweighted Least 
Squares followed by Adaptive Sparse PLS regression for binary response, by K-fold cross-validation}
\usage{
	rirls.spls.tune(X, Y, lambda.ridge.range, lambda.l1.range, ncomp.range, 
                    adapt=TRUE, maxIter=100, svd.decompose=TRUE, return.grid=FALSE, 
                    ncores=1, nfolds=10)
}
\arguments{
	\item{X}{a (n x p) data matrix of predictors. \code{X} must be a matrix. Each row corresponds 
	to an observation and each column to a predictor variable.}
	\item{Y}{a ntrain vector of responses. \code{Y} must be a vector or a one column matrix. 
    	\code{Y} is a \{0,1\}-valued vector and contains the response variable for each observation.}
	\item{lambda.ridge.range}{a vector of positive real values. \code{lambda.ridge} is the ridge 
	regularization parameter for the RIRLS algorithm (see details), the optimal values will be 
	chosen among \code{lambda.ridge.range}.}
	\item{lambda.l1.range}{a vecor of positive real values, in [0,1]. \code{lambda.l1} is the 
	sparse penalty parameter for the dimension reduction step by sparse PLS (see details), the 
	optimal values will be chosen among \code{lambda.l1.range}.}
	\item{ncomp.range}{a vector of positive integers. \code{ncomp} is the number of PLS 
	components. If \code{ncomp}=0,then the Ridge regression is performed without dimension 
	reduction. The optimal values will be chosen among \code{ncomp.range}.}
	\item{adapt}{a boolean value, indicating whether the sparse PLS selection step sould be 
	adaptive or nor.}
	\item{maxIter}{a positive integer. \code{maxIter} is the maximal number of iterations in the 
    	Newton-Raphson parts in the RIRLS algorithm (see details).}
	\item{svd.decompose}{a boolean parameter. \code{svd.decompose} indicates wether or not should 
	the design matrix X be decomposed by SVD (singular values decomposition) for the RIRLS step 
	(see details).}
    	\item{return.grid}{a boolean values indicating whether the grid of hyper-parameters values 
    	with corresponding mean prediction error rate over the folds should be returned or not.}
    	\item{ncores}{a positve integer, indicating if the cross-validation procedure should be 
    	parallelized over the folds (ncores > nfolds would lead to the generation of unused child 
    	process). If ncores>1, the procedure generates ncores child process over the cpu corresponding 
    	number of cpu cores (see details).}
    	\item{nfolds}{a positive integer indicating the number of folds in K-folds cross-validation 
    	procedure, nfolds=n corresponds to leave-one-out cross-validation.}
}
\description{
	The function \code{rirls.spls.tune} tuns the hyper-parameter values used in the 
	\code{rirls.spls} procedure, by minimizing the prediction error rate over the hyper-parameter 
	grid, using Durif et al. (2015) RIRLS-SPLS algorithm.
	
}

\details{
	The columns of the data matrices \code{X} may not be standardized, 
	since standardizing is performed by the function \code{rirls.spls} as a preliminary step
	before the algorithm is run. 
	
	The procedure is described in Durif et al. (2015). The K-fold cross-validation can be 
	summarize as follow: the train set is partitioned into K folds, for each value of hyper-
	parameters the model is fit K times, using each fold to compute the prediction error rate, and 
	fitting the model on the remaining observations. The cross-validation procedure returns the 
	optimal hyper-parameters values, meaning the one that minimize the prediction error rate 
	averaged over all the folds.
	
	This procedures uses the \code{mclapply} from the \code{parallel} package, available on 
	GNU/Linux and MacOS. Users of Microsoft Windows can refer to the README file in the source to 
	be able to use a mclapply type function.
	
}

\value{
	A list with the following components:
	\item{lambda.ridge.opt}{the optimal value in \code{lambda.ridge.range}.}
	\item{lambda.l1.opt}{the optimal value in \code{lambda.l1.range}.}
	\item{ncomp.opt}{the optimal value in \code{ncomp.range}.}
	\item{conv.per}{the overall percentage of models that converge during the cross-validation 
	procedure.}
	\item{cv.grid}{the grid of hyper-parameters and corresponding prediction error rate over the 
	nfolds.
	\code{cv.grid} is NULL if \code{return.grid} is set to FALSE.}
}

\references{
	
	G. Durif, F. Picard, S. Lambert-Lacroix (2015). Adaptive sparse PLS for logistic regression, 
	(in prep), available on (\url{http://arxiv.org/abs/1502.05933}). 
	
}

\author{
	Ghislain Durif (\url{http://lbbe.univ-lyon1.fr/-Durif-Ghislain-.html}). 
	
}
\seealso{\code{\link{rirls.spls}}.}

\examples{
\dontrun{
## between 5~15 seconds
### load plsgenomics library
library(plsgenomics)

### generating data
n <- 50
p <- 100
sample1 <- sample.bin(n=n, p=p, kstar=20, lstar=2, beta.min=0.25, beta.max=0.75, mean.H=0.2, 
                    sigma.H=10, sigma.F=5)

X <- sample1$X
Y <- sample1$Y

### hyper-parameters values to test
lambda.l1.range <- seq(0.05,0.95,by=0.3) # between 0 and 1
ncomp.range <- 1:2

# log-linear range between 0.01 a,d 1000 for lambda.ridge.range
logspace <- function( d1, d2, n) exp(log(10)*seq(d1, d2, length.out=n)) 
lambda.ridge.range <- signif(logspace(d1 <- -2, d2 <- 3, n=6), digits=3)

### tuning the hyper-parameters
cv1 <- rirls.spls.tune(X=X, Y=Y, lambda.ridge.range=lambda.ridge.range, 
                         lambda.l1.range=lambda.l1.range, ncomp.range=ncomp.range, 
                         adapt=TRUE, maxIter=100, svd.decompose=TRUE, 
                         return.grid=TRUE, ncores=1, nfolds=10)
str(cv1)
	
}
}
