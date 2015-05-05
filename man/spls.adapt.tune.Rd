\name{spls.adapt.tune}
\alias{spls.adapt.tune}

\title{Tuning parameters (ncomp, lambda.l1) for Adaptive Sparse PLS regression for continuous 
response, by K-fold cross-validation}
\usage{
spls.adapt.tune(X, Y, lambda.l1.range, ncomp.range, weight.mat=NULL, adapt=TRUE, 
                    center.X=TRUE, center.Y=TRUE, scale.X=TRUE, scale.Y=TRUE, 
                    weighted.center=FALSE, return.grid=FALSE, ncores=1, nfolds=10)
}
\arguments{
	\item{X}{a (n x p) data matrix of predictors. \code{Xtrain} must be a matrix. 
	Each row corresponds to an observation and each column to a predictor variable.}
	\item{Y}{a (n) vector of (continuous) responses. \code{Ytrain} must be a vector or a one 
	column matrix. and contains the response variable for each observation.}
	\item{lambda.l1.range}{a vecor of positive real values, in [0,1]. \code{lambda.l1} is the 
	sparse penalty parameter for the dimension reduction step by sparse PLS (see details), the 
	optimal values will be chosen among \code{lambda.l1.range}.}
	\item{ncomp.range}{a vector of positive integers. \code{ncomp} is the number of PLS 
	components. If \code{ncomp}=0,then the Ridge regression is performed without dimension 
	reduction. The optimal values will be chosen among \code{ncomp.range}.}
	\item{weight.mat}{a (n x n) matrix used to weight the l2 metric in observation spase 
	if necessary, especially the covariance inverse of the Ytrain observations in heteroskedastic 
	context. If NULL, the l2 metric is the standard one, corresponding to homoskedastic model.}
	\item{adapt}{a boolean value, indicating whether the sparse PLS selection step sould be 
	adaptive or not.}
	\item{center.X}{a boolean value indicating whether the design matrices Xtrain in train set 
	and Xtest in test set if non NULL should be centered or not}
	\item{scale.X}{a boolean value indicating whether the design matrices Xtrain in train set 
	and Xtest in test set if non NULL should be scaled or not, scale.X=TRUE implies center.X=TRUE}
	\item{center.Y}{a boolean value indicating whether the response Ytrain in train set should be 
	centered or not}
	\item{scale.Y}{a boolean value indicating whether the response Ytrain should be scaled or not, 
	scale.Y=TRUE implies center.Y=TRUE}
	\item{weighted.center}{a boolean value indicating whether should the centering take into 
	account the weighted l2 metric or not (if TRUE, it implies that weighted.mat is non NULL).}
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
	\code{spls.adapt} procedure, by minimizing the mean squared error of prediction over the 
	hyper-parameter grid, using Durif et al. (2015) adaptive SPLS algorithm.

}

\details{
	The columns of the data matrices \code{Xtrain} and \code{Xtest} may not be standardized, 
	since standardizing is can be performed by the function \code{spls.adapt.tune} as a 
	preliminary step before the algorithm is run. 
	
	The procedure is described in Durif et al. (2015). The K-fold cross-validation can be 
	summarize as follow: the train set is partitioned into K folds, for each value of hyper-
	parameters the model is fit K times, using each fold to compute the prediction error rate, and 
	fitting the model on the remaining observations. The cross-validation procedure returns the 
	optimal hyper-parameters values, meaning the one that minimize the mean squared error of 
	prediction averaged over all the folds.
	
	This procedures uses the \code{mclapply} from the \code{parallel} package, available on 
	GNU/Linux and MacOS. Users of Microsoft Windows can refer to the README file in the source to 
	be able to use a mclapply type function.
 
}

\value{
	A list with the following components:
	\item{lambda.l1.opt}{the optimal value in \code{lambda.l1.range}.}
	\item{ncomp.opt}{the optimal value in \code{ncomp.range}.}
	\item{cv.grid}{the grid of hyper-parameters and corresponding prediction error rate over 
	the nfolds. \code{cv.grid} is NULL if \code{return.grid} is set to FALSE.}
}

\references{
 
     G. Durif, F. Picard, S. Lambert-Lacroix (2015). Adaptive sparse PLS for logistic regression, 
     (in prep), available on (\url{http://arxiv.org/abs/1502.05933}).
     
     Chun, H., & Keles, S. (2010). Sparse partial least squares regression for simultaneous 
     dimension reduction and variable selection.  Journal of the Royal Statistical Society.
     Series B (Methodological), 72(1), 3-25. doi:10.1111/j.1467-9868.2009.00723.x
     
     Chung, D., & Keles, S. (2010). Sparse partial least squares classification for high 
     dimensional data. Statistical Applications in Genetics and Molecular Biology, 9, Article17. 
     doi:10.2202/1544-6115.1492

}

\author{
	Ghislain Durif (\url{http://lbbe.univ-lyon1.fr/-Durif-Ghislain-.html}). 

}
\seealso{\code{\link{spls.adapt}}.}

\examples{
### load plsgenomics library
library(plsgenomics)

### generating data
n <- 100
p <- 1000
sample1 <- sample.cont(n=100, p=1000, kstar=20, lstar=2, beta.min=0.25, beta.max=0.75, mean.H=0.2, 
					sigma.H=10, sigma.F=5, sigma.E=5)

X <- sample1$X
Y <- sample1$Y

### tuning the hyper-parameters
cv1 <- spls.adapt.tune(X=X, Y=Y, lambda.l1.range=seq(0.05, 0.95, by=0.3), ncomp.range=1:2, 
                         weight.mat=NULL, adapt=TRUE, center.X=TRUE, 
                         center.Y=TRUE, scale.X=TRUE, scale.Y=TRUE, weighted.center=FALSE, 
                         return.grid=TRUE, ncores=1, nfolds=10)
str(cv1)

### otpimal values
cv1$lambda.l1.opt
cv1$ncomp.opt

}
