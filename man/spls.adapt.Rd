\name{spls.adapt}
\alias{spls.adapt}

\title{Classification by Ridge Iteratively Reweighted Least Squares followed by 
Adaptive Sparse PLS regression for binary response}
\usage{
spls.adapt(Xtrain, Ytrain, lambda.l1, ncomp, weight.mat=NULL, Xtest=NULL, adapt=TRUE, 
               center.X=TRUE, center.Y=TRUE, scale.X=TRUE, 
               scale.Y=TRUE, weighted.center=FALSE)
}
\arguments{
	\item{Xtrain}{a (ntrain x p) data matrix of predictors. \code{Xtrain} must be a matrix. 
	Each row corresponds to an observation and each column to a predictor variable.}
	\item{Ytrain}{a (ntrain) vector of (continuous) responses. \code{Ytrain} must be a vector or a 
	one column matrix, and contains the response variable for each observation.}
	\item{lambda.l1}{a positive real value, in [0,1]. \code{lambda.l1} is the sparse penalty 
	parameter for the dimension reduction step by sparse PLS (see details).}
	\item{ncomp}{a positive integer. \code{ncomp} is the number of PLS components. 
	If \code{ncomp}=0,then the Ridge regression is performed without dimension reduction.}
	\item{weight.mat}{a (ntrain x ntrain) matrix used to weight the l2 metric in observation 
	space if necessary, especially the covariance inverse of the Ytrain observations in 
	heteroskedastic context. If NULL, the l2 metric is the standard one, corresponding to 
	homoskedastic model.}
	\item{Xtest}{a (ntest x p) matrix containing the predictors for the test data set. 
	\code{Xtest} may also be a vector of length p (corresponding to only one test observation). 
	If \code{Xtest} is not equal to NULL, then the prediction step is made for these new predictor 
	variables.}
	\item{adapt}{a boolean value, indicating whether the sparse PLS selection step sould be 
	adaptive or not.}
	\item{center.X}{a boolean value indicating whether the design matrices Xtrain in train set and 
	Xtest in test set if non NULL should be centered or not}
	\item{scale.X}{a boolean value indicating whether the design matrices Xtrain in train set and 
	Xtest in test set if non NULL should be scaled or not, scale.X=TRUE implies center.X=TRUE}
	\item{center.Y}{a boolean value indicating whether the response Ytrain in train set should 
	be centered or not}
	\item{scale.Y}{a boolean value indicating whether the response Ytrain should be scaled or not, 
	scale.Y=TRUE implies center.Y=TRUE}
	\item{weighted.center}{a boolean value indicating whether should the centering take into 
	account the weighted l2 metric or not (if TRUE, it implies that weighted.mat is non NULL).}
}
\description{
	The function \code{spls.adapt} performs compression, variable selection in regression context 
	(with possible prediction) using Durif et al. (2015) adaptive SPLS algorithm, based on sparse 
	PLS developped by Chun and Keles (2010).
  
}

\details{
	The columns of the data matrices \code{Xtrain} and \code{Xtest} may not be standardized, 
	since standardizing is can be performed by the function \code{spls.adapt} as a preliminary 
	step before the algorithm is run. 
	
	The procedure described in Durif et al. (2015) is used to determine
	latent sparse components to be used for regression and when \code{Xtest} 
	is not equal to NULL, the procedure predicts the response for these new 
	predictor variables.  
 
}

\value{
	A list with the following components:
	\item{Xtrain}{the design matrix.}
	\item{Ytrain}{the response observations.}
	\item{sXtrain}{the centered if so and scaled if so design matrix.}
	\item{sYtrain}{the centered if so and scaled if so response.}
	\item{betahat}{the linear coefficients in model 
	\code{sYtrain = sXtrain \%*\% betahat + residuals}.}
	\item{betahat.nc}{the (p+1) vector containing the coefficients and intercept for the non 
	centered and non scaled model 
	\code{Ytrain = cbind(rep(1,ntrain),Xtrain) \%*\% betahat.nc + residuals.nc}.}
	\item{meanXtrain}{the (p) vector of Xtrain column mean, used for centering if so.}
	\item{sigmaXtrain}{the (p) vector of Xtrain column standard deviation, used for scaling if so.}
	\item{meanYtrain}{the mean of Ytrain, used for centering if so.}
	\item{sigmaYtrain}{the standard deviation of Ytrain, used for centering if so.}
	\item{X.score}{a (n x ncomp) matrix being the observations coordinates or scores in the 
	new component basis produced by the compression step (sparse PLS). Each column t.k of 
	\code{X.score} is a new component.}
	\item{X.score.low}{a (n x ncomp) matrix being the observations coordinates in the subspace of 
	selected variables.}
	\item{X.loading}{the (ncomp x p) matrix of coefficients in regression of Xtrain over the new 
	components \code{X.score}.}
	\item{Y.loading}{the (ncomp) vector of coefficients in regression of Ytrain over the new 
	components \code{X.score}.}
	\item{X.weight}{a (p x ncomp) matrix being the coefficients of predictors in each components 
	produced by sparse PLS. Each column w.k of \code{X.weight} verifies t.k = Xtrain x w.k 
	(as a matrix product).}
	\item{residuals}{the (ntrain) vector of residuals in the model 
	\code{sYtrain = sXtrain \%*\% betahat + residuals}.}
	\item{residuals.nc}{the (ntrain) vector of residuals in the non centered non scaled model 
	\code{Ytrain = cbind(rep(1,ntrain),Xtrain) \%*\% betahat.nc + residuals.nc}.}
	\item{hatY}{the (ntrain) vector containing the estimated reponse value on the train set of 
	centered and scaled (if so) predictors \code{sXtrain}, \code{hatY = sXtrain \%*\% betahat}.}
	\item{hatY.nc}{the (ntrain) vector containing the estimated reponse value on the train set of 
	non centered and non scaled predictors \code{Xtrain}, 
	\code{hatY.nc = cbind(rep(1,ntrain),Xtrain) \%*\% betahat.nc}.}
	\item{hatYtest}{the (ntest) vector containing the predicted values for the response on the 
	centered and scaled test set \code{sXtest} if \code{Xtest} non NULL, 
	\code{hatYtest = sXtest \%*\% betahat}.}
	\item{hatYtest.nc}{the (ntest) vector containing the predicted values for the response on the 
	non centered and non scaled test set \code{Xtest} if non NULL, 
	\code{hatYtest.nc = cbind(rep(1,ntest),Xtest) \%*\% betahat.nc}.}	
	\item{A}{the active set of predictors selected by the procedures. \code{A} is a subset 
	of \code{1:p}.}
	\item{betamat}{a (ncomp) list of coefficient vector betahat in the model with \code{k} 
	components, for \code{k=1,...,ncomp}.}
	\item{new2As}{a (ncomp) list of subset of \code{(1:p)} indicating the variables selected when 
	constructing the components \code{k}, for \code{k=1,...,ncomp}.}
	\item{lambda.l1}{the sparse hyper-parameter used to fit the model.}
	\item{ncomp}{the number of components used to fit the model.}
	\item{V}{the (ntrain x ntrain) matrix used to weight the metric in the sparse PLS step.}
	\item{adapt}{a boolean value, indicating whether the sparse PLS selection step was adaptive 
	or not.}
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
     
     Adapted in part from spls code by H. Chun, D. Chung and S.Keles 
     (\url{http://cran.r-project.org/web/packages/spls/index.html}).

}
\seealso{\code{\link{spls.adapt.tune}}.}

\examples{
### load plsgenomics library
library(plsgenomics)

### generating data
n <- 100
p <- 1000
sample1 <- sample.cont(n=n, p=p, kstar=20, lstar=2, beta.min=0.25, beta.max=0.75, 
                         mean.H=0.2, sigma.H=10, sigma.F=5, sigma.E=5)

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
model1 <- spls.adapt(Xtrain=Xtrain, Ytrain=Ytrain, lambda.l1=0.5, ncomp=2, weight.mat=NULL, 
                         Xtest=Xtest, adapt=TRUE, center.X=TRUE, center.Y=TRUE, scale.X=TRUE, 
                         scale.Y=TRUE, weighted.center=FALSE)
str(model1)

### plotting the estimation versus real values for the non centered response
plot(model1$Ytrain, model1$hatY.nc, xlab="real Ytrain", ylab="Ytrain estimates")
points(-1000:1000,-1000:1000, type="l")

### plotting residuals versus centered response values
plot(model1$sYtrain, model1$residuals, xlab="sYtrain", ylab="residuals")

### plotting the predictor coefficients
plot(model1$betahat.nc, xlab="variable index", ylab="coeff")

### mean squares error of prediction on test sample
sYtest <- as.matrix(scale(Ytest, center=model1$meanYtrain, scale=model1$sigmaYtrain))
sum((model1$hatYtest - sYtest)^2) / length(index.test)

### plotting predicted values versus non centered real response values on the test set
plot(model1$hatYtest, sYtest, xlab="real Ytest", ylab="predicted values")
points(-1000:1000,-1000:1000, type="l")

}
