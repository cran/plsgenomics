### spls.R  (2014-10)
###
###    Adaptive Sparse PLS regression for continuous response
###
### Copyright 2014-10 Ghislain DURIF
###
### Adapted from R package "spls"
### Reference: Chun H and Keles S (2010)
### "Sparse partial least squares for simultaneous dimension reduction and variable selection",
### Journal of the Royal Statistical Society - Series B, Vol. 72, pp. 3--25.
###
### This file is part of the `plsgenomics' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA

#' @title 
#' Adaptive Sparse Partial Least Squares (SPLS) regression
#' @aliases spls
#' 
#' @description 
#' The function \code{spls.adapt} performs compression and variable selection 
#' in the context of linear regression (with possible prediction) 
#' using Durif et al. (2017) adaptive SPLS algorithm.
#' 
#' @details 
#' The columns of the data matrices \code{Xtrain} and \code{Xtest} may 
#' not be standardized, since standardizing can be performed by the function 
#' \code{spls} as a preliminary step.
#' 
#' The procedure described in Durif et al. (2017) is used to compute 
#' latent sparse components that are used in a regression model.
#' In addition, when a matrix \code{Xtest} is supplied, the procedure 
#' predicts the response associated to these new values of the predictors.
#' 
#' @param Xtrain a (ntrain x p) data matrix of predictor values. 
#' \code{Xtrain} must be a matrix. Each row corresponds to an observation 
#' and each column to a predictor variable.
#' @param Ytrain a (ntrain) vector of (continuous) responses. \code{Ytrain} 
#' must be a vector or a one column matrix, and contains the response variable 
#' for each observation.
#' @param lambda.l1 a positive real value, in [0,1]. \code{lambda.l1} is the 
#' sparse penalty parameter for the dimension reduction step by sparse PLS 
#' (see details).
#' @param ncomp a positive integer. \code{ncomp} is the number of PLS 
#' components.
#' @param weight.mat a (ntrain x ntrain) matrix used to weight the l2 metric 
#' in the observation space, it can be the covariance inverse of the Ytrain 
#' observations in a heteroskedastic context. If NULL, the l2 metric is the 
#' standard one, corresponding to homoskedastic model (\code{weight.mat} is the 
#' identity matrix).
#' @param Xtest a (ntest x p) matrix containing the predictor values for the 
#' test data set. \code{Xtest} may also be a vector of length p 
#' (corresponding to only one test observation). Default value is NULL, 
#' meaning that no prediction is performed.
#' @param adapt a boolean value, indicating whether the sparse PLS selection 
#' step sould be adaptive or not (see details).
#' @param center.X a boolean value indicating whether the data matrices 
#' \code{Xtrain} and \code{Xtest} (if provided) should be centered or not.
#' @param scale.X a boolean value indicating whether the data matrices 
#' \code{Xtrain} and \code{Xtest} (if provided) should be scaled or not 
#' (\code{scale.X=TRUE} implies \code{center.X=TRUE}).
#' @param center.Y a boolean value indicating whether the response values 
#' \code{Ytrain} set should be centered or not.
#' @param scale.Y a boolean value indicating whether the response values 
#' \code{Ytrain} should be scaled or not (\code{scale.Y=TRUE} implies 
#' \code{center.Y=TRUE}).
#' @param weighted.center a boolean value indicating whether the centering 
#' should take into account the weighted l2 metric or not 
#' (if TRUE, it requires that weighted.mat is non NULL).
#' 
#' @return An object of class \code{spls} with the following attributes
#' \item{Xtrain}{the ntrain x p predictor matrix.}
#' \item{Ytrain}{the response observations.}
#' \item{sXtrain}{the centered if so and scaled if so predictor matrix.}
#' \item{sYtrain}{the centered if so and scaled if so response.}
#' \item{betahat}{the linear coefficients in model 
#' \code{sYtrain = sXtrain \%*\% betahat + residuals}.}
#' \item{betahat.nc}{the (p+1) vector containing the coefficients and intercept 
#' for the non centered and non scaled model 
#' \code{Ytrain = cbind(rep(1,ntrain),Xtrain) \%*\% betahat.nc + residuals.nc}.}
#' \item{meanXtrain}{the (p) vector of Xtrain column mean, 
#' used for centering if so.}
#' \item{sigmaXtrain}{the (p) vector of Xtrain column standard deviation, 
#' used for scaling if so.}
#' \item{meanYtrain}{the mean of Ytrain, used for centering if so.}
#' \item{sigmaYtrain}{the standard deviation of Ytrain, used for centering 
#' if so.}
#' \item{X.score}{a (n x ncomp) matrix being the observations coordinates or 
#' scores in the new component basis produced by the compression step 
#' (sparse PLS). Each column t.k of \code{X.score} is a SPLS component.}
#' \item{X.score.low}{a (n x ncomp) matrix being the PLS components only 
#' computed with the selected predictors.}
#' \item{X.loading}{the (ncomp x p) matrix of coefficients in regression of 
#' Xtrain over the new components \code{X.score}.}
#' \item{Y.loading}{the (ncomp) vector of coefficients in regression of Ytrain 
#' over the SPLS components \code{X.score}.}
#' \item{X.weight}{a (p x ncomp) matrix being the coefficients of predictors 
#' in each components produced by sparse PLS. Each column w.k of 
#' \code{X.weight} verifies t.k = Xtrain x w.k (as a matrix product).}
#' \item{residuals}{the (ntrain) vector of residuals in the model 
#' \code{sYtrain = sXtrain \%*\% betahat + residuals}.}
#' \item{residuals.nc}{the (ntrain) vector of residuals in the non centered 
#' and non scaled model 
#' \code{Ytrain = cbind(rep(1,ntrain),Xtrain) \%*\% betahat.nc + residuals.nc}.}
#' \item{hatY}{the (ntrain) vector containing the estimated reponse values 
#' on the train set of centered and scaled (if so) predictors 
#' \code{sXtrain}, \code{hatY = sXtrain \%*\% betahat}.}
#' \item{hatY.nc}{the (ntrain) vector containing the estimated reponse value 
#' on the train set of non centered and non scaled predictors \code{Xtrain}, 
#' \code{hatY.nc = cbind(rep(1,ntrain),Xtrain) \%*\% betahat.nc}.}
#' \item{hatYtest}{the (ntest) vector containing the predicted values 
#' for the response on the centered and scaled test set \code{sXtest}
#' (if provided), \code{hatYtest = sXtest \%*\% betahat}.}
#' \item{hatYtest.nc}{the (ntest) vector containing the predicted values 
#' for the response on the non centered and non scaled test set \code{Xtest} 
#' (if provided), 
#' \code{hatYtest.nc = cbind(rep(1,ntest),Xtest) \%*\% betahat.nc}.}	
#' \item{A}{the active set of predictors selected by the procedures. \code{A} 
#' is a subset of \code{1:p}.}
#' \item{betamat}{a (ncomp) list of coefficient vector betahat in the model 
#' with \code{k} components, for \code{k=1,...,ncomp}.}
#' \item{new2As}{a (ncomp) list of subset of \code{(1:p)} indicating the 
#' variables that are selected when constructing the 
#' components \code{k}, for \code{k=1,...,ncomp}.}
#' \item{lambda.l1}{the sparse hyper-parameter used to fit the model.}
#' \item{ncomp}{the number of components used to fit the model.}
#' \item{V}{the (ntrain x ntrain) matrix used to weight the metric in 
#' the sparse PLS step.}
#' \item{adapt}{a boolean value, indicating whether the sparse PLS selection 
#' step was adaptive or not.}
#' 
#' @references 
#' Durif G., Modolo L., Michaelsson J., Mold J. E., Lambert-Lacroix S., 
#' Picard F. (2017). High Dimensional Classification with combined Adaptive 
#' Sparse PLS and Logistic Regression, (in prep), 
#' available on (\url{http://arxiv.org/abs/1502.05933}).
#' 
#' Chun, H., & Keles, S. (2010). Sparse partial least squares regression for 
#' simultaneous dimension reduction and variable selection.  Journal of the 
#' Royal Statistical Society. Series B (Methodological), 72(1), 3-25. 
#' doi:10.1111/j.1467-9868.2009.00723.x
#' 
#' @author
#' Ghislain Durif (\url{http://thoth.inrialpes.fr/people/gdurif/}). 
#' 
#' Adapted in part from spls code by H. Chun, D. Chung and S.Keles 
#' (\url{https://CRAN.R-project.org/package=spls}).
#' 
#' @seealso \code{\link{spls.cv}}
#' 
#' @examples
#' ### load plsgenomics library
#' library(plsgenomics)
#' 
#' ### generating data
#' n <- 100
#' p <- 100
#' sample1 <- sample.cont(n=n, p=p, kstar=10, lstar=2, beta.min=0.25, 
#'                        beta.max=0.75, mean.H=0.2, sigma.H=10, 
#'                        sigma.F=5, sigma.E=5)
#' X <- sample1$X
#' Y <- sample1$Y
#' ### splitting between learning and testing set
#' index.train <- sort(sample(1:n, size=round(0.7*n)))
#' index.test <- (1:n)[-index.train]
#' Xtrain <- X[index.train,]
#' Ytrain <- Y[index.train,]
#' Xtest <- X[index.test,]
#' Ytest <- Y[index.test,]
#' 
#' ### fitting the model, and predicting new observations
#' model1 <- spls(Xtrain=Xtrain, Ytrain=Ytrain, lambda.l1=0.5, ncomp=2, 
#'                weight.mat=NULL, Xtest=Xtest, adapt=TRUE, center.X=TRUE, 
#'                center.Y=TRUE, scale.X=TRUE, scale.Y=TRUE, 
#'                weighted.center=FALSE)
#' 
#' str(model1)
#' 
#' ### plotting the estimation versus real values for the non centered response
#' plot(model1$Ytrain, model1$hatY.nc, 
#'      xlab="real Ytrain", ylab="Ytrain estimates")
#' points(-1000:1000,-1000:1000, type="l")
#' 
#' ### plotting residuals versus centered response values
#' plot(model1$sYtrain, model1$residuals, xlab="sYtrain", ylab="residuals")
#' 
#' ### plotting the predictor coefficients
#' plot(model1$betahat.nc, xlab="variable index", ylab="coeff")
#' 
#' ### mean squares error of prediction on test sample
#' sYtest <- as.matrix(scale(Ytest, center=model1$meanYtrain, scale=model1$sigmaYtrain))
#' sum((model1$hatYtest - sYtest)^2) / length(index.test)
#' 
#' ### plotting predicted values versus non centered real response values 
#' ## on the test set
#' plot(model1$hatYtest, sYtest, xlab="real Ytest", ylab="predicted values")
#' points(-1000:1000,-1000:1000, type="l")
#' 
#' @export
spls <- function(Xtrain, Ytrain, lambda.l1, ncomp, weight.mat=NULL, Xtest=NULL, 
                 adapt=TRUE, center.X=TRUE, center.Y=TRUE, 
                 scale.X=TRUE, scale.Y=TRUE, weighted.center=FALSE) {
	
	#####################################################################
	#### Initialisation
	#####################################################################
	Xtrain <- as.matrix(Xtrain)
	ntrain <- nrow(Xtrain) # nb observations
	p <- ncol(Xtrain) # nb covariates
	index.p <- c(1:p)
	Ytrain <- as.matrix(Ytrain)
	q <- ncol(Ytrain)
	
	if(!is.null(Xtest)) {
		ntest <- nrow(Xtest)
	}
	
	cnames <- NULL
	if(!is.null(colnames(Xtrain))) {
	     cnames <- colnames(Xtrain)
	} else {
	     cnames <- paste0(1:p)
	}
	
	
	#####################################################################
	#### Tests on type input
	#####################################################################
	# On Xtrain
	if ((!is.matrix(Xtrain)) || (!is.numeric(Xtrain))) {
		stop("Message from spls: Xtrain is not of valid type")
	}
	
	if (p==1) {
		stop("Message from spls: p=1 is not valid")}
	
	# On Xtest if necessary
	if (!is.null(Xtest)) {
		
		if (is.vector(Xtest)==TRUE) {
			Xtest <- matrix(Xtest,nrow=1)
		}
		
		Xtest <- as.matrix(Xtest)
		ntest <- nrow(Xtest) 
		
		if ((!is.matrix(Xtest)) || (!is.numeric(Xtest))) {
			stop("Message from spls: Xtest is not of valid type")}
		
		if (p != ncol(Xtest)) {
			stop("Message from spls: columns of Xtest and columns of Xtrain must be equal")
		}	
	}
	
	# On Ytrain
	if ((!is.matrix(Ytrain)) || (!is.numeric(Ytrain))) {
		stop("Message from spls: Ytrain is not of valid type")
	}
	
	if (q != 1) {
		stop("Message from spls: Ytrain must be univariate")
	}
	
	if (nrow(Ytrain)!=ntrain) {
		stop("Message from spls: the number of observations in Ytrain is not equal to the Xtrain row number")
	}
	
	# On weighting matrix V
	if(!is.null(weight.mat)) { # weighting in scalar product (in observation space of dimension n)
		V <- as.matrix(weight.mat) 
		
		if ((!is.matrix(V)) || (!is.numeric(V))) {
			stop("Message from spls: V is not of valid type")}
		
		if ((ntrain != ncol(V)) || (ntrain != nrow(V))) {
			stop("Message from spls: wrong dimension for V, must be a square matrix of size the number of observations in Xtrain")
		}
	} else { # no weighting in scalar product
		V <- diag(rep(1, ntrain), nrow=ntrain, ncol=ntrain)
	}
	
	# On hyper parameter: lambda.ridge, lambda.l1
	if ((!is.numeric(lambda.l1)) || (lambda.l1<0) || (lambda.l1>1)) {
		stop("Message from spls: lambda is not of valid type")
	}
	
	# ncomp type
	if ((!is.numeric(ncomp)) || (round(ncomp)-ncomp!=0) || (ncomp<1) || (ncomp>p)) {
		stop("Message from spls: ncomp is not of valid type")
	}
	
	# On weighted.center
	if ( (weighted.center) && (is.null(weight.mat))) {
		stop("Message from spls: if the centering is weighted, the weighting matrix V should be provided")
	}
	
	
	#####################################################################
	#### centering and scaling
	#####################################################################
	if (!weighted.center) {
		
		# Xtrain mean
		meanXtrain <- apply(Xtrain, 2, mean)
		
		# Xtrain sd
		sigmaXtrain <- apply(Xtrain, 2, sd)
		# test if predictors with null variance
		if ( any( sigmaXtrain < .Machine$double.eps )) {
			stop("Some of the columns of the predictor matrix have zero variance.")
		}
		
		# centering & eventually scaling X
		if(center.X && scale.X) {
			sXtrain <- scale( Xtrain, center=meanXtrain, scale=sigmaXtrain)
		} else if(center.X && !scale.X) {
			sXtrain <- scale( Xtrain, center=meanXtrain, scale=FALSE)
		} else {
			sXtrain <- Xtrain
		}
		
		# Y mean
		meanYtrain <- apply(Ytrain, 2, mean)
		
		# Y sd
		sigmaYtrain <- apply(Ytrain, 2, sd)
		# test if predictors with null variance
		if ( any( sigmaYtrain < .Machine$double.eps )) {
			stop("The response matrix has zero variance.")
		}
		# centering & eventually scaling Y
		if(center.Y && scale.Y) {
			sYtrain <- scale( Ytrain, center=meanYtrain, scale=sigmaYtrain )
		} else if(center.Y && !scale.Y) {
			sYtrain <- scale( Ytrain, center=meanYtrain, scale=FALSE )
		} else {
			sYtrain <- Ytrain
		}
		
		# Xtest
		if (!is.null(Xtest)) {
			
			## centering and scaling depend on Xtest
			if(center.X && scale.X) {
				sXtest <- scale( Xtest, center=meanXtrain, scale=sigmaXtrain)
			} else if(center.X && !scale.X) {
				sXtest <- scale( Xtest, center=meanXtrain, scale=FALSE)
			} else {
				sXtest <- Xtest
			}
			
		}
		
	} else { # weighted scaling
		
		sumV <- sum(diag(V))
		
		# X mean
		meanXtrain <- matrix(diag(V), nrow=1) %*% Xtrain / sumV
		
		# X sd
		sigmaXtrain <- apply(Xtrain, 2, sd)
		# test if predictors with null variance
		if ( any( sigmaXtrain < .Machine$double.eps ) ) {
			stop("Some of the columns of the predictor matrix have zero variance.")
		}
		# centering & eventually scaling X
		sXtrain <- scale( Xtrain, center=meanXtrain, scale=FALSE )
		
		# Y mean
		meanYtrain <- matrix(diag(V), nrow=1) %*% Ytrain / sumV
		
		# Y sd
		sigmaYtrain <- apply(Ytrain, 2, sd)
		# test if predictors with null variance
		if ( any( sigmaYtrain < .Machine$double.eps ) ) {
			stop("The response matrix have zero variance.")
		}
		# centering & eventually scaling Y
		sYtrain <- scale( Ytrain, center=meanYtrain, scale=FALSE )
		
		# Xtest
		if (!is.null(Xtest)) {
			sXtest <- scale( Xtest, center=meanXtrain, scale=FALSE )
		}
		
	}
	
	
	#####################################################################
	#### Result objects
	#####################################################################
	betahat <- matrix(0, nrow=p, ncol=1)
	betamat <- list()
	X1 <- sXtrain
	Y1 <- sYtrain
	
	W <- matrix(data=NA, nrow=p, ncol=ncomp) # spls weight over each component
	T <- matrix(data=NA, nrow=ntrain, ncol=ncomp) # spls components
	P <- matrix(data=NA, nrow=ncomp, ncol=p) # regression of X over T
	Q <- matrix(data=NA, nrow=ncomp, ncol=q) # regression of Y over T
	
	
	#####################################################################
	#### Main iteration
	#####################################################################
	
	if ( is.null(colnames(Xtrain)) ) {
		Xnames <- index.p
	} else { 
		Xnames <- colnames(Xtrain)
	}
	
	new2As <- list()
	
	## SPLS
	for (k in 1:ncomp) {
		
		## define M
		M <- t(X1) %*% (V %*% Y1)
		
		#### soft threshold
		Mnorm1 <- median( abs(M) )
		
		M <- M / Mnorm1
		
		## adpative version
		if (adapt) {
			wi <- 1/abs(M)
			
			what <- ust.adapt(M, lambda.l1, wi)
			
		} else {
		## non adaptive version
			what <- ust(M, lambda.l1)
		}
		
		#### construct active set A
		A <- unique( index.p[ what!=0 | betahat[,1]!=0 ] )
		new2A <- index.p[ what!=0 & betahat[,1]==0 ]
		
		#### fit pls with selected predictors (meaning in A)
		X.A <- sXtrain[ , A, drop=FALSE ]           
		plsfit <- wpls( Xtrain=X.A, Ytrain=sYtrain, weight.mat=V, ncomp=min(k,length(A)), type="pls1", center.X=FALSE, scale.X=FALSE, center.Y=FALSE, scale.Y=FALSE, weighted.center=FALSE )
		
		
		#### output storage
		
		# weights
		w.k <- matrix(data=what, ncol=1)
		w.k <- w.k / sqrt(as.numeric(t(w.k) %*% w.k))
		W[,k] <- w.k
		
		# components on total observation space
		t.k <- (X1 %*% w.k) / as.numeric(t(w.k) %*% w.k)
		T[,k] <- t.k
		
		# regression of X over T
		p.k <- (t(X1) %*% (V %*% t.k)) / as.numeric(t(t.k) %*% (V %*% t.k))
		P[k,] <- t(p.k)
		
		# regression of Y over T
		q.k <- (t(Y1) %*% (V %*% t.k)) / as.numeric(t(t.k) %*% (V %*% t.k))
		Q[k,] <- t(q.k)
		
		
		## update
		
		Y1 <- sYtrain - plsfit$T %*% plsfit$Q
		X1 <- sXtrain
		X1[,A] <- sXtrain[,A] - plsfit$T %*% plsfit$P
		
		
		betahat <- matrix( 0, p, q )
		betahat[A,] <- matrix( plsfit$coeff, length(A), q )
		betamat[[k]] <- betahat # for cv.spls
		
		# variables that join the active set
		new2As[[k]] <- new2A
		
		
	}
	
	##### return objects
	
	hatY <- numeric(ntrain)
	hatY.nc <- numeric(ntrain)
	
	## components in lower subspace of selected variables
	T.low <- plsfit$T
	
	## estimations
	hatY <- sXtrain %*% betahat
	
	## residuals
	residuals <- sYtrain - hatY
	
	#### betahat for non centered and non scaled data
	if((!scale.X) || (weighted.center)) { # if X non scaled, betahat don't have to be corrected regards sd.x
		sd.X <- rep(1, p)
	} else { # if X is scaled, it has to
		sd.X <- sigmaXtrain
	}
	if((!scale.Y) || (weighted.center)) {
		sd.Y <- 1
	} else {
		sd.Y <- sigmaYtrain
	}
	
	betahat.nc <- sd.Y * betahat / sd.X
	intercept <- meanYtrain - ( sd.Y * (drop( (meanXtrain / sd.X) %*% betahat)) )
	betahat.nc <- as.matrix(c(intercept, betahat.nc))
	
	#### non centered non scaled version of estimation and residuals
	hatY.nc <- cbind(rep(1,ntrain),Xtrain) %*% betahat.nc
	residuals.nc <- Ytrain - hatY.nc
	
	
	## predictions
	if(!is.null(Xtest)) {
		hatYtest <- sXtest %*% betahat
		hatYtest.nc <- cbind(rep(1,ntest),Xtest) %*% betahat.nc
	} else {
		hatYtest <- NULL
		hatYtest.nc <- NULL
	}
	
	rownames(betahat) <- cnames
	
	Anames <- cnames[A]
	
	
	#### return object
	result <- list( Xtrain=Xtrain, Ytrain=Ytrain, sXtrain=sXtrain, sYtrain=sYtrain,
				 betahat=betahat, betahat.nc=betahat.nc,
				 meanXtrain=meanXtrain, meanYtrain=meanYtrain, sigmaXtrain=sigmaXtrain, sigmaYtrain=sigmaYtrain,
				 X.score=T, X.score.low=T.low, X.loading=P, Y.loading=Q, X.weight=W, 
				 residuals=residuals, residuals.nc=residuals.nc,
				 hatY=hatY, hatY.nc=hatY.nc,
				 hatYtest=hatYtest, hatYtest.nc=hatYtest.nc,
				 A=A, Anames=Anames, betamat=betamat, new2As=new2As,
				 lambda.l1=lambda.l1, ncomp=ncomp,
				 V=V, adapt=adapt)
	
	class(result) <- "spls"
	return(result)
	
	
}
