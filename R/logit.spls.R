### logit.spls.R  (2014-10)
###
###    Ridge Iteratively Reweighted Least Squares followed by Adaptive Sparse PLS regression for binary response
###
### Copyright 2014-10 Ghislain DURIF
###
### Adapted from rpls function in plsgenomics package, copyright 2006-01 Sophie Lambert-Lacroix
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
#' Classification procedure for binary response based on a logistic model, 
#' solved by a combination of the Ridge Iteratively Reweighted Least Squares 
#' (RIRLS) algorithm and the Adaptive Sparse PLS (SPLS) regression
#' @aliases logit.spls
#' 
#' @description 
#' The function \code{logit.spls} performs compression and variable selection 
#' in the context of binary classification (with possible prediction) 
#' using Durif et al. (2018) algorithm based on Ridge IRLS and sparse PLS.
#' 
#' @details 
#' The columns of the data matrices \code{Xtrain} and \code{Xtest} may 
#' not be standardized, since standardizing can be performed by the function 
#' \code{logit.spls} as a preliminary step.
#' 
#' The procedure described in Durif et al. (2018) is used to compute
#' latent sparse components that are used in a logistic regression model.
#' In addition, when a matrix \code{Xtest} is supplied, the procedure 
#' predicts the response associated to these new values of the predictors.
#' 
#' @param Xtrain a (ntrain x p) data matrix of predictor values. 
#' \code{Xtrain} must be a matrix. Each row corresponds to an observation 
#' and each column to a predictor variable.
#' @param Ytrain a (ntrain) vector of (continuous) responses. \code{Ytrain} 
#' must be a vector or a one column matrix, and contains the response variable 
#' for each observation. \code{Ytrain} should take values in \{0,1\}.
#' @param lambda.ridge a positive real value. \code{lambda.ridge} is the Ridge 
#' regularization parameter for the RIRLS algorithm (see details).
#' @param lambda.l1 a positive real value, in [0,1]. \code{lambda.l1} is the 
#' sparse penalty parameter for the dimension reduction step by sparse PLS 
#' (see details).
#' @param ncomp a positive integer. \code{ncomp} is the number of 
#' PLS components. If \code{ncomp=0},then the Ridge regression is performed 
#' without any dimension reduction (no SPLS step).
#' @param Xtest a (ntest x p) matrix containing the predictor values for the 
#' test data set. \code{Xtest} may also be a vector of length p 
#' (corresponding to only one test observation). Default value is NULL, 
#' meaning that no prediction is performed.
#' @param adapt a boolean value, indicating whether the sparse PLS selection 
#' step sould be adaptive or not (see details).
#' @param maxIter a positive integer. \code{maxIter} is the maximal number of 
#' iterations in the Newton-Raphson parts in the RIRLS algorithm (see details).
#' @param svd.decompose a boolean parameter. \code{svd.decompose} indicates 
#' wether or not the predictor matrix \code{Xtrain} should be decomposed by 
#' SVD (singular values decomposition) for the RIRLS step (see details).
#' @param center.X a boolean value indicating whether the data matrices 
#' \code{Xtrain} and \code{Xtest} (if provided) should be centered or not.
#' @param scale.X a boolean value indicating whether the data matrices 
#' \code{Xtrain} and \code{Xtest} (if provided) should be scaled or not 
#' (\code{scale.X=TRUE} implies \code{center.X=TRUE}) in the spls step.
#' @param weighted.center a boolean value indicating whether the centering 
#' should take into account the weighted l2 metric or not in the SPLS step.
#' 
#' @return An object of class \code{logit.spls} with the following attributes
#' \item{Coefficients}{the (p+1) vector containing the linear coefficients 
#' associated to the predictors and intercept in the logistic model 
#' explaining the response Y.}
#' \item{hatY}{the (ntrain) vector containing the estimated response value on 
#' the train set \code{Xtrain}.}
#' \item{hatYtest}{the (ntest) vector containing the predicted labels 
#' for the observations from \code{Xtest} (if provided).}
#' \item{DeletedCol}{the vector containing the indexes of columns with null 
#' variance in \code{Xtrain} that were skipped in the procedure.}
#' \item{A}{the active set of predictors selected by the procedures. 
#' \code{A} is a subset of 1:p.}
#' \item{Anames}{Vector of selected predictor names, i.e. the names of the 
#' columns from \code{Xtrain} that are in \code{A}.}
#' \item{converged}{a \{0,1\} value indicating whether the RIRLS algorithm did
#' converge in less than \code{maxIter} iterations or not.}
#' \item{X.score}{a (n x ncomp) matrix being the observations coordinates or 
#' scores in the new component basis produced by the SPLS step (sparse PLS). 
#' Each column t.k of \code{X.score} is a SPLS component.}
#' \item{X.weight}{a (p x ncomp) matrix being the coefficients of predictors 
#' in each components produced by sparse PLS. Each column w.k of 
#' \code{X.weight} verifies t.k = Xtrain x w.k (as a matrix product).}
#' \item{Xtrain}{the design matrix.}
#' \item{sXtrain}{the scaled predictor matrix.}
#' \item{Ytrain}{the response observations.}
#' \item{sPseudoVar}{the scaled pseudo-response produced by the RIRLS
#' algorithm.}
#' \item{lambda.ridge}{the Ridge hyper-parameter used to fit the model.}
#' \item{lambda.l1}{the sparse hyper-parameter used to fit the model.}
#' \item{ncomp}{the number of components used to fit the model.}
#' \item{V}{the (ntrain x ntrain) matrix used to weight the metric in the 
#' sparse PLS step. \code{V} is the inverse of the covariance matrix of the 
#' pseudo-response produced by the RIRLS step.}
#' \item{proba}{the (ntrain) vector of estimated probabilities for the 
#' observations in code \code{Xtrain}, that are used to estimate the 
#' \code{hatY} labels.}
#' \item{proba.test}{the (ntest) vector of predicted probabilities for the 
#' new observations in \code{Xtest}, that are used to predict the 
#' \code{hatYtest} labels.}
#' 
#' @references 
#' Durif, G., Modolo, L., Michaelsson, J., Mold, J.E., Lambert-Lacroix, S., 
#' Picard, F., 2018. High dimensional classification with combined 
#' adaptive sparse PLS and logistic regression. Bioinformatics 34, 
#' 485--493. \doi{10.1093/bioinformatics/btx571}.
#' Available at \url{http://arxiv.org/abs/1502.05933}.
#' 
#' @author
#' Ghislain Durif (\url{https://gdurif.perso.math.cnrs.fr/}).
#' 
#' @seealso \code{\link{spls}}, \code{\link{logit.spls.cv}}
#' 
#' @examples
#' \dontrun{
#' ### load plsgenomics library
#' library(plsgenomics)
#' 
#' ### generating data
#' n <- 100
#' p <- 100
#' sample1 <- sample.bin(n=n, p=p, kstar=20, lstar=2, 
#'                       beta.min=0.25, beta.max=0.75, 
#'                       mean.H=0.2, sigma.H=10, sigma.F=5)
#' X <- sample1$X
#' Y <- sample1$Y
#' 
#' ### splitting between learning and testing set
#' index.train <- sort(sample(1:n, size=round(0.7*n)))
#' index.test <- (1:n)[-index.train]
#' 
#' Xtrain <- X[index.train,]
#' Ytrain <- Y[index.train,]
#' Xtest <- X[index.test,]
#' Ytest <- Y[index.test,]
#' 
#' ### fitting the model, and predicting new observations
#' model1 <- logit.spls(Xtrain=Xtrain, Ytrain=Ytrain, lambda.ridge=2, 
#'                      lambda.l1=0.5, ncomp=2, Xtest=Xtest, adapt=TRUE, 
#'                      maxIter=100, svd.decompose=TRUE)
#'                      
#' str(model1)
#' 
#' ### prediction error rate
#' sum(model1$hatYtest!=Ytest) / length(index.test)
#' }
#' 
#' @export
logit.spls <- function(Xtrain, Ytrain, lambda.ridge, lambda.l1, ncomp, 
                       Xtest=NULL, adapt=TRUE, maxIter=100, svd.decompose=TRUE, 
                       center.X=TRUE, scale.X=FALSE, weighted.center=TRUE) {
	
	
	#####################################################################
	#### Initialisation
	#####################################################################
	Xtrain <- as.matrix(Xtrain)
	ntrain <- nrow(Xtrain) # nb observations
	p <- ncol(Xtrain) # nb covariates
	index.p <- c(1:p)
	if(is.factor(Ytrain)) {
	     Ytrain <- as.numeric(levels(Ytrain))[Ytrain]
	}
	Ytrain <- as.integer(Ytrain)
	Ytrain <- as.matrix(Ytrain)
	q <- ncol(Ytrain)
	one <- matrix(1,nrow=1,ncol=ntrain)
	
	cnames <- NULL
	if(!is.null(colnames(Xtrain))) {
	     cnames <- colnames(Xtrain)
	} else {
	     cnames <- paste0(1:p)
	}
	
	
	#####################################################################
	#### Tests on type input
	#####################################################################
	
	# if multicategorical response
	if(length(table(Ytrain)) > 2) {
	     warning("message from logit.spls: multicategorical response")
	     results = multinom.spls(Xtrain=Xtrain, Ytrain=Ytrain, lambda.ridge=lambda.ridge, 
	                             lambda.l1=lambda.l1, ncomp=ncomp, Xtest=Xtest, adapt=adapt, 
	                             maxIter=maxIter, svd.decompose=svd.decompose, 
	                             center.X=center.X, scale.X=scale.X, weighted.center=weighted.center)
	     return(results)
	}
	
	# On Xtrain
	if ((!is.matrix(Xtrain)) || (!is.numeric(Xtrain))) {
		stop("Message from logit.spls: Xtrain is not of valid type")
	}
	
	if (ncomp > p) {
	     warning("Message from logit.spls: ncomp>p is not valid, ncomp is set to p")
	     ncomp <- p
	}
	
	if (p==1) {
	     # stop("Message from logit.spls: p=1 is not valid")
	     warning("Message from logit.spls: p=1 is not valid, ncomp is set to 0")
	     ncomp <- 0
	}
	
	# On Xtest if necessary
	if (!is.null(Xtest)) {
		
		if (is.vector(Xtest)==TRUE) {
			Xtest <- matrix(Xtest,ncol=p)
		}
          
          Xtest <- as.matrix(Xtest)
		ntest <- nrow(Xtest) 
		
		if ((!is.matrix(Xtest)) || (!is.numeric(Xtest))) {
			stop("Message from logit.spls: Xtest is not of valid type")
		}
		
		if (p != ncol(Xtest)) {
			stop("Message from logit.spls: columns of Xtest and columns of Xtrain must be equal")
		}	
	}
	
	# On Ytrain
	if ((!is.matrix(Ytrain)) || (!is.numeric(Ytrain))) {
		stop("Message from logit.spls: Ytrain is not of valid type")
	}
	
	if (q != 1) {
		stop("Message from logit.spls: Ytrain must be univariate")
	}
	
	if (nrow(Ytrain)!=ntrain) {
		stop("Message from logit.spls: the number of observations in Ytrain is not equal to the Xtrain row number")
	}
	
	# On Ytrain value
	if (sum(is.na(Ytrain))!=0) {
		stop("Message from logit.spls: NA values in Ytrain")
	}
	
	if (sum(!(Ytrain %in% c(0,1)))!=0) {
		stop("Message from logit.spls: Ytrain is not of valid type")
	}
	
	if (sum(as.numeric(table(Ytrain))==0)!=0) {
		stop("Message from logit.spls: there are empty classes")
	}
	
	# On hyper parameter: lambda.ridge, lambda.l1
	if ((!is.numeric(lambda.ridge)) || (lambda.ridge<0) || (!is.numeric(lambda.l1)) || (lambda.l1<0)) {
		stop("Message from logit.spls: lambda is not of valid type")
	}
	
	# ncomp type
	if ((!is.numeric(ncomp)) || (round(ncomp)-ncomp!=0) || (ncomp<0) || (ncomp>p)) {
		stop("Message from logit.spls: ncomp is not of valid type")
	}
	
	# maxIter
	if ((!is.numeric(maxIter)) || (round(maxIter)-maxIter!=0) || (maxIter<1)) {
		stop("Message from logit.spls: maxIter is not of valid type")
	}
	
	
	#####################################################################
	#### Move into the reduced space
	#####################################################################
	
	r <- min(p, ntrain)
	DeletedCol <- NULL
	
	### Standardize the Xtrain matrix
	# standard deviation (biased one) of Xtrain
	sigma2train <- apply(Xtrain, 2, var) * (ntrain-1)/(ntrain)
	
	# test on sigma2train
	# predictor with null variance ?
	if (sum(sigma2train < .Machine$double.eps)!=0){
		
		# predicteur with non null variance < 2 ?
		if (sum(sigma2train < .Machine$double.eps)>(p-2)){
			stop("Message from logit.spls: the procedure stops because number of predictor variables with no null variance is less than 1.")
		}
		
		warning("There are covariables with nul variance")
		
		# remove predictor with null variance
		Xtrain <- Xtrain[,which(sigma2train >= .Machine$double.eps)]
		if (!is.null(Xtest)) {
			Xtest <- Xtest[,which(sigma2train>= .Machine$double.eps)]
		}
		
		# list of removed predictors
		DeletedCol <- index.p[which(sigma2train < .Machine$double.eps)]
		
		# removed null standard deviation
		sigma2train <- sigma2train[which(sigma2train >= .Machine$double.eps)]
		
		# new number of predictor
		p <- ncol(Xtrain)
		r <- min(p,ntrain)
	}
	
	if(!is.null(DeletedCol)) {
	     cnames <- cnames[-DeletedCol]
	}
	
	# mean of Xtrain
	meanXtrain <- apply(Xtrain,2,mean)
	
	# center and scale Xtrain
     if(center.X && scale.X) {
          sXtrain <- scale(Xtrain, center=meanXtrain, scale=sqrt(sigma2train))
     } else if(center.X && !scale.X) {
          sXtrain <- scale(Xtrain, center=meanXtrain, scale=FALSE)
     } else {
          sXtrain <- Xtrain
     }
	
	sXtrain.nosvd <- sXtrain # keep in memory if svd decomposition
	
	# Compute the svd when necessary -> case p > ntrain (high dim)
	if ((p > ntrain) && (svd.decompose) && (p>1)) {
		# svd de sXtrain
		svd.sXtrain <- svd(t(sXtrain))
		# number of singular value non null
		r <- length(svd.sXtrain$d[abs(svd.sXtrain$d)>10^(-13)])
		V <- svd.sXtrain$u[,1:r]
		D <- diag(c(svd.sXtrain$d[1:r]))
		U <- svd.sXtrain$v[,1:r]
		sXtrain <- U %*% D
		rm(D)
		rm(U)
		rm(svd.sXtrain)
	}
	
	# center and scale Xtest if necessary
	if (!is.null(Xtest)) {
		
		meanXtest <- apply(Xtest,2,mean)
		sigma2test <- apply(Xtest,2,var)
		
		if(center.X && scale.X) {
		     sXtest <- scale(Xtest, center=meanXtrain, scale=sqrt(sigma2train))
		} else if(center.X && !scale.X) {
		     sXtest <- scale(Xtest, center=meanXtrain, scale=FALSE)
		} else {
		     sXtest <- Xtest
		}
		
		sXtest.nosvd <- sXtest # keep in memory if svd decomposition
		
		# if svd decomposition
		if ((p > ntrain) && (svd.decompose)) {
			sXtest <- sXtest%*%V
		}
	}
	
	
	#####################################################################
	#### Ridge IRLS step
	#####################################################################
	
	fit <- wirrls(Y=Ytrain, Z=cbind(rep(1,ntrain),sXtrain), Lambda=lambda.ridge, NbrIterMax=maxIter, WKernel=diag(rep(1,ntrain)))
	
	converged=fit$Cvg
	
	#  Check WIRRLS convergence
	if (converged==0) {
		warning("Message from logit.spls : Ridge IRLS did not converge; try another lambda.ridge value")
	}
	
	# if ncomp == 0 then wirrls without spls step
	if (ncomp==0) {
		BETA <- fit$Coefficients
	}
	
	
	#####################################################################
	#### weighted SPLS step
	#####################################################################
	
	# if ncomp > 0
	if (ncomp!=0) {
		
		#Compute ponderation matrix V and pseudo variable z
		#Pseudovar = Eta + W^-1 Psi
		
		# Eta = X * betahat (covariate summary)
		Eta <- cbind(rep(1, ntrain), sXtrain) %*% fit$Coefficients
		
		## Run SPLS on Xtrain without svd decomposition
		sXtrain = sXtrain.nosvd
		if (!is.null(Xtest)) {
			sXtest = sXtest.nosvd
		}
		
		# mu = h(Eta)
		mu = 1 / (1 + exp(-Eta))
		
		# ponderation matrix : V
		diagV <- mu * (1-mu)
		V <- diag(c(diagV))
		
		# inv de V
		diagVinv = 1/ifelse(diagV!=0, diagV, diagV+0.00000001)
		Vinv = diag(c(diagVinv))
		
		Psi <- Ytrain-mu
		
		pseudoVar = Eta + Vinv %*% Psi
		
		# V-Center the sXtrain and pseudo variable
		sumV=sum(diagV)
		
		# Weighted centering of Pseudo variable
		VmeanPseudoVar <- sum(V %*% Eta + Psi ) / sumV
		
		# Weighted centering of sXtrain
		VmeansXtrain <- t(diagV)%*%sXtrain/sumV
		
		# SPLS(X, pseudo-var, weighting = V)
		resSPLS = spls(Xtrain=sXtrain, Ytrain=pseudoVar, ncomp=ncomp, 
		               weight.mat=V, lambda.l1=lambda.l1, adapt=adapt, 
		               center.X=center.X, scale.X=scale.X, 
		               center.Y=TRUE, scale.Y=FALSE, 
		               weighted.center=weighted.center)
		
		BETA = resSPLS$betahat.nc
		
	} else {
	     resSPLS = list(A=NULL, X.score=NULL, X.weight=NULL, sXtrain=sXtrain, sYtrain=NULL, V=NULL)
	}
	
	
	#####################################################################
	#### classification step
	#####################################################################
	
	hatY <- numeric(ntrain)
	hatY <- cbind(rep(1,ntrain),sXtrain) %*% BETA
	hatY <- as.numeric(hatY>0)
	proba <- numeric(ntrain)
	proba <- inv.logit( cbind(rep(1,ntrain),sXtrain) %*% BETA )
	
	if (!is.null(Xtest)) {
		
		hatYtest <- cbind(rep(1,ntest),sXtest) %*% BETA
		hatYtest <- as.numeric(hatYtest>0)
		proba.test = inv.logit( cbind(rep(1,ntest),sXtest) %*% BETA )
		
	} else {
		
		hatYtest <- NULL
		proba.test <- NULL
		
	}
	
	
	#####################################################################
	#### Conclude
	#####################################################################
	
	Coefficients=BETA
	
	if(p > 1) {
	     Coefficients[-1] <- diag(c(1/sqrt(sigma2train)))%*%BETA[-1]
	} else {
	     Coefficients[-1] <- (1/sqrt(sigma2train))%*%BETA[-1]
	}
	
	Coefficients[1] <- BETA[1] - meanXtrain %*% Coefficients[-1]
	
	
	#### RETURN
	
	Anames <- cnames[resSPLS$A]
	
	result <- list(Coefficients=Coefficients, hatY=hatY, hatYtest=hatYtest, DeletedCol=DeletedCol, 
	               A=resSPLS$A, Anames=Anames, converged=converged, 
	               X.score=resSPLS$X.score, X.weight=resSPLS$X.weight, 
	               sXtrain=resSPLS$sXtrain, sPseudoVar=resSPLS$sYtrain, 
	               lambda.ridge=lambda.ridge, lambda.l1=lambda.l1, ncomp=ncomp, 
	               V=resSPLS$V, proba=proba, proba.test=proba.test, Xtrain=Xtrain, Ytrain=Ytrain)
	class(result) <- "logit.spls"
	return(result)
	
	
	
	
}
