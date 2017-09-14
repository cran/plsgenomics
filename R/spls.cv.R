### spls.cv.R  (2014-10)
###
###    Tuning parameters (ncomp, lambda.l1) for adaptive sparse PLS regression 
###    for continuous response, by K-fold cross-validation
###
### Copyright 2014-10 Ghislain DURIF
###
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
#' Cross-validation procedure to calibrate the parameters (ncomp, lambda.l1) 
#' of the Adaptive Sparse PLS regression
#' @aliases spls.cv
#' 
#' @description 
#' The function \code{spls.cv} chooses the optimal values for the 
#' hyper-parameter of the \code{spls} procedure, by minimizing the mean 
#' squared error of prediction over the hyper-parameter grid, 
#' using Durif et al. (2017) adaptive SPLS algorithm.
#' 
#' @details 
#' The columns of the data matrices \code{Xtrain} and \code{Xtest} may not 
#' be standardized, since standardizing can be performed by the function 
#' \code{spls.cv} as a preliminary step.
#' 
#' The procedure is described in Durif et al. (2017). The K-fold 
#' cross-validation can be summarize as follow: the train set is partitioned 
#' into K folds, for each value of hyper-parameters the model is fit K times, 
#' using each fold to compute the prediction error rate, and fitting the 
#' model on the remaining observations. The cross-validation procedure returns 
#' the optimal hyper-parameters values, meaning the one that minimize 
#' the mean squared error of prediction averaged over all the folds.
#' 
#' This procedures uses the \code{mclapply} from the \code{parallel} package, 
#' available on GNU/Linux and MacOS. Users of Microsoft Windows can refer to 
#' the README file in the source to be able to use a mclapply type function.
#' 
#' @param X a (n x p) data matrix of predictors. \code{X} must be a matrix. 
#' Each row corresponds to an observation and each column to a 
#' predictor variable.
#' @param Y a (n) vector of (continuous) responses. \code{Y} must be a 
#' vector or a one column matrix. It contains the response variable for 
#' each observation.
#' @param lambda.l1.range a vecor of positive real values, in [0,1]. 
#' \code{lambda.l1} is the sparse penalty parameter for the dimension 
#' reduction step by sparse PLS (see details), the optimal value will be 
#' chosen among \code{lambda.l1.range}.
#' @param ncomp.range a vector of positive integers. \code{ncomp} is the 
#' number of PLS components. The optimal value will be chosen 
#' among \code{ncomp.range}.
#' @param weight.mat a (ntrain x ntrain) matrix used to weight the l2 metric 
#' in the observation space, it can be the covariance inverse of the Ytrain 
#' observations in a heteroskedastic context. If NULL, the l2 metric is the 
#' standard one, corresponding to homoskedastic model (\code{weight.mat} is the 
#' identity matrix).
#' @param adapt a boolean value, indicating whether the sparse PLS selection 
#' step sould be adaptive or not (see details).
#' @param center.X a boolean value indicating whether the data matrices 
#' \code{Xtrain} and \code{Xtest} (if provided) should be centered or not.
#' @param scale.X}{aa boolean value indicating whether the data matrices 
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
#' @param return.grid a boolean values indicating whether the grid of 
#' hyper-parameters values with corresponding mean prediction error rate over 
#' the folds should be returned or not.
#' @param ncores a positve integer, indicating the number of cores that the 
#' cross-validation is allowed to use for parallel computation (see details).
#' @param nfolds a positive integer indicating the number of folds in the 
#' K-folds cross-validation procedure, \code{nfolds=n} corresponds 
#' to the leave-one-out cross-validation, default is 10.
#' @param nrun a positive integer indicating how many times the K-folds cross-
#' validation procedure should be repeated, default is 1.
#' @param verbose a boolean value indicating verbosity.
#' 
#' @return An object with the following attributes
#' \item{lambda.l1.opt}{the optimal value in \code{lambda.l1.range}.}
#' \item{ncomp.opt}{the optimal value in \code{ncomp.range}.}
#' \item{cv.grid}{the grid of hyper-parameters and corresponding prediction 
#' error rate over the folds. 
#' \code{cv.grid} is NULL if \code{return.grid} is set to FALSE.}
#' 
#' @references 
#' Durif G., Modolo L., Michaelsson J., Mold J. E., Lambert-Lacroix S., 
#' Picard F. (2017). High Dimensional Classification with combined Adaptive 
#' Sparse PLS and Logistic Regression, (in prep), 
#' available on (\url{http://arxiv.org/abs/1502.05933}).
#' 
#' @author
#' Ghislain Durif (\url{http://thoth.inrialpes.fr/people/gdurif/}).
#' 
#' @seealso \code{\link{spls}}
#' 
#' @examples
#' \dontrun{
#' ### load plsgenomics library
#' library(plsgenomics)
#' 
#' ### generating data
#' n <- 100
#' p <- 100
#' sample1 <- sample.cont(n=n, p=p, kstar=10, lstar=2, 
#'                        beta.min=0.25, beta.max=0.75, mean.H=0.2, 
#'                        sigma.H=10, sigma.F=5, sigma.E=5)
#'                        
#' X <- sample1$X
#' Y <- sample1$Y
#' 
#' ### hyper-parameters values to test
#' lambda.l1.range <- seq(0.05,0.95,by=0.1) # between 0 and 1
#' ncomp.range <- 1:10
#' 
#' ### tuning the hyper-parameters
#' cv1 <- spls.cv(X=X, Y=Y, lambda.l1.range=lambda.l1.range, 
#'                ncomp.range=ncomp.range, weight.mat=NULL, adapt=TRUE, 
#'                center.X=TRUE, center.Y=TRUE, 
#'                scale.X=TRUE, scale.Y=TRUE, weighted.center=FALSE, 
#'                return.grid=TRUE, ncores=1, nfolds=10, nrun=1)
#' str(cv1)
#' 
#' ### otpimal values
#' cv1$lambda.l1.opt
#' cv1$ncomp.opt
#' }
#' 
#' @export
spls.cv <- function(X, Y, lambda.l1.range, ncomp.range, weight.mat=NULL, 
                    adapt=TRUE, center.X=TRUE, center.Y=TRUE, 
                    scale.X=TRUE, scale.Y=TRUE, weighted.center=FALSE, 
                    return.grid=FALSE, ncores=1, nfolds=10, nrun=1,
                    verbose=FALSE) {
	
	#####################################################################
	#### Initialisation
	#####################################################################
	X <- as.matrix(X)
	n <- nrow(X) # nb observations
	p <- ncol(X) # nb covariates
	index.p <- c(1:p)
	Y <- as.matrix(Y)
	q <- ncol(Y)
	one <- matrix(1,nrow=1,ncol=n)
	
	#####################################################################
	#### Tests on type input
	#####################################################################
	
	# On X
	if ((!is.matrix(X)) || (!is.numeric(X))) {
	     stop("Message from spls.cv: X is not of valid type")
	}
	
	if (p==1) {
	     stop("Message from spls.cv: p=1 is not valid")
	}
	
	# On Y
	if ((!is.matrix(Y)) || (!is.numeric(Y))) {
	     stop("Message from spls.cv: Y is not of valid type")
	}
	
	if (q != 1) {
	     stop("Message from spls.cv: Y must be univariate")
	}
	
	if (nrow(Y)!=n) {
	     stop("Message from spls.cv: the number of observations in Y is not equal to the number of row in X")
	}
	
	# On Y value
	if (sum(is.na(Y))!=0) {
	     stop("Message from spls.cv: NA values in Ytrain")
	}
	
	# On weighting matrix V
	if(!is.null(weight.mat)) { # weighting in scalar product (in observation space of dimension n)
		Vfull <- as.matrix(weight.mat) 
		
		if ((!is.matrix(Vfull)) || (!is.numeric(Vfull))) {
			stop("Message from spls.adapt: Vfull is not of valid type")}
		
		if ((n != ncol(Vfull)) || (n != nrow(Vfull))) {
			stop("Message from spls.adapt: wrong dimension for Vfull, must be a square matrix of size the number of observations in Xtrain")
		}
	} else { # no weighting in sclar product
		Vfull <- diag(rep(1,n), nrow=n, ncol=n)
	}
	
	# On weighted.center
	if ( (weighted.center) && (is.null(weight.mat))) {
		stop("Message from spls.adapt: if the centering is weighted, the weighting matrix V should be provided")
	}
	
	# ncores
	if ((!is.numeric(ncores)) || (round(ncores)-ncores!=0) || (ncores<1)) {
		stop("Message from spls.cv: ncores is not of valid type")
	}
	
	# nfolds
	if ((!is.numeric(nfolds)) || (round(nfolds)-nfolds!=0) || (nfolds<1)) {
		stop("Message from spls.cv: nfolds is not of valid type")
	}
	
	# nrun
	if ((!is.numeric(nrun)) || (round(nrun)-nrun!=0) || (nrun<1)) {
	     stop("Message from spls.cv: nfolds is not of valid type")
	}
	
	# On hyper parameter: lambda.l1
	if ( (sum(!is.numeric(lambda.l1.range))) || (sum(lambda.l1.range<0)) || (sum(lambda.l1.range>1)) ) {
	     stop("Message from spls.adapt: lambda is not of valid type")
	}
	
	# ncomp type
	if ( (sum(!is.numeric(ncomp.range))) || (sum(round(ncomp.range)-ncomp.range!=0)) || (sum(ncomp.range<1)) || (sum(ncomp.range>p)) ) {
	     stop("Message from spls.adapt: ncomp is not of valid type")
	}
	
	#####################################################################
	#### Cross-validation: computation on each fold over the entire grid
	#####################################################################
	
	## the train set is partitioned into nfolds part, each observation is assigned into a fold
	## draw nrun nfolds partition
	folds.obs = sapply(1:nrun, function(run) {
	     return(sample(x=rep(1:nfolds, length.out = n), size=n, replace=FALSE))
	})
	
	## folds x run grid
	folds.grid = expand.grid(k=1:nfolds, run=1:nrun, KEEP.OUT.ATTRS=FALSE)
	
	### normalization of train and test set by fold
	ntrain_values <- matrix(NA, nrow=nfolds, ncol=nrun)
	ntest_values <- matrix(NA, nrow=nfolds, ncol=nrun)
	meanXtrain_values <- list()
	sigmaXtrain_values <- list()
	meanYtrain_values <- list()
	sigmaYtrain_values <- list()
	
	for(index in 1:nrow(folds.grid)) {
	     
	     k = folds.grid$k[index]
	     run = folds.grid$run[index]
	     
	     #### train and test variable
	     Xtrain <- subset(X, folds.obs[,run] != k)
	     Ytrain <- subset(Y, folds.obs[,run] != k)
	     
	     ntrain <- nrow(Xtrain)
	     
	     Xtest <- subset(X, folds.obs[,run] == k)
	     Ytest <- subset(Y, folds.obs[,run] == k)
	     
	     ntest <- nrow(Xtest)
	     
	     V <- Vfull[folds.obs != k, folds.obs != k]
	     
	     if (is.vector(Xtest)==TRUE) {
	          Xtest <- matrix(Xtest,nrow=1)
	     }
	     
	     Xtest <- as.matrix(Xtest)
	     ntest <- nrow(Xtest)
	     
	     DeletedCol <- NULL
	     
	     #####################################################################
	     #### centering and scaling
	     #####################################################################
	     if (!weighted.center) {
	          
	          # Xtrain sd
	          sigmaXtrain <- apply(Xtrain, 2, sd)
	          # predictor with null variance ?
	          if (sum(sigmaXtrain < .Machine$double.eps) !=0) {
                    # predicteur with non null variance < 2 ?
                    if (sum(sigmaXtrain < .Machine$double.eps)>(p-2)){
                         stop("Message from spls.cv: the procedure stops because number of predictor variables with no null variance is less than 1.")
                    }
                    
                    warning("Message from spls.cv: There are covariables with null variance in the current sub-sampling, they will be ignored.")
                    
                    # remove predictor with null variance
                    Xtrain <- Xtrain[,which(sigmaXtrain >= .Machine$double.eps)]
                    Xtest <- Xtest[,which(sigmaXtrain>= .Machine$double.eps)]
                    
                    # list of removed predictors
                    DeletedCol <- index.p[which(sigmaXtrain < .Machine$double.eps)]
                    
                    sigmaXtrain <- sigmaXtrain[-DeletedCol]
	          
	          }
	          
	          # Xtrain mean
	          meanXtrain <- apply(Xtrain, 2, mean)
	          
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
	               sYtrain <- scale( Ytrain, center=meanYtrain, scale=sigmaYtrain)
	          } else if(center.Y && !scale.Y) {
	               sYtrain <- scale( Ytrain, center=meanYtrain, scale=FALSE)
	          } else {
	               sYtrain <- Ytrain
	          }
	          
	          if(center.Y && scale.Y) {
	               sYtest <- scale( Ytest, center=meanYtrain, scale=sigmaYtrain)
	          } else if(center.Y && !scale.Y) {
	               sYtest <- scale( Ytest, center=meanYtrain, scale=FALSE)
	          } else {
	               sYtest <- Ytest
	          }
	          
	          # Xtest	
	          ## centering and scaling depend on Xtest
	          if(center.X && scale.X) {
	               sXtest <- scale( Xtest, center=meanXtrain, scale=sigmaXtrain )
	          } else if(center.X && !scale.X) {
	               sXtest <- scale( Xtest, center=meanXtrain, scale=FALSE )
	          } else {
	               sXtest <- Xtest
	          }
	          
	     } else { # weighted scaling
	          
	          sumV <- sum(diag(V))
	          
	          # X sd
	          # predictor with null variance ?
	          if (sum(sigmaXtrain < .Machine$double.eps) !=0) {
	               # predicteur with non null variance < 2 ?
	               if (sum(sigmaXtrain < .Machine$double.eps)>(p-2)){
	                    stop("Message from spls.cv: the procedure stops because number of predictor variables with no null variance is less than 1.")
	               }
	               
	               warning("Message from spls.cv: There are covariables with nul variance in the current sub-sampling, they will be ignored.")
	               
	               # remove predictor with null variance
	               Xtrain <- Xtrain[,which(sigmaXtrain >= .Machine$double.eps)]
	               Xtest <- Xtest[,which(sigmaXtrain>= .Machine$double.eps)]
	               
	               # list of removed predictors
	               DeletedCol <- index.p[which(sigmaXtrain < .Machine$double.eps)]
	               
	               sigmaXtrain <- sigmaXtrain[-DeletedCol]
	               
	          }
	          
	          # X mean
	          meanXtrain <- matrix(diag(V), nrow=1) %*% Xtrain / sumV
	          
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
	          sYtest <- scale( Ytest, center=meanYtrain, scale=FALSE )
	          
	          # Xtest
	          sXtest <- scale( Xtest, center=meanXtrain, scale=FALSE )
	          
	     }
	     
	     ## store scaled matrices
	     assign(paste0("sXtrain_", k, "_", run), sXtrain)
	     assign(paste0("sXtest_", k, "_", run), sXtest)
	     
	     assign(paste0("sYtrain_", k, "_", run), sYtrain)
	     assign(paste0("sYtest_", k, "_", run), sYtest)
	     
	     assign(paste0("DeletedCol_", k, "_", run), DeletedCol)
	     
	     ntrain_values[k,run] <- ntrain
	     ntest_values[k,run] <- ntest
	     
	     meanXtrain_values[[k + (run-1)*nfolds]] <- meanXtrain
	     sigmaXtrain_values[[k + (run-1)*nfolds]] <- sigmaXtrain
	     
	     meanYtrain_values[[k + (run-1)*nfolds]] <- meanYtrain
	     sigmaYtrain_values[[k + (run-1)*nfolds]] <- sigmaYtrain
	
	}
	
	### K-fold cross-validation grid (fold x run x parameters)
	paramGrid <- expand.grid(fold=1:nfolds, run=1:nrun, lambdaL1=lambda.l1.range, ncomp=ncomp.range, KEEP.OUT.ATTRS=FALSE)
	
	## computations
	res_cv <- Reduce("rbind", mclapply(1:nrow(paramGrid), function(gridRow) {
	     
	     ## GRID: fold, run, lambdaL1, lambdaL2, ncomp
	     k <- paramGrid$fold[gridRow]
	     run <- paramGrid$run[gridRow]
	     
	     ntrain <- ntrain_values[k,run]
	     ntest <- ntest_values[k,run]
	     
	     #### train and test variable
	     Xtrain <- subset(X, folds.obs[,run] != k)
	     Ytrain <- subset(Y, folds.obs[,run] != k)
	     
	     Xtest <- subset(X, folds.obs[,run] == k)
	     Ytest <- subset(Y, folds.obs[,run] == k)
	     
	     V <- Vfull[folds.obs != k, folds.obs != k]
	     
	     if(!is.null(get(paste0("DeletedCol_", k, "_", run)))) {
	          Xtrain <- Xtrain[,-get(paste0("DeletedCol_", k, "_", run))]
	          Xtest <- Xtest[,-get(paste0("DeletedCol_", k, "_", run))]
	     }
          
	     ### computations
	     model <- tryCatch( spls.aux(Xtrain=Xtrain,
	                                 sXtrain=get(paste0("sXtrain_", k, "_", run)),
	                                 Ytrain=Ytrain,
	                                 sYtrain=get(paste0("sYtrain_", k, "_", run)),
	                                 lambda.l1=paramGrid$lambdaL1[gridRow],
	                                 ncomp=paramGrid$ncomp[gridRow],
	                                 weight.mat=V,
	                                 Xtest=Xtest,
	                                 sXtest=get(paste0("sXtest_", k, "_", run)),
	                                 adapt=adapt,
	                                 meanXtrain=meanXtrain_values[[k + (run-1)*nfolds]],
	                                 meanYtrain=meanYtrain_values[[k + (run-1)*nfolds]],
	                                 sigmaXtrain=sigmaXtrain_values[[k + (run-1)*nfolds]],
	                                 sigmaYtrain=sigmaYtrain_values[[k + (run-1)*nfolds]],
	                                 center.X=center.X, center.Y=center.Y,
	                                 scale.X=scale.X, scale.Y=scale.Y,
	                                 weighted.center=weighted.center),
	                        error = function(e) { warnings("Message from spls.cv: error when fitting a model in crossvalidation"); return(NULL);} )
	     
	     ## results
	     res = numeric(6)
	     
	     if(!is.null(model)) {
	          res = c(k, run, paramGrid$lambdaL1[gridRow], paramGrid$ncomp[gridRow], 
	                  model$lenA, sum((model$hatYtest - get(paste0("sYtest_", k, "_", run)))^2) / ntest)
	     } else {
	          res = c(k, run, paramGrid$lambdaL1[gridRow], paramGrid$ncomp[gridRow], 0, NA)
	     }
	     
	     return(res)
	     
	}, mc.cores=ncores, mc.silent=!verbose))
	rownames(res_cv) <- paste(1:nrow(res_cv))
	res_cv = data.frame(res_cv)
	colnames(res_cv) = c("nfold", "nrun", "lambda.l1", "ncomp", "lenA", "error")
	
	
	#####################################################################
	#### Find the optimal point in the grid
	#####################################################################
	
	## compute the number of NAs (=fails)
	cv.grid.fails <- data.frame( as.matrix( with( res_cv, aggregate(error, list(lambda.l1, ncomp), function(x) {sum(is.na(x))}))))
	colnames(cv.grid.fails) <- c("lambda.l1", "ncomp", "nb.fail")
	
	## check number of NAs (=fails)
	if(sum(cv.grid.fails$nb.fail>0.8*nrun*nfolds) > (0.8*nrow(paramGrid))) {
	     warnings("Message from spls.cv: too many errors during the cross-validation process, the grid is not enough filled")
	}     
	
	## compute the mean error over the folds for each point of the grid
	cv.grid.error <- data.frame( as.matrix( with( res_cv, aggregate(error, list(lambda.l1, ncomp), function(x) c(mean(x, na.rm=TRUE), sd(x, na.rm=TRUE)) ))))
	colnames(cv.grid.error) <- c("lambda.l1", "ncomp", "error", "error.sd")
	
	## merge the three tables
	cv.grid <- merge(cv.grid.error, cv.grid.fails, by = c("lambda.l1", "ncomp"))
	
	index_min <- which.min(cv.grid$error)
	lambdaL1.opt <- cv.grid$lambda.l1[index_min]
	ncomp.opt <- cv.grid$ncomp[index_min]
	
	##### return
	if(return.grid) {
	     
	     return( list(lambda.l1.opt=lambdaL1.opt, ncomp.opt=ncomp.opt, cv.grid=cv.grid) )
	     
	} else {
	     
	     return( list(lambda.l1.opt=lambdaL1.opt, ncomp.opt=ncomp.opt, cv.grid=NULL) )
	     
	}
	
}
