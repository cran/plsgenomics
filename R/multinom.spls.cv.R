### multinom.spls.cv.R  (2015-10)
###
###    Tuning parameters (ncomp, lambda.l1, lambda.ridge) for Ridge Iteratively Reweighted 
###    Least Squares followed by Adaptive Sparse PLS regression for multicategorial response, 
###    by K-fold cross-validation
###
### Copyright 2015-10 Ghislain DURIF
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
#' Cross-validation procedure to calibrate the parameters (ncomp, lambda.l1, 
#' lambda.ridge) for the multinomial-SPLS method
#' @aliases multinom.spls.cv
#' 
#' @description 
#' The function \code{multinom.spls.cv} chooses the optimal values for the 
#' hyper-parameter of the \code{multinom.spls} procedure, by minimizing the 
#' averaged error of prediction over the hyper-parameter grid, 
#' using Durif et al. (2018) multinomial-SPLS algorithm.
#' 
#' @details
#' The columns of the data matrices \code{X} may not be standardized, 
#' since standardizing is performed by the function \code{multinom.spls.cv} 
#' as a preliminary step. 
#' 
#' The procedure is described in Durif et al. (2018). The K-fold 
#' cross-validation can be summarize as follow: the train set is partitioned 
#' into K folds, for each value of hyper-parameters the model is fit K times, 
#' using each fold to compute the prediction error rate, and fitting the 
#' model on the remaining observations. The cross-validation procedure returns 
#' the optimal hyper-parameters values, meaning the one that minimize 
#' the averaged error of prediction averaged over all the folds.
#' 
#' This procedures uses \code{mclapply} from the \code{parallel} package, 
#' available on GNU/Linux and MacOS. Users of Microsoft Windows can refer to 
#' the README file in the source to be able to use a mclapply type function.
#' 
#' @param X a (n x p) data matrix of predictors. \code{X} must be a matrix. 
#' Each row corresponds to an observation and each column to a 
#' predictor variable.
#' @param Y a (n) vector of (continuous) responses. \code{Y} must be a 
#' vector or a one column matrix. It contains the response variable for 
#' each observation. \code{Y} should take values in \{0,...,nclass-1\},
#' where nclass is the number of class.
#' @param lambda.ridge.range a vector of positive real values. 
#' \code{lambda.ridge} is the Ridge regularization parameter for the 
#' RIRLS algorithm (see details), the optimal value will be chosen among
#' \code{lambda.ridge.range}.
#' @param lambda.l1.range a vecor of positive real values, in [0,1]. 
#' \code{lambda.l1} is the sparse penalty parameter for the dimension 
#' reduction step by sparse PLS (see details), the optimal value will be 
#' chosen among \code{lambda.l1.range}.
#' @param ncomp.range a vector of positive integers. \code{ncomp} is the 
#' number of PLS components. The optimal value will be chosen 
#' among \code{ncomp.range}.
#' @param adapt a boolean value, indicating whether the sparse PLS selection 
#' step sould be adaptive or not (see details).
#' @param maxIter a positive integer, the maximal number of iterations in the 
#' RIRLS algorithm (see details).
#' @param svd.decompose a boolean parameter. \code{svd.decompose} indicates 
#' wether or not the predictor matrix \code{Xtrain} should be decomposed by 
#' SVD (singular values decomposition) for the RIRLS step (see details).
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
#' @param center.X a boolean value indicating whether the data matrices 
#' \code{Xtrain} and \code{Xtest} (if provided) should be centered or not.
#' @param scale.X a boolean value indicating whether the data matrices 
#' \code{Xtrain} and \code{Xtest} (if provided) should be scaled or not 
#' (\code{scale.X=TRUE} implies \code{center.X=TRUE}) in the spls step.
#' @param weighted.center a boolean value indicating whether the centering 
#' should take into account the weighted l2 metric or not in the SPLS step.
#' @param seed a positive integer value (default is NULL). If non NULL, 
#' the seed for pseudo-random number generation is set accordingly.
#' @param verbose a boolean parameter indicating the verbosity.
#' 
#' @return An object of class \code{multinom.spls} with the following attributes
#' \item{lambda.ridge.opt}{the optimal value in \code{lambda.ridge.range}.}
#' \item{lambda.l1.opt}{the optimal value in \code{lambda.l1.range}.}
#' \item{ncomp.opt}{the optimal value in \code{ncomp.range}.}
#' \item{conv.per}{the overall percentage of models that converge during the 
#' cross-validation procedure.}
#' \item{cv.grid}{the grid of hyper-parameters and corresponding prediction 
#' error rate averaged over the folds. \code{cv.grid} is NULL if 
#' \code{return.grid} is set to FALSE.}
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
#' @seealso \code{\link{multinom.spls}}, \code{\link{multinom.spls.stab}}
#' 
#' @examples
#' \dontrun{
#' ### load plsgenomics library
#' library(plsgenomics)
#' 
#' ### generating data
#' n <- 100
#' p <- 100
#' nclass <- 3
#' sample1 <- sample.multinom(n=n, p=p, nb.class=nclass, kstar=10, lstar=2, 
#'                            beta.min=0.25, beta.max=0.75, mean.H=0.2, 
#'                            sigma.H=10, sigma.F=5)
#' 
#' X <- sample1$X
#' Y <- sample1$Y
#' 
#' ### hyper-parameters values to test
#' lambda.l1.range <- seq(0.05,0.95,by=0.1) # between 0 and 1
#' ncomp.range <- 1:10
#' # log-linear range between 0.01 a,d 1000 for lambda.ridge.range
#' logspace <- function( d1, d2, n) exp(log(10)*seq(d1, d2, length.out=n))
#' lambda.ridge.range <- signif(logspace(d1 <- -2, d2 <- 3, n=21), digits=3)
#' 
#' ### tuning the hyper-parameters
#' cv1 <- multinom.spls.cv(X=X, Y=Y, lambda.ridge.range=lambda.ridge.range, 
#'                         lambda.l1.range=lambda.l1.range, 
#'                         ncomp.range=ncomp.range, 
#'                         adapt=TRUE, maxIter=100, svd.decompose=TRUE, 
#'                         return.grid=TRUE, ncores=1, nfolds=10)
#'                        
#' str(cv1)
#' }
#' 
#' @export
multinom.spls.cv <- function(X, Y, lambda.ridge.range, lambda.l1.range, ncomp.range, adapt=TRUE, maxIter=100, svd.decompose=TRUE, 
                             return.grid=FALSE, ncores=1, nfolds=10, nrun=1, 
                             center.X=TRUE, scale.X=FALSE, weighted.center=TRUE, 
                             seed=NULL, verbose=TRUE) {
     
     #####################################################################
     #### Initialisation
     #####################################################################
     X <- as.matrix(X)
     n <- nrow(X) # nb observations
     p <- ncol(X) # nb covariates
     index.p <- c(1:p)
     if(is.factor(Y)) {
          Y <- as.numeric(levels(Y))[Y]
     }
     Y <- as.integer(Y)
     Y <- as.matrix(Y)
     nclass <- length(unique(Y))
     q <- ncol(Y)
     one <- matrix(1,nrow=1,ncol=n)
     
     if(!is.null(seed)) {
          set.seed(seed)
     }
     
     
     #####################################################################
     #### Tests on type input
     #####################################################################
     
     # if Binary response
     if(length(table(Y)) == 2) {
          warning("message from multinom.spls.cv: binary response")
          results = multinom.spls.cv(X=X, Y=Y, lambda.ridge.range=lambda.ridge.range, 
                                  lambda.l1.range=lambda.l1.range, ncomp.range=ncomp.range, 
                                  adapt=adapt, maxIter=maxIter, svd.decompose=svd.decompose, 
                                  return.grid=return.grid, ncores=ncores, 
                                  nfolds=nfolds, nrun=nrun, 
                                  center.X=center.X, scale.X=scale.X, 
                                  weighted.center=weighted.center, 
                                  seed=seed, verbose=verbose)
          return(results)
     }
     
     # On X
     if ((!is.matrix(X)) || (!is.numeric(X))) {
          stop("Message from multinom.spls.cv: X is not of valid type")
     }
     
     if (p==1) {
          # stop("Message from multinom.spls.cv: p=1 is not valid")
          warning("Message from multinom.spls.cv: p=1 is not valid, ncomp.range is set to 0")
          ncomp.range <- 0
     }
     
     # On Y
     if ((!is.matrix(Y)) || (!is.numeric(Y))) {
          stop("Message from multinom.spls.cv: Y is not of valid type")
     }
     
     if (q != 1) {
          stop("Message from multinom.spls.cv: Y must be univariate")
     }
     
     if (nrow(Y)!=n) {
          stop("Message from multinom.spls.cv: the number of observations in Y is not equal to the Xtrain row number")
     }
     
     # On Y value
     if (sum(is.na(Y))!=0) {
          stop("Message from multinom.spls.cv: NA values in Y")
     }
     
     if((sum(floor(Y)-Y)!=0)||(sum(Y<0)>0)) {
          stop("Message from multinom.spls.cv: Y is not of valid type")
     }
     
     if(any(!Y %in% c(0:(nclass-1)))) {
          stop("Message from multinom.spls.cv: Y should be in {0,...,nclass-1}")
     }
     
     if (sum(as.numeric(table(Y))==0)!=0) {
          stop("Message from multinom.spls.cv: there are empty classes")
     }
     
     # On hyper parameter: lambda.ridge, lambda.l1
     if ( any(!is.numeric(lambda.ridge.range)) || any(lambda.ridge.range<0) || any(!is.numeric(lambda.l1.range)) || any(lambda.l1.range<0)) {
          stop("Message from multinom.spls: lambda is not of valid type")
     }
     
     # ncomp type
     if ( any(!is.numeric(ncomp.range)) || any(round(ncomp.range)-ncomp.range!=0) || any(ncomp.range<0) || any(ncomp.range>p)) {
          stop("Message from multinom.spls: ncomp.range is not of valid type")
     }
     
     
     # maxIter
     if ((!is.numeric(maxIter)) || (round(maxIter)-maxIter!=0) || (maxIter<1)) {
          stop("Message from multinom.spls.cv: maxIter is not of valid type")
     }
     
     # ncores
     if ((!is.numeric(ncores)) || (round(ncores)-ncores!=0) || (ncores<1)) {
          stop("Message from multinom.spls.cv: ncores is not of valid type")
     }
     
     # nfolds
     if ((!is.numeric(nfolds)) || (round(nfolds)-nfolds!=0) || (nfolds<1)) {
          stop("Message from multinom.spls.cv: nfolds is not of valid type")
     }
     
     # necessary to insure that both classes are represented in each folds
     if( any(as.vector(table(Y))<nfolds) ) {
          stop("Message from multinom.spls.cv: there is a class defined by Y that has less members than the number of folds nfold")
     }
     
     ## fold size
     fold.size = n %/% nfolds
     if( length(unique(Y))>fold.size ) {
          stop("Message from multinom.spls.cv: there are too many classes compared to the size of folds")
     }
     
     
     #####################################################################
     #### Cross-validation: computation on each fold over the entire grid
     #####################################################################
     
     ## draw nrun nfolds partition
     folds.obs = sapply(1:nrun, function(run) {
          ## the train set is partitioned into nfolds part, each observation is assigned into a fold
          
          # observations from the different classes
          classes <- sort(unique(Y))
          index <- lapply(classes, function(g) {
               return((1:n)[Y==g])
          })
          
          # each folds contains at least one elements of each class
          # choose the represent of each class in each fold
          folds <- lapply(classes+1, function(g) {
               return(c(1:nfolds, rep(0, length(index[[g]]) - nfolds))[sample(x=1:length(index[[g]]), size=length(index[[g]]), replace=FALSE)])
          })
          
          ind.folds <- lapply(classes+1, function(g) {
               return((index[[g]])[folds[[g]]!=0])
          })
          
          # group the remaining observation and equally dispatch them into the folds
          rest <- unlist(lapply(classes+1, function(g) {
               return((index[[g]])[folds[[g]]==0])
          }))
          
          folds.rest = rep(1:nfolds, length.out = length(rest))[sample(x=1:length(rest), size=length(rest), replace=FALSE)]
          
          folds.res = rep(NA, n)
          
          for(g in (classes+1)) {
               folds.res[ind.folds[[g]]] <- (folds[[g]])[folds[[g]]!=0]
          }
          folds.res[rest] = folds.rest
          
          #table(folds.res)
          #sapply(1:nfolds, function(x) table(Y[folds.res==x]))
          
          return(folds.res)
     })
     #for(i in 1:nfolds) print(table(Y[folds.obs[,1]==i,]))
     
     
     
     ## folds x run grid
     folds.grid = expand.grid(k=1:nfolds, run=1:nrun, KEEP.OUT.ATTRS=FALSE)
     
     ### normalization of train and test set by fold
     ntrain_values <- matrix(NA, nrow=nfolds, ncol=nrun)
     ntest_values <- matrix(NA, nrow=nfolds, ncol=nrun)
     meanXtrain_values <- list()
     sigma2train_values <- list()
     
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
          
          #### Move into the reduced space
          
          r <- p #min(p, ntrain)
          DeletedCol <- NULL
          
          ### Standardize the Xtrain matrix
          # standard deviation (biased one) of Xtrain
          sigma2train <- apply(Xtrain, 2, var) * (ntrain-1)/(ntrain)
          
          # test on sigma2train
          # predictor with null variance ?
          if (sum(sigma2train < .Machine$double.eps)!=0){
               
               # predicteur with non null variance < 2 ?
               if (sum(sigma2train < .Machine$double.eps)>(p-2)){
                    stop("Message from multinom.spls.cv: the procedure stops because number of predictor variables with no null variance is less than 1.")
               }
               
               warning("Message from multinom.spls.cv: There are covariables with nul variance in the current sub-sampling, they will be ignored.")
               
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
               r <- p #min(p,ntrain)
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
          
          sXtrain.nosvd = sXtrain # keep in memory if svd decomposition
          
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
          
          # center and scale Xtest	
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
          
          ## store scaled matrices
          assign(paste0("sXtrain_", k, "_", run), sXtrain)
          assign(paste0("sXtest_", k, "_", run), sXtest)
          
          assign(paste0("sXtrain.nosvd_", k, "_", run), sXtrain.nosvd)
          assign(paste0("sXtest.nosvd_", k, "_", run), sXtest.nosvd)
          
          ntrain_values[k,run] <- ntrain
          ntest_values[k,run] <- ntest
          
          meanXtrain_values[[k + (run-1)*nfolds]] <- meanXtrain
          sigma2train_values[[k + (run-1)*nfolds]] <- sigma2train
          
          
     }
     
     ### K-fold cross-validation grid (fold x run x parameters)
     paramGrid <- expand.grid(fold=1:nfolds, run=1:nrun, lambdaL1=lambda.l1.range, lambdaL2=lambda.ridge.range, ncomp=ncomp.range, KEEP.OUT.ATTRS=FALSE)
     
     ## computations
     res_cv <- Reduce("rbind", mclapply(split(paramGrid, f=row.names(paramGrid)), function(gridRow) {
          
          ## GRID: fold, run, lambdaL1, lambdaL2, ncomp
          k <- gridRow$fold
          run <- gridRow$run
          
          ntrain <- ntrain_values[k,run]
          ntest <- ntest_values[k,run]
          
          Ytrain <- subset(Y, folds.obs[,run] != k)
          Ytest <- subset(Y, folds.obs[,run] == k)
          
          ### computations
          model <- tryCatch( multinom.spls.aux(sXtrain=get(paste0("sXtrain_", k, "_", run)), 
                                              sXtrain.nosvd=get(paste0("sXtrain.nosvd_", k, "_", run)), 
                                              Ytrain=Ytrain, lambda.ridge=gridRow$lambdaL2, 
                                              lambda.l1=gridRow$lambdaL1, ncomp=gridRow$ncomp, 
                                              sXtest=get(paste0("sXtest_", k, "_", run)), 
                                              sXtest.nosvd=get(paste0("sXtest.nosvd_", k, "_", run)), 
                                              adapt=adapt, maxIter=maxIter, svd.decompose=svd.decompose, 
                                              meanXtrain=meanXtrain_values[[k + (run-1)*nfolds]], 
                                              sigma2train=sigma2train_values[[k + (run-1)*nfolds]], 
                                              center.X=center.X, scale.X=scale.X, weighted.center=weighted.center), 
                             error = function(e) { print(e); warnings("Message from multinom.spls.cv: error when fitting a model in crossvalidation"); return(NULL);} )
          
          
          ## results
          res = numeric(8)
          
          if(!is.null(model)) {
               res = c(gridRow$fold, gridRow$run, gridRow$lambdaL1, gridRow$lambdaL2, gridRow$ncomp, model$lenA, model$converged, sum(model$hatYtest != Ytest) / ntest)
          } else {
               res = c(gridRow$fold, gridRow$run, gridRow$lambdaL1, gridRow$lambdaL2, gridRow$ncomp, 0, NA, NA)
          }
          
          return(res)
          
     }, mc.cores=ncores, mc.silent=!verbose))
     rownames(res_cv) <- paste(1:nrow(res_cv))
     res_cv = data.frame(res_cv)
     colnames(res_cv) = c("nfold", "nrun", "lambda.l1", "lambda.ridge", "ncomp", "lenA", "converged", "error")
     
     
     #####################################################################
     #### Find the optimal point in the grid
     #####################################################################
     
     ## compute the number of NAs (=fails)
     cv.grid.fails <- data.frame( as.matrix( with( res_cv, aggregate(error, list(lambda.ridge, lambda.l1, ncomp), function(x) {sum(is.na(x))}))))
     colnames(cv.grid.fails) <- c("lambda.ridge", "lambda.l1", "ncomp", "nb.fail")
     
     ## check number of NAs (=fails)
     if(sum(cv.grid.fails$nb.fail>0.8*nrun*nfolds) > (0.8*nrow(paramGrid))) {
          warnings("Message from multinom.spls.cv: too many errors during the cross-validation process, the grid is not enough filled")
     }     
     
     ## compute the mean error over the folds for each point of the grid
     cv.grid.error <- data.frame( as.matrix( with( res_cv, aggregate(error, list(lambda.ridge, lambda.l1, ncomp), function(x) c(mean(x, na.rm=TRUE), sd(x, na.rm=TRUE)) ))))
     colnames(cv.grid.error) <- c("lambda.ridge", "lambda.l1", "ncomp", "error", "error.sd")
     
     ## compute the percentage of convergence over the folds for each point of the grid
     cv.grid.conv <- data.frame( as.matrix( with( res_cv, aggregate(converged, list(lambda.ridge, lambda.l1, ncomp), mean, na.rm=TRUE))))
     colnames(cv.grid.conv) <- c("lambda.ridge", "lambda.l1", "ncomp", "converged")
     
     ## % of convergence over all cross-validation process
     conv.per <- mean(cv.grid.conv$converged)
     
     ## merge the three tables
     cv.grid1 <- merge(cv.grid.error, cv.grid.conv, by = c("lambda.ridge", "lambda.l1", "ncomp"))
     cv.grid <- merge(cv.grid1, cv.grid.fails, by = c("lambda.ridge", "lambda.l1", "ncomp"))
     
     index_min <- which.min(cv.grid$error)
     lambdaL1.opt <- cv.grid$lambda.l1[index_min]
     lambdaL2.opt <- cv.grid$lambda.ridge[index_min]
     ncomp.opt <- cv.grid$ncomp[index_min]
     
     ##### return
     if(return.grid) {
          
          return( list(lambda.ridge.opt=lambdaL2.opt, lambda.l1.opt=lambdaL1.opt, ncomp.opt=ncomp.opt, conv.per=conv.per, cv.grid=cv.grid) )
          
     } else {
          
          return( list(lambda.ridge.opt=lambdaL2.opt, lambda.l1.opt=lambdaL1.opt, ncomp.opt=ncomp.opt, conv.per=conv.per, cv.grid=NULL) )
          
     }
     
     
}

