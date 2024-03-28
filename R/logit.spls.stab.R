### logit.spls.stab.R  (2015-10)
###
###    Stability selection procedure for logit.spls
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
#' Stability selection procedure to estimate probabilities of selection of 
#' covariates for the LOGIT-SPLS method
#' @aliases logit.spls.stab
#' 
#' @description 
#' The function \code{logit.spls.stab} train a logit-spls model for each 
#' candidate values \code{(ncomp, lambda.l1, lambda.ridge)} of hyper-parameters 
#' on multiple sub-samplings in the data. The stability selection procedure 
#' selects the covariates that are selected by most of the models among the 
#' grid of hyper-parameters, following the procedure described in 
#' Durif et al. (2018). Candidates values for \code{ncomp}, \code{lambda.l1} 
#' and \code{lambda.l2} are respectively given by 
#' the input arguments \code{ncomp.range}, \code{lambda.l1.range} 
#' and \code{lambda.l2.range}.
#' 
#' 
#' @details
#' The columns of the data matrices \code{X} may not be standardized, 
#' since standardizing is performed by the function \code{logit.spls.stab} 
#' as a preliminary step. 
#' 
#' The procedure is described in Durif et al. (2018). The stability selection 
#' procedure can be summarize as follow (c.f. Meinshausen and Buhlmann, 2010).
#' 
#' (i) For each candidate values \code{(ncomp, lambda.l1, lambda.ridge)} of 
#' hyper-parameters, a logit-SPLS is trained on \code{nresamp} resamplings 
#' of the data. Then, for each triplet \code{(ncomp, lambda.l1, lambda.ridge)}, 
#' the probability that a covariate (i.e. a column in \code{X}) is selected is 
#' computed among the resamplings.
#' 
#' (ii) Eventually, the set of "stable selected" variables corresponds to the 
#' set of covariates that were selected by most of the training among the 
#' grid of hyper-parameters candidate values.
#' 
#' This function achieves the first step (i) of the stability selection 
#' procedure. The second step (ii) is achieved by the function 
#' \code{\link{stability.selection}}.
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
#' each observation. \code{Y} should take values in \{0,1\}.
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
#' @param ncores a positve integer, indicating the number of cores that the 
#' cross-validation is allowed to use for parallel computation (see details).
#' @param nresamp number of resamplings of the data to estimate the probility 
#' of selection for each covariate, default is 100.
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
#' @return An object with the following attributes
#' \item{q.Lambda}{A table with values of q.Lambda (c.f. Durif 
#' et al. (2018) for the notation), being the averaged number of covariates
#' selected among the entire grid of hyper-parameters candidates values,
#' for increasing size of hyper-parameter grid.}
#' \item{probs.lambda}{A table with estimated probability of selection for each 
#' covariates depending on the candidates values for hyper-parameters.}
#' \item{p}{An integer values indicating the number of covariates in the 
#' model.}
#' 
#' @references 
#' Durif, G., Modolo, L., Michaelsson, J., Mold, J.E., Lambert-Lacroix, S., 
#' Picard, F., 2018. High dimensional classification with combined 
#' adaptive sparse PLS and logistic regression. Bioinformatics 34, 
#' 485--493. \doi{10.1093/bioinformatics/btx571}.
#' Available at \url{http://arxiv.org/abs/1502.05933}.
#' 
#' Meinshausen, N., Buhlmann P. (2010). Stability Selection. Journal of the 
#' Royal Statistical Society: Series B (Statistical Methodology) 
#' 72, no. 4, 417-473.
#' 
#' @author
#' Ghislain Durif (\url{https://gdurif.perso.math.cnrs.fr/}).
#' 
#' @seealso \code{\link{logit.spls}}, \code{\link{stability.selection}}, 
#' \code{\link{stability.selection.heatmap}}
#' 
#' @examples
#' \dontrun{
#' ### load plsgenomics library
#' library(plsgenomics)
#' 
#' ### generating data
#' n <- 100
#' p <- 100
#' sample1 <- sample.bin(n=n, p=p, kstar=10, lstar=2, 
#'                       beta.min=0.25, beta.max=0.75, mean.H=0.2, 
#'                       sigma.H=10, sigma.F=5)
#' 
#' X <- sample1$X
#' Y <- sample1$Y
#' 
#' ### pertinent covariates id
#' sample1$sel
#' 
#' ### hyper-parameters values to test
#' lambda.l1.range <- seq(0.05,0.95,by=0.1) # between 0 and 1
#' ncomp.range <- 1:10
#' # log-linear range between 0.01 a,d 1000 for lambda.ridge.range
#' logspace <- function( d1, d2, n) exp(log(10)*seq(d1, d2, length.out=n))
#' lambda.ridge.range <- signif(logspace(d1 <- -2, d2 <- 3, n=21), digits=3)
#' 
#' ### tuning the hyper-parameters
#' stab1 <- logit.spls.stab(X=X, Y=Y, lambda.ridge.range=lambda.ridge.range, 
#'                          lambda.l1.range=lambda.l1.range, 
#'                          ncomp.range=ncomp.range, 
#'                          adapt=TRUE, maxIter=100, svd.decompose=TRUE, 
#'                          ncores=1, nresamp=100)
#'                        
#' str(stab1)
#' 
#' ### heatmap of estimated probabilities
#' stability.selection.heatmap(stab1)
#' 
#' ### selected covariates
#' stability.selection(stab1, piThreshold=0.6, rhoError=10)
#' }
#' 
#' @export
logit.spls.stab <- function(X, Y, lambda.ridge.range, lambda.l1.range, 
                            ncomp.range, 
                            adapt=TRUE, maxIter=100, svd.decompose=TRUE, 
                            ncores=1, nresamp=100, 
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
     q <- ncol(Y)
     one <- matrix(1,nrow=1,ncol=n)
     
     if(!is.null(seed)) {
          set.seed(seed)
     }
     
     cnames <- NULL
     if(!is.null(colnames(X))) {
          cnames <- colnames(X)
     } else {
          cnames <- paste0(1:p)
     }
     
     
     #####################################################################
     #### Tests on type input
     #####################################################################
     
     # if multicategorical response
     if(length(table(Y)) > 2) {
          warning("message from logit.spls.stab: multicategorical response")
          results <- multinom.spls.stab(X=X, Y=Y, lambda.ridge.range=lambda.ridge.range, lambda.l1.range=lambda.l1.range,
                                        ncomp.range=ncomp.range,
                                        adapt=adapt, maxIter=maxIter, svd.decompose=svd.decompose,
                                        ncores=ncores, nresamp=nresamp,
                                        center.X=center.X, scale.X=scale.X, weighted.center=weighted.center,
                                        seed=seed, verbose=verbose)
          return(results)
     }
     
     # On X
     if ((!is.matrix(X)) || (!is.numeric(X))) {
          stop("message from logit.spls.stab: X is not of valid type")
     }
     
     if (p==1) {
          stop("message from logit.spls.stab: p=1 is not valid")
     }
     
     # On Y
     if ((!is.matrix(Y)) || (!is.numeric(Y))) {
          stop("message from logit.spls.stab: Y is not of valid type")
     }
     
     if (q != 1) {
          stop("message from logit.spls.stab: Y must be univariate")
     }
     
     if (nrow(Y)!=n) {
          stop("message from logit.spls.stab: the number of observations in Y is not equal to the number of row in X")
     }
     
     # On Y value
     if (sum(is.na(Y))!=0) {
          stop("message from logit.spls.stab: NA values in Ytrain")
     }
     
     if (sum(!(Y %in% c(0,1)))!=0) {
          stop("message from logit.spls.stab: Y is not of valid type")
     }
     
     if (sum(as.numeric(table(Y))==0)!=0) {
          stop("message from logit.spls.stab: there are empty classes")
     }
     
     # On hyper parameter: lambda.ridge, lambda.l1
     if (any(!is.numeric(lambda.ridge.range)) || any(lambda.ridge.range<0) 
         || any(!is.numeric(lambda.l1.range)) || any(lambda.l1.range<0) || any(lambda.l1.range>1)) {
          stop("Message from logit.spls.stab: lambda is not of valid type")
     }
     
     # ncomp type
     if (any(!is.numeric(ncomp.range)) || any(round(ncomp.range)-ncomp.range!=0) || any(ncomp.range<0) || any(ncomp.range>p)) {
          stop("Message from logit.spls.stab: ncomp is not of valid type")
     }
     
     # maxIter
     if ((!is.numeric(maxIter)) || (round(maxIter)-maxIter!=0) || (maxIter<1)) {
          stop("message from logit.spls.stab: maxIter is not of valid type")
     }
     
     # ncores
     if ((!is.numeric(ncores)) || (round(ncores)-ncores!=0) || (ncores<1)) {
          stop("message from logit.spls.stab: ncores is not of valid type")
     }
     
     # nresamp
     if ((!is.numeric(nresamp)) || (round(nresamp)-nresamp!=0) || (nresamp<1)) {
          stop("message from logit.spls.stab: nresamp is not of valid type")
     }
     
     
     #####################################################################
     #### Stability selection procedure
     #####################################################################
     
     ## computation on the folds x run grid
     grid.resampling <- as.matrix( Reduce("rbind", mclapply(1:nresamp, function(id.samp) {
          
          #### train and test variable
          ntrain = floor(0.5*n)
          
          index.train = sort(sample(1:n, size=ntrain))
          
          Xtrain = X[index.train,]
          Ytrain = Y[index.train]
          
          condition = length(table(Ytrain))<2
          test = 0
          while(condition & test<100) {
               index.train = sort(sample(1:n, size=ntrain))
               
               Xtrain = X[index.train,]
               Ytrain = Y[index.train]
               
               condition = length(table(Ytrain))<2
               test = test+1
          }
          
          if(test==100) {
               ## at least one observation of each class in the train set (if random splitting did not work)
               ind0 = which(Y==0)
               ind1 = which(Y==1)
               index.train = sample(ind0, size=1)
               index.train = c(index.train, sample(ind1, size=1))
               index.train = c(index.train, sample((1:n)[which(!(1:n) %in% index.train)], size=ntrain-2))
               
               Xtrain = X[index.train,]
               Ytrain = Y[index.train]
          }
          
          #### hyper-parameter grid
          paramGrid <- expand.grid(lambdaL1=lambda.l1.range, 
                                   lambdaL2=lambda.ridge.range, 
                                   ncomp=ncomp.range, KEEP.OUT.ATTRS=FALSE)
          
          #### fit the model for the different lambda.l1
          grid_out <- as.matrix( Reduce("rbind", lapply(1:nrow(paramGrid), function(gridRow) {
               
               lambdaL1 <- paramGrid$lambdaL1[gridRow]
               lambdaL2 <- paramGrid$lambdaL2[gridRow]
               ncomp <- paramGrid$ncomp[gridRow]
               
               ## fit the model for the chosen lambda ridge
               fit_out <- logit.spls(Xtrain=Xtrain, Ytrain=Ytrain, 
                                     lambda.ridge=lambdaL2, 
                                     lambda.l1=lambdaL1, 
                                     ncomp=ncomp, Xtest=NULL, adapt=adapt, 
                                     maxIter=maxIter, svd.decompose=svd.decompose, 
                                     center.X=center.X, scale.X=scale.X, 
                                     weighted.center=weighted.center)
               
               ## selected variables ?
               sel_var <- fit_out$Anames
               status_var <- rep(0, length(cnames))
               status_var[which(cnames %in% sel_var)] <- rep(1, length(which(cnames %in% sel_var)))
               
               tmp <- c(lambdaL1, lambdaL2, ncomp, id.samp, sum(status_var), status_var)
               return(tmp)
               
          })))
          
          rownames(grid_out) <- NULL
          return(grid_out)
          
     }, mc.cores=ncores, mc.silent=!verbose)))
     
     grid.resampling <- data.frame(grid.resampling)
     colnames(grid.resampling) <- c("lambdaL1", "lambdaL2", "ncomp", 
                                    "id", "nbVar", cnames)
     
     grid.resampling$point <- paste0(grid.resampling$lambdaL1, "_",
                                     grid.resampling$lambdaL2, "_",
                                     grid.resampling$ncomp)
     
     o.grid <- order(grid.resampling$nbVar)
     grid.resampling <- grid.resampling[o.grid,]
     
     if(any(table(grid.resampling$point)<nresamp)) {
          warning("message from logit.spls.stab: empty classe in a resampling")
          print(table(grid.resampling$point))
     }
     
     #####################################################################
     #### Compute q_lambda
     #####################################################################
     
     ## increasing value of q_lambda
     tmp_qLambda <- as.matrix( Reduce("rbind", mclapply(1:nresamp, function(id.samp) {
          
          tmp1 <- subset(grid.resampling, grid.resampling$id==id.samp)
          tmp2 <- apply(t(tmp1)[-c(1:5,tail(1:ncol(grid.resampling), 1)),],1,cumsum) # cumsum by genes
          tmp3 <- apply(t(tmp2), 2, function(x) return(sum(x!=0)))
          
          tmp4 <- cbind(tmp1[,c(1:4)], unname(tmp3))
          
          return(tmp4)
          
     }, mc.cores=ncores, mc.silent=!verbose)))
     
     tmp_qLambda <- data.frame(tmp_qLambda)
     colnames(tmp_qLambda) <- c("lambdaL1", "lambdaL2", "ncomp", "id.samp", "qLambda")
     
     qLambda <- ddply(tmp_qLambda, c("lambdaL1", "lambdaL2", "ncomp"), 
                      function(x) colMeans(x[c("qLambda")], na.rm=TRUE))
     o.qLambda <- order(qLambda$qLambda)
     qLambda <- qLambda[o.qLambda,]
     
     #####################################################################
     #### Compute p_j(lambda)
     #####################################################################
     
     probs_lambda <- ddply(grid.resampling, c("lambdaL1", "lambdaL2", "ncomp"), 
                           function(x) colMeans(x[cnames], na.rm=TRUE))
     probs_lambda <- probs_lambda[o.qLambda,]
     
     # select variables in another function
     
     
     ##### return
     return(list(q.Lambda=qLambda, probs.lambda=probs_lambda, p=p))
     
}


