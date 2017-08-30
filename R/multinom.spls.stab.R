### multinom.spls.stab.R  (2015-10)
###
###    Stability selection procedure for multinomial-spls
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
#' Stability selection procedure  to estimate probabilities of selection of
#' covariates for the multinomial-SPLS method
#' @aliases multinom.spls.stab
#' 
#' @description 
#' The function \code{multinom.spls.stab} train a multinomial-spls model for 
#' each candidate values \code{(ncomp, lambda.l1, lambda.ridge)} of 
#' hyper-parameters on multiple sub-samplings in the data. The stability 
#' selection procedure selects the covariates that are selected by most of the 
#' models among the grid of hyper-parameters, following the procedure 
#' described in Durif et al.  (2017). Candidates values for \code{ncomp}, 
#' \code{lambda.l1} and \code{lambda.l2} are respectively given by 
#' the input arguments \code{ncomp.range}, \code{lambda.l1.range} 
#' and \code{lambda.l2.range}.
#' 
#' 
#' @details
#' The columns of the data matrices \code{X} may not be standardized, 
#' since standardizing is performed by the function \code{multinom.spls.stab} 
#' as a preliminary step. 
#' 
#' The procedure is described in Durif et al. (2017). The stability selection 
#' procedure can be summarize as follow (c.f. Meinshausen and Buhlmann, 2010).
#' 
#' (i) For each candidate values \code{(ncomp, lambda.l1, lambda.ridge)} of 
#' hyper-parameters, a multinomial-spls is trained on \code{nresamp} 
#' resamplings of the data. Then, for each triplet 
#' \code{(ncomp, lambda.l1, lambda.ridge)}, the probability that a covariate 
#' (i.e. a column in \code{X}) is selected is computed among the resamplings.
#' 
#' The estimated probabilities can be visualized as a heatmap with the 
#' function \code{\link{stability.selection.heatmap}}.
#' 
#' (ii) Eventually, the set of "stable selected" variables corresponds to the 
#' set of covariates that were selected by most of the training among the 
#' grid of hyper-parameters candidate values.
#' 
#' This function achieves the first step (i) of the stability selection 
#' procedure. The second step (ii) is achieved by the function 
#' \code{\link{stability.selection}}
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
#' each observation.  \code{Y} should take values in \{0,...,nclass-1\}, 
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
#' et al. (2017) for the notation), being the averaged number of covariates
#' selected among the entire grid of hyper-parameters candidates values,
#' for increasing size of hyper-parameter grid.}
#' \item{probs.lambda}{A table with estimated probability of selection for each 
#' covariates depending on the candidates values for hyper-parameters.}
#' \item{p}{An integer values indicating the number of covariates in the 
#' model.}
#' 
#' @references 
#' Durif G., Modolo L., Michaelsson J., Mold J. E., Lambert-Lacroix S., 
#' Picard F. (2017). High Dimensional Classification with combined Adaptive 
#' Sparse PLS and Logistic Regression, (in prep), 
#' available on (\url{http://arxiv.org/abs/1502.05933}).
#' 
#' Meinshausen, N., Buhlmann P. (2010). Stability Selection. Journal of the 
#' Royal Statistical Society: Series B (Statistical Methodology) 
#' 72, no. 4, 417-473.
#' 
#' @author
#' Ghislain Durif (\url{http://thoth.inrialpes.fr/people/gdurif/}).
#' 
#' @seealso \code{\link{multinom.spls}}, \code{\link{stability.selection}}, 
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
#' nclass <- 3
#' sample1 <- sample.multinom(n, p, nb.class=nclass, kstar=20, lstar=2, 
#'                            beta.min=0.25, beta.max=0.75, 
#'                            mean.H=0.2, sigma.H=10, sigma.F=5)
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
#' stab1 <- multinom.spls.stab(X=X, Y=Y, lambda.ridge.range=lambda.ridge.range, 
#'                             lambda.l1.range=lambda.l1.range, 
#'                             ncomp.range=ncomp.range, 
#'                             adapt=TRUE, maxIter=100, svd.decompose=TRUE, 
#'                             ncores=1, nresamp=100)
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
multinom.spls.stab <- function(X, Y, lambda.ridge.range, lambda.l1.range, 
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
     nclass <- length(unique(Y))
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
     
     # if binary response
     if(length(table(Y)) < 3) {
          warning("message from multinom.spls.stab: multicategorical response")
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
          stop("message from multinom.spls.stab: X is not of valid type")
     }
     
     if (p==1) {
          stop("message from multinom.spls.stab: p=1 is not valid")
     }
     
     # On Y
     if ((!is.matrix(Y)) || (!is.numeric(Y))) {
          stop("Message from multinom.spls.stab: Y is not of valid type")
     }
     
     if (q != 1) {
          stop("Message from multinom.spls.stab: Y must be univariate")
     }
     
     if (nrow(Y)!=n) {
          stop("Message from multinom.spls.stab: the number of observations in Y is not equal to the Xtrain row number")
     }
     
     # On Y value
     if (sum(is.na(Y))!=0) {
          stop("Message from multinom.spls.stab: NA values in Y")
     }
     
     if((sum(floor(Y)-Y)!=0)||(sum(Y<0)>0)) {
          stop("Message from multinom.spls.stab: Y is not of valid type")
     }
     
     if(any(!Y %in% c(0:(nclass-1)))) {
          stop("Message from multinom.spls.cv: Y should be in {0,...,nclass-1}")
     }
     
     if (sum(as.numeric(table(Y))==0)!=0) {
          stop("Message from multinom.spls.stab: there are empty classes")
     }
     
     
     # On hyper parameter: lambda.ridge, lambda.l1
     if (any(!is.numeric(lambda.ridge.range)) || any(lambda.ridge.range<0) 
         || any(!is.numeric(lambda.l1.range)) || any(lambda.l1.range<0) || any(lambda.l1.range>1)) {
          stop("Message from multinom.spls.stab: lambda is not of valid type")
     }
     
     # ncomp type
     if (any(!is.numeric(ncomp.range)) || any(round(ncomp.range)-ncomp.range!=0) || any(ncomp.range<0) || any(ncomp.range>p)) {
          stop("Message from multinom.spls.stab: ncomp is not of valid type")
     }
     
     # maxIter
     if ((!is.numeric(maxIter)) || (round(maxIter)-maxIter!=0) || (maxIter<1)) {
          stop("message from multinom.spls.stab: maxIter is not of valid type")
     }
     
     # ncores
     if ((!is.numeric(ncores)) || (round(ncores)-ncores!=0) || (ncores<1)) {
          stop("message from multinom.spls.stab: ncores is not of valid type")
     }
     
     # nresamp
     if ((!is.numeric(nresamp)) || (round(nresamp)-nresamp!=0) || (nresamp<1)) {
          stop("message from multinom.spls.stab: nresamp is not of valid type")
     }
     
     
     #####################################################################
     #### Stability selection procedure
     #####################################################################
     
     print("test 1")
     
     ## computation on the folds x run grid
     grid.resampling <- as.matrix( Reduce("rbind", mclapply(1:nresamp, function(id.samp) {
          
          #### train and test variable
          ntrain = floor(0.5*n)
          
          index.train = sort(sample(1:n, size=ntrain))
          
          Xtrain = X[index.train,]
          Ytrain = Y[index.train]
          
          condition = length(table(Ytrain))<nclass
          test = 0
          while(condition & test<100) {
               index.train = sort(sample(1:n, size=ntrain))
               
               Xtrain = X[index.train,]
               Ytrain = Y[index.train]
               
               condition = length(table(Ytrain))<nclass
               test = test+1
          }
          
          if(test==100) {
               ## at least one observation of each class in the train set (if random splitting did not work)
               index.train <- sapply(0:(nclass-1), function(ind_class) {
                    return(sample(as.vector(which(Y==ind_class)), size=1, replace=FALSE))
               })
               index.train = c(index.train, sample((1:n)[which(!(1:n) %in% index.train)], size=ntrain-nclass))
               
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
               fit_out <- multinom.spls(Xtrain=Xtrain, Ytrain=Ytrain, 
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
     
     print("test 2")
     
     grid.resampling <- data.frame(grid.resampling)
     colnames(grid.resampling) <- c("lambdaL1", "lambdaL2", "ncomp", 
                                    "id", "nbVar", cnames)
     
     grid.resampling$point = paste0(grid.resampling$lambdaL1, "_",
                                    grid.resampling$lambdaL2, "_",
                                    grid.resampling$ncomp)
     
     o.grid <- order(grid.resampling$nbVar)
     grid.resampling <- grid.resampling[o.grid,]
     
     if(any(table(grid.resampling$point)<nresamp)) {
          warning("message from multinom.spls.stab: empty classe in a resampling")
          print(table(grid.resampling$point))
     }
     
     #####################################################################
     #### Compute q_lambda
     #####################################################################
     
     print("test 3")
     
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
