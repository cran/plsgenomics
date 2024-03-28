### sample.bin.R  (2014-10)
###
###    Generates covariate matrix X with correlated block of covariates and a binary random reponse depening on X through logit model
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

#' @title Generates covariate matrix X with correlated block of covariates and 
#' a binary random reponse depening on X through a logistic model
#' @aliases sample.bin
#' 
#' @description
#' The function \code{sample.bin} generates a random sample of n observations, 
#' composed of p predictors, collected in the n x p matrix X, and a binary 
#' response, in a vector Y of length n, thanks to a logistic model, where the 
#' response Y is generated as a Bernoulli random variable of parameter 
#' logit^\{-1\}(XB), the coefficients B are sparse. In addition, the covariate 
#' matrix X is composed of correlated blocks of predictors.
#' 
#' @details
#' The set (1:p) of predictors is partitioned into kstar block. 
#' Each block k (k=1,...,kstar) depends on a latent variable H.k which are 
#' independent and identically distributed following a Gaussian distribution 
#' N(mean.H, sigma.H^2). Each columns X.j of the matrix X is generated 
#' as H.k + F.j for j in the block k, where F.j is independent and identically 
#' distributed gaussian noise N(0,sigma.F^2).
#' 
#' The coefficients B are generated as random between beta.min and beta.max 
#' on lstar blocks, randomly chosen, and null otherwise. The variables with 
#' non null coefficients are then relevant to explain the response, whereas 
#' the ones with null coefficients are not.
#' 
#' The response is generated as by drawing one observation of n different 
#' Bernoulli random variables of parameters logit^\{-1\}(XB).
#' 
#' The details of the procedure are developped by Durif et al. (2018).
#' 
#' @param n the number of observations in the sample.
#' @param p the number of covariates in the sample.
#' @param kstar the number of underlying latent variables used to generates 
#' the covariate matrix \code{X}, \code{kstar <= p}. \code{kstar} is also the 
#' number of blocks in the covariate matrix (see details).
#' @param lstar the number of blocks in the covariate matrix \code{X} that are 
#' used to generates the response \code{Y}, i.e. with non null coefficients 
#' in vector \code{B}, \code{lstar <= kstar}.
#' @param beta.min the inf bound for non null coefficients (see details).
#' @param beta.max the sup bound for non null coefficients (see details).
#' @param mean.H the mean of latent variables used to generates \code{X}.
#' @param sigma.H the standard deviation of latent variables used to 
#' generates \code{X}.
#' @param sigma.F the standard deviation of the noise added to latent 
#' variables used to generates \code{X}.
#' @param seed an positive integer, if non NULL it fix the seed (with the 
#' command \code{set.seed}) used for random number generation.
#' 
#' @return An object with the following attributes:
#' \item{X}{the (n x p) covariate matrix, containing the \code{n} observations 
#' for each of the \code{p} predictors.}
#' \item{Y}{the (n) vector of Y observations.}
#' \item{proba}{the n vector of Bernoulli parameters used to generate the 
#' response, in particular \code{logit^{-1}(X \%*\% B)}.}
#' \item{sel}{the index in (1:p) of covariates with non null coefficients 
#' in \code{B}.}
#' \item{nosel}{the index in (1:p) of covariates with null coefficients 
#' in \code{B}.}
#' \item{B}{the (n) vector of coefficients.}
#' \item{block.partition}{a (p) vector indicating the block of each predictors 
#' in (1:kstar).}
#' \item{p}{the number of covariates in the sample.}
#' \item{kstar}{the number of underlying latent variables used to generates 
#' the covariate matrix \code{X}, \code{kstar <= p}. \code{kstar} is also 
#' the number of blocks in the covariate matrix (see details).}
#' \item{lstar}{the number of blocks in the covariate matrix \code{X} that 
#' are used to generates the response \code{Y}, i.e. with non null 
#' coefficients in vector \code{B}, \code{lstar <= kstar}.}
#' \item{p0}{the number of predictors with non null coefficients in \code{B}.}
#' \item{block.sel}{a (lstar) vector indicating the index in (1:kstar) of 
#' blocks with predictors having non null coefficient in \code{B}.}
#' \item{beta.min}{the inf bound for non null coefficients (see details).}
#' \item{beta.max}{the sup bound for non null coefficients (see details).}
#' \item{mean.H}{the mean of latent variables used to generates \code{X}.}
#' \item{sigma.H}{the standard deviation of latent variables used to 
#' generates \code{X}.}
#' \item{sigma.F}{the standard deviation of the noise added to latent 
#' variables used to generates \code{X}.}
#' \item{seed}{an positive integer, if non NULL it fix the seed 
#' (with the command \code{set.seed}) used for random number generation.}
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
#' @seealso \code{\link{sample.cont}}
#' 
#' @examples
#' ### load plsgenomics library
#' library(plsgenomics)
#' 
#' ### generating data
#' n <- 100
#' p <- 1000
#' sample1 <- sample.bin(n=n, p=p, kstar=20, lstar=2, beta.min=0.25, 
#'                       beta.max=0.75, mean.H=0.2, 
#'                       sigma.H=10, sigma.F=5)
#' 
#' str(sample1)
#' 
#' @export
sample.bin = function(n, p, kstar, lstar, beta.min, beta.max, 
                      mean.H=0, sigma.H=1, sigma.F=1, seed=NULL) {
	
	### input
	# n : sample size
	# p : number of covariates
	# kstar : number of latent variables
	# lstar : number of blocks of covariates used to generate the response Y
	# seed : seed for the random generator
	
	### tests on input
	
	if((!is.null(seed)) && (!is.numeric(seed)) && (round(seed)-seed!=0)) {
		stop("Message from sample.cont: seed must be integer")
	}
	
	if((!is.numeric(n)) || (round(n)-n!=0) || (!is.numeric(p)) || (round(p)-p!=0) || (!is.numeric(kstar)) || (round(kstar)-kstar!=0) || (!is.numeric(lstar)) || (round(lstar)-lstar!=0) ) {
		stop("Message from sample.cont: n, p, kstar, lstar must be integer")
	}
	
	if((!is.numeric(mean.H)) || (!is.numeric(sigma.H)) || (!is.numeric(sigma.F)) ) {
		stop("Message from sample.cont: mean.H, sigma.H, sigma.F are not of valid type")
	}
	
	if((sigma.H<0) || (sigma.F<0)) {
		stop("Message from sample.cont: sigma.H, sigma.F are not of valid type")
	}
	
	if(n<1) {
		stop("Message from sample.cont: n<1, must be strict positive integer")
	}
	
	if(p<1) {
		stop("Message from sample.cont: p<1, must be strict positive integer")
	}
	
	if(kstar<1) {
		stop("Message from sample.cont: kstar<1, must be strict positive integer")
	}
	
	if(lstar>kstar) {
		stop("Message from sample.cont: kstar<lstar, try to use more blocks in X than it actually exists")
	}
	
	if(p<kstar) {
		stop("Message from sample.cont: p<kstar, more blocks than actual covariates")
	}
	
	### random generation
	if(!is.null(seed)) {
		set.seed(seed)
	}
	
	# block size
	block.size <- p %/% kstar
	last.block.size <- p %/% kstar + p %% kstar
	
	# latent variables
	for(i in 1:kstar) {
		assign(paste0("H", i), rnorm(n, mean=mean.H, sd=sigma.H))
	}
	
	### generation of X
	# split into k blocks
	if(last.block.size==0) {
		block.partition <- rep(1:kstar, each=block.size)
	} else {
		block.partition <- c( rep(1:(kstar-1), each=block.size), rep(kstar, each=last.block.size) ) # the last block has more covariates
	}
	
	# index of the jth column in X determines which latent variables to use to generate it
	X <- matrix(data=NA, nrow=n, ncol=p)
	
	X <- sapply(1:p, function(j) {
		
		F <- rnorm(n, mean=0, sd=sigma.F) # noise of column j
		
		return(get(paste0("H", block.partition[j]), inherits=TRUE) + F) # on utilise Hj suivant l'intervalle trouve
		
	})
	
	### generation of Y (binary)
	Y <- numeric(n)
	proba <- numeric(n)
	
	# selected blocks
	block.sel <- sort(sample.int(kstar, size=lstar))
	
	# linear coefficients: non null in the lstar selected blocks
	sel <- (1:p)[block.partition %in% block.sel]
	
	nosel <- (1:p)[-sel]
	
	p0 <- length(sel)
	
	B <- numeric(p)
	B[sel] <- signif(runif(n=p0, min=beta.min, max=beta.max), digits=2)
	B[nosel] <- rep(0, length.out=p-p0)
	
	eta <- X %*% B
	proba <- inv.logit( drop( eta))
	Y <- rbinom(1, size=n, prob = proba)
	
	for(i in 1:n) {
		
		proba[i] <- inv.logit( drop( X[i,] %*% B))
		Y[i] <- rbinom(1,1, prob = proba[i])
		
	}
	
	Y <- as.matrix(Y)
	
	### output:
	# sel: index of variables used to generate Y
	# nosel: index of unused variables
	
	return(list(X=X, Y=Y, proba=proba, sel=sel, nosel=nosel, B=B, block.partition=block.partition, n=n, p=p, kstar=kstar, lstar=lstar, p0=p0, block.sel=block.sel, beta.min=beta.min, beta.max=beta.max, mean.H=mean.H, sigma.H=sigma.H, sigma.F=sigma.F, seed=seed))

}
