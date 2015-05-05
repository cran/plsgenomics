\name{sample.bin}
\alias{sample.bin}

\title{Generates design matrix X with correlated block of covariates and a binary random reponse 
depening on X through logit model}
\usage{
sample.bin(n, p, kstar, lstar, beta.min, beta.max, mean.H=0, sigma.H, sigma.F, seed=NULL)
}
\arguments{
	\item{n}{the number of observations in the sample.}
	\item{p}{the number of covariates in the sample.}
	\item{kstar}{the number of underlying latent variables used to generates the design matrix 
	\code{X}, \code{kstar <= p}. \code{kstar} is also the number of blocks in the design matrix 
	(see details).}
	\item{lstar}{the number of blocks in the design matrix \code{X} that are used to generates the 
	response \code{Y}, i.e. with non null coefficients in vector \code{B}, \code{lstar <= kstar}.}
	\item{beta.min}{the inf bound for non null coefficients (see details).}
	\item{beta.max}{the sup bound for non null coefficients (see details).}
	\item{mean.H}{the mean of latent variables used to generates \code{X}.}
	\item{sigma.H}{the standard deviation of latent variables used to generates \code{X}.}
	\item{sigma.F}{the standard deviation of the noise added to latent variables used to generates 
	\code{X}.}
	\item{seed}{an positive integer, if non NULL it fix the seed (with the command 
	\code{set.seed}) used for random number generation.}
	
}
\description{
	
	The function \code{sample.bin} generates a random sample with p predictors X, a binary 
	response Y, and n observations, through a logistic model, where the response Y is generated as 
	a Bernoulli random variable of parameter logit^\{-1\}(XB), the coefficients B are sparse, and 
	the design matrix X is composed of correlated blocks of predictors.
	
}

\details{
	The set (1:p) of predictors is partitioned into kstar block. Each block k (k=1,...,kstar) 
	depends on a latent variable H.k which are independent and identically distributed following a 
	distribution N(mean.H, sigma.H^2). Each columns X.j of the matrix X is generated as H.k + F.j 
	for j in the block k, where F.j is independent and identically distributed gaussian noise 
	N(0,sigma.F^2).
	
	The coefficients B are generated as random between beta.min and beta.max on lstar blocks, 
	randomly chosen, and null otherwise. The variables with non null coefficients are then 
	relevant to explain the response, whereas the ones with null coefficients are not.
	
	The response is generated as by a drawing a Bernoulli random variable of parameter 
	logit^\{-1\}(XB).
	
	The details of the procedure are developped by Durif et al. (2015).
 
}

\value{
	A list with the following components:
	\item{X}{the (n x p) design matrix, containing the \code{n} observations for each of the 
	\code{p} predictors.}
	\item{Y}{the (n) vector of Y observations.}
	\item{proba}{the n vector of Bernoulli parameters used to generate the response, in particular 
	\code{logit^{-1}(X \%*\% B)}. }
	\item{sel}{the index in (1:p) of covariates with non null coefficients in \code{B}.}
	\item{nosel}{the index in (1:p) of covariates with null coefficients in \code{B}.}
	\item{B}{the (n) vector of coefficients.}
	\item{block.partition}{a (p) vector indicating the block of each predictors in (1:kstar).}
	\item{p}{the number of covariates in the sample.}
	\item{kstar}{the number of underlying latent variables used to generates the design matrix 
	\code{X}, \code{kstar <= p}. \code{kstar} is also the number of blocks in the design matrix 
	(see details).}
	\item{lstar}{the number of blocks in the design matrix \code{X} that are used to generates the 
	response \code{Y}, i.e. with non null coefficients in vector \code{B}, \code{lstar <= kstar}.}
	\item{p0}{the number of predictors with non null coefficients in \code{B}.}
	\item{block.sel}{a (lstar) vector indicating the index in (1:kstar) of blocks with predictors 
	having non null coefficient in \code{B}.}
	\item{beta.min}{the inf bound for non null coefficients (see details).}
	\item{beta.max}{the sup bound for non null coefficients (see details).}
	\item{mean.H}{the mean of latent variables used to generates \code{X}.}
	\item{sigma.H}{the standard deviation of latent variables used to generates \code{X}.}
	\item{sigma.F}{the standard deviation of the noise added to latent variables used to 
	generates \code{X}.}
	\item{seed}{an positive integer, if non NULL it fix the seed (with the command 
	\code{set.seed}) used for random number generation.}
}

\references{
 
	G. Durif, F. Picard, S. Lambert-Lacroix (2015). Adaptive sparse PLS for logistic regression, 
	(in prep), available on (\url{http://arxiv.org/}).

}

\author{
	Ghislain Durif (\url{http://lbbe.univ-lyon1.fr/-Durif-Ghislain-.html}). 

}
\seealso{\code{\link{sample.cont}}.}

\examples{
### load plsgenomics library
library(plsgenomics)

### generating data
n <- 100
p <- 1000
sample1 <- sample.bin(n=n, p=p, kstar=20, lstar=2, beta.min=0.25, beta.max=0.75, mean.H=0.2, 
				sigma.H=10, sigma.F=5)
str(sample1)

}
