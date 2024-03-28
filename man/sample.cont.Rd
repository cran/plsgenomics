\name{sample.cont}
\alias{sample.cont}

\title{Generates design matrix X with correlated block of covariates and a continuous random 
reponse Y depening on X through gaussian linear model Y=XB+E}
\usage{
sample.cont(n, p, kstar, lstar, beta.min, beta.max, mean.H=0, sigma.H, 
               sigma.F, sigma.E, seed=NULL)
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
	\item{sigma.F}{the standard deviation of the noise added to latent variables used to 
	generates \code{X}.}
	\item{sigma.E}{the standard deviation of the noise in the linear model 
	\code{Y = X \%*\%* B + E} used to generates \code{Y}.}
	\item{seed}{an positive integer, if non NULL it fix the seed (with the command 
	\code{set.seed}) used for random number generation.}
	
}
\description{
	
	The function \code{sample.cont} generates a random sample with p predictors X, a response Y, 
	and n observations, through a linear model Y=XB+E, where the noise E is gaussian, 
	the coefficients B are sparse, and the design matrix X is composed of correlated blocks of
	predictors.
	
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
	
	The response is generated as Y = X \%*\% B + E, where E is some gaussian noise N(0,sigma.E^2).
	
	The details of the procedure are developped by Durif et al. (2018).
 
}

\value{
	A list with the following components:
	\item{X}{the (n x p) design matrix, containing the \code{n} observations for each of the 
	\code{p} predictors.}
	\item{Y}{the (n) vector of Y observations.}
	\item{residuals}{the (n) vector corresponding to the noise \code{E} in the model 
	\code{Y = X \%*\% B + E}.}
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
	\item{sigma.E}{the standard deviation of the noise in the linear model.}
	\item{seed}{an positive integer, if non NULL it fix the seed (with the command 
	\code{set.seed}) used for random number generation.}
}

\references{
 
	Durif, G., Modolo, L., Michaelsson, J., Mold, J.E., Lambert-Lacroix, S., 
     Picard, F., 2018. High dimensional classification with combined 
     adaptive sparse PLS and logistic regression. Bioinformatics 34, 
     485–493. \doi{10.1093/bioinformatics/btx571}.
     Available at \url{http://arxiv.org/abs/1502.05933}.

}

\author{
	Ghislain Durif (\url{https://gdurif.perso.math.cnrs.fr/}). 

}
\seealso{\code{\link{sample.bin}}.}

\examples{
### load plsgenomics library
library(plsgenomics)

### generating data
n <- 100
p <- 1000
sample1 <- sample.cont(n=n, p=p, kstar=20, lstar=2, beta.min=0.25, beta.max=0.75, mean.H=0.2, 
					sigma.H=10, sigma.F=5, sigma.E=5)
str(sample1)

}
