\name{pls.regression}
\alias{pls.regression}

\title{Multivariate Partial Least Squares Regression}
\usage{
pls.regression(Xtrain, Ytrain, Xtest=NULL, ncomp=NULL,  unit.weights=TRUE)
}
\arguments{
  \item{Xtrain}{a (ntrain x p) data matrix of predictors. \code{Xtrain} may be a matrix or a
  data frame. Each row corresponds to an observation and each column to a predictor variable.}
  \item{Ytrain}{a (ntrain x m) data matrix of responses. \code{Ytrain} may be a vector (if m=1), a matrix or a data frame. 
  If \code{Ytrain} is a matrix or a data
  frame, each row corresponds to an observation and each column to a response variable. If
  \code{Ytrain} is a vector, it contains the unique response variable for each
  observation.}
  \item{Xtest}{a (ntest x p) matrix containing the predictors for the test data
  set. \code{Xtest} may also be a
  vector of length p (corresponding to only one test observation).}
  \item{ncomp}{the number of latent components to be used for regression. If
  \code{ncomp} is a vector of integers, the regression model is built
  successively with each number of components. If \code{ncomp=NULL}, the maximal
  number of components min(ntrain,p) is chosen.}
   \item{unit.weights}{if \code{TRUE} then the latent components
    will be constructed from weight vectors that are standardized to length 1, 
    otherwise the weight vectors do not have length 1 but the latent components have
    norm 1.}
}
\description{
 The function \code{pls.regression} performs pls multivariate regression (with several response variables and 
  several predictor variables) using de Jong's SIMPLS algorithm. This function
  is an adaptation of R. Wehrens' code from the package pls.pcr. 
  
}

\details{
  The columns of the data matrices \code{Xtrain} and \code{Ytrain} must not be centered to have mean
  zero, since centering is performed by the function \code{pls.regression} as a preliminary step
  before the SIMPLS algorithm is run. 

  In the original definition of SIMPLS by de Jong (1993), the weight vectors have
  length 1. If the weight vectors are standardized to have length 1, 
  they satisfy a simple optimality criterion (de Jong, 1993). However, it is
  also usual (and computationally efficient) to standardize the latent
  components to have length 1. 
 
  In contrast to the original version found in the package \code{pls.pcr}, 
  the prediction for the observations from \code{Xtest} is performed after
  centering the columns of \code{Xtest} by substracting the columns means
  calculated from \code{Xtrain}.
}

\value{
  A list with the following components:
  \item{B}{the (p x m x length(\code{ncomp})) matrix containing the regression coefficients. Each row corresponds
  to a predictor variable and each column to a response variable. The third
  dimension of the matrix \code{B} corresponds to the number of PLS components used
  to compute the regression coefficients. If \code{ncomp} has length 1, \code{B}
  is just a (p x m) matrix.}
  \item{Ypred}{the (ntest x m x length(\code{ncomp})) containing the predicted
  values of the response variables for the observations from \code{Xtest}. The
  third dimension of the matrix \code{Ypred} corresponds to the number of PLS
  components used to compute the regression coefficients. } 
  \item{P}{the (p x max(\code{ncomp})) matrix containing the X-loadings.}
  \item{Q}{the (m x max(\code{ncomp})) matrix containing the Y-loadings.}
  \item{T}{the (ntrain x max(\code{ncomp})) matrix containing the X-scores (latent components)}
  \item{R}{the (p x max(\code{ncomp})) matrix containing the weights used to construct the
  latent components.}
  \item{meanX}{the p-vector containing the means of the columns of \code{Xtrain}.}
  }

\references{
 
 S. de Jong (1993). SIMPLS: an alternative approach to partial least squares 
 regression, Chemometrics Intell. Lab. Syst. \bold{18}, 251--263.

 C. J. F. ter Braak and S. de Jong (1993). The objective function of partial
 least squares regression, Journal of Chemometrics \bold{12}, 41--54.

}

\author{
  Anne-Laure Boulesteix (\url{http://www.ibe.med.uni-muenchen.de/organisation/mitarbeiter/020_professuren/boulesteix/index.html}) 
   and 
  Korbinian Strimmer (\url{http://strimmerlab.org/}).

 Adapted in part from 
  pls.pcr code by R. Wehrens (\url{http://cran.r-project.org/src/contrib/Descriptions/pls.html}). 


}
\seealso{\code{\link{pls.lda}}, \code{\link{TFA.estimate}},
\code{\link{pls.regression.cv}}.}

\examples{
# load plsgenomics library
library(plsgenomics)

# load the Ecoli data
data(Ecoli)

# perform pls regression
# with unit latent components
pls.regression(Xtrain=Ecoli$CONNECdata,Ytrain=Ecoli$GEdata,Xtest=Ecoli$CONNECdata,ncomp=1:3,unit.weights=FALSE)

# with unit weight vectors
pls.regression(Xtrain=Ecoli$CONNECdata,Ytrain=Ecoli$GEdata,Xtest=Ecoli$CONNECdata,ncomp=1:3,unit.weights=TRUE)


}
\keyword{multivariate}
