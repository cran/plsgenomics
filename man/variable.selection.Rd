\name{variable.selection}
\alias{variable.selection}

\title{Variable selection using the PLS weights}
\usage{
variable.selection(X, Y, nvar=NULL)
}
\arguments{
  \item{X}{a (n x p) data matrix of predictors. X may be a matrix or a
  data frame. Each row corresponds to an observation and each column corresponds
  to a predictor variable.}
  \item{Y}{a vector of length n giving the classes of the n observations.  The two classes must be coded as 1,2.}
  \item{nvar}{the number of variables to be returned. If \code{nvar=NULL}, all the variables are
   returned.}
}
\description{
The function \code{variable.selection}  performs variable selection for binary classification.
   
}

\details{
  
  The function \code{variable.selection} orders the variables according to 
  the absolute value of the weight defining the first PLS
  component. This ordering is equivalent to the ordering obtained with the 
  F-statistic and t-test with equal variances (Boulesteix, 2004).
 

  For computational reasons, the function \code{variable.selection} does not use 
  the pls algorithm, but the obtained ordering of the variables is exactly 
  equivalent to the ordering obtained using the PLS weights output by
  \code{\link{pls.regression}}.
  
}

\value{
  A vector of length \code{nvar} (or of length p if \code{nvar=NULL}) containing the indices of the variables to be selected. 
  The variables are ordered from the best to the worst variable. 
}

\references{

 A. L. Boulesteix (2004). PLS dimension reduction for classification with microarray data,
 Statistical Applications in Genetics and Molecular Biology \bold{3}, Issue 1, Article 33.
 
 S. de Jong (1993). SIMPLS: an alternative approach to partial least squares 
 regression, Chemometrics Intell. Lab. Syst. \bold{18}, 251--263.

}

\author{
  Anne-Laure Boulesteix (\url{http://www.statistik.lmu.de/~socher/}) 
   
 

}
\seealso{\code{\link{pls.regression}}.}

\examples{
# load plsgenomics library
library(plsgenomics)

# generate X and Y (4 observations and 3 variables)
X<-matrix(c(4,3,3,4,1,0,6,7,3,5,5,9),4,3,byrow=FALSE)
Y<-c(1,1,2,2)

# select the 2 best variables
variable.selection(X,Y,nvar=2)
# order the 3 variables
variable.selection(X,Y)

# load the leukemia data 
data(leukemia)

# select the 50 best variables from the leukemia data
variable.selection(leukemia$X,leukemia$Y,nvar=50)

}

\keyword{multivariate}
