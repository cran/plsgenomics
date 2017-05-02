\name{mrpls.cv}
\alias{mrpls.cv}

\title{Determination of the ridge regularization parameter and the number of PLS 
components to be used for classification with RPLS for categorical data}
\usage{
mrpls.cv(Ytrain,Xtrain,LambdaRange,ncompMax,NbIterMax=50)
}
\arguments{
  \item{Xtrain}{a (ntrain x p) data matrix of predictors. \code{Xtrain} must be a matrix. 
  Each row corresponds to an observation and each column to a predictor variable.}
  \item{Ytrain}{a ntrain vector of responses. \code{Ytrain} must be a vector. 
  \code{Ytrain} is a \{1,...,c+1\}-valued vector and contains the response variable for each
  observation. c+1 is the number of classes.}
  \item{LambdaRange}{the vector of positive real value from which the best ridge regularization 
  parameter has to be chosen by cross-validation.}
  \item{ncompMax}{a positive integer. The best number of components is chosen from  
  1,...,\code{ncompMax}. If \code{ncompMax}=0,then the Ridge regression is performed without 
  reduction dimension. }
  \item{NbIterMax}{a positive integer. \code{NbIterMax} is the maximal number of iterations in the 
  Newton-Rapson parts.}
}
\description{
 The function \code{mrpls.cv} determines the best ridge regularization parameter and the best 
 number of PLS components to be used for classification for Fort et al. (2005) MRPLS algorithm.

}

\details{
  A cross-validation procedure is used to determine the best ridge regularization parameter and 
  number of PLS components to be used for classification with MRPLS for categorical data 
  (for binary data see \code{\link{rpls}} and \code{\link{rpls.cv}}).
  At each cross-validation run, \code{Xtrain} is split into a pseudo training
  set (ntrain-1 samples) and a pseudo test set (1 sample) and the classification error rate is 
  determined for each value of ridge regularization parameter and number of components. Finally, 
  the function \code{mrpls.cv} returns the values of the ridge regularization parameter and 
  bandwidth for which the mean classification error rate is minimal. 
 
}


\value{
  A list with the following components:
  \item{Lambda}{the optimal regularization parameter.} 
  \item{ncomp}{the optimal number of PLS components.}
  }

\references{
 
 G. Fort, S. Lambert-Lacroix and Julie Peyre (2005). Reduction de dimension dans les modeles 
 lineaires generalises : application a la classification supervisee de donnees issues des biopuces.
 Journal de la SFDS, tome 146, n1-2, 117-152. 

}

\author{
  Sophie Lambert-Lacroix 
  (\url{http://membres-timc.imag.fr/Sophie.Lambert/}). 

}
\seealso{\code{\link{mrpls}}, \code{\link{rpls}}, \code{\link{rpls.cv}}.}

\examples{
\dontrun{
## between 5~15 seconds
# load plsgenomics library
library(plsgenomics)

# load SRBCT data
data(SRBCT)
IndexLearn <- c(sample(which(SRBCT$Y==1),10),sample(which(SRBCT$Y==2),4),
			sample(which(SRBCT$Y==3),7),sample(which(SRBCT$Y==4),9))

# Determine optimum ncomp and Lambda
nl <- mrpls.cv(Ytrain=SRBCT$Y[IndexLearn],Xtrain=SRBCT$X[IndexLearn,],
			LambdaRange=c(0.1,1),ncompMax=3)

# perform prediction by MRPLS
res <- mrpls(Ytrain=SRBCT$Y[IndexLearn],Xtrain=SRBCT$X[IndexLearn,],Lambda=nl$Lambda,
			ncomp=nl$ncomp,Xtest=SRBCT$X[-IndexLearn,])
sum(res$Ytest!=SRBCT$Y[-IndexLearn])
}
}

