### spls.in.R  (2014-10)
###
###    Sparse PLS regression called in m.rirls.spls
###
### Copyright 2015-10 Ghislain DURIF
###
### Adapted from R package "spls"
### Reference: Chun H and Keles S (2010)
### "Sparse partial least squares for simultaneous dimension reduction and variable selection",
### Journal of the Royal Statistical Society - Series B, Vol. 72, pp. 3--25.
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

spls.in <- function(Xtrain, Ytrain, lambda.l1, ncomp, weight.mat=NULL, 
                    adapt=TRUE, center.X=TRUE, center.Y=TRUE, 
                    scale.X=TRUE, scale.Y=TRUE, weighted.center=FALSE) {
     
     #####################################################################
     #### Initialisation
     #####################################################################
     ntrain <- nrow(Xtrain) # nb observations
     p <- ncol(Xtrain) # nb covariates
     index.p <- c(1:p)
     Ytrain <- as.matrix(Ytrain)
     q <- ncol(Ytrain)
     
     # weighting matrix V
     if(!is.null(weight.mat)) { # weighting in scalar product (in observation space of dimension n)
          V <- as.matrix(weight.mat) 
     } else { # no weighting in scalar product
          V <- diag(rep(1, ntrain), nrow=ntrain, ncol=ntrain)
     }
     
     
     #####################################################################
     #### centering and scaling
     #####################################################################
     if (!weighted.center) {
          
          # Xtrain mean
          meanXtrain <- apply(Xtrain, 2, mean)
          
          # Xtrain sd
          sigmaXtrain <- apply(Xtrain, 2, sd)
          # test if predictors with null variance
          if ( any( sigmaXtrain < .Machine$double.eps )) {
               stop("Some of the columns of the predictor matrix have zero variance.")
          }
          
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
               sYtrain <- scale( Ytrain, center=meanYtrain, scale=sigmaYtrain )
          } else if(center.Y && !scale.Y) {
               sYtrain <- scale( Ytrain, center=meanYtrain, scale=FALSE )
          } else {
               sYtrain <- Ytrain
          }
          
          
     } else { # weighted scaling
          
          sumV <- sum(diag(V))
          
          # X mean
          meanXtrain <- matrix(diag(V), nrow=1) %*% Xtrain / sumV
          
          # X sd
          sigmaXtrain <- apply(Xtrain, 2, sd)
          # test if predictors with null variance
          if ( any( sigmaXtrain < .Machine$double.eps ) ) {
               stop("Some of the columns of the predictor matrix have zero variance.")
          }
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
          
     }
     
     
     #####################################################################
     #### Result objects
     #####################################################################
     betahat <- matrix(0, nrow=p, ncol=1)
     betamat <- list()
     X1 <- sXtrain
     Y1 <- sYtrain
     
     W <- matrix(data=NA, nrow=p, ncol=ncomp) # spls weight over each component
     T <- matrix(data=NA, nrow=ntrain, ncol=ncomp) # spls components
     P <- matrix(data=NA, nrow=ncomp, ncol=p) # regression of X over T
     Q <- matrix(data=NA, nrow=ncomp, ncol=q) # regression of Y over T
     
     
     #####################################################################
     #### Main iteration
     #####################################################################
     
     if ( is.null(colnames(Xtrain)) ) {
          Xnames <- index.p
     } else { 
          Xnames <- colnames(Xtrain)
     }
     
     new2As <- list()
     
     ## SPLS
     for (k in 1:ncomp) {
          
          ## define M
          M <- t(X1) %*% (V %*% Y1)
          
          #### soft threshold
          Mnorm1 <- median( abs(M) )
          
          M <- M / Mnorm1
          
          ## adpative version
          if (adapt) {
               wi <- 1/abs(M)
               
               what <- ust.adapt(M, lambda.l1, wi)
               
          } else {
               ## non adaptive version
               what <- ust(M, lambda.l1)
          }
          
          #### construct active set A
          A <- unique( index.p[ what!=0 | betahat[,1]!=0 ] )
          new2A <- index.p[ what!=0 & betahat[,1]==0 ]
          
          #### fit pls with selected predictors (meaning in A)
          X.A <- sXtrain[ , A, drop=FALSE ]           
          plsfit <- wpls( Xtrain=X.A, Ytrain=sYtrain, weight.mat=V, ncomp=min(k,length(A)), type="pls1", center.X=FALSE, scale.X=FALSE, center.Y=FALSE, scale.Y=FALSE, weighted.center=FALSE )
          
          
          #### output storage
          
          # weights
          w.k <- matrix(data=what, ncol=1)
          w.k <- w.k / sqrt(as.numeric(t(w.k) %*% w.k))
          W[,k] <- w.k
          
          # components on total observation space
          t.k <- (X1 %*% w.k) / as.numeric(t(w.k) %*% w.k)
          T[,k] <- t.k
          
          # regression of X over T
          p.k <- (t(X1) %*% (V %*% t.k)) / as.numeric(t(t.k) %*% (V %*% t.k))
          P[k,] <- t(p.k)
          
          # regression of Y over T
          q.k <- (t(Y1) %*% (V %*% t.k)) / as.numeric(t(t.k) %*% (V %*% t.k))
          Q[k,] <- t(q.k)
          
          
          ## update
          
          Y1 <- sYtrain - plsfit$T %*% plsfit$Q
          X1 <- sXtrain
          X1[,A] <- sXtrain[,A] - plsfit$T %*% plsfit$P
          
          
          betahat <- matrix( 0, p, q )
          betahat[A,] <- matrix( plsfit$coeff, length(A), q )
          betamat[[k]] <- betahat # for cv.spls
          
          # variables that join the active set
          new2As[[k]] <- new2A
          
          
     }
     
     ##### return objects
     
     ## components in lower subspace of selected variables
     T.low <- plsfit$T
     
     #### betahat for non centered and non scaled data
     if((!scale.X) || (weighted.center)) { # if X non scaled, betahat don't have to be corrected regards sd.x
          sd.X <- rep(1, p)
     } else { # if X is scaled, it has to
          sd.X <- sigmaXtrain
     }
     if((!scale.Y) || (weighted.center)) {
          sd.Y <- 1
     } else {
          sd.Y <- sigmaYtrain
     }
     
     betahat.nc <- sd.Y * betahat / sd.X
     intercept <- meanYtrain - ( sd.Y * (drop( (meanXtrain / sd.X) %*% betahat)) )
     betahat.nc <- as.matrix(c(intercept, betahat.nc))
     
     #### return object
     result <- list( betahat=betahat, betahat.nc=betahat.nc, 
                     X.score=T, X.score.low=T.low, X.loading=P, Y.loading=Q, X.weight=W, 
                     A=A, lenA=length(A))
     
     return(result)
     
     
}
