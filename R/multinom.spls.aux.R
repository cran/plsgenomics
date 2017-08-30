### multinom.spls.aux.R  (2015-10)
###
###    Ridge Iteratively Reweighted Least Squares followed by Adaptive Sparse PLS regression for 
###    multicategorial response
###    Short version for multiple call in cross-validation procedure
###
### Copyright 2015-10 Ghislain DURIF
###
### Adapted from rpls function in plsgenomics package, copyright 2006-01 Sophie Lambert-Lacroix
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


multinom.spls.aux <- function(sXtrain, sXtrain.nosvd=NULL, Ytrain, lambda.ridge, lambda.l1, ncomp, sXtest, sXtest.nosvd=NULL, 
                              adapt=TRUE, maxIter=100, svd.decompose=TRUE, 
                              meanXtrain, sigma2train, 
                              center.X=TRUE, scale.X=FALSE, weighted.center=TRUE) {
     
     
     #####################################################################
     #### Initialisation
     #####################################################################
     sXtrain <- as.matrix(sXtrain)
     ntrain <- nrow(sXtrain) # nb observations
     p <- ncol(sXtrain) # nb covariates
     index.p <- c(1:p)
     Ytrain <- as.matrix(Ytrain)
     q <- ncol(Ytrain)
     one <- matrix(1,nrow=1,ncol=ntrain)
     ntest <- nrow(sXtest)
     
     r <- p #min(p, ntrain)
     G <- max(Ytrain)
     
     #Compute Zblock 
     Z <- cbind(rep(1,ntrain),sXtrain)
     Zbloc <- matrix(0,nrow=ntrain*G,ncol=G*(r+1))
     
     Zt <- cbind(rep(1,ntest),sXtest)
     Ztestbloc <- matrix(0,nrow=ntest*G,ncol=G*(r+1))
     
     for (g in 1:G) {
          row <- (0:(ntrain-1))*G+g
          col <- (r+1)*(g-1)+1:(r+1)
          Zbloc[row,col] <- Z
          row <- (0:(ntest-1))*G+g
          Ztestbloc[row,col] <- Zt
     }
     rm(Z)
     Zt <- NULL
     
     
     #####################################################################
     #### Ridge IRLS step
     #####################################################################
     
     fit <- mwirrls(Y=Ytrain, Z=Zbloc, Lambda=lambda.ridge, NbrIterMax=maxIter, WKernel=diag(rep(1,ntrain*G)))
     
     converged=fit$Cvg
     
     #  Check WIRRLS convergence
     if (converged==0) {
          warning("Message from multinom.spls.aux : Ridge IRLS did not converge; try another lambda.ridge value")
     }
     
     # if ncomp == 0 then wirrls without spls step
     if (ncomp==0) {
          BETA <- fit$Coefficients
     }
     
     
     #####################################################################
     #### weighted SPLS step
     #####################################################################
     
     # if ncomp > 0
     if (ncomp!=0) {
          
          #Compute ponderation matrix V and pseudo variable z
          #Pseudovar = Eta + W^-1 Psi
          
          # Eta = X * betahat (covariate summary)
          Eta <- Zbloc %*% fit$Coefficients
          
          ## Run SPLS on Xtrain without svd decomposition
          if(svd.decompose) {
               
               p <- ncol(sXtrain.nosvd)
               r <- p
               sXtrain = sXtrain.nosvd
               sXtest = sXtest.nosvd
               
               #Compute Zblock (for X without svd)
               Z <- cbind(rep(1,ntrain),sXtrain)
               Zbloc <- matrix(0,nrow=ntrain*G,ncol=G*(r+1))
               
               Zt <- cbind(rep(1,ntest),sXtest)
               Ztestbloc <- matrix(0,nrow=ntest*G,ncol=G*(r+1))
               
               for (g in 1:G) {
                    row <- (0:(ntrain-1))*G+g
                    col <- (r+1)*(g-1)+1:(r+1)
                    Zbloc[row,col] <- Z
                    row <- (0:(ntest-1))*G+g
                    Ztestbloc[row,col] <- Zt
               }
               rm(Z)
               Zt <- NULL
          }
          
          # compute ponderation matrix V and mean parameter mu
          mu <- rep(0, length(Eta))
          V <- matrix(0, length(mu), length(mu))
          Vinv <- matrix(0, length(mu), length(mu))
          for (kk in 1:ntrain) {
               mu[G*(kk-1)+(1:G)] <- exp(Eta[G*(kk-1)+(1:G)])/(1+sum(exp(Eta[G*(kk-1)+(1:G)])))
               Blocmu <- mu[G*(kk-1)+(1:G)]
               BlocV <- -Blocmu %*% t(Blocmu)
               BlocV <- BlocV + diag(Blocmu)
               V[G*(kk-1)+(1:G), G*(kk-1)+(1:G)] <- BlocV
               Vinv[G*(kk-1)+(1:G), G*(kk-1)+(1:G)] <- tryCatch(solve(BlocV), error = function(e) return(ginv(BlocV)))
          }
          Psi <- fit$Ybloc - mu
          
          # V-Center the Xtrain and pseudo variable
          col.intercept <- seq(from=1, to=G*(p+1), by=(p+1)) # column corresponding to intercept
          index <- 1:(G*(p+1))
          Xbloc <- Zbloc[,-col.intercept]
          Cte <- Zbloc[,col.intercept]
          # Weighted centering of Pseudo variable
          H <- t(Cte) %*% V %*% Cte
          VMeanPseudoVar <- solve(H, t(Cte) %*% (V %*% Eta + Psi))
          VCtrPsi <- Psi
          VCtrEta <- Eta - Cte %*% VMeanPseudoVar
          # Weighted centering of sXtrain
          VMeansXtrain <- solve(H, t(Cte) %*% V %*% Xbloc)
          VCtrsXtrain <- Xbloc - Cte %*% VMeansXtrain
          rm(H)
          
          pseudoVar = Eta + Vinv %*% Psi
          pseudoVar = pseudoVar - Cte %*% VMeanPseudoVar
          
          if(center.X && weighted.center) {
               sXtrain <- VCtrsXtrain
          } else {
               sXtrain <- Xbloc
          }
          
          # SPLS(X, pseudo-var, weighting = V)
          resSPLS = spls.in(Xtrain=sXtrain, Ytrain=pseudoVar, ncomp=ncomp, weight.mat=V, lambda.l1=lambda.l1, adapt=adapt,
                            center.X=FALSE, center.Y=FALSE, scale.X=FALSE, scale.Y=FALSE, weighted.center=FALSE)
          
          #Express regression coefficients w.r.t. the columns of [1 sX] for ncomp
          BETA <- matrix(0, nrow=G*(r+1), ncol=1)
          BETA[-col.intercept,] <- resSPLS$betahat
          BETA[col.intercept,] <- VMeanPseudoVar - VMeansXtrain %*% BETA[-col.intercept,]
          
     }
     
     
     #####################################################################
     #### classification step
     #####################################################################
     
     hatYtest <- numeric(ntest)
     Eta.test <- matrix(0, nrow=G+1, ncol=1)
     proba.test <- matrix(0, nrow=ntest, ncol=G+1)
     
     Eta.test <- cbind(rep(0,ntest),matrix(Ztestbloc%*%BETA,nrow=ntest,byrow=TRUE))
     proba.test <- t(apply(exp(Eta.test), 1, function(x) x/sum(x)))
     hatYtest <- as.matrix(apply(proba.test,1,which.max)-1)
     
     
     #####################################################################
     #### Conclude
     #####################################################################
     
     ##Compute the coefficients w.r.t. [1 X]
     Beta <- t(matrix(BETA,nrow=G,byrow=TRUE))
     Coefficients <- t(matrix(0, nrow=G, ncol=(p+1)))
     if(p > 1) {
          Coefficients[-1,] <- diag(c(1/sqrt(sigma2train))) %*% Beta[-1,]
     } else {
          Coefficients[-1,] <- (1/sqrt(sigma2train))%*%Beta[-1,]
     }
     Coefficients[1,] <- Beta[1,] - meanXtrain %*% Coefficients[-1,]
     
     
     
     #### RETURN
     
     result <- list(Coefficients=Coefficients, hatYtest=hatYtest, converged=converged, lenA=resSPLS$lenA)
     class(result) <- "multinom.spls.aux"
     return(result)
     
}
