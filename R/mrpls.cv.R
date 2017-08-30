### mrpls.cv.R  (2015-10)
###
###    Determination by Cross-validation of Ridge Partial Least square
###                 hyper-parameters for categorical data
###
### Copyright 2015-01 Sophie Lambert-Lacroix, Ghislain Durif
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

mrpls.cv <- function (Ytrain, Xtrain, LambdaRange, ncompMax, NbIterMax=50, ncores=1) {

     ##    INPUT VARIABLES
     #########################
     ##  Xtrain   : matrix ntrain x p
     ##      train data matrix
     ##  Ytrain   : vector ntrain
     ##      response variable {0,...c}-valued vector
     ##  LambdaRange : vector nLambda
     ##      possible values for the regularization parameter Lambda
     ##  NbIterMax : positive integer
     ##      max number of iteration in the MWIRRLS part
     ##  ncompMax : positive integer
     ##      maximal number of PLS components
     ##      if ncompMax=0 -> Ridge
     ## ncores : positive integer
     ##      number of cores to be used in parallel computing
     
     
     ##    OUTPUT VARIABLES
     ##########################
     ## Lambda : optimal regularization parameter Lambda
     ## ncomp : optimal number of PLS components
     
     
     ##  TEST ON INPUT VARIABLES
     ##############################
     #On Xtrain
     Xtrain <- as.matrix(Xtrain)
     if ((is.matrix(Xtrain)==FALSE)||(is.numeric(Xtrain)==FALSE)) {
          stop("Message from mrpls.cv.R: Xtrain is not of valid type")
     }
     
     if (ncompMax > dim(Xtrain)[2]) {
          warning("Message from mrpls.cv: ncomp>p is not valid, ncompMax is set to p")
          ncomp <- dim(Xtrain)[2]
     }
     
     if (dim(Xtrain)[2]==1) {
          # stop("Message from mrpls.cv: p=1 is not valid")
          warning("Message from mrpls.cv: p=1 is not valid, ncompMax is set to 0")
          ncompMax <- 0
     }
     
     ntrain <- dim(Xtrain)[1]
     
     #On Ytrain
     Ytrain <- as.matrix(Ytrain)
     if ((ncol(Ytrain)>1)||(is.numeric(Ytrain)==FALSE)) {
          stop("Message from mrpls.cv.R: Ytrain is not of valid type")
     }
     
     if (length(Ytrain)!=ntrain) {
          stop("Message from mrpls.cv.R: the length of Ytrain is not equal to the Xtrain row number")
     }
     
     ##Ytrain <- Ytrain-1
     
     if ((sum(floor(Ytrain)-Ytrain)!=0)||(sum(Ytrain<0)>0)) {
          stop("Message from mrpls.cv.R: Ytrain is not of valid type")
     }
     
     c <- max(Ytrain)
     eff<-rep(0,(c+1))
     for (i in 0:c) {
          eff[(i+1)]<-sum(Ytrain==i)
     }
     if (sum(eff<=1)>0) {
          stop("Message from mrpls.cv.R: there are not enough samples for each class")
     }
     
     if (c==1) {
          stop("Message from mrpls.cv.R: Ytrain is a binary vector, use rpls.cv.R")
     }
     
     #On hyper parameters range
     
     if ((is.numeric(LambdaRange)==FALSE)||(is.vector(LambdaRange)==FALSE)||(sum(LambdaRange<0)>0)) {
          stop("Message from mrpls.cv.R: LambdaRange is not of valid type")
     }
     
     if ((is.numeric(ncompMax)==FALSE)||(round(ncompMax)-ncompMax!=0)||(ncompMax<0)) {
          stop("Message from mrpls.cv.R: ncompMax is not of valid type")
     }
     
     if ((is.numeric(NbIterMax)==FALSE)||(round(NbIterMax)-NbIterMax!=0)||(NbIterMax<1)) {
          stop("Message from mrpls.cv.R: NbIterMax is not of valid type")
     }
     
     if ((is.numeric(ncores)==FALSE)||(round(ncores)-ncores!=0)||(ncores<1)) {
          stop("Message from mrpls.cv.R: ncores is not of valid type")
     }
     
     
     ## CV LOOP
     ############
     
     #Some initializations
     LambdaRange <- sort(LambdaRange)
     ntrainCV <- dim(Xtrain)[1]-1
     
     LambdaIndex <- 1:length(LambdaRange)
     if (ncompMax <=1) {
          nc <- 1
     }
     if (ncompMax > 1) {
          nc <- ncompMax
     }
     ResCV <- matrix(0,nrow=length(LambdaRange),ncol=nc)
     ResCVbis <- matrix(0,nrow=length(LambdaRange),ncol=nc)
     
     # for (ncv in 1:ntrain) {
     res_cv <- Reduce("rbind", mclapply(1:ntrain, function(ncv) {
          # Determine the data matrix
          
          cvXtrain <- as.matrix(Xtrain[-ncv,])
          cvXtest <- matrix(Xtrain[ncv,],nrow=1)
          p <- dim(cvXtrain)[2]
          r <- min(p,ntrainCV)
          #  Standardize the cvXtrain matrix
          
          Sigma2train <- apply(cvXtrain,2,var)*(ntrainCV-1)/ntrainCV
          if (sum(Sigma2train==0)!=0) {
               if (sum(Sigma2train==0)>(p-2)) {
                    stop("Message from mrpls.cv.R: the procedure stops because, after leaving one sample, number of predictor variables with no null variance is less than 1.")
               }
               cvXtrain <- as.matrix(cvXtrain[,which(Sigma2train!=0)])
               cvXtest <- matrix(cvXtest[,which(Sigma2train!=0)],nrow=1)
               Sigma2train <-Sigma2train[which(Sigma2train!=0)]
               p <- dim(cvXtrain)[2]
               r <- min(p,ntrainCV)
          }
          MeancvXtrain <- apply(cvXtrain,2,mean)
          sXtrain <- sweep(cvXtrain,2,MeancvXtrain,FUN="-")
          sXtrain <- sweep(sXtrain,2,sqrt(Sigma2train),FUN="/")
          
          # Move in the reduced space when necessary
          if (p>ntrainCV) {
               svd.sXtrain <- svd(t(sXtrain))
               r<-length(svd.sXtrain$d[abs(svd.sXtrain$d)>10^(-13)])
               V <- svd.sXtrain$u[,1:r]
               D <- diag(c(svd.sXtrain$d[1:r]))
               U <- svd.sXtrain$v[,1:r]
               sXtrain <- U%*%D
               rm(D)
               rm(U)
               rm(svd.sXtrain)
          }
          sXtest <- sweep(cvXtest,2,MeancvXtrain,FUN="-")
          sXtest <- sweep(sXtest,2,sqrt(Sigma2train),FUN="/")
          if (p>ntrainCV) {
               sXtest <- sXtest%*%V
               rm(V)
          }
          rm(cvXtrain)
          
          #Compute Zblock
          Z <- cbind(rep(1,ntrainCV),sXtrain)
          Zt <- cbind(rep(1,1),sXtest)
          Zbloc <- matrix(0,nrow=ntrainCV*c,ncol=c*(r+1))
          Ztestbloc <- matrix(0,nrow=c,ncol=c*(r+1))
          for (cc in 1:c) {
               row <- (0:(ntrainCV-1))*c+cc
               col <- (r+1)*(cc-1)+1:(r+1)
               Zbloc[row,col] <- Z
               row <- cc
               Ztestbloc[row,col] <- Zt
          }
          rm(Z)
          rm(Zt)
          
          
          # LIndexaux <- LambdaIndex
          # 
          # for (i in LambdaIndex) {
          #      res <- mrplsaux(Ytrain[-ncv],Zbloc,LambdaRange[i],ncompMax,Ztestbloc,NbIterMax=NbIterMax)
          #      if (res$Convergence==0) {
          #           LIndexaux <- LIndexaux[LIndexaux!=i]
          #      }
          #      if (res$Convergence==1) {
          #           ResCV[i,] <- ResCV[i,]+abs(res$hatY-Ytrain[ncv])
          #           ResCVbis[i,] <- ResCVbis[i,]+as.numeric(res$hatY!=Ytrain[ncv])
          #      }
          # }
          # 
          # LambdaIndex <- LIndexaux
          
          resInter <- sapply(LambdaIndex, function(i) {
               res <- mrplsaux(Ytrain[-ncv],Zbloc,LambdaRange[i],ncompMax,Ztestbloc,NbIterMax=NbIterMax)
               return(c(ncv, i, res$Convergence, abs(res$hatY-Ytrain[ncv]), as.numeric(res$hatY!=Ytrain[ncv])))
          })
          return(t(resInter))
     }, mc.cores=ncores))
     #}
     
     ## formatting results
     res_cv <- data.frame(res_cv)
     nc <- ncompMax
     index_nc <- 1:nc
     if(ncompMax==0) {
          nc <- 1
          index_nc <- 0
     }
     colnames(res_cv) <- c("ncv", "LambdaIndex", "convergence", paste0(index_nc), paste0(index_nc))
     
     nconv <- which(res_cv$convergence == 0)
     badLambdaIndex <- unique(res_cv$LambdaIndex[nconv])
     keepIndex <- which(! res_cv$LambdaIndex %in% badLambdaIndex)
     
     if(length(keepIndex)==0) {
          stop("No optimal Lambda for the given LambdaRange")
     }
     
     res_cv <- res_cv[keepIndex,]
     
     res_cv1 <- res_cv[,c(1:3, (1:nc)+3)]
     res_cv2 <- res_cv[,c(1:3, (1:nc)+nc+3)]
     
     resCV <- with(res_cv1, aggregate(res_cv1[,(1:nc)+3], data.frame(LambdaIndex), sum))
     resCVbis <- with(res_cv2, aggregate(res_cv2[,(1:nc)+3], data.frame(LambdaIndex), mean))
     
     resCV <- resCV[order(resCV$LambdaIndex),]
     resCVbis <- resCVbis[order(resCVbis$LambdaIndex),]
     
     ResCV = as.matrix(resCV[,-1])
     ResCVbis = as.matrix(resCVbis[,-1])
     
     
     ## CONCLUDE
     ##############
     # if (length(LambdaIndex)==0) {
     #      stop("No optimal Lambda for the given LambdaRange")
     # }
     
     # ResCVbis = ResCVbis/ntrain
     
     #else
     # ResCV <- ResCV[LambdaIndex,]
     # if (length(LambdaIndex)==1)  {
     #      ResCV <- matrix(ResCV,nrow=1)
     # }
     # Determine optimal Lambda and ncomp
     aux <- which.min(ResCV)
     if (ncompMax <=1) {
          ncomp <- ncompMax
          Lambda <- LambdaRange[LambdaIndex[aux]]
     }
     if (ncompMax > 1) {
          nl <- dim(ResCV)[1]
          ncomp <- (aux-1)%/%nl+1
          Lambda <- LambdaRange[LambdaIndex[(aux-1)%%nl+1]]
     }
     
     ncompRange <- 1:ncompMax
     if(ncompMax==0) {
          ncompRange <- 0
     }
     cv.grid <- melt(ResCVbis, id=c(LambdaRange, ncompRange))
     cv.grid <- data.frame(t(sapply(1:nrow(cv.grid), function(ind) {
          return(c(LambdaRange[cv.grid[ind,1]], ncompRange[cv.grid[ind,2]], cv.grid[ind,3]))
     })))
     colnames(cv.grid) <- c("lambda", "ncomp", "error")
     
     return(list(ncomp=ncomp, Lambda=Lambda, cv.grid=cv.grid))

}

