#OUTLINE-----------------------------------------------
#Function 1: Data generation for simulation Sudy
#Function 2: Estimators based on Penalty functions
#Function 3: Calculation of evaluation metrics
#------------------------------------------------------
#Requaired packages------------------------------------
library(glmnet)
library(ncvreg)
library(Rlab)
library(pracma)
library(psych)
library(pracma)
library(ggplot2)
library(gridExtra)
library(reshape)
library(corrplot)
library(MASS)
library(tidyverse)
library(GGally)
library(stats)
library(cowplot)
library(latticeExtra) 
library(lattice)
library(viridis)
library(plotly)
library(heatmaply)
#FUNCTIONS----------------------------------------------------------------------
simdata <- function(n,pn,sc,rho){
  #n : Sample size
  #pn : No covariates in parametric component should be pn=(300, 1000, 5000)
  #sc: Selection of scenario, sc=1; \beta=c(5,-3), sc=2; \beta=c(1,-0.5)
  #rho: Correlation coefficient rho=(0.80,0.95)
  
  # Generation of correlated parametric covariates X----------------------------
  # create the variance covariance matrix
  sigma1<-diag(5)
  for (i in 1:5){
    for (j in 1:5){
      if (i==j){sigma1[i,j] <- 1}
      else{
        sigma1[i,j] <- rho 
      }
    }
  }
  mu<-c(runif(5))
  df<-as.data.frame(mvrnorm(n=n, mu=mu, Sigma=sigma1))
  p1 <- pn-5
  Xa <-  matrix(rnorm(n*p1), ncol=p1)
  X <- matrix(c(df$V1,df$V2,df$V3,df$V4,df$V5,Xa),n,pn)
  #----------------------------------------------------------------------------
  #Generation of Regression coefficients---------------------------------------
  beta <- 0
  if(sc==1){
    for (j in 1:pn){
      beta[j] <- 5*I(j>=1 && j<=10)-0.3*I(j>15 && j<=25)
    }
  }
  if(sc==2){
    for (j in 1:pn){
      beta[j] <- 1*I(j>=1 && j<=10)-0.5*I(j>15 && j<=25)
    }
  } 
  #Generation of nonparametric covariate and function---------------------------
  t <- 0
  for (j in 1:n){
    t[j] <- 2.4*((j-0.5)/n)
  }
  f <- -t*sin(-t^2)
  #-----------------------------------------------------------------------------
  #Generation of High-dimensional PLM in equation (4)---------------------------
  e <- rnorm(n,sd=0.5)
  y <- X%*%beta+f+e
  
  #X <- (X-mean(y))/(max(y)-min(y))
  #y <- (y-mean(y))/(max(y)-min(y))
  #t <- (t-mean(y))/(max(y)-min(y))
  #f <- (f-mean(y))/(max(y)-min(y))
  #-----------------------------------------------------------------------------
  out      <- new.env()
  out$X    <- X
  out$beta <- beta
  out$t    <- t
  out$f    <- f
  out$y    <- y
  out$pn   <- pn
  out$n    <- n
  out$e    <- e 
  return(out)
}
#Fig6 Figure for generated data-------------------------------------------------
aicfunc <- function(y,obj,lam,f){
  x       <- obj$X_S1 
  pn      <-length(obj$S1) 
  S       <- Smatrix(t,2,lam)
  xtil    <- (diag(n)-S)%*%x
  ytil    <- (diag(n)-S)%*%y
  betahat <- solve(t(xtil)%*%xtil+0.1*diag(pn),tol=1e-100)%*%t(xtil)%*%ytil
  fhat    <- S%*%(y-x%*%betahat)
  S1      <- obj$S1
  var     <- var((f-fhat))
  score   <- log(var)+1+((2*pn+1)/(n-pn-2))
  return(score)
}
#-------------------------------------------------------------------------------
fig6 <- function(mydata){
  x     <- mydata$X
  y     <- mydata$y
  t     <- mydata$t
  f     <- mydata$f 
  
  fe   <- f+mydata$e
  df1  <- data.frame(x[,1:7])
  dffe <- data.frame(t,fe)
  dff  <- data.frame(t,f) 
  ggpairs(df1)        
  pheat <- heatmaply(x[,1:50],dendrogram = "none", xlab = "", ylab = "", main = "",
                     scale = "column", margins = c(20,10,20,20), grid_color = "white",
                     grid_width = 0.00001, titleX = FALSE,hide_colorbar = TRUE,
                     branches_lwd = 0.1, label_names = c("Country", "Feature:", "Value"),
                     fontsize_row = 5, fontsize_col = 5,heatmap_layers = theme(axis.line=element_blank()))
  
  pf <- ggplot()+geom_point(data=dffe,aes(x=t,y=fe))+geom_line(data=dff,aes(x=t,y=f),lwd=1)+xlab("t")+ylab("f(t)")+ theme(plot.title = element_text(hjust = 0.5))
  par(mar=c(5,5,5,5))
  beta <- mydata$beta[1:25]
  names(beta) <- c("1","2","3","4","5","6","7","8","9","10",".",".",".",".",".","16","17","18","19","20","21","22","23","24","25")
  barplot(beta,ylab="Beta",xlab="No. Coefficient")
  
  print(ggpairs(df1))
  print(pheat)
  print(pf)
}
#Adaptive-Lasso-----------------------------------------------------------------
adalasso<-function(x,y){
  #ymean <- mean(y)
  #  y <- y-mean(y)  
  #  xmean <- colMeans(x)
  #  xnorm <- sqrt(n-1)*apply(x,2,sd)
  #  x <- scale(x, center = xmean, scale = xnorm)
  n <- length(y)
  pn <- nrow(x)
  beta.init <- solve(t(x)%*%x,tol=1e-100)%*%t(x)%*%y
  
  # calculate weights
  w  <- abs(beta.init)  
  x2 <- scale(x, center=FALSE, scale=1/w)  
  
  la_eq <- cv.glmnet(x2,y,nfolds=20,intercept=F)
  coefs <- coef(la_eq,c(la_eq$lambda.min))
  S1a  <- coefs@i+1                              #Subset of Strong Signals (index)
  S1r  <- coefs@i
  if ((length(S1a))>n){
    S1_beta <- coefs@x[1:n]
    S1 <- S1r[1:n]
  } else {
    S1_beta <- coefs@x
    S1 <- S1r
  }
  S1_beta <- coefs@x
  X_S1 <- x2[,c(S1)]                               #Subset of Strong signals
  X_S1c <- x[,-S1]                                #Subset of complement
  alas <- new.env()
  
  alas$X_S1    <- X_S1
  alas$X_S1c   <- X_S1c
  alas$S1_beta <- S1_beta
  alas$S1      <- S1
  alas$coefs   <- coefs
  return(alas)
}
#Lasso--------------------------------------------------------------------------
lasso   <- function(x,y){
  #ymean <- mean(y)
  #y <- y-mean(y)  
  #xmean <- colMeans(x)
  #xnorm <- sqrt(n-1)*apply(x,2,sd)
  #x <- scale(x, center = xmean, scale = xnorm)
  
  la_eq <- cv.glmnet(x,y,nfolds=10,intercept=F)
  coefs <- coef(la_eq,c(la_eq$lambda.min))
  S1a  <- coefs@i+1                   #Subset of Strong Signals (index)
  S1r  <- coefs@i
  if ((length(S1a))>n){
    S1_beta <- coefs@x[1:n]
    S1 <- S1r[1:n]
  } else {
    S1_beta <- coefs@x
    S1 <- S1r
  }
  
  
  #S1_beta <- coefs@x
  X_S1 <- x[,S1]                               #Subset of Strong signals
  X_S1c <- x[,-S1]                                #Subset of complement
  las <- new.env()
  
  las$X_S1    <- X_S1
  las$X_S1c   <- X_S1c
  las$S1_beta <- S1_beta
  las$S1      <- S1
  las$coefs   <-coefs 
  
  return(las)
}
#SCAD---------------------------------------------------------------------------
scad    <- function(x,y){
  # n <- length(y)
  #dimx <- dim(x)
  #pn <- dimx[2]
  #ymean <- mean(y)
  #y <- y-mean(y)  
  #xmean <- colMeans(x)
  #xnorm <- sqrt(n-1)*apply(x,2,sd)
  #x <- scale(x, center = xmean, scale = xnorm)
  lambdaseq <- seq(1,0.1,length.out=100)
  scadm <- ncvreg(x,y,family="gaussian",penalty="SCAD",lambda=lambdaseq,intercept=F)
  S1_beta <- as.numeric(scadm$beta[,80])
  coefs <- S1_beta
  S1_beta2 <- S1_beta[S1_beta!=0]
  lobeta <- length(S1_beta2)
  
  if (lobeta>=n){
    S1_beta2 <- S1_beta2[1:n]
  } else {
    S1_beta2 <- S1_beta2
  }
  S1 <- 0
  for (i in 1:pn){
    if (S1_beta[i]!=0){
      S1[i] <- i
    }
  }
  S1a   <- na.omit(S1)
  ls1a  <-length(S1a) 
  S1_beta2 <- S1_beta2[1:ls1a]
  X_S1 <- x[,S1a] 
  X_S1c <- x[,-S1a]
  
  scad <- new.env()
  scad$S1_beta <- S1_beta2
  scad$X_S1    <- X_S1
  scad$X_S1c   <- X_S1c
  scad$S1      <- S1a
  scad$coefs   <- coefs 
  return(scad)
} 
#SMOOTHING MATRIX FOR SMOOTHING SPLINES-----------------------------------------
Smatrix = function(x, df,lam){
  n = length(x);
  A = matrix(0, n, n);
  for(i in 1:n){
    y = rep(0, n); y[i]=1;
    yi = smooth.spline(x, y, df=df,spar=lam)$y;
    A[,i]= yi;
  }
  return(A)
}    #FROM Feng Liang https://publish.illinois.edu/liangf/teaching/stat-424/rcode-smoothing-splines/
#WEIGTHTED RIDGE REGRESSION-----------------------------------------------------
WR <- function(obj,S,x,y,lambdar){
  #obj: object obtained from penalty function
  #thr: ridge threshold 
  #STEP 1. Restricted Estimator (beta^RE)
  S2index <- 0
  dimx    <-dim(x)
  pn      <-dimx[2] 
  n       <-length(y)
  
  X_S1    <- (diag(n)-S)%*%obj$X_S1
  X_S1c   <- (diag(n)-S)%*%obj$X_S1c
  S1      <- length(obj$S1)
  beta_S1 <- obj$S1_beta[1:S1]
  Y       <- (diag(n)-S)%*%y
  ps1c    <-pn-length(obj$S1) 
  
  betaRE <-solve(t(X_S1)%*%X_S1,tol=1e-100)%*%t(X_S1)%*%Y 
  #STEP 2. WR estimators (beta^WR)
  betaWR   <- solve(t(X_S1c)%*%X_S1c+lambdar*diag(ps1c),tol=1e-100)%*%t(X_S1c)%*%Y
  betaWR2  <- solve(t(X_S1)%*%X_S1+lambdar*diag(S1),tol=1e-100)%*%t(X_S1)%*%Y
  
  thr_seq  <- seq(min(abs(betaWR))*3,max(abs(betaWR))-0.1,length.out=50)                       #involves 50 elements
  k <- 5
  MPE  <- matrix(0,k,50)
  S2p2 <- matrix(0,k,50)
  AMPE <- 0
  #Selection of ridge threshold with k-fold
  for (j in 1:50){
    k <- 5
    nk <- n/k
    r <- 0
    for (j2 in 1:k){
      testXS1  <-X_S1[(j2+(r-1)*I(j2>1)):(nk+r),] 
      trainXS1 <-X_S1[-c((j2+(r-1)*I(j2>1)):(nk+r)),] 
      xtr      <-x[-c((j2+(r-1)*I(j2>1)):(nk+r)),] 
      xtest    <-x[(j2+(r-1)*I(j2>1)):(nk+r),] 
      Ytrain   <-y[-c((j2+(r-1)*I(j2>1)):(nk+r))]
      Ytest    <-y[c((j2+(r-1)*I(j2>1)):(nk+r))]
      
      
      S2indexp <- 0
      beta_S2a <- as.numeric(betaWR*I(abs(betaWR)>thr_seq[j]))
      for (i in 1:ps1c){
        if (beta_S2a[i]!=0)
          S2indexp[i] <- i
      }
      S2indexp <-na.omit(S2indexp)   
      
      beta_S2p   <- beta_S2a[beta_S2a!=0]       #Weak signals 
      S2p        <- length(beta_S2p)            #no. S2 elements
      S2p2[j2,j] <- S2p
      X_S2p      <- X_S1c[-c((j2+(r-1)*I(j2>1)):(nk+r)),S2indexp]            #X matrix of weak signals
      #-------------------------------------------------------------------------------
      #STEP 3 PSE estimator
      dimt     <- dim(trainXS1)
      ntrain   <- dimt[1]
      Mp       <- diag(ntrain)-(trainXS1%*%solve(t(trainXS1)%*%trainXS1,tol=1e-100)%*%t(trainXS1))
      varp     <- var(Ytrain-trainXS1%*%betaRE)
      Trp      <- (t(beta_S2p)%*%(t(X_S2p)%*%Mp%*%X_S2p)%*%beta_S2p)/varp
      betaPSEp <- (beta_S1-(min(((S2p-2)/Trp),1)))*(betaWR2-betaRE)
      
      fit <- testXS1%*%betaPSEp
      MPE[j2,j] <- mean((Ytest-fit)^2)
      r        <- r+nk
    }
    AMPE[j] <-mean(MPE[,j]) 
  }
  
  for (j3 in 1:50){
    if (AMPE[j3]==min(AMPE)){
      pnt <- thr_seq[j3]
    }
  }
  #plot(thr_seq,AMPE,type="l",lwd=2,main="Selection of optimal ridge-threshold with k-fold CV",ylab="Prediction error (TEST)",xlab="Ridge-threshold")
  #par(new=TRUE)
  #grid()
  #abline(v=pnt,col="red",lty=2,lwd=2)
  opt.thr <- pnt
  #-------------------------------------------------------------------------------
  X_S1    <- (diag(n)-S)%*%obj$X_S1
  X_S1c   <- (diag(n)-S)%*%obj$X_S1c
  S1      <- length(obj$S1)
  beta_S1 <- obj$S1_beta
  Y       <- (diag(n)-S)%*%y
  S2index <- 0
  beta_S2a <- as.numeric(betaWR*I(abs(betaWR)>opt.thr))
  for (i in 1:ps1c){
    if (beta_S2a[i]!=0){
      S2index[i] <- i
    }
  }
  S2index <-na.omit(S2index)   
  
  beta_S2 <- beta_S2a[beta_S2a!=0]       #Weak signals 
  S2      <- length(beta_S2)                  #no. S2 elements
  X_S2    <- X_S1c[,S2index]             #X matrix of weak signals
  totalsig <- S1+S2
  #if (totalsig>=n) stop('Threshold is too small, please give another value')
  #-------------------------------------------------------------------------------
  #STEP 3 PSE estimator
  M       <- diag(n)-(X_S1%*%solve(t(X_S1)%*%X_S1,tol=1e-100)%*%t(X_S1))
  var     <- var(y-obj$X_S1%*%betaRE)
  Tr      <- (t(beta_S2)%*%(t(X_S2)%*%M%*%X_S2)%*%beta_S2)/var
  betaPSE <- betaWR2-((min(((S2-2)/Tr),1)))*(abs(betaWR2-betaRE))
  
  WR <- new.env()
  
  WR$betaRE    <- betaRE
  WR$betaWR    <- betaWR
  WR$betaWR2   <- betaWR2
  WR$X_S2      <- X_S2
  WR$beta_S2   <- beta_S2
  WR$betaPSE   <- betaPSE 
  WR$S2        <- S2 
  WR$opt.thr   <- opt.thr 
  return(WR)
}
#EVALUATION METRICS------------------------------------------------------------
rmsebeta <- function(realbeta,estbeta){
  k <- length(realbeta)
    rmse_beta <-sqrt((1/k)*t(realbeta-estbeta)%*%(realbeta-estbeta)) #RMSE for BETA
return(rmse_beta)
    }
rsq <- function(y,yhat,p1){
  n <- length(y)
  rsq       <- 1-((sum((y-yhat)^2)/(n-p1))/sum(((y-mean(y))^2)/(n-1)))  #R-Square for Model  
return(rsq)
  }  
recall <- function(realbeta,estbeta){
  #Values of Confussion Matrix----------------------------------------------------
  k <- length(realbeta)
  a <- 0
  b <- 0
  c <- 0
  d <- 0
  for (i in 1:k){
    if (estbeta[i]!=0 & realbeta[i]!=0){
      d <- d+1
    }
    if (estbeta[i]==0 & realbeta[i]==0){
      a <- a+1
    }
    if (estbeta[i]!=0 & realbeta[i]==0){
      c <- c+1
    }
    if (estbeta[i]==0 & realbeta[i]!=0){
      b <- b+1
    }
  }
  acc <- (a+d)/(a+b+c+d)
  sens <- d/(c+d)
  spec <- a/(a+b)
  if (spec=="NaN"){
    spec <- 0.01
  }
  Gscore <- sqrt(sens*spec)
res <- new.env()
res$accuracy    <- acc
res$sensitivity <- sens
res$specifcity  <- spec 
res$Gscore      <- Gscore 

return(res)
}
msef <- function(f,fhat){
  mse <- mean((f-fhat)^2)
  return(mse)
}  #Evaluation of Nonparametric compoenent-----------------------------------------
#--------------------------------------------------------------------------------------------------------------
ni <- c(50,100,120)
pni <- c(300,350,400)
rhoi <- c(0.80,0.95)
sim <- 100                        #Number of simulation
for (s1 in 1:3){                  #Sample size index
  for (s2 in 1:3){                #pn index
        if (s1==1){n=50} 
        if (s1==2){n=100}  
        if (s1==3){n=200}
        if (s2==1){pn=300} 
        if (s2==2){pn=600} 
        if (s2==3){pn=1200}
        rho <- 0.80
        sc  <- 1
        #ZERO MATRICES----------------------------------------------------------
        rmsebeta_las    <- 0
        rmsebeta_alas   <- 0
        rmsebeta_scad   <- 0
        
        Rsq_las         <- 0
        Rsq_alas        <- 0
        Rsq_scad        <- 0
        
        msef_las         <- 0
        msef_alas        <- 0
        msef_scad        <- 0
        
        sens_las <- 0
        spec_las <- 0
        acc_las  <- 0
        G_las    <- 0
        sens_alas <- 0
        spec_alas <- 0
        acc_alas  <- 0
        G_alas    <- 0
        sens_scad <- 0
        spec_scad <- 0
        acc_scad  <- 0
        G_scad    <- 0
        
        strong_las  <- 0
        weak_las    <- 0
        strong_alas <- 0
        weak_alas   <- 0 
        strong_scad <- 0 
        weak_scad   <- 0
        
        Est.f.LassoWR  <- matrix(0,n,sim)
        Est.f.aLassoWR <- matrix(0,n,sim)
        Est.f.SCADWR   <- matrix(0,n,sim)
        
        Est.f.Lasso  <- matrix(0,n,sim)
        Est.f.aLasso <- matrix(0,n,sim)
        Est.f.SCAD   <- matrix(0,n,sim)
        
        Est.y.LassoWR  <- matrix(0,n,sim)
        Est.y.aLassoWR <- matrix(0,n,sim)
        Est.y.SCADWR   <- matrix(0,n,sim)
        
        Est.y.Lasso  <- matrix(0,n,sim)
        Est.y.aLasso <- matrix(0,n,sim)
        Est.y.SCAD   <- matrix(0,n,sim)
        #------------------------------------------------------------------------------
        no.seed <- 19345
        aa <- 1
        for (s in 1:sim){
          #Set seed number needs to be changed in every repeat----------------------------
          set.seed(no.seed+aa)
          mydata <- simdata(n,pn,sc,rho)
          set.seed(no.seed+s)
          #TEST AREA (To be Deleted)------------------------------------------------------
          x     <- mydata$X
          y     <- mydata$y
          t     <- mydata$t
          f     <- mydata$f 
          nol   <- 25
          aiclas  <- 0
          aicalas <- 0
          aicscad <- 0
          
          las     <- lasso(x,y)
          alas    <- adalasso(x,y)
          scadm   <- scad(x,y)
          
          lams    <- seq(0.1,1.5,length.out=nol)
          for (j in 1:nol){
            aiclas[j]  <- aicfunc(y,las,lams[j],f)
            aicalas[j]  <- aicfunc(y,alas,lams[j],f)
            aicscad[j]  <- aicfunc(y,scadm,lams[j],f)
          }
          #PLOTS for smoothing parameter selection (AICc)---------------------------------
          dflas  <- data.frame(lams,aiclas)
          dfalas <- data.frame(lams,aicalas)
          dfscad <- data.frame(lams,aicscad)
          #ggplot()+geom_line(data=dflas,aes(x=lams,aiclas,col="Lasso_WR"),lwd=1)+geom_line(data=dfalas,aes(x=lams,aicalas,col="Adapt.-Lasso_WR"),lwd=1)+geom_line(data=dfscad,aes(x=lams,aicscad,col="SCAD_WR"),lwd=1)+xlab("lambda_spline")+ylab("Improved_AIC score")+ggtitle("(e) n=100, pn=3000, rho=0.95")+theme(legend.title = element_blank(),legend.position = "bottom")
          #-------------------------------------------------------------------------------
          for (j2 in 1:nol){
            if (aiclas[j2]==min(aiclas)){
              lam_las <-  lams[j2]
            }
            if (aicalas[j2]==min(aicalas)){
              lam_alas <-  lams[j2]
            }
            if (aicscad[j2]==min(aicscad)){
              lam_scad <-  lams[j2]
            }
          }
          S_las  <- Smatrix(t,2,lam_las)
          S_alas <- Smatrix(t,2,lam_alas)
          S_scad <- Smatrix(t,2,lam_scad)
          
          
          ytil_las  <- (diag(n)-S_las)%*%y
          ytil_alas <- (diag(n)-S_alas)%*%y
          ytil_scad <- (diag(n)-S_scad)%*%y
          
          
          xtil_las  <- (diag(n)-S_las)%*%x
          xtil_alas <- (diag(n)-S_alas)%*%x
          xtil_scad <- (diag(n)-S_scad)%*%x
          
          newX_las  <- (diag(n)-S_las)%*%las$X_S1
          newX_alas <- (diag(n)-S_alas)%*%alas$X_S1
          newX_scad <- (diag(n)-S_scad)%*%scadm$X_S1
          #Estimation based on Penalty functions
          las_beta      <- solve(t(newX_las)%*%newX_las+0.01*diag(las$S1))%*%t(newX_las)%*%ytil_las
          fhat_las      <- S_las%*%(y-las$X_S1%*%las_beta)
          alas_beta     <- solve(t(newX_alas)%*%newX_alas+0.01*diag(alas$S1))%*%t(newX_alas)%*%ytil_alas
          fhat_alas     <- S_alas%*%(y-alas$X_S1%*%alas_beta)
          scad_beta     <- solve(t(newX_scad)%*%newX_scad+0.01*diag(scadm$S1))%*%t(newX_scad)%*%ytil_scad
          fhat_scad     <- S_scad%*%(y-scadm$X_S1%*%scad_beta)
          
          yhat_las  <- las$X_S1%*%las_beta+fhat_las
          yhat_alas <- alas$X_S1%*%alas_beta+fhat_alas
          yhat_scad <- scadm$X_S1%*%scad_beta+fhat_scad
          
          WRobj_las  <- WR(las,S_las,x,y,0.01)
          WRobj_alas <- WR(alas,S_alas,x,y,0.01)
          WRobj_scad <- WR(scadm,S_scad,x,y,0.01)
          
          WR_fhat_las  <- S_las%*%(y-las$X_S1%*%WRobj_las$betaPSE)
          WR_fhat_alas <- S_alas%*%(y-alas$X_S1%*%WRobj_alas$betaPSE)
          WR_fhat_scad <- S_scad%*%(y-scadm$X_S1%*%WRobj_scad$betaPSE)
          
          WR_yhat_las  <- las$X_S1%*%WRobj_las$betaPSE+fhat_las
          WR_yhat_alas <- alas$X_S1%*%WRobj_alas$betaPSE+fhat_alas
          WR_yhat_scad <- scadm$X_S1%*%WRobj_scad$betaPSE
          #---------------------------------------------------------------------
          las_index       <- las$S1
          real_beta_las   <- mydata$beta[las_index]
          alas_index      <- alas$S1
          real_beta_alas  <- mydata$beta[alas_index]
          scad_index      <- scadm$S1
          real_beta_scad  <- mydata$beta[scad_index]
          
          rmsebeta_las[s]  <- rmsebeta(real_beta_las,WRobj_las$betaPSE)
          rmsebeta_alas[s] <- rmsebeta(real_beta_alas,WRobj_alas$betaPSE)
          rmsebeta_scad[s] <- rmsebeta(real_beta_scad,WRobj_scad$betaPSE)
          
          Rsq_las[s]  <- rsq(mydata$y,WR_yhat_las,length(las_index))
          Rsq_alas[s] <- rsq(mydata$y,WR_yhat_alas,length(alas_index))
          Rsq_scad[s] <- rsq(mydata$y,WR_yhat_scad,length(scad_index))
          
          msef_las[s]  <- msef(mydata$f,WR_fhat_las)
          msef_alas[s] <- msef(mydata$f,WR_fhat_alas)
          msef_scad[s] <- msef(mydata$f,WR_fhat_scad)
          
          recall_las  <- recall(real_beta_las,WRobj_las$betaPSE)
          sens_las[s] <- recall_las$sensitivity
          spec_las[s] <- recall_las$specifcity
          acc_las[s]  <- recall_las$accuracy
          G_las[s]    <- recall_las$Gscore
          
          recall_alas  <- recall(real_beta_alas,WRobj_alas$betaPSE)
          sens_alas[s] <- recall_alas$sensitivity
          spec_alas[s] <- recall_alas$specifcity
          acc_alas[s]  <- recall_alas$accuracy
          G_alas[s]    <- recall_alas$Gscore
          
          recall_scad  <- recall(real_beta_scad,WRobj_scad$betaPSE)
          sens_scad[s] <- recall_scad$sensitivity
          spec_scad[s] <- recall_scad$specifcity
          acc_scad[s]  <- recall_scad$accuracy
          G_scad[s]    <- recall_scad$Gscore
          
          strong_las[s]<- length(WRobj_las$betaPSE)
          weak_las[s]  <- WRobj_las$S2
          
          strong_alas[s]<- length(WRobj_alas$betaPSE)
          weak_alas[s]  <- WRobj_alas$S2
          
          strong_scad[s]<- length(WRobj_scad$betaPSE)
          weak_scad[s]  <- WRobj_scad$S2
          #---------------------------------------------------------------------
          Est.f.LassoWR[,s]  <- WR_fhat_las
          Est.f.aLassoWR[,s] <- WR_fhat_alas
          Est.f.SCADWR[,s]   <- WR_fhat_scad
          
          Est.f.Lasso[,s]  <- fhat_las
          Est.f.aLasso[,s] <- fhat_alas
          Est.f.SCAD[,s]   <- fhat_scad
          
          Est.y.Lasso[,s]  <- yhat_las
          Est.y.aLasso[,s] <- yhat_alas
          Est.y.SCAD[,s]   <- yhat_scad
          
          Est.y.LassoWR[,s]  <- WR_yhat_las
          Est.y.aLassoWR[,s] <- WR_yhat_alas
          Est.y.SCADWR[,s]   <- WR_yhat_scad
          
          #plot(f,type="l",ylim=c(min(f),max(f)))
          #par(new=TRUE)
          #plot(WR_fhat_las,type="l",ylim=c(min(f),max(f)),col="red")
          #par(new=TRUE)
          #plot(WR_fhat_alas,type="l",ylim=c(min(f),max(f)),col="blue")
          #par(new=TRUE)
          #plot(WR_fhat_scad,type="l",ylim=c(min(f),max(f)),col="green")
          aa <- aa+1
          message(s, "th Simulation ends for n=",n, " pn=",pn, " and rho=", rho, " in Scenario ",sc)
        }
        if (s2==1 & s1==1){
          WR.lasso.fhat.p300.n50  <-rowMeans(Est.f.LassoWR) 
          WR.alasso.fhat.p300.n50 <-rowMeans(Est.f.aLassoWR)
          WR.scad.fhat.p300.n50   <-rowMeans(Est.f.SCADWR)
          
          rmsebeta_las.p300.n50  <- mean(rmsebeta_las)
          rmsebeta_alas.p300.n50 <- mean(rmsebeta_alas)
          rmsebeta_scad.p300.n50 <- mean(rmsebeta_scad)
          
          Rsq_las.p300.n50  <- mean(Rsq_las)
          Rsq_alas.p300.n50 <- mean(Rsq_alas)
          Rsq_scad.p300.n50 <- mean(Rsq_scad)
          
          msef_las.p300.n50  <- mean(msef_las)
          msef_alas.p300.n50 <- mean(msef_alas)
          msef_scad.p300.n50 <- mean(msef_scad)
          
          sens_las.p300.n50 <- mean(sens_las)
          spec_las.p300.n50 <- mean(spec_las) 
          acc_las.p300.n50  <- mean(acc_las)
          G_las.p300.n50    <- mean(G_las)
          
          sens_alas.p300.n50 <- mean(sens_alas)
          spec_alas.p300.n50 <- mean(spec_alas) 
          acc_alas.p300.n50  <- mean(acc_alas)
          G_alas.p300.n50    <- mean(G_alas)
          
          sens_scad.p300.n50 <- mean(sens_scad)
          spec_scad.p300.n50 <- mean(spec_scad) 
          acc_scad.p300.n50  <- mean(acc_scad)
          G_scad.p300.n50    <- mean(G_scad)
          
          stong_signals.p300.n50 <- data.frame((strong_las),(strong_alas),(strong_scad))
          weak_signals.p300.n50 <- data.frame((weak_las),(weak_alas),(weak_scad))
          
          
          RMSE_beta.p300.n50             <- data.frame(rmsebeta_las.p300.n50,rmsebeta_alas.p300.n50,rmsebeta_scad.p300.n50)
          colnames(RMSE_beta.p300.n50)   <- c("Lasso","Adaptive-Lasso","SCAD")
          
          RSQ.p300.n50                   <- data.frame(Rsq_las.p300.n50,Rsq_alas.p300.n50,Rsq_scad.p300.n50)
          colnames(RSQ.p300.n50)         <- c("Lasso","Adaptive-Lasso","SCAD")
          
          MSEF.p300.n50                  <- data.frame(msef_las.p300.n50,msef_alas.p300.n50,msef_scad.p300.n50)
          colnames(MSEF.p300.n50)        <- c("Lasso","Adaptive-Lasso","SCAD")
          
          metrics_all.p300.n50           <- rbind(RMSE_beta.p300.n50,RSQ.p300.n50,MSEF.p300.n50)
          rownames(metrics_all.p300.n50) <- c("REMSE","R-square","MSEf") 
          
          RECALL.las.p300.n50                     <- data.frame(sens_las.p300.n50,spec_las.p300.n50,acc_las.p300.n50,G_las.p300.n50)
          colnames(RECALL.las.p300.n50)           <- c("Sensitivity","Specificity","Accuracy","G-score")
          RECALL.alas.p300.n50                    <- data.frame(sens_alas.p300.n50,spec_alas.p300.n50,acc_alas.p300.n50,G_alas.p300.n50)
          colnames(RECALL.alas.p300.n50)          <- c("Sensitivity","Specificity","Accuracy","G-score")
          RECALL.scad.p300.n50                    <- data.frame(sens_scad.p300.n50,spec_scad.p300.n50,acc_scad.p300.n50,G_scad.p300.n50)
          colnames(RECALL.scad.p300.n50)          <- c("Sensitivity","Specificity","Accuracy","G-score")
          RECALL_all.p300.n50                     <- rbind(RECALL.las.p300.n50,RECALL.alas.p300.n50,RECALL.scad.p300.n50)
          rownames(RECALL_all.p300.n50)           <- c("Lasso","Adaptive-Lasso","SCAD")
          #-PLOTS of fs---------------------------------------------------------
          plot(f,type="l",ylim=c(min(f),max(f)))
          par(new=TRUE)
          plot(WR.lasso.fhat.p300.n50,type="l",ylim=c(min(f),max(f)),col="red")
          par(new=TRUE)
          plot(WR.alasso.fhat.p300.n50,type="l",ylim=c(min(f),max(f)),col="blue")
          par(new=TRUE)
          plot(WR.scad.fhat.p300.n50,type="l",ylim=c(min(f),max(f)),col="green")
        }
        if (s2==2 & s1==1){
          WR.lasso.fhat.p1000.n50  <-rowMeans(Est.f.LassoWR) 
          WR.alasso.fhat.p1000.n50 <-rowMeans(Est.f.aLassoWR)
          WR.scad.fhat.p1000.n50   <-rowMeans(Est.f.SCADWR)
          
          rmsebeta_las.p1000.n50  <- mean(rmsebeta_las)
          rmsebeta_alas.p1000.n50 <- mean(rmsebeta_alas)
          rmsebeta_scad.p1000.n50 <- mean(rmsebeta_scad)
          
          Rsq_las.p1000.n50  <- mean(Rsq_las)
          Rsq_alas.p1000.n50 <- mean(Rsq_alas)
          Rsq_scad.p1000.n50 <- mean(Rsq_scad)
          
          msef_las.p1000.n50  <- mean(msef_las)
          msef_alas.p1000.n50 <- mean(msef_alas)
          msef_scad.p1000.n50 <- mean(msef_scad)
          
          sens_las.p1000.n50 <- mean(sens_las)
          spec_las.p1000.n50 <- mean(spec_las) 
          acc_las.p1000.n50  <- mean(acc_las)
          G_las.p1000.n50    <- mean(G_las)
          
          sens_alas.p1000.n50 <- mean(sens_alas)
          spec_alas.p1000.n50 <- mean(spec_alas) 
          acc_alas.p1000.n50  <- mean(acc_alas)
          G_alas.p1000.n50    <- mean(G_alas)
          
          sens_scad.p1000.n50 <- mean(sens_scad)
          spec_scad.p1000.n50 <- mean(spec_scad) 
          acc_scad.p1000.n50  <- mean(acc_scad)
          G_scad.p1000.n50    <- mean(G_scad)
          
          RMSE_beta.p1000.n50             <- data.frame(rmsebeta_las.p1000.n50,rmsebeta_alas.p1000.n50,rmsebeta_scad.p1000.n50)
          colnames(RMSE_beta.p1000.n50)   <- c("Lasso","Adaptive-Lasso","SCAD")
          
          RSQ.p1000.n50                   <- data.frame(Rsq_las.p1000.n50,Rsq_alas.p1000.n50,Rsq_scad.p1000.n50)
          colnames(RSQ.p1000.n50)         <- c("Lasso","Adaptive-Lasso","SCAD")
          
          MSEF.p1000.n50                  <- data.frame(msef_las.p1000.n50,msef_alas.p1000.n50,msef_scad.p1000.n50)
          colnames(MSEF.p1000.n50)        <- c("Lasso","Adaptive-Lasso","SCAD")
          
          metrics_all.p1000.n50           <- rbind(RMSE_beta.p1000.n50,RSQ.p1000.n50,MSEF.p1000.n50)
          rownames(metrics_all.p1000.n50) <- c("REMSE","R-square","MSEf") 
          
          stong_signals.p1000.n50 <- data.frame((strong_las),(strong_alas),(strong_scad))
          weak_signals.p1000.n50 <- data.frame((weak_las),(weak_alas),(weak_scad))
          
          RECALL.las.p1000.n50                     <- data.frame(sens_las.p1000.n50,spec_las.p1000.n50,acc_las.p1000.n50,G_las.p1000.n50)
          colnames(RECALL.las.p1000.n50)           <- c("Sensitivity","Specificity","Accuracy","G-score")
          RECALL.alas.p1000.n50                    <- data.frame(sens_alas.p1000.n50,spec_alas.p1000.n50,acc_alas.p1000.n50,G_alas.p1000.n50)
          colnames(RECALL.alas.p1000.n50)          <- c("Sensitivity","Specificity","Accuracy","G-score")
          RECALL.scad.p1000.n50                    <- data.frame(sens_scad.p1000.n50,spec_scad.p1000.n50,acc_scad.p1000.n50,G_scad.p1000.n50)
          colnames(RECALL.scad.p1000.n50)          <- c("Sensitivity","Specificity","Accuracy","G-score")
          RECALL_all.p1000.n50                     <- rbind(RECALL.las.p1000.n50,RECALL.alas.p1000.n50,RECALL.scad.p1000.n50)
          rownames(RECALL_all.p1000.n50)           <- c("Lasso","Adaptive-Lasso","SCAD")
        }
        if (s2==3 & s1==1){
          WR.lasso.fhat.p3000.n50  <-rowMeans(Est.f.LassoWR) 
          WR.alasso.fhat.p3000.n50 <-rowMeans(Est.f.aLassoWR)
          WR.scad.fhat.p3000.n50   <-rowMeans(Est.f.SCADWR)
          
          rmsebeta_las.p3000.n50  <- mean(rmsebeta_las)
          rmsebeta_alas.p3000.n50 <- mean(rmsebeta_alas)
          rmsebeta_scad.p3000.n50 <- mean(rmsebeta_scad)
          
          Rsq_las.p3000.n50  <- mean(Rsq_las)
          Rsq_alas.p3000.n50 <- mean(Rsq_alas)
          Rsq_scad.p3000.n50 <- mean(Rsq_scad)
          
          msef_las.p3000.n50  <- mean(msef_las)
          msef_alas.p3000.n50 <- mean(msef_alas)
          msef_scad.p3000.n50 <- mean(msef_scad)
          
          sens_las.p3000.n50 <- mean(sens_las)
          spec_las.p3000.n50 <- mean(spec_las) 
          acc_las.p3000.n50  <- mean(acc_las)
          G_las.p3000.n50    <- mean(G_las)
          
          sens_alas.p3000.n50 <- mean(sens_alas)
          spec_alas.p3000.n50 <- mean(spec_alas) 
          acc_alas.p3000.n50  <- mean(acc_alas)
          G_alas.p3000.n50    <- mean(G_alas)
          
          sens_scad.p3000.n50 <- mean(sens_scad)
          spec_scad.p3000.n50 <- mean(spec_scad) 
          acc_scad.p3000.n50  <- mean(acc_scad)
          G_scad.p3000.n50    <- mean(G_scad)
          
          RMSE_beta.p3000.n50             <- data.frame(rmsebeta_las.p3000.n50,rmsebeta_alas.p3000.n50,rmsebeta_scad.p3000.n50)
          colnames(RMSE_beta.p3000.n50)   <- c("Lasso","Adaptive-Lasso","SCAD")
          
          RSQ.p3000.n50                   <- data.frame(Rsq_las.p3000.n50,Rsq_alas.p3000.n50,Rsq_scad.p3000.n50)
          colnames(RSQ.p3000.n50)         <- c("Lasso","Adaptive-Lasso","SCAD")
          
          MSEF.p3000.n50                  <- data.frame(msef_las.p3000.n50,msef_alas.p3000.n50,msef_scad.p3000.n50)
          colnames(MSEF.p3000.n50)        <- c("Lasso","Adaptive-Lasso","SCAD")
          
          metrics_all.p3000.n50           <- rbind(RMSE_beta.p3000.n50,RSQ.p3000.n50,MSEF.p3000.n50)
          rownames(metrics_all.p3000.n50) <- c("REMSE","R-square","MSEf") 
          
          stong_signals.p3000.n50 <- data.frame((strong_las),(strong_alas),(strong_scad))
          weak_signals.p3000.n50 <- data.frame((weak_las),(weak_alas),(weak_scad))
          
          RECALL.las.p3000.n50                     <- data.frame(sens_las.p3000.n50,spec_las.p3000.n50,acc_las.p3000.n50,G_las.p3000.n50)
          colnames(RECALL.las.p3000.n50)           <- c("Sensitivity","Specificity","Accuracy","G-score")
          RECALL.alas.p3000.n50                    <- data.frame(sens_alas.p3000.n50,spec_alas.p3000.n50,acc_alas.p3000.n50,G_alas.p3000.n50)
          colnames(RECALL.alas.p3000.n50)          <- c("Sensitivity","Specificity","Accuracy","G-score")
          RECALL.scad.p3000.n50                    <- data.frame(sens_scad.p3000.n50,spec_scad.p3000.n50,acc_scad.p3000.n50,G_scad.p3000.n50)
          colnames(RECALL.scad.p3000.n50)          <- c("Sensitivity","Specificity","Accuracy","G-score")
          RECALL_all.p3000.n50                     <- rbind(RECALL.las.p3000.n50,RECALL.alas.p3000.n50,RECALL.scad.p3000.n50)
          rownames(RECALL_all.p3000.n50)           <- c("Lasso","Adaptive-Lasso","SCAD")
        }
        if (s2==1 & s1==2){
          WR.lasso.fhat.p300.n100  <-rowMeans(Est.f.LassoWR) 
          WR.alasso.fhat.p300.n100 <-rowMeans(Est.f.aLassoWR)
          WR.scad.fhat.p300.n100   <-rowMeans(Est.f.SCADWR)
          
          rmsebeta_las.p300.n100   <- mean(rmsebeta_las)
          rmsebeta_alas.p300.n100 <- mean(rmsebeta_alas)
          rmsebeta_scad.p300.n100 <- mean(rmsebeta_scad)
          
          Rsq_las.p300.n100  <- mean(Rsq_las)
          Rsq_alas.p300.n100 <- mean(Rsq_alas)
          Rsq_scad.p300.n100 <- mean(Rsq_scad)
          
          msef_las.p300.n100  <- mean(msef_las)
          msef_alas.p300.n100 <- mean(msef_alas)
          msef_scad.p300.n100 <- mean(msef_scad)
          
          sens_las.p300.n100 <- mean(sens_las)
          spec_las.p300.n100 <- mean(spec_las) 
          acc_las.p300.n100  <- mean(acc_las)
          G_las.p300.n100    <- mean(G_las)
          
          sens_alas.p300.n100 <- mean(sens_alas)
          spec_alas.p300.n100 <- mean(spec_alas) 
          acc_alas.p300.n100  <- mean(acc_alas)
          G_alas.p300.n100    <- mean(G_alas)
          
          sens_scad.p300.n100 <- mean(sens_scad)
          spec_scad.p300.n100 <- mean(spec_scad) 
          acc_scad.p300.n100  <- mean(acc_scad)
          G_scad.p300.n100    <- mean(G_scad)
          
          RMSE_beta.p300.n100             <- data.frame(rmsebeta_las.p300.n100,rmsebeta_alas.p300.n100,rmsebeta_scad.p300.n100)
          colnames(RMSE_beta.p300.n100)   <- c("Lasso","Adaptive-Lasso","SCAD")
          
          RSQ.p300.n100                   <- data.frame(Rsq_las.p300.n100,Rsq_alas.p300.n100,Rsq_scad.p300.n100)
          colnames(RSQ.p300.n100)         <- c("Lasso","Adaptive-Lasso","SCAD")
          
          MSEF.p300.n100                  <- data.frame(msef_las.p300.n100,msef_alas.p300.n100,msef_scad.p300.n100)
          colnames(MSEF.p300.n100)        <- c("Lasso","Adaptive-Lasso","SCAD")
          
          metrics_all.p300.n100           <- rbind(RMSE_beta.p300.n100,RSQ.p300.n100,MSEF.p300.n100)
          rownames(metrics_all.p300.n100) <- c("REMSE","R-square","MSEf") 
          
          stong_signals.p300.n100 <- data.frame((strong_las),(strong_alas),(strong_scad))
          weak_signals.p300.n100 <- data.frame((weak_las),(weak_alas),(weak_scad))
          
          RECALL.las.p300.n100                     <- data.frame(sens_las.p300.n100,spec_las.p300.n100,acc_las.p300.n100,G_las.p300.n100)
          colnames(RECALL.las.p300.n100)           <- c("Sensitivity","Specificity","Accuracy","G-score")
          RECALL.alas.p300.n100                    <- data.frame(sens_alas.p300.n100,spec_alas.p300.n100,acc_alas.p300.n100,G_alas.p300.n100)
          colnames(RECALL.alas.p300.n100)          <- c("Sensitivity","Specificity","Accuracy","G-score")
          RECALL.scad.p300.n100                    <- data.frame(sens_scad.p300.n100,spec_scad.p300.n100,acc_scad.p300.n100,G_scad.p300.n100)
          colnames(RECALL.scad.p300.n100)          <- c("Sensitivity","Specificity","Accuracy","G-score")
          RECALL_all.p300.n100                     <- rbind(RECALL.las.p300.n100,RECALL.alas.p300.n100,RECALL.scad.p300.n100)
          rownames(RECALL_all.p300.n100)           <- c("Lasso","Adaptive-Lasso","SCAD")
          
        }
        if (s2==2 & s1==2){
          WR.lasso.fhat.p1000.n100  <-rowMeans(Est.f.LassoWR) 
          WR.alasso.fhat.p1000.n100 <-rowMeans(Est.f.aLassoWR)
          WR.scad.fhat.p1000.n100   <-rowMeans(Est.f.SCADWR)
          
          rmsebeta_las.p1000.n100   <- mean(rmsebeta_las)
          rmsebeta_alas.p1000.n100  <- mean(rmsebeta_alas)
          rmsebeta_scad.p1000.n100  <- mean(rmsebeta_scad)
          
          Rsq_las.p1000.n100  <- mean(Rsq_las)
          Rsq_alas.p1000.n100 <- mean(Rsq_alas)
          Rsq_scad.p1000.n100 <- mean(Rsq_scad)
          
          msef_las.p1000.n100  <- mean(msef_las)
          msef_alas.p1000.n100 <- mean(msef_alas)
          msef_scad.p1000.n100 <- mean(msef_scad)
          
          sens_las.p1000.n100 <- mean(sens_las)
          spec_las.p1000.n100 <- mean(spec_las) 
          acc_las.p1000.n100  <- mean(acc_las)
          G_las.p1000.n100    <- mean(G_las)
          
          sens_alas.p1000.n100 <- mean(sens_alas)
          spec_alas.p1000.n100 <- mean(spec_alas) 
          acc_alas.p1000.n100  <- mean(acc_alas)
          G_alas.p1000.n100    <- mean(G_alas)
          
          sens_scad.p1000.n100 <- mean(sens_scad)
          spec_scad.p1000.n100 <- mean(spec_scad) 
          acc_scad.p1000.n100  <- mean(acc_scad)
          G_scad.p1000.n100    <- mean(G_scad)
          
          RMSE_beta.p1000.n100             <- data.frame(rmsebeta_las.p1000.n100,rmsebeta_alas.p1000.n100,rmsebeta_scad.p1000.n100)
          colnames(RMSE_beta.p1000.n100)   <- c("Lasso","Adaptive-Lasso","SCAD")
          
          RSQ.p1000.n100                   <- data.frame(Rsq_las.p1000.n100,Rsq_alas.p1000.n100,Rsq_scad.p1000.n100)
          colnames(RSQ.p1000.n100)         <- c("Lasso","Adaptive-Lasso","SCAD")
          
          MSEF.p1000.n100                  <- data.frame(msef_las.p1000.n100,msef_alas.p1000.n100,msef_scad.p1000.n100)
          colnames(MSEF.p1000.n100)        <- c("Lasso","Adaptive-Lasso","SCAD")
          
          metrics_all.p1000.n100           <- rbind(RMSE_beta.p1000.n100,RSQ.p1000.n100,MSEF.p1000.n100)
          rownames(metrics_all.p1000.n100) <- c("REMSE","R-square","MSEf") 
          
          stong_signals.p1000.n100 <- data.frame((strong_las),(strong_alas),(strong_scad))
          weak_signals.p1000.n100 <- data.frame((weak_las),(weak_alas),(weak_scad))
          
          RECALL.las.p1000.n100                     <- data.frame(sens_las.p1000.n100,spec_las.p1000.n100,acc_las.p1000.n100,G_las.p1000.n100)
          colnames(RECALL.las.p1000.n100)           <- c("Sensitivity","Specificity","Accuracy","G-score")
          RECALL.alas.p1000.n100                    <- data.frame(sens_alas.p1000.n100,spec_alas.p1000.n100,acc_alas.p1000.n100,G_alas.p1000.n100)
          colnames(RECALL.alas.p1000.n100)          <- c("Sensitivity","Specificity","Accuracy","G-score")
          RECALL.scad.p1000.n100                    <- data.frame(sens_scad.p1000.n100,spec_scad.p1000.n100,acc_scad.p1000.n100,G_scad.p1000.n100)
          colnames(RECALL.scad.p1000.n100)          <- c("Sensitivity","Specificity","Accuracy","G-score")
          RECALL_all.p1000.n100                     <- rbind(RECALL.las.p1000.n100,RECALL.alas.p1000.n100,RECALL.scad.p1000.n100)
          rownames(RECALL_all.p1000.n100)           <- c("Lasso","Adaptive-Lasso","SCAD")
          
        }
        if (s2==3 & s1==2){
          WR.lasso.fhat.p3000.n100  <-rowMeans(Est.f.LassoWR) 
          WR.alasso.fhat.p3000.n100 <-rowMeans(Est.f.aLassoWR)
          WR.scad.fhat.p3000.n100   <-rowMeans(Est.f.SCADWR)
          
          rmsebeta_las.p3000.n100   <- mean(rmsebeta_las)
          rmsebeta_alas.p3000.n100  <- mean(rmsebeta_alas)
          rmsebeta_scad.p3000.n100  <- mean(rmsebeta_scad)
          
          Rsq_las.p3000.n100  <- mean(Rsq_las)
          Rsq_alas.p3000.n100 <- mean(Rsq_alas)
          Rsq_scad.p3000.n100 <- mean(Rsq_scad)
          
          msef_las.p3000.n100  <- mean(msef_las)
          msef_alas.p3000.n100 <- mean(msef_alas)
          msef_scad.p3000.n100 <- mean(msef_scad)
          
          sens_las.p3000.n100 <- mean(sens_las)
          spec_las.p3000.n100 <- mean(spec_las) 
          acc_las.p3000.n100  <- mean(acc_las)
          G_las.p3000.n100    <- mean(G_las)
          
          sens_alas.p3000.n100 <- mean(sens_alas)
          spec_alas.p3000.n100 <- mean(spec_alas) 
          acc_alas.p3000.n100  <- mean(acc_alas)
          G_alas.p3000.n100    <- mean(G_alas)
          
          sens_scad.p3000.n100 <- mean(sens_scad)
          spec_scad.p3000.n100 <- mean(spec_scad) 
          acc_scad.p3000.n100  <- mean(acc_scad)
          G_scad.p3000.n100    <- mean(G_scad)
          
          RMSE_beta.p3000.n100             <- data.frame(rmsebeta_las.p3000.n100,rmsebeta_alas.p3000.n100,rmsebeta_scad.p3000.n100)
          colnames(RMSE_beta.p3000.n100)   <- c("Lasso","Adaptive-Lasso","SCAD")
          
          RSQ.p3000.n100                   <- data.frame(Rsq_las.p3000.n100,Rsq_alas.p3000.n100,Rsq_scad.p3000.n100)
          colnames(RSQ.p3000.n100)         <- c("Lasso","Adaptive-Lasso","SCAD")
          
          MSEF.p3000.n100                  <- data.frame(msef_las.p3000.n100,msef_alas.p3000.n100,msef_scad.p3000.n100)
          colnames(MSEF.p3000.n100)        <- c("Lasso","Adaptive-Lasso","SCAD")
          
          metrics_all.p3000.n100           <- rbind(RMSE_beta.p3000.n100,RSQ.p3000.n100,MSEF.p3000.n100)
          rownames(metrics_all.p3000.n100) <- c("REMSE","R-square","MSEf") 
          
          stong_signals.p3000.n100 <- data.frame((strong_las),(strong_alas),(strong_scad))
          weak_signals.p3000.n100 <- data.frame((weak_las),(weak_alas),(weak_scad))
          
          RECALL.las.p3000.n100                     <- data.frame(sens_las.p3000.n100,spec_las.p3000.n100,acc_las.p3000.n100,G_las.p3000.n100)
          colnames(RECALL.las.p3000.n100)           <- c("Sensitivity","Specificity","Accuracy","G-score")
          RECALL.alas.p3000.n100                    <- data.frame(sens_alas.p3000.n100,spec_alas.p3000.n100,acc_alas.p3000.n100,G_alas.p3000.n100)
          colnames(RECALL.alas.p3000.n100)          <- c("Sensitivity","Specificity","Accuracy","G-score")
          RECALL.scad.p3000.n100                    <- data.frame(sens_scad.p3000.n100,spec_scad.p3000.n100,acc_scad.p3000.n100,G_scad.p3000.n100)
          colnames(RECALL.scad.p3000.n100)          <- c("Sensitivity","Specificity","Accuracy","G-score")
          RECALL_all.p3000.n100                     <- rbind(RECALL.las.p3000.n100,RECALL.alas.p3000.n100,RECALL.scad.p3000.n100)
          rownames(RECALL_all.p3000.n100)           <- c("Lasso","Adaptive-Lasso","SCAD")
        }
        if (s2==1 & s1==3){
          WR.lasso.fhat.p300.n200  <-rowMeans(Est.f.LassoWR) 
          WR.alasso.fhat.p300.n200 <-rowMeans(Est.f.aLassoWR)
          WR.scad.fhat.p300.n200   <-rowMeans(Est.f.SCADWR)
          
          rmsebeta_las.p300.n200   <- mean(rmsebeta_las)
          rmsebeta_alas.p300.n200 <- mean(rmsebeta_alas)
          rmsebeta_scad.p300.n200 <- mean(rmsebeta_scad)
          
          Rsq_las.p300.n200  <- mean(Rsq_las)
          Rsq_alas.p300.n200 <- mean(Rsq_alas)
          Rsq_scad.p300.n200 <- mean(Rsq_scad)
          
          msef_las.p300.n200  <- mean(msef_las)
          msef_alas.p300.n200 <- mean(msef_alas)
          msef_scad.p300.n200 <- mean(msef_scad)
          
          sens_las.p300.n200 <- mean(sens_las)
          spec_las.p300.n200 <- mean(spec_las) 
          acc_las.p300.n200  <- mean(acc_las)
          G_las.p300.n200    <- mean(G_las)
          
          sens_alas.p300.n200 <- mean(sens_alas)
          spec_alas.p300.n200 <- mean(spec_alas) 
          acc_alas.p300.n200  <- mean(acc_alas)
          G_alas.p300.n200    <- mean(G_alas)
          
          sens_scad.p300.n200 <- mean(sens_scad)
          spec_scad.p300.n200 <- mean(spec_scad) 
          acc_scad.p300.n200  <- mean(acc_scad)
          G_scad.p300.n200    <- mean(G_scad)
          
          RMSE_beta.p300.n200             <- data.frame(rmsebeta_las.p300.n200,rmsebeta_alas.p300.n200,rmsebeta_scad.p300.n200)
          colnames(RMSE_beta.p300.n200)   <- c("Lasso","Adaptive-Lasso","SCAD")
          
          RSQ.p300.n200                   <- data.frame(Rsq_las.p300.n200,Rsq_alas.p300.n200,Rsq_scad.p300.n200)
          colnames(RSQ.p300.n200)         <- c("Lasso","Adaptive-Lasso","SCAD")
          
          MSEF.p300.n200                  <- data.frame(msef_las.p300.n200,msef_alas.p300.n200,msef_scad.p300.n200)
          colnames(MSEF.p300.n200)        <- c("Lasso","Adaptive-Lasso","SCAD")
          
          metrics_all.p300.n200           <- rbind(RMSE_beta.p300.n200,RSQ.p300.n200,MSEF.p300.n200)
          rownames(metrics_all.p300.n200) <- c("REMSE","R-square","MSEf") 
          
          stong_signals.p300.n200 <- data.frame((strong_las),(strong_alas),(strong_scad))
          weak_signals.p300.n200  <- data.frame((weak_las),(weak_alas),(weak_scad))
          
          RECALL.las.p300.n200                     <- data.frame(sens_las.p300.n200,spec_las.p300.n200,acc_las.p300.n200,G_las.p300.n200)
          colnames(RECALL.las.p300.n200)           <- c("Sensitivity","Specificity","Accuracy","G-score")
          RECALL.alas.p300.n200                    <- data.frame(sens_alas.p300.n200,spec_alas.p300.n200,acc_alas.p300.n200,G_alas.p300.n200)
          colnames(RECALL.alas.p300.n200)          <- c("Sensitivity","Specificity","Accuracy","G-score")
          RECALL.scad.p300.n200                    <- data.frame(sens_scad.p300.n200,spec_scad.p300.n200,acc_scad.p300.n200,G_scad.p300.n200)
          colnames(RECALL.scad.p300.n200)          <- c("Sensitivity","Specificity","Accuracy","G-score")
          RECALL_all.p300.n200                     <- rbind(RECALL.las.p300.n200,RECALL.alas.p300.n200,RECALL.scad.p300.n200)
          rownames(RECALL_all.p300.n200)           <- c("Lasso","Adaptive-Lasso","SCAD")
        
        }
        if (s2==2 & s1==3){
          WR.lasso.fhat.p1000.n200  <-rowMeans(Est.f.LassoWR) 
          WR.alasso.fhat.p1000.n200 <-rowMeans(Est.f.aLassoWR)
          WR.scad.fhat.p1000.n200   <-rowMeans(Est.f.SCADWR)
          
          rmsebeta_las.p1000.n200   <- mean(rmsebeta_las)
          rmsebeta_alas.p1000.n200 <- mean(rmsebeta_alas)
          rmsebeta_scad.p1000.n200 <- mean(rmsebeta_scad)
          
          Rsq_las.p1000.n200  <- mean(Rsq_las)
          Rsq_alas.p1000.n200 <- mean(Rsq_alas)
          Rsq_scad.p1000.n200 <- mean(Rsq_scad)
          
          msef_las.p1000.n200  <- mean(msef_las)
          msef_alas.p1000.n200 <- mean(msef_alas)
          msef_scad.p1000.n200 <- mean(msef_scad)
          
          sens_las.p1000.n200 <- mean(sens_las)
          spec_las.p1000.n200 <- mean(spec_las) 
          acc_las.p1000.n200  <- mean(acc_las)
          G_las.p1000.n200    <- mean(G_las)
          
          sens_alas.p1000.n200 <- mean(sens_alas)
          spec_alas.p1000.n200 <- mean(spec_alas) 
          acc_alas.p1000.n200  <- mean(acc_alas)
          G_alas.p1000.n200    <- mean(G_alas)
          
          sens_scad.p1000.n200 <- mean(sens_scad)
          spec_scad.p1000.n200 <- mean(spec_scad) 
          acc_scad.p1000.n200  <- mean(acc_scad)
          G_scad.p1000.n200    <- mean(G_scad)
          
          stong_signals.p1000.n200 <- data.frame((strong_las),(strong_alas),(strong_scad))
          weak_signals.p1000.n200  <- data.frame((weak_las),(weak_alas),(weak_scad))
          
          RMSE_beta.p1000.n200             <- data.frame(rmsebeta_las.p1000.n200,rmsebeta_alas.p1000.n200,rmsebeta_scad.p1000.n200)
          colnames(RMSE_beta.p1000.n200)   <- c("Lasso","Adaptive-Lasso","SCAD")
          
          RSQ.p1000.n200                   <- data.frame(Rsq_las.p1000.n200,Rsq_alas.p1000.n200,Rsq_scad.p1000.n200)
          colnames(RSQ.p1000.n200)         <- c("Lasso","Adaptive-Lasso","SCAD")
          
          MSEF.p1000.n200                  <- data.frame(msef_las.p1000.n200,msef_alas.p1000.n200,msef_scad.p1000.n200)
          colnames(MSEF.p1000.n200)        <- c("Lasso","Adaptive-Lasso","SCAD")
          
          metrics_all.p1000.n200           <- rbind(RMSE_beta.p1000.n200,RSQ.p1000.n200,MSEF.p1000.n200)
          rownames(metrics_all.p1000.n200) <- c("REMSE","R-square","MSEf") 
          
          RECALL.las.p1000.n200                     <- data.frame(sens_las.p1000.n200,spec_las.p1000.n200,acc_las.p1000.n200,G_las.p1000.n200)
          colnames(RECALL.las.p1000.n200)           <- c("Sensitivity","Specificity","Accuracy","G-score")
          RECALL.alas.p1000.n200                    <- data.frame(sens_alas.p1000.n200,spec_alas.p1000.n200,acc_alas.p1000.n200,G_alas.p1000.n200)
          colnames(RECALL.alas.p1000.n200)          <- c("Sensitivity","Specificity","Accuracy","G-score")
          RECALL.scad.p1000.n200                    <- data.frame(sens_scad.p1000.n200,spec_scad.p1000.n200,acc_scad.p1000.n200,G_scad.p1000.n200)
          colnames(RECALL.scad.p1000.n200)          <- c("Sensitivity","Specificity","Accuracy","G-score")
          RECALL_all.p1000.n200                     <- rbind(RECALL.las.p1000.n200,RECALL.alas.p1000.n200,RECALL.scad.p1000.n200)
          rownames(RECALL_all.p1000.n200)           <- c("Lasso","Adaptive-Lasso","SCAD")
          
        }
        if (s2==3 & s1==3){
          WR.lasso.fhat.p3000.n200  <-rowMeans(Est.f.LassoWR) 
          WR.alasso.fhat.p3000.n200 <-rowMeans(Est.f.aLassoWR)
          WR.scad.fhat.p3000.n200   <-rowMeans(Est.f.SCADWR)
          
          rmsebeta_las.p3000.n200   <- mean(rmsebeta_las)
          rmsebeta_alas.p3000.n200 <- mean(rmsebeta_alas)
          rmsebeta_scad.p3000.n200 <- mean(rmsebeta_scad)
          
          Rsq_las.p3000.n200  <- mean(Rsq_las)
          Rsq_alas.p3000.n200 <- mean(Rsq_alas)
          Rsq_scad.p3000.n200 <- mean(Rsq_scad)
          
          msef_las.p3000.n200  <- mean(msef_las)
          msef_alas.p3000.n200 <- mean(msef_alas)
          msef_scad.p3000.n200 <- mean(msef_scad)
          
          sens_las.p3000.n200 <- mean(sens_las)
          spec_las.p3000.n200 <- mean(spec_las) 
          acc_las.p3000.n200  <- mean(acc_las)
          G_las.p3000.n200    <- mean(G_las)
          
          sens_alas.p3000.n200 <- mean(sens_alas)
          spec_alas.p3000.n200 <- mean(spec_alas) 
          acc_alas.p3000.n200  <- mean(acc_alas)
          G_alas.p3000.n200    <- mean(G_alas)
          
          sens_scad.p3000.n200 <- mean(sens_scad)
          spec_scad.p3000.n200 <- mean(spec_scad) 
          acc_scad.p3000.n200  <- mean(acc_scad)
          G_scad.p3000.n200    <- mean(G_scad)
          
          RMSE_beta.p3000.n200             <- data.frame(rmsebeta_las.p3000.n200,rmsebeta_alas.p3000.n200,rmsebeta_scad.p3000.n200)
          colnames(RMSE_beta.p3000.n200)   <- c("Lasso","Adaptive-Lasso","SCAD")
          
          RSQ.p3000.n200                   <- data.frame(Rsq_las.p3000.n200,Rsq_alas.p3000.n200,Rsq_scad.p3000.n200)
          colnames(RSQ.p3000.n200)         <- c("Lasso","Adaptive-Lasso","SCAD")
          
          MSEF.p3000.n200                  <- data.frame(msef_las.p3000.n200,msef_alas.p3000.n200,msef_scad.p3000.n200)
          colnames(MSEF.p3000.n200)        <- c("Lasso","Adaptive-Lasso","SCAD")
          
          metrics_all.p3000.n200           <- rbind(RMSE_beta.p3000.n200,RSQ.p3000.n200,MSEF.p3000.n200)
          rownames(metrics_all.p3000.n200) <- c("REMSE","R-square","MSEf") 
          
          stong_signals.p3000.n200 <- data.frame((strong_las),(strong_alas),(strong_scad))
          weak_signals.p3000.n200  <- data.frame((weak_las),(weak_alas),(weak_scad))
          
          RECALL.las.p3000.n200                     <- data.frame(sens_las.p3000.n200,spec_las.p3000.n200,acc_las.p3000.n200,G_las.p3000.n200)
          colnames(RECALL.las.p3000.n200)           <- c("Sensitivity","Specificity","Accuracy","G-score")
          RECALL.alas.p3000.n200                    <- data.frame(sens_alas.p3000.n200,spec_alas.p3000.n200,acc_alas.p3000.n200,G_alas.p3000.n200)
          colnames(RECALL.alas.p3000.n200)          <- c("Sensitivity","Specificity","Accuracy","G-score")
          RECALL.scad.p3000.n200                    <- data.frame(sens_scad.p3000.n200,spec_scad.p3000.n200,acc_scad.p3000.n200,G_scad.p3000.n200)
          colnames(RECALL.scad.p3000.n200)          <- c("Sensitivity","Specificity","Accuracy","G-score")
          RECALL_all.p3000.n200                     <- rbind(RECALL.las.p3000.n200,RECALL.alas.p3000.n200,RECALL.scad.p3000.n200)
          rownames(RECALL_all.p3000.n200)           <- c("Lasso","Adaptive-Lasso","SCAD")
          
        }
    }
  }
df_list.p300  <- rbind(metrics_all.p300.n50,metrics_all.p300.n100,metrics_all.p300.n200)
df_list.p1000 <- rbind(metrics_all.p1000.n50,metrics_all.p1000.n100,metrics_all.p1000.n200)
df_list.p3000 <- rbind(metrics_all.p3000.n50,metrics_all.p3000.n100,metrics_all.p3000.n200)


df_list.p300
df_list.p1000 
df_list.p3000

strong.signals.p300 <-cbind(stong_signals.p300.n50,stong_signals.p300.n100,stong_signals.p300.n200) 
weak.signals.p300   <-cbind(weak_signals.p300.n50,weak_signals.p300.n100,weak_signals.p300.n200)

strong.signals.p1000 <-cbind(stong_signals.p1000.n50,stong_signals.p1000.n100,stong_signals.p1000.n200) 
weak.signals.p1000   <-cbind(weak_signals.p1000.n50,weak_signals.p1000.n100,weak_signals.p1000.n200)

strong.signals.p3000 <-cbind(stong_signals.p300.n50,stong_signals.p300.n100,stong_signals.p300.n200) 
weak.signals.p3000   <-cbind(weak_signals.p300.n50,weak_signals.p300.n100,weak_signals.p300.n200)

par(mfrow=c(3,1))
boxplot(strong.signals.p300,main="Strong signals for pn=300")
boxplot(strong.signals.p1000,main="Strong signals for pn=1000")
boxplot(strong.signals.p3000,main="Strong signals for pn=3000")

par(mfrow=c(3,1))
boxplot(weak.signals.p300,main="Weak signals for pn=300")
boxplot(weak.signals.p1000,main="Weak signals for pn=1000")
boxplot(weak.signals.p3000,main="Weak signals for pn=3000")

recall.p300 <- data.frame(sens_las.p300.n50,sens_alas.p300.n50,sens_scad.p300.n50,acc_las.p300.n50,acc_alas.p300.n50,acc_scad.p300.n50)
recall.p1000 <- data.frame(sens_las.p1000.n50,sens_alas.p1000.n50,sens_scad.p1000.n50,acc_las.p1000.n50,acc_alas.p1000.n50,acc_scad.p1000.n50)
recall.p3000 <- data.frame(sens_las.p3000.n50,sens_alas.p3000.n50,sens_scad.p3000.n50,acc_las.p3000.n50,acc_alas.p3000.n50,acc_scad.p3000.n50)

recall.p300
recall.p1000
recall.p3000
###################################################################################################
plot(f,type="l",ylim=c(min(f),max(f)),col=1,main="n=200, pn=3000, Scenario 1, rho=0.80",xlab="t",ylab="f(t)")
par(new=TRUE)
plot(WR.lasso.fhat.p3000.n200,type="l",ylim=c(min(f),max(f)),col=2,xlab="t",ylab="f(t)")
par(new=TRUE)
plot(WR.alasso.fhat.p3000.n200+0.1,type="l",ylim=c(min(f),max(f)),col=3,xlab="t",ylab="f(t)")
par(new=TRUE)
plot(WR.scad.fhat.p3000.n200,type="l",ylim=c(min(f),max(f)),col=4,xlab="t",ylab="f(t)")
legend("bottomleft",legend=c("Real f", "Lasso-PSE","aLasso-PSE","SCAD-PSE"),
       col=c(1, 2,3,4),lty=c(1,1,1,1))
grid()

write.csv(df_list.p300,"p300sc1rho1.csv")
write.csv(df_list.p1000,"p800sc1rho1.csv")
write.csv(df_list.p3000,"p1500sc1rho1.csv")

write.csv(recall.p300,"RECALL.p300sc1rho1.csv")
write.csv(recall.p1000,"RECALL.p800sc1rho1.csv")
write.csv(recall.p3000,"RECALL.p1500sc1rho1.csv")


############################################################################################################
#===========================================================================================================
############################################################################################################
#SECOND SIMULATION for SCENARIO 2
#OUTLINE-----------------------------------------------
#Function 1: Data generation for simulation Sudy
#Function 2: Estimators based on Penalty functions
#Function 3: Calculation of evaluation metrics
#------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------
ni <- c(50,100,120)
pni <- c(300,350,400)
rhoi <- c(0.80,0.95)              #Number of simulation
for (s1 in 1:3){                  #Sample size index
  for (s2 in 1:3){                #pn index
    if (s1==1){n=50} 
    if (s1==2){n=100}  
    if (s1==3){n=200}
    if (s2==1){pn=300} 
    if (s2==2){pn=600} 
    if (s2==3){pn=1200}
    rho <- 0.80
    sc  <- 2
    #ZERO MATRICES----------------------------------------------------------
    rmsebeta_las    <- 0
    rmsebeta_alas   <- 0
    rmsebeta_scad   <- 0
    
    Rsq_las         <- 0
    Rsq_alas        <- 0
    Rsq_scad        <- 0
    
    msef_las         <- 0
    msef_alas        <- 0
    msef_scad        <- 0
    
    sens_las <- 0
    spec_las <- 0
    acc_las  <- 0
    G_las    <- 0
    sens_alas <- 0
    spec_alas <- 0
    acc_alas  <- 0
    G_alas    <- 0
    sens_scad <- 0
    spec_scad <- 0
    acc_scad  <- 0
    G_scad    <- 0
    
    strong_las  <- 0
    weak_las    <- 0
    strong_alas <- 0
    weak_alas   <- 0 
    strong_scad <- 0 
    weak_scad   <- 0
    
    Est.f.LassoWR  <- matrix(0,n,sim)
    Est.f.aLassoWR <- matrix(0,n,sim)
    Est.f.SCADWR   <- matrix(0,n,sim)
    
    Est.f.Lasso  <- matrix(0,n,sim)
    Est.f.aLasso <- matrix(0,n,sim)
    Est.f.SCAD   <- matrix(0,n,sim)
    
    Est.y.LassoWR  <- matrix(0,n,sim)
    Est.y.aLassoWR <- matrix(0,n,sim)
    Est.y.SCADWR   <- matrix(0,n,sim)
    
    Est.y.Lasso  <- matrix(0,n,sim)
    Est.y.aLasso <- matrix(0,n,sim)
    Est.y.SCAD   <- matrix(0,n,sim)
    #------------------------------------------------------------------------------
    no.seed <- 19345
    aa <- 1
    for (s in 1:sim){
      #Set seed number needs to be changed in every repeat----------------------------
      set.seed(no.seed+aa)
      mydata <- simdata(n,pn,sc,rho)
      set.seed(no.seed+s)
      #TEST AREA (To be Deleted)------------------------------------------------------
      x     <- mydata$X
      y     <- mydata$y
      t     <- mydata$t
      f     <- mydata$f 
      nol   <- 25
      aiclas  <- 0
      aicalas <- 0
      aicscad <- 0
      
      las     <- lasso(x,y)
      alas    <- adalasso(x,y)
      scadm   <- scad(x,y)
      
      lams    <- seq(0.1,1.5,length.out=nol)
      for (j in 1:nol){
        aiclas[j]  <- aicfunc(y,las,lams[j],f)
        aicalas[j]  <- aicfunc(y,alas,lams[j],f)
        aicscad[j]  <- aicfunc(y,scadm,lams[j],f)
      }
      #PLOTS for smoothing parameter selection (AICc)---------------------------------
      dflas  <- data.frame(lams,aiclas)
      dfalas <- data.frame(lams,aicalas)
      dfscad <- data.frame(lams,aicscad)
      #ggplot()+geom_line(data=dflas,aes(x=lams,aiclas,col="Lasso_WR"),lwd=1)+geom_line(data=dfalas,aes(x=lams,aicalas,col="Adapt.-Lasso_WR"),lwd=1)+geom_line(data=dfscad,aes(x=lams,aicscad,col="SCAD_WR"),lwd=1)+xlab("lambda_spline")+ylab("Improved_AIC score")+ggtitle("(e) n=100, pn=3000, rho=0.95")+theme(legend.title = element_blank(),legend.position = "bottom")
      #-------------------------------------------------------------------------------
      for (j2 in 1:nol){
        if (aiclas[j2]==min(aiclas)){
          lam_las <-  lams[j2]
        }
        if (aicalas[j2]==min(aicalas)){
          lam_alas <-  lams[j2]
        }
        if (aicscad[j2]==min(aicscad)){
          lam_scad <-  lams[j2]
        }
      }
      S_las  <- Smatrix(t,2,lam_las)
      S_alas <- Smatrix(t,2,lam_alas)
      S_scad <- Smatrix(t,2,lam_scad)
      
      
      ytil_las  <- (diag(n)-S_las)%*%y
      ytil_alas <- (diag(n)-S_alas)%*%y
      ytil_scad <- (diag(n)-S_scad)%*%y
      
      
      xtil_las  <- (diag(n)-S_las)%*%x
      xtil_alas <- (diag(n)-S_alas)%*%x
      xtil_scad <- (diag(n)-S_scad)%*%x
      
      newX_las  <- (diag(n)-S_las)%*%las$X_S1
      newX_alas <- (diag(n)-S_alas)%*%alas$X_S1
      newX_scad <- (diag(n)-S_scad)%*%scadm$X_S1
      #Estimation based on Penalty functions
      las_beta      <- solve(t(newX_las)%*%newX_las+0.01*diag(las$S1))%*%t(newX_las)%*%ytil_las
      fhat_las      <- S_las%*%(y-las$X_S1%*%las_beta)
      alas_beta     <- solve(t(newX_alas)%*%newX_alas+0.01*diag(alas$S1))%*%t(newX_alas)%*%ytil_alas
      fhat_alas     <- S_alas%*%(y-alas$X_S1%*%alas_beta)
      scad_beta     <- solve(t(newX_scad)%*%newX_scad+0.01*diag(scadm$S1))%*%t(newX_scad)%*%ytil_scad
      fhat_scad     <- S_scad%*%(y-scadm$X_S1%*%scad_beta)
      
      yhat_las  <- las$X_S1%*%las_beta+fhat_las
      yhat_alas <- alas$X_S1%*%alas_beta+fhat_alas
      yhat_scad <- scadm$X_S1%*%scad_beta+fhat_scad
      
      WRobj_las  <- WR(las,S_las,x,y,0.01)
      WRobj_alas <- WR(alas,S_alas,x,y,0.01)
      WRobj_scad <- WR(scadm,S_scad,x,y,0.01)
      
      WR_fhat_las  <- S_las%*%(y-las$X_S1%*%WRobj_las$betaPSE)
      WR_fhat_alas <- S_alas%*%(y-alas$X_S1%*%WRobj_alas$betaPSE)
      WR_fhat_scad <- S_scad%*%(y-scadm$X_S1%*%WRobj_scad$betaPSE)
      
      WR_yhat_las  <- las$X_S1%*%WRobj_las$betaPSE+fhat_las
      WR_yhat_alas <- alas$X_S1%*%WRobj_alas$betaPSE+fhat_alas
      WR_yhat_scad <- scadm$X_S1%*%WRobj_scad$betaPSE
      #---------------------------------------------------------------------
      las_index       <- las$S1
      real_beta_las   <- mydata$beta[las_index]
      alas_index      <- alas$S1
      real_beta_alas  <- mydata$beta[alas_index]
      scad_index      <- scadm$S1
      real_beta_scad  <- mydata$beta[scad_index]
      
      rmsebeta_las[s]  <- rmsebeta(real_beta_las,WRobj_las$betaPSE)
      rmsebeta_alas[s] <- rmsebeta(real_beta_alas,WRobj_alas$betaPSE)
      rmsebeta_scad[s] <- rmsebeta(real_beta_scad,WRobj_scad$betaPSE)
      
      Rsq_las[s]  <- rsq(mydata$y,WR_yhat_las,length(las_index))
      Rsq_alas[s] <- rsq(mydata$y,WR_yhat_alas,length(alas_index))
      Rsq_scad[s] <- rsq(mydata$y,WR_yhat_scad,length(scad_index))
      
      msef_las[s]  <- msef(mydata$f,WR_fhat_las)
      msef_alas[s] <- msef(mydata$f,WR_fhat_alas)
      msef_scad[s] <- msef(mydata$f,WR_fhat_scad)
      
      recall_las  <- recall(real_beta_las,WRobj_las$betaPSE)
      sens_las[s] <- recall_las$sensitivity
      spec_las[s] <- recall_las$specifcity
      acc_las[s]  <- recall_las$accuracy
      G_las[s]    <- recall_las$Gscore
      
      recall_alas  <- recall(real_beta_alas,WRobj_alas$betaPSE)
      sens_alas[s] <- recall_alas$sensitivity
      spec_alas[s] <- recall_alas$specifcity
      acc_alas[s]  <- recall_alas$accuracy
      G_alas[s]    <- recall_alas$Gscore
      
      recall_scad  <- recall(real_beta_scad,WRobj_scad$betaPSE)
      sens_scad[s] <- recall_scad$sensitivity
      spec_scad[s] <- recall_scad$specifcity
      acc_scad[s]  <- recall_scad$accuracy
      G_scad[s]    <- recall_scad$Gscore
      
      strong_las[s]<- length(WRobj_las$betaPSE)
      weak_las[s]  <- WRobj_las$S2
      
      strong_alas[s]<- length(WRobj_alas$betaPSE)
      weak_alas[s]  <- WRobj_alas$S2
      
      strong_scad[s]<- length(WRobj_scad$betaPSE)
      weak_scad[s]  <- WRobj_scad$S2
      #---------------------------------------------------------------------
      Est.f.LassoWR[,s]  <- WR_fhat_las
      Est.f.aLassoWR[,s] <- WR_fhat_alas
      Est.f.SCADWR[,s]   <- WR_fhat_scad
      
      Est.f.Lasso[,s]  <- fhat_las
      Est.f.aLasso[,s] <- fhat_alas
      Est.f.SCAD[,s]   <- fhat_scad
      
      Est.y.Lasso[,s]  <- yhat_las
      Est.y.aLasso[,s] <- yhat_alas
      Est.y.SCAD[,s]   <- yhat_scad
      
      Est.y.LassoWR[,s]  <- WR_yhat_las
      Est.y.aLassoWR[,s] <- WR_yhat_alas
      Est.y.SCADWR[,s]   <- WR_yhat_scad
      
      #plot(f,type="l",ylim=c(min(f),max(f)))
      #par(new=TRUE)
      #plot(WR_fhat_las,type="l",ylim=c(min(f),max(f)),col="red")
      #par(new=TRUE)
      #plot(WR_fhat_alas,type="l",ylim=c(min(f),max(f)),col="blue")
      #par(new=TRUE)
      #plot(WR_fhat_scad,type="l",ylim=c(min(f),max(f)),col="green")
      aa <- aa+1
      message(s, "th Simulation ends for n=",n, " pn=",pn, " and rho=", rho, " in Scenario ",sc)
    }
    if (s2==1 & s1==1){
      WR.lasso.fhat.p300.n50  <-rowMeans(Est.f.LassoWR) 
      WR.alasso.fhat.p300.n50 <-rowMeans(Est.f.aLassoWR)
      WR.scad.fhat.p300.n50   <-rowMeans(Est.f.SCADWR)
      
      rmsebeta_las.p300.n50  <- mean(rmsebeta_las)
      rmsebeta_alas.p300.n50 <- mean(rmsebeta_alas)
      rmsebeta_scad.p300.n50 <- mean(rmsebeta_scad)
      
      Rsq_las.p300.n50  <- mean(Rsq_las)
      Rsq_alas.p300.n50 <- mean(Rsq_alas)
      Rsq_scad.p300.n50 <- mean(Rsq_scad)
      
      msef_las.p300.n50  <- mean(msef_las)
      msef_alas.p300.n50 <- mean(msef_alas)
      msef_scad.p300.n50 <- mean(msef_scad)
      
      sens_las.p300.n50 <- mean(sens_las)
      spec_las.p300.n50 <- mean(spec_las) 
      acc_las.p300.n50  <- mean(acc_las)
      G_las.p300.n50    <- mean(G_las)
      
      sens_alas.p300.n50 <- mean(sens_alas)
      spec_alas.p300.n50 <- mean(spec_alas) 
      acc_alas.p300.n50  <- mean(acc_alas)
      G_alas.p300.n50    <- mean(G_alas)
      
      sens_scad.p300.n50 <- mean(sens_scad)
      spec_scad.p300.n50 <- mean(spec_scad) 
      acc_scad.p300.n50  <- mean(acc_scad)
      G_scad.p300.n50    <- mean(G_scad)
      
      stong_signals.p300.n50 <- data.frame((strong_las),(strong_alas),(strong_scad))
      weak_signals.p300.n50 <- data.frame((weak_las),(weak_alas),(weak_scad))
      
      
      RMSE_beta.p300.n50             <- data.frame(rmsebeta_las.p300.n50,rmsebeta_alas.p300.n50,rmsebeta_scad.p300.n50)
      colnames(RMSE_beta.p300.n50)   <- c("Lasso","Adaptive-Lasso","SCAD")
      
      RSQ.p300.n50                   <- data.frame(Rsq_las.p300.n50,Rsq_alas.p300.n50,Rsq_scad.p300.n50)
      colnames(RSQ.p300.n50)         <- c("Lasso","Adaptive-Lasso","SCAD")
      
      MSEF.p300.n50                  <- data.frame(msef_las.p300.n50,msef_alas.p300.n50,msef_scad.p300.n50)
      colnames(MSEF.p300.n50)        <- c("Lasso","Adaptive-Lasso","SCAD")
      
      metrics_all.p300.n50           <- rbind(RMSE_beta.p300.n50,RSQ.p300.n50,MSEF.p300.n50)
      rownames(metrics_all.p300.n50) <- c("REMSE","R-square","MSEf") 
      
      RECALL.las.p300.n50                     <- data.frame(sens_las.p300.n50,spec_las.p300.n50,acc_las.p300.n50,G_las.p300.n50)
      colnames(RECALL.las.p300.n50)           <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.alas.p300.n50                    <- data.frame(sens_alas.p300.n50,spec_alas.p300.n50,acc_alas.p300.n50,G_alas.p300.n50)
      colnames(RECALL.alas.p300.n50)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.scad.p300.n50                    <- data.frame(sens_scad.p300.n50,spec_scad.p300.n50,acc_scad.p300.n50,G_scad.p300.n50)
      colnames(RECALL.scad.p300.n50)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL_all.p300.n50                     <- rbind(RECALL.las.p300.n50,RECALL.alas.p300.n50,RECALL.scad.p300.n50)
      rownames(RECALL_all.p300.n50)           <- c("Lasso","Adaptive-Lasso","SCAD")
      #-PLOTS of fs---------------------------------------------------------
      plot(f,type="l",ylim=c(min(f),max(f)))
      par(new=TRUE)
      plot(WR.lasso.fhat.p300.n50,type="l",ylim=c(min(f),max(f)),col="red")
      par(new=TRUE)
      plot(WR.alasso.fhat.p300.n50,type="l",ylim=c(min(f),max(f)),col="blue")
      par(new=TRUE)
      plot(WR.scad.fhat.p300.n50,type="l",ylim=c(min(f),max(f)),col="green")
    }
    if (s2==2 & s1==1){
      WR.lasso.fhat.p1000.n50  <-rowMeans(Est.f.LassoWR) 
      WR.alasso.fhat.p1000.n50 <-rowMeans(Est.f.aLassoWR)
      WR.scad.fhat.p1000.n50   <-rowMeans(Est.f.SCADWR)
      
      rmsebeta_las.p1000.n50  <- mean(rmsebeta_las)
      rmsebeta_alas.p1000.n50 <- mean(rmsebeta_alas)
      rmsebeta_scad.p1000.n50 <- mean(rmsebeta_scad)
      
      Rsq_las.p1000.n50  <- mean(Rsq_las)
      Rsq_alas.p1000.n50 <- mean(Rsq_alas)
      Rsq_scad.p1000.n50 <- mean(Rsq_scad)
      
      msef_las.p1000.n50  <- mean(msef_las)
      msef_alas.p1000.n50 <- mean(msef_alas)
      msef_scad.p1000.n50 <- mean(msef_scad)
      
      sens_las.p1000.n50 <- mean(sens_las)
      spec_las.p1000.n50 <- mean(spec_las) 
      acc_las.p1000.n50  <- mean(acc_las)
      G_las.p1000.n50    <- mean(G_las)
      
      sens_alas.p1000.n50 <- mean(sens_alas)
      spec_alas.p1000.n50 <- mean(spec_alas) 
      acc_alas.p1000.n50  <- mean(acc_alas)
      G_alas.p1000.n50    <- mean(G_alas)
      
      sens_scad.p1000.n50 <- mean(sens_scad)
      spec_scad.p1000.n50 <- mean(spec_scad) 
      acc_scad.p1000.n50  <- mean(acc_scad)
      G_scad.p1000.n50    <- mean(G_scad)
      
      RMSE_beta.p1000.n50             <- data.frame(rmsebeta_las.p1000.n50,rmsebeta_alas.p1000.n50,rmsebeta_scad.p1000.n50)
      colnames(RMSE_beta.p1000.n50)   <- c("Lasso","Adaptive-Lasso","SCAD")
      
      RSQ.p1000.n50                   <- data.frame(Rsq_las.p1000.n50,Rsq_alas.p1000.n50,Rsq_scad.p1000.n50)
      colnames(RSQ.p1000.n50)         <- c("Lasso","Adaptive-Lasso","SCAD")
      
      MSEF.p1000.n50                  <- data.frame(msef_las.p1000.n50,msef_alas.p1000.n50,msef_scad.p1000.n50)
      colnames(MSEF.p1000.n50)        <- c("Lasso","Adaptive-Lasso","SCAD")
      
      metrics_all.p1000.n50           <- rbind(RMSE_beta.p1000.n50,RSQ.p1000.n50,MSEF.p1000.n50)
      rownames(metrics_all.p1000.n50) <- c("REMSE","R-square","MSEf") 
      
      stong_signals.p1000.n50 <- data.frame((strong_las),(strong_alas),(strong_scad))
      weak_signals.p1000.n50 <- data.frame((weak_las),(weak_alas),(weak_scad))
      
      RECALL.las.p1000.n50                     <- data.frame(sens_las.p1000.n50,spec_las.p1000.n50,acc_las.p1000.n50,G_las.p1000.n50)
      colnames(RECALL.las.p1000.n50)           <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.alas.p1000.n50                    <- data.frame(sens_alas.p1000.n50,spec_alas.p1000.n50,acc_alas.p1000.n50,G_alas.p1000.n50)
      colnames(RECALL.alas.p1000.n50)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.scad.p1000.n50                    <- data.frame(sens_scad.p1000.n50,spec_scad.p1000.n50,acc_scad.p1000.n50,G_scad.p1000.n50)
      colnames(RECALL.scad.p1000.n50)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL_all.p1000.n50                     <- rbind(RECALL.las.p1000.n50,RECALL.alas.p1000.n50,RECALL.scad.p1000.n50)
      rownames(RECALL_all.p1000.n50)           <- c("Lasso","Adaptive-Lasso","SCAD")
    }
    if (s2==3 & s1==1){
      WR.lasso.fhat.p3000.n50  <-rowMeans(Est.f.LassoWR) 
      WR.alasso.fhat.p3000.n50 <-rowMeans(Est.f.aLassoWR)
      WR.scad.fhat.p3000.n50   <-rowMeans(Est.f.SCADWR)
      
      rmsebeta_las.p3000.n50  <- mean(rmsebeta_las)
      rmsebeta_alas.p3000.n50 <- mean(rmsebeta_alas)
      rmsebeta_scad.p3000.n50 <- mean(rmsebeta_scad)
      
      Rsq_las.p3000.n50  <- mean(Rsq_las)
      Rsq_alas.p3000.n50 <- mean(Rsq_alas)
      Rsq_scad.p3000.n50 <- mean(Rsq_scad)
      
      msef_las.p3000.n50  <- mean(msef_las)
      msef_alas.p3000.n50 <- mean(msef_alas)
      msef_scad.p3000.n50 <- mean(msef_scad)
      
      sens_las.p3000.n50 <- mean(sens_las)
      spec_las.p3000.n50 <- mean(spec_las) 
      acc_las.p3000.n50  <- mean(acc_las)
      G_las.p3000.n50    <- mean(G_las)
      
      sens_alas.p3000.n50 <- mean(sens_alas)
      spec_alas.p3000.n50 <- mean(spec_alas) 
      acc_alas.p3000.n50  <- mean(acc_alas)
      G_alas.p3000.n50    <- mean(G_alas)
      
      sens_scad.p3000.n50 <- mean(sens_scad)
      spec_scad.p3000.n50 <- mean(spec_scad) 
      acc_scad.p3000.n50  <- mean(acc_scad)
      G_scad.p3000.n50    <- mean(G_scad)
      
      RMSE_beta.p3000.n50             <- data.frame(rmsebeta_las.p3000.n50,rmsebeta_alas.p3000.n50,rmsebeta_scad.p3000.n50)
      colnames(RMSE_beta.p3000.n50)   <- c("Lasso","Adaptive-Lasso","SCAD")
      
      RSQ.p3000.n50                   <- data.frame(Rsq_las.p3000.n50,Rsq_alas.p3000.n50,Rsq_scad.p3000.n50)
      colnames(RSQ.p3000.n50)         <- c("Lasso","Adaptive-Lasso","SCAD")
      
      MSEF.p3000.n50                  <- data.frame(msef_las.p3000.n50,msef_alas.p3000.n50,msef_scad.p3000.n50)
      colnames(MSEF.p3000.n50)        <- c("Lasso","Adaptive-Lasso","SCAD")
      
      metrics_all.p3000.n50           <- rbind(RMSE_beta.p3000.n50,RSQ.p3000.n50,MSEF.p3000.n50)
      rownames(metrics_all.p3000.n50) <- c("REMSE","R-square","MSEf") 
      
      stong_signals.p3000.n50 <- data.frame((strong_las),(strong_alas),(strong_scad))
      weak_signals.p3000.n50 <- data.frame((weak_las),(weak_alas),(weak_scad))
      
      RECALL.las.p3000.n50                     <- data.frame(sens_las.p3000.n50,spec_las.p3000.n50,acc_las.p3000.n50,G_las.p3000.n50)
      colnames(RECALL.las.p3000.n50)           <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.alas.p3000.n50                    <- data.frame(sens_alas.p3000.n50,spec_alas.p3000.n50,acc_alas.p3000.n50,G_alas.p3000.n50)
      colnames(RECALL.alas.p3000.n50)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.scad.p3000.n50                    <- data.frame(sens_scad.p3000.n50,spec_scad.p3000.n50,acc_scad.p3000.n50,G_scad.p3000.n50)
      colnames(RECALL.scad.p3000.n50)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL_all.p3000.n50                     <- rbind(RECALL.las.p3000.n50,RECALL.alas.p3000.n50,RECALL.scad.p3000.n50)
      rownames(RECALL_all.p3000.n50)           <- c("Lasso","Adaptive-Lasso","SCAD")
    }
    if (s2==1 & s1==2){
      WR.lasso.fhat.p300.n100  <-rowMeans(Est.f.LassoWR) 
      WR.alasso.fhat.p300.n100 <-rowMeans(Est.f.aLassoWR)
      WR.scad.fhat.p300.n100   <-rowMeans(Est.f.SCADWR)
      
      rmsebeta_las.p300.n100   <- mean(rmsebeta_las)
      rmsebeta_alas.p300.n100 <- mean(rmsebeta_alas)
      rmsebeta_scad.p300.n100 <- mean(rmsebeta_scad)
      
      Rsq_las.p300.n100  <- mean(Rsq_las)
      Rsq_alas.p300.n100 <- mean(Rsq_alas)
      Rsq_scad.p300.n100 <- mean(Rsq_scad)
      
      msef_las.p300.n100  <- mean(msef_las)
      msef_alas.p300.n100 <- mean(msef_alas)
      msef_scad.p300.n100 <- mean(msef_scad)
      
      sens_las.p300.n100 <- mean(sens_las)
      spec_las.p300.n100 <- mean(spec_las) 
      acc_las.p300.n100  <- mean(acc_las)
      G_las.p300.n100    <- mean(G_las)
      
      sens_alas.p300.n100 <- mean(sens_alas)
      spec_alas.p300.n100 <- mean(spec_alas) 
      acc_alas.p300.n100  <- mean(acc_alas)
      G_alas.p300.n100    <- mean(G_alas)
      
      sens_scad.p300.n100 <- mean(sens_scad)
      spec_scad.p300.n100 <- mean(spec_scad) 
      acc_scad.p300.n100  <- mean(acc_scad)
      G_scad.p300.n100    <- mean(G_scad)
      
      RMSE_beta.p300.n100             <- data.frame(rmsebeta_las.p300.n100,rmsebeta_alas.p300.n100,rmsebeta_scad.p300.n100)
      colnames(RMSE_beta.p300.n100)   <- c("Lasso","Adaptive-Lasso","SCAD")
      
      RSQ.p300.n100                   <- data.frame(Rsq_las.p300.n100,Rsq_alas.p300.n100,Rsq_scad.p300.n100)
      colnames(RSQ.p300.n100)         <- c("Lasso","Adaptive-Lasso","SCAD")
      
      MSEF.p300.n100                  <- data.frame(msef_las.p300.n100,msef_alas.p300.n100,msef_scad.p300.n100)
      colnames(MSEF.p300.n100)        <- c("Lasso","Adaptive-Lasso","SCAD")
      
      metrics_all.p300.n100           <- rbind(RMSE_beta.p300.n100,RSQ.p300.n100,MSEF.p300.n100)
      rownames(metrics_all.p300.n100) <- c("REMSE","R-square","MSEf") 
      
      stong_signals.p300.n100 <- data.frame((strong_las),(strong_alas),(strong_scad))
      weak_signals.p300.n100 <- data.frame((weak_las),(weak_alas),(weak_scad))
      
      RECALL.las.p300.n100                     <- data.frame(sens_las.p300.n100,spec_las.p300.n100,acc_las.p300.n100,G_las.p300.n100)
      colnames(RECALL.las.p300.n100)           <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.alas.p300.n100                    <- data.frame(sens_alas.p300.n100,spec_alas.p300.n100,acc_alas.p300.n100,G_alas.p300.n100)
      colnames(RECALL.alas.p300.n100)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.scad.p300.n100                    <- data.frame(sens_scad.p300.n100,spec_scad.p300.n100,acc_scad.p300.n100,G_scad.p300.n100)
      colnames(RECALL.scad.p300.n100)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL_all.p300.n100                     <- rbind(RECALL.las.p300.n100,RECALL.alas.p300.n100,RECALL.scad.p300.n100)
      rownames(RECALL_all.p300.n100)           <- c("Lasso","Adaptive-Lasso","SCAD")
      
    }
    if (s2==2 & s1==2){
      WR.lasso.fhat.p1000.n100  <-rowMeans(Est.f.LassoWR) 
      WR.alasso.fhat.p1000.n100 <-rowMeans(Est.f.aLassoWR)
      WR.scad.fhat.p1000.n100   <-rowMeans(Est.f.SCADWR)
      
      rmsebeta_las.p1000.n100   <- mean(rmsebeta_las)
      rmsebeta_alas.p1000.n100  <- mean(rmsebeta_alas)
      rmsebeta_scad.p1000.n100  <- mean(rmsebeta_scad)
      
      Rsq_las.p1000.n100  <- mean(Rsq_las)
      Rsq_alas.p1000.n100 <- mean(Rsq_alas)
      Rsq_scad.p1000.n100 <- mean(Rsq_scad)
      
      msef_las.p1000.n100  <- mean(msef_las)
      msef_alas.p1000.n100 <- mean(msef_alas)
      msef_scad.p1000.n100 <- mean(msef_scad)
      
      sens_las.p1000.n100 <- mean(sens_las)
      spec_las.p1000.n100 <- mean(spec_las) 
      acc_las.p1000.n100  <- mean(acc_las)
      G_las.p1000.n100    <- mean(G_las)
      
      sens_alas.p1000.n100 <- mean(sens_alas)
      spec_alas.p1000.n100 <- mean(spec_alas) 
      acc_alas.p1000.n100  <- mean(acc_alas)
      G_alas.p1000.n100    <- mean(G_alas)
      
      sens_scad.p1000.n100 <- mean(sens_scad)
      spec_scad.p1000.n100 <- mean(spec_scad) 
      acc_scad.p1000.n100  <- mean(acc_scad)
      G_scad.p1000.n100    <- mean(G_scad)
      
      RMSE_beta.p1000.n100             <- data.frame(rmsebeta_las.p1000.n100,rmsebeta_alas.p1000.n100,rmsebeta_scad.p1000.n100)
      colnames(RMSE_beta.p1000.n100)   <- c("Lasso","Adaptive-Lasso","SCAD")
      
      RSQ.p1000.n100                   <- data.frame(Rsq_las.p1000.n100,Rsq_alas.p1000.n100,Rsq_scad.p1000.n100)
      colnames(RSQ.p1000.n100)         <- c("Lasso","Adaptive-Lasso","SCAD")
      
      MSEF.p1000.n100                  <- data.frame(msef_las.p1000.n100,msef_alas.p1000.n100,msef_scad.p1000.n100)
      colnames(MSEF.p1000.n100)        <- c("Lasso","Adaptive-Lasso","SCAD")
      
      metrics_all.p1000.n100           <- rbind(RMSE_beta.p1000.n100,RSQ.p1000.n100,MSEF.p1000.n100)
      rownames(metrics_all.p1000.n100) <- c("REMSE","R-square","MSEf") 
      
      stong_signals.p1000.n100 <- data.frame((strong_las),(strong_alas),(strong_scad))
      weak_signals.p1000.n100 <- data.frame((weak_las),(weak_alas),(weak_scad))
      
      RECALL.las.p1000.n100                     <- data.frame(sens_las.p1000.n100,spec_las.p1000.n100,acc_las.p1000.n100,G_las.p1000.n100)
      colnames(RECALL.las.p1000.n100)           <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.alas.p1000.n100                    <- data.frame(sens_alas.p1000.n100,spec_alas.p1000.n100,acc_alas.p1000.n100,G_alas.p1000.n100)
      colnames(RECALL.alas.p1000.n100)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.scad.p1000.n100                    <- data.frame(sens_scad.p1000.n100,spec_scad.p1000.n100,acc_scad.p1000.n100,G_scad.p1000.n100)
      colnames(RECALL.scad.p1000.n100)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL_all.p1000.n100                     <- rbind(RECALL.las.p1000.n100,RECALL.alas.p1000.n100,RECALL.scad.p1000.n100)
      rownames(RECALL_all.p1000.n100)           <- c("Lasso","Adaptive-Lasso","SCAD")
      
    }
    if (s2==3 & s1==2){
      WR.lasso.fhat.p3000.n100  <-rowMeans(Est.f.LassoWR) 
      WR.alasso.fhat.p3000.n100 <-rowMeans(Est.f.aLassoWR)
      WR.scad.fhat.p3000.n100   <-rowMeans(Est.f.SCADWR)
      
      rmsebeta_las.p3000.n100   <- mean(rmsebeta_las)
      rmsebeta_alas.p3000.n100  <- mean(rmsebeta_alas)
      rmsebeta_scad.p3000.n100  <- mean(rmsebeta_scad)
      
      Rsq_las.p3000.n100  <- mean(Rsq_las)
      Rsq_alas.p3000.n100 <- mean(Rsq_alas)
      Rsq_scad.p3000.n100 <- mean(Rsq_scad)
      
      msef_las.p3000.n100  <- mean(msef_las)
      msef_alas.p3000.n100 <- mean(msef_alas)
      msef_scad.p3000.n100 <- mean(msef_scad)
      
      sens_las.p3000.n100 <- mean(sens_las)
      spec_las.p3000.n100 <- mean(spec_las) 
      acc_las.p3000.n100  <- mean(acc_las)
      G_las.p3000.n100    <- mean(G_las)
      
      sens_alas.p3000.n100 <- mean(sens_alas)
      spec_alas.p3000.n100 <- mean(spec_alas) 
      acc_alas.p3000.n100  <- mean(acc_alas)
      G_alas.p3000.n100    <- mean(G_alas)
      
      sens_scad.p3000.n100 <- mean(sens_scad)
      spec_scad.p3000.n100 <- mean(spec_scad) 
      acc_scad.p3000.n100  <- mean(acc_scad)
      G_scad.p3000.n100    <- mean(G_scad)
      
      RMSE_beta.p3000.n100             <- data.frame(rmsebeta_las.p3000.n100,rmsebeta_alas.p3000.n100,rmsebeta_scad.p3000.n100)
      colnames(RMSE_beta.p3000.n100)   <- c("Lasso","Adaptive-Lasso","SCAD")
      
      RSQ.p3000.n100                   <- data.frame(Rsq_las.p3000.n100,Rsq_alas.p3000.n100,Rsq_scad.p3000.n100)
      colnames(RSQ.p3000.n100)         <- c("Lasso","Adaptive-Lasso","SCAD")
      
      MSEF.p3000.n100                  <- data.frame(msef_las.p3000.n100,msef_alas.p3000.n100,msef_scad.p3000.n100)
      colnames(MSEF.p3000.n100)        <- c("Lasso","Adaptive-Lasso","SCAD")
      
      metrics_all.p3000.n100           <- rbind(RMSE_beta.p3000.n100,RSQ.p3000.n100,MSEF.p3000.n100)
      rownames(metrics_all.p3000.n100) <- c("REMSE","R-square","MSEf") 
      
      stong_signals.p3000.n100 <- data.frame((strong_las),(strong_alas),(strong_scad))
      weak_signals.p3000.n100 <- data.frame((weak_las),(weak_alas),(weak_scad))
      
      RECALL.las.p3000.n100                     <- data.frame(sens_las.p3000.n100,spec_las.p3000.n100,acc_las.p3000.n100,G_las.p3000.n100)
      colnames(RECALL.las.p3000.n100)           <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.alas.p3000.n100                    <- data.frame(sens_alas.p3000.n100,spec_alas.p3000.n100,acc_alas.p3000.n100,G_alas.p3000.n100)
      colnames(RECALL.alas.p3000.n100)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.scad.p3000.n100                    <- data.frame(sens_scad.p3000.n100,spec_scad.p3000.n100,acc_scad.p3000.n100,G_scad.p3000.n100)
      colnames(RECALL.scad.p3000.n100)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL_all.p3000.n100                     <- rbind(RECALL.las.p3000.n100,RECALL.alas.p3000.n100,RECALL.scad.p3000.n100)
      rownames(RECALL_all.p3000.n100)           <- c("Lasso","Adaptive-Lasso","SCAD")
    }
    if (s2==1 & s1==3){
      WR.lasso.fhat.p300.n200  <-rowMeans(Est.f.LassoWR) 
      WR.alasso.fhat.p300.n200 <-rowMeans(Est.f.aLassoWR)
      WR.scad.fhat.p300.n200   <-rowMeans(Est.f.SCADWR)
      
      rmsebeta_las.p300.n200   <- mean(rmsebeta_las)
      rmsebeta_alas.p300.n200 <- mean(rmsebeta_alas)
      rmsebeta_scad.p300.n200 <- mean(rmsebeta_scad)
      
      Rsq_las.p300.n200  <- mean(Rsq_las)
      Rsq_alas.p300.n200 <- mean(Rsq_alas)
      Rsq_scad.p300.n200 <- mean(Rsq_scad)
      
      msef_las.p300.n200  <- mean(msef_las)
      msef_alas.p300.n200 <- mean(msef_alas)
      msef_scad.p300.n200 <- mean(msef_scad)
      
      sens_las.p300.n200 <- mean(sens_las)
      spec_las.p300.n200 <- mean(spec_las) 
      acc_las.p300.n200  <- mean(acc_las)
      G_las.p300.n200    <- mean(G_las)
      
      sens_alas.p300.n200 <- mean(sens_alas)
      spec_alas.p300.n200 <- mean(spec_alas) 
      acc_alas.p300.n200  <- mean(acc_alas)
      G_alas.p300.n200    <- mean(G_alas)
      
      sens_scad.p300.n200 <- mean(sens_scad)
      spec_scad.p300.n200 <- mean(spec_scad) 
      acc_scad.p300.n200  <- mean(acc_scad)
      G_scad.p300.n200    <- mean(G_scad)
      
      RMSE_beta.p300.n200             <- data.frame(rmsebeta_las.p300.n200,rmsebeta_alas.p300.n200,rmsebeta_scad.p300.n200)
      colnames(RMSE_beta.p300.n200)   <- c("Lasso","Adaptive-Lasso","SCAD")
      
      RSQ.p300.n200                   <- data.frame(Rsq_las.p300.n200,Rsq_alas.p300.n200,Rsq_scad.p300.n200)
      colnames(RSQ.p300.n200)         <- c("Lasso","Adaptive-Lasso","SCAD")
      
      MSEF.p300.n200                  <- data.frame(msef_las.p300.n200,msef_alas.p300.n200,msef_scad.p300.n200)
      colnames(MSEF.p300.n200)        <- c("Lasso","Adaptive-Lasso","SCAD")
      
      metrics_all.p300.n200           <- rbind(RMSE_beta.p300.n200,RSQ.p300.n200,MSEF.p300.n200)
      rownames(metrics_all.p300.n200) <- c("REMSE","R-square","MSEf") 
      
      stong_signals.p300.n200 <- data.frame((strong_las),(strong_alas),(strong_scad))
      weak_signals.p300.n200  <- data.frame((weak_las),(weak_alas),(weak_scad))
      
      RECALL.las.p300.n200                     <- data.frame(sens_las.p300.n200,spec_las.p300.n200,acc_las.p300.n200,G_las.p300.n200)
      colnames(RECALL.las.p300.n200)           <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.alas.p300.n200                    <- data.frame(sens_alas.p300.n200,spec_alas.p300.n200,acc_alas.p300.n200,G_alas.p300.n200)
      colnames(RECALL.alas.p300.n200)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.scad.p300.n200                    <- data.frame(sens_scad.p300.n200,spec_scad.p300.n200,acc_scad.p300.n200,G_scad.p300.n200)
      colnames(RECALL.scad.p300.n200)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL_all.p300.n200                     <- rbind(RECALL.las.p300.n200,RECALL.alas.p300.n200,RECALL.scad.p300.n200)
      rownames(RECALL_all.p300.n200)           <- c("Lasso","Adaptive-Lasso","SCAD")
      
    }
    if (s2==2 & s1==3){
      WR.lasso.fhat.p1000.n200  <-rowMeans(Est.f.LassoWR) 
      WR.alasso.fhat.p1000.n200 <-rowMeans(Est.f.aLassoWR)
      WR.scad.fhat.p1000.n200   <-rowMeans(Est.f.SCADWR)
      
      rmsebeta_las.p1000.n200   <- mean(rmsebeta_las)
      rmsebeta_alas.p1000.n200 <- mean(rmsebeta_alas)
      rmsebeta_scad.p1000.n200 <- mean(rmsebeta_scad)
      
      Rsq_las.p1000.n200  <- mean(Rsq_las)
      Rsq_alas.p1000.n200 <- mean(Rsq_alas)
      Rsq_scad.p1000.n200 <- mean(Rsq_scad)
      
      msef_las.p1000.n200  <- mean(msef_las)
      msef_alas.p1000.n200 <- mean(msef_alas)
      msef_scad.p1000.n200 <- mean(msef_scad)
      
      sens_las.p1000.n200 <- mean(sens_las)
      spec_las.p1000.n200 <- mean(spec_las) 
      acc_las.p1000.n200  <- mean(acc_las)
      G_las.p1000.n200    <- mean(G_las)
      
      sens_alas.p1000.n200 <- mean(sens_alas)
      spec_alas.p1000.n200 <- mean(spec_alas) 
      acc_alas.p1000.n200  <- mean(acc_alas)
      G_alas.p1000.n200    <- mean(G_alas)
      
      sens_scad.p1000.n200 <- mean(sens_scad)
      spec_scad.p1000.n200 <- mean(spec_scad) 
      acc_scad.p1000.n200  <- mean(acc_scad)
      G_scad.p1000.n200    <- mean(G_scad)
      
      stong_signals.p1000.n200 <- data.frame((strong_las),(strong_alas),(strong_scad))
      weak_signals.p1000.n200  <- data.frame((weak_las),(weak_alas),(weak_scad))
      
      RMSE_beta.p1000.n200             <- data.frame(rmsebeta_las.p1000.n200,rmsebeta_alas.p1000.n200,rmsebeta_scad.p1000.n200)
      colnames(RMSE_beta.p1000.n200)   <- c("Lasso","Adaptive-Lasso","SCAD")
      
      RSQ.p1000.n200                   <- data.frame(Rsq_las.p1000.n200,Rsq_alas.p1000.n200,Rsq_scad.p1000.n200)
      colnames(RSQ.p1000.n200)         <- c("Lasso","Adaptive-Lasso","SCAD")
      
      MSEF.p1000.n200                  <- data.frame(msef_las.p1000.n200,msef_alas.p1000.n200,msef_scad.p1000.n200)
      colnames(MSEF.p1000.n200)        <- c("Lasso","Adaptive-Lasso","SCAD")
      
      metrics_all.p1000.n200           <- rbind(RMSE_beta.p1000.n200,RSQ.p1000.n200,MSEF.p1000.n200)
      rownames(metrics_all.p1000.n200) <- c("REMSE","R-square","MSEf") 
      
      RECALL.las.p1000.n200                     <- data.frame(sens_las.p1000.n200,spec_las.p1000.n200,acc_las.p1000.n200,G_las.p1000.n200)
      colnames(RECALL.las.p1000.n200)           <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.alas.p1000.n200                    <- data.frame(sens_alas.p1000.n200,spec_alas.p1000.n200,acc_alas.p1000.n200,G_alas.p1000.n200)
      colnames(RECALL.alas.p1000.n200)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.scad.p1000.n200                    <- data.frame(sens_scad.p1000.n200,spec_scad.p1000.n200,acc_scad.p1000.n200,G_scad.p1000.n200)
      colnames(RECALL.scad.p1000.n200)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL_all.p1000.n200                     <- rbind(RECALL.las.p1000.n200,RECALL.alas.p1000.n200,RECALL.scad.p1000.n200)
      rownames(RECALL_all.p1000.n200)           <- c("Lasso","Adaptive-Lasso","SCAD")
      
    }
    if (s2==3 & s1==3){
      WR.lasso.fhat.p3000.n200  <-rowMeans(Est.f.LassoWR) 
      WR.alasso.fhat.p3000.n200 <-rowMeans(Est.f.aLassoWR)
      WR.scad.fhat.p3000.n200   <-rowMeans(Est.f.SCADWR)
      
      rmsebeta_las.p3000.n200   <- mean(rmsebeta_las)
      rmsebeta_alas.p3000.n200 <- mean(rmsebeta_alas)
      rmsebeta_scad.p3000.n200 <- mean(rmsebeta_scad)
      
      Rsq_las.p3000.n200  <- mean(Rsq_las)
      Rsq_alas.p3000.n200 <- mean(Rsq_alas)
      Rsq_scad.p3000.n200 <- mean(Rsq_scad)
      
      msef_las.p3000.n200  <- mean(msef_las)
      msef_alas.p3000.n200 <- mean(msef_alas)
      msef_scad.p3000.n200 <- mean(msef_scad)
      
      sens_las.p3000.n200 <- mean(sens_las)
      spec_las.p3000.n200 <- mean(spec_las) 
      acc_las.p3000.n200  <- mean(acc_las)
      G_las.p3000.n200    <- mean(G_las)
      
      sens_alas.p3000.n200 <- mean(sens_alas)
      spec_alas.p3000.n200 <- mean(spec_alas) 
      acc_alas.p3000.n200  <- mean(acc_alas)
      G_alas.p3000.n200    <- mean(G_alas)
      
      sens_scad.p3000.n200 <- mean(sens_scad)
      spec_scad.p3000.n200 <- mean(spec_scad) 
      acc_scad.p3000.n200  <- mean(acc_scad)
      G_scad.p3000.n200    <- mean(G_scad)
      
      RMSE_beta.p3000.n200             <- data.frame(rmsebeta_las.p3000.n200,rmsebeta_alas.p3000.n200,rmsebeta_scad.p3000.n200)
      colnames(RMSE_beta.p3000.n200)   <- c("Lasso","Adaptive-Lasso","SCAD")
      
      RSQ.p3000.n200                   <- data.frame(Rsq_las.p3000.n200,Rsq_alas.p3000.n200,Rsq_scad.p3000.n200)
      colnames(RSQ.p3000.n200)         <- c("Lasso","Adaptive-Lasso","SCAD")
      
      MSEF.p3000.n200                  <- data.frame(msef_las.p3000.n200,msef_alas.p3000.n200,msef_scad.p3000.n200)
      colnames(MSEF.p3000.n200)        <- c("Lasso","Adaptive-Lasso","SCAD")
      
      metrics_all.p3000.n200           <- rbind(RMSE_beta.p3000.n200,RSQ.p3000.n200,MSEF.p3000.n200)
      rownames(metrics_all.p3000.n200) <- c("REMSE","R-square","MSEf") 
      
      stong_signals.p3000.n200 <- data.frame((strong_las),(strong_alas),(strong_scad))
      weak_signals.p3000.n200  <- data.frame((weak_las),(weak_alas),(weak_scad))
      
      RECALL.las.p3000.n200                     <- data.frame(sens_las.p3000.n200,spec_las.p3000.n200,acc_las.p3000.n200,G_las.p3000.n200)
      colnames(RECALL.las.p3000.n200)           <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.alas.p3000.n200                    <- data.frame(sens_alas.p3000.n200,spec_alas.p3000.n200,acc_alas.p3000.n200,G_alas.p3000.n200)
      colnames(RECALL.alas.p3000.n200)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.scad.p3000.n200                    <- data.frame(sens_scad.p3000.n200,spec_scad.p3000.n200,acc_scad.p3000.n200,G_scad.p3000.n200)
      colnames(RECALL.scad.p3000.n200)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL_all.p3000.n200                     <- rbind(RECALL.las.p3000.n200,RECALL.alas.p3000.n200,RECALL.scad.p3000.n200)
      rownames(RECALL_all.p3000.n200)           <- c("Lasso","Adaptive-Lasso","SCAD")
      
    }
  }
}
df_list.p300  <- rbind(metrics_all.p300.n50,metrics_all.p300.n100,metrics_all.p300.n200)
df_list.p1000 <- rbind(metrics_all.p1000.n50,metrics_all.p1000.n100,metrics_all.p1000.n200)
df_list.p3000 <- rbind(metrics_all.p3000.n50,metrics_all.p3000.n100,metrics_all.p3000.n200)


df_list.p300
df_list.p1000 
df_list.p3000

strong.signals.p300 <-cbind(stong_signals.p300.n50,stong_signals.p300.n100,stong_signals.p300.n200) 
weak.signals.p300   <-cbind(weak_signals.p300.n50,weak_signals.p300.n100,weak_signals.p300.n200)

strong.signals.p1000 <-cbind(stong_signals.p1000.n50,stong_signals.p1000.n100,stong_signals.p1000.n200) 
weak.signals.p1000   <-cbind(weak_signals.p1000.n50,weak_signals.p1000.n100,weak_signals.p1000.n200)

strong.signals.p3000 <-cbind(stong_signals.p300.n50,stong_signals.p300.n100,stong_signals.p300.n200) 
weak.signals.p3000   <-cbind(weak_signals.p300.n50,weak_signals.p300.n100,weak_signals.p300.n200)

par(mfrow=c(3,1))
boxplot(strong.signals.p300,main="Strong signals for pn=300")
boxplot(strong.signals.p1000,main="Strong signals for pn=1000")
boxplot(strong.signals.p3000,main="Strong signals for pn=3000")

par(mfrow=c(3,1))
boxplot(weak.signals.p300,main="Weak signals for pn=300")
boxplot(weak.signals.p1000,main="Weak signals for pn=1000")
boxplot(weak.signals.p3000,main="Weak signals for pn=3000")

recall.p300 <- data.frame(sens_las.p300.n50,sens_alas.p300.n50,sens_scad.p300.n50,acc_las.p300.n50,acc_alas.p300.n50,acc_scad.p300.n50)
recall.p1000 <- data.frame(sens_las.p1000.n50,sens_alas.p1000.n50,sens_scad.p1000.n50,acc_las.p1000.n50,acc_alas.p1000.n50,acc_scad.p1000.n50)
recall.p3000 <- data.frame(sens_las.p3000.n50,sens_alas.p3000.n50,sens_scad.p3000.n50,acc_las.p3000.n50,acc_alas.p3000.n50,acc_scad.p3000.n50)

recall.p300
recall.p1000
recall.p3000
###################################################################################################
plot(f,type="l",ylim=c(min(f),max(f)),col=1,main="n=200, pn=3000, Scenario 1, rho=0.80",xlab="t",ylab="f(t)")
par(new=TRUE)
plot(WR.lasso.fhat.p3000.n200,type="l",ylim=c(min(f),max(f)),col=2,xlab="t",ylab="f(t)")
par(new=TRUE)
plot(WR.alasso.fhat.p3000.n200+0.1,type="l",ylim=c(min(f),max(f)),col=3,xlab="t",ylab="f(t)")
par(new=TRUE)
plot(WR.scad.fhat.p3000.n200,type="l",ylim=c(min(f),max(f)),col=4,xlab="t",ylab="f(t)")
legend("bottomleft",legend=c("Real f", "Lasso-PSE","aLasso-PSE","SCAD-PSE"),
       col=c(1, 2,3,4),lty=c(1,1,1,1))
grid()

write.csv(df_list.p300,"p300sc2rho1.csv")
write.csv(df_list.p1000,"p800sc2rho1.csv")
write.csv(df_list.p3000,"p1500sc2rho1.csv")

write.csv(recall.p300,"RECALL.p300sc2rho1.csv")
write.csv(recall.p1000,"RECALL.p800sc2rho1.csv")
write.csv(recall.p3000,"RECALL.p1500sc2rho1.csv")

############################################################################################################
#===========================================================================================================
############################################################################################################
#SCENARIO 1 RHO 0.95

#OUTLINE-----------------------------------------------
#Function 1: Data generation for simulation Sudy
#Function 2: Estimators based on Penalty functions
#Function 3: Calculation of evaluation metrics
#------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------
ni <- c(50,100,120)
pni <- c(300,350,400)
rhoi <- c(0.80,0.95)              #Number of simulation
for (s1 in 1:3){                  #Sample size index
  for (s2 in 1:3){                #pn index
    if (s1==1){n=50} 
    if (s1==2){n=100}  
    if (s1==3){n=200}
    if (s2==1){pn=300} 
    if (s2==2){pn=600} 
    if (s2==3){pn=1200}
    rho <- 0.95
    sc  <- 1
    #ZERO MATRICES----------------------------------------------------------
    rmsebeta_las    <- 0
    rmsebeta_alas   <- 0
    rmsebeta_scad   <- 0
    
    Rsq_las         <- 0
    Rsq_alas        <- 0
    Rsq_scad        <- 0
    
    msef_las         <- 0
    msef_alas        <- 0
    msef_scad        <- 0
    
    sens_las <- 0
    spec_las <- 0
    acc_las  <- 0
    G_las    <- 0
    sens_alas <- 0
    spec_alas <- 0
    acc_alas  <- 0
    G_alas    <- 0
    sens_scad <- 0
    spec_scad <- 0
    acc_scad  <- 0
    G_scad    <- 0
    
    strong_las  <- 0
    weak_las    <- 0
    strong_alas <- 0
    weak_alas   <- 0 
    strong_scad <- 0 
    weak_scad   <- 0
    
    Est.f.LassoWR  <- matrix(0,n,sim)
    Est.f.aLassoWR <- matrix(0,n,sim)
    Est.f.SCADWR   <- matrix(0,n,sim)
    
    Est.f.Lasso  <- matrix(0,n,sim)
    Est.f.aLasso <- matrix(0,n,sim)
    Est.f.SCAD   <- matrix(0,n,sim)
    
    Est.y.LassoWR  <- matrix(0,n,sim)
    Est.y.aLassoWR <- matrix(0,n,sim)
    Est.y.SCADWR   <- matrix(0,n,sim)
    
    Est.y.Lasso  <- matrix(0,n,sim)
    Est.y.aLasso <- matrix(0,n,sim)
    Est.y.SCAD   <- matrix(0,n,sim)
    #------------------------------------------------------------------------------
    no.seed <- 19345
    aa <- 1
    for (s in 1:sim){
      #Set seed number needs to be changed in every repeat----------------------------
      set.seed(no.seed+aa)
      mydata <- simdata(n,pn,sc,rho)
      set.seed(no.seed+s)
      #TEST AREA (To be Deleted)------------------------------------------------------
      x     <- mydata$X
      y     <- mydata$y
      t     <- mydata$t
      f     <- mydata$f 
      nol   <- 25
      aiclas  <- 0
      aicalas <- 0
      aicscad <- 0
      
      las     <- lasso(x,y)
      alas    <- adalasso(x,y)
      scadm   <- scad(x,y)
      
      lams    <- seq(0.1,1.5,length.out=nol)
      for (j in 1:nol){
        aiclas[j]  <- aicfunc(y,las,lams[j],f)
        aicalas[j]  <- aicfunc(y,alas,lams[j],f)
        aicscad[j]  <- aicfunc(y,scadm,lams[j],f)
      }
      #PLOTS for smoothing parameter selection (AICc)---------------------------------
      dflas  <- data.frame(lams,aiclas)
      dfalas <- data.frame(lams,aicalas)
      dfscad <- data.frame(lams,aicscad)
      #ggplot()+geom_line(data=dflas,aes(x=lams,aiclas,col="Lasso_WR"),lwd=1)+geom_line(data=dfalas,aes(x=lams,aicalas,col="Adapt.-Lasso_WR"),lwd=1)+geom_line(data=dfscad,aes(x=lams,aicscad,col="SCAD_WR"),lwd=1)+xlab("lambda_spline")+ylab("Improved_AIC score")+ggtitle("(e) n=100, pn=3000, rho=0.95")+theme(legend.title = element_blank(),legend.position = "bottom")
      #-------------------------------------------------------------------------------
      for (j2 in 1:nol){
        if (aiclas[j2]==min(aiclas)){
          lam_las <-  lams[j2]
        }
        if (aicalas[j2]==min(aicalas)){
          lam_alas <-  lams[j2]
        }
        if (aicscad[j2]==min(aicscad)){
          lam_scad <-  lams[j2]
        }
      }
      S_las  <- Smatrix(t,2,lam_las)
      S_alas <- Smatrix(t,2,lam_alas)
      S_scad <- Smatrix(t,2,lam_scad)
      
      
      ytil_las  <- (diag(n)-S_las)%*%y
      ytil_alas <- (diag(n)-S_alas)%*%y
      ytil_scad <- (diag(n)-S_scad)%*%y
      
      
      xtil_las  <- (diag(n)-S_las)%*%x
      xtil_alas <- (diag(n)-S_alas)%*%x
      xtil_scad <- (diag(n)-S_scad)%*%x
      
      newX_las  <- (diag(n)-S_las)%*%las$X_S1
      newX_alas <- (diag(n)-S_alas)%*%alas$X_S1
      newX_scad <- (diag(n)-S_scad)%*%scadm$X_S1
      #Estimation based on Penalty functions
      las_beta      <- solve(t(newX_las)%*%newX_las+0.01*diag(las$S1))%*%t(newX_las)%*%ytil_las
      fhat_las      <- S_las%*%(y-las$X_S1%*%las_beta)
      alas_beta     <- solve(t(newX_alas)%*%newX_alas+0.01*diag(alas$S1))%*%t(newX_alas)%*%ytil_alas
      fhat_alas     <- S_alas%*%(y-alas$X_S1%*%alas_beta)
      scad_beta     <- solve(t(newX_scad)%*%newX_scad+0.01*diag(scadm$S1))%*%t(newX_scad)%*%ytil_scad
      fhat_scad     <- S_scad%*%(y-scadm$X_S1%*%scad_beta)
      
      yhat_las  <- las$X_S1%*%las_beta+fhat_las
      yhat_alas <- alas$X_S1%*%alas_beta+fhat_alas
      yhat_scad <- scadm$X_S1%*%scad_beta+fhat_scad
      
      WRobj_las  <- WR(las,S_las,x,y,0.01)
      WRobj_alas <- WR(alas,S_alas,x,y,0.01)
      WRobj_scad <- WR(scadm,S_scad,x,y,0.01)
      
      WR_fhat_las  <- S_las%*%(y-las$X_S1%*%WRobj_las$betaPSE)
      WR_fhat_alas <- S_alas%*%(y-alas$X_S1%*%WRobj_alas$betaPSE)
      WR_fhat_scad <- S_scad%*%(y-scadm$X_S1%*%WRobj_scad$betaPSE)
      
      WR_yhat_las  <- las$X_S1%*%WRobj_las$betaPSE+fhat_las
      WR_yhat_alas <- alas$X_S1%*%WRobj_alas$betaPSE+fhat_alas
      WR_yhat_scad <- scadm$X_S1%*%WRobj_scad$betaPSE
      #---------------------------------------------------------------------
      las_index       <- las$S1
      real_beta_las   <- mydata$beta[las_index]
      alas_index      <- alas$S1
      real_beta_alas  <- mydata$beta[alas_index]
      scad_index      <- scadm$S1
      real_beta_scad  <- mydata$beta[scad_index]
      
      rmsebeta_las[s]  <- rmsebeta(real_beta_las,WRobj_las$betaPSE)
      rmsebeta_alas[s] <- rmsebeta(real_beta_alas,WRobj_alas$betaPSE)
      rmsebeta_scad[s] <- rmsebeta(real_beta_scad,WRobj_scad$betaPSE)
      
      Rsq_las[s]  <- rsq(mydata$y,WR_yhat_las,length(las_index))
      Rsq_alas[s] <- rsq(mydata$y,WR_yhat_alas,length(alas_index))
      Rsq_scad[s] <- rsq(mydata$y,WR_yhat_scad,length(scad_index))
      
      msef_las[s]  <- msef(mydata$f,WR_fhat_las)
      msef_alas[s] <- msef(mydata$f,WR_fhat_alas)
      msef_scad[s] <- msef(mydata$f,WR_fhat_scad)
      
      recall_las  <- recall(real_beta_las,WRobj_las$betaPSE)
      sens_las[s] <- recall_las$sensitivity
      spec_las[s] <- recall_las$specifcity
      acc_las[s]  <- recall_las$accuracy
      G_las[s]    <- recall_las$Gscore
      
      recall_alas  <- recall(real_beta_alas,WRobj_alas$betaPSE)
      sens_alas[s] <- recall_alas$sensitivity
      spec_alas[s] <- recall_alas$specifcity
      acc_alas[s]  <- recall_alas$accuracy
      G_alas[s]    <- recall_alas$Gscore
      
      recall_scad  <- recall(real_beta_scad,WRobj_scad$betaPSE)
      sens_scad[s] <- recall_scad$sensitivity
      spec_scad[s] <- recall_scad$specifcity
      acc_scad[s]  <- recall_scad$accuracy
      G_scad[s]    <- recall_scad$Gscore
      
      strong_las[s]<- length(WRobj_las$betaPSE)
      weak_las[s]  <- WRobj_las$S2
      
      strong_alas[s]<- length(WRobj_alas$betaPSE)
      weak_alas[s]  <- WRobj_alas$S2
      
      strong_scad[s]<- length(WRobj_scad$betaPSE)
      weak_scad[s]  <- WRobj_scad$S2
      #---------------------------------------------------------------------
      Est.f.LassoWR[,s]  <- WR_fhat_las
      Est.f.aLassoWR[,s] <- WR_fhat_alas
      Est.f.SCADWR[,s]   <- WR_fhat_scad
      
      Est.f.Lasso[,s]  <- fhat_las
      Est.f.aLasso[,s] <- fhat_alas
      Est.f.SCAD[,s]   <- fhat_scad
      
      Est.y.Lasso[,s]  <- yhat_las
      Est.y.aLasso[,s] <- yhat_alas
      Est.y.SCAD[,s]   <- yhat_scad
      
      Est.y.LassoWR[,s]  <- WR_yhat_las
      Est.y.aLassoWR[,s] <- WR_yhat_alas
      Est.y.SCADWR[,s]   <- WR_yhat_scad
      
      #plot(f,type="l",ylim=c(min(f),max(f)))
      #par(new=TRUE)
      #plot(WR_fhat_las,type="l",ylim=c(min(f),max(f)),col="red")
      #par(new=TRUE)
      #plot(WR_fhat_alas,type="l",ylim=c(min(f),max(f)),col="blue")
      #par(new=TRUE)
      #plot(WR_fhat_scad,type="l",ylim=c(min(f),max(f)),col="green")
      aa <- aa+1
      message(s, "th Simulation ends for n=",n, " pn=",pn, " and rho=", rho, " in Scenario ",sc)
    }
    if (s2==1 & s1==1){
      WR.lasso.fhat.p300.n50  <-rowMeans(Est.f.LassoWR) 
      WR.alasso.fhat.p300.n50 <-rowMeans(Est.f.aLassoWR)
      WR.scad.fhat.p300.n50   <-rowMeans(Est.f.SCADWR)
      
      rmsebeta_las.p300.n50  <- mean(rmsebeta_las)
      rmsebeta_alas.p300.n50 <- mean(rmsebeta_alas)
      rmsebeta_scad.p300.n50 <- mean(rmsebeta_scad)
      
      Rsq_las.p300.n50  <- mean(Rsq_las)
      Rsq_alas.p300.n50 <- mean(Rsq_alas)
      Rsq_scad.p300.n50 <- mean(Rsq_scad)
      
      msef_las.p300.n50  <- mean(msef_las)
      msef_alas.p300.n50 <- mean(msef_alas)
      msef_scad.p300.n50 <- mean(msef_scad)
      
      sens_las.p300.n50 <- mean(sens_las)
      spec_las.p300.n50 <- mean(spec_las) 
      acc_las.p300.n50  <- mean(acc_las)
      G_las.p300.n50    <- mean(G_las)
      
      sens_alas.p300.n50 <- mean(sens_alas)
      spec_alas.p300.n50 <- mean(spec_alas) 
      acc_alas.p300.n50  <- mean(acc_alas)
      G_alas.p300.n50    <- mean(G_alas)
      
      sens_scad.p300.n50 <- mean(sens_scad)
      spec_scad.p300.n50 <- mean(spec_scad) 
      acc_scad.p300.n50  <- mean(acc_scad)
      G_scad.p300.n50    <- mean(G_scad)
      
      stong_signals.p300.n50 <- data.frame((strong_las),(strong_alas),(strong_scad))
      weak_signals.p300.n50 <- data.frame((weak_las),(weak_alas),(weak_scad))
      
      
      RMSE_beta.p300.n50             <- data.frame(rmsebeta_las.p300.n50,rmsebeta_alas.p300.n50,rmsebeta_scad.p300.n50)
      colnames(RMSE_beta.p300.n50)   <- c("Lasso","Adaptive-Lasso","SCAD")
      
      RSQ.p300.n50                   <- data.frame(Rsq_las.p300.n50,Rsq_alas.p300.n50,Rsq_scad.p300.n50)
      colnames(RSQ.p300.n50)         <- c("Lasso","Adaptive-Lasso","SCAD")
      
      MSEF.p300.n50                  <- data.frame(msef_las.p300.n50,msef_alas.p300.n50,msef_scad.p300.n50)
      colnames(MSEF.p300.n50)        <- c("Lasso","Adaptive-Lasso","SCAD")
      
      metrics_all.p300.n50           <- rbind(RMSE_beta.p300.n50,RSQ.p300.n50,MSEF.p300.n50)
      rownames(metrics_all.p300.n50) <- c("REMSE","R-square","MSEf") 
      
      RECALL.las.p300.n50                     <- data.frame(sens_las.p300.n50,spec_las.p300.n50,acc_las.p300.n50,G_las.p300.n50)
      colnames(RECALL.las.p300.n50)           <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.alas.p300.n50                    <- data.frame(sens_alas.p300.n50,spec_alas.p300.n50,acc_alas.p300.n50,G_alas.p300.n50)
      colnames(RECALL.alas.p300.n50)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.scad.p300.n50                    <- data.frame(sens_scad.p300.n50,spec_scad.p300.n50,acc_scad.p300.n50,G_scad.p300.n50)
      colnames(RECALL.scad.p300.n50)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL_all.p300.n50                     <- rbind(RECALL.las.p300.n50,RECALL.alas.p300.n50,RECALL.scad.p300.n50)
      rownames(RECALL_all.p300.n50)           <- c("Lasso","Adaptive-Lasso","SCAD")
      #-PLOTS of fs---------------------------------------------------------
      plot(f,type="l",ylim=c(min(f),max(f)))
      par(new=TRUE)
      plot(WR.lasso.fhat.p300.n50,type="l",ylim=c(min(f),max(f)),col="red")
      par(new=TRUE)
      plot(WR.alasso.fhat.p300.n50,type="l",ylim=c(min(f),max(f)),col="blue")
      par(new=TRUE)
      plot(WR.scad.fhat.p300.n50,type="l",ylim=c(min(f),max(f)),col="green")
    }
    if (s2==2 & s1==1){
      WR.lasso.fhat.p1000.n50  <-rowMeans(Est.f.LassoWR) 
      WR.alasso.fhat.p1000.n50 <-rowMeans(Est.f.aLassoWR)
      WR.scad.fhat.p1000.n50   <-rowMeans(Est.f.SCADWR)
      
      rmsebeta_las.p1000.n50  <- mean(rmsebeta_las)
      rmsebeta_alas.p1000.n50 <- mean(rmsebeta_alas)
      rmsebeta_scad.p1000.n50 <- mean(rmsebeta_scad)
      
      Rsq_las.p1000.n50  <- mean(Rsq_las)
      Rsq_alas.p1000.n50 <- mean(Rsq_alas)
      Rsq_scad.p1000.n50 <- mean(Rsq_scad)
      
      msef_las.p1000.n50  <- mean(msef_las)
      msef_alas.p1000.n50 <- mean(msef_alas)
      msef_scad.p1000.n50 <- mean(msef_scad)
      
      sens_las.p1000.n50 <- mean(sens_las)
      spec_las.p1000.n50 <- mean(spec_las) 
      acc_las.p1000.n50  <- mean(acc_las)
      G_las.p1000.n50    <- mean(G_las)
      
      sens_alas.p1000.n50 <- mean(sens_alas)
      spec_alas.p1000.n50 <- mean(spec_alas) 
      acc_alas.p1000.n50  <- mean(acc_alas)
      G_alas.p1000.n50    <- mean(G_alas)
      
      sens_scad.p1000.n50 <- mean(sens_scad)
      spec_scad.p1000.n50 <- mean(spec_scad) 
      acc_scad.p1000.n50  <- mean(acc_scad)
      G_scad.p1000.n50    <- mean(G_scad)
      
      RMSE_beta.p1000.n50             <- data.frame(rmsebeta_las.p1000.n50,rmsebeta_alas.p1000.n50,rmsebeta_scad.p1000.n50)
      colnames(RMSE_beta.p1000.n50)   <- c("Lasso","Adaptive-Lasso","SCAD")
      
      RSQ.p1000.n50                   <- data.frame(Rsq_las.p1000.n50,Rsq_alas.p1000.n50,Rsq_scad.p1000.n50)
      colnames(RSQ.p1000.n50)         <- c("Lasso","Adaptive-Lasso","SCAD")
      
      MSEF.p1000.n50                  <- data.frame(msef_las.p1000.n50,msef_alas.p1000.n50,msef_scad.p1000.n50)
      colnames(MSEF.p1000.n50)        <- c("Lasso","Adaptive-Lasso","SCAD")
      
      metrics_all.p1000.n50           <- rbind(RMSE_beta.p1000.n50,RSQ.p1000.n50,MSEF.p1000.n50)
      rownames(metrics_all.p1000.n50) <- c("REMSE","R-square","MSEf") 
      
      stong_signals.p1000.n50 <- data.frame((strong_las),(strong_alas),(strong_scad))
      weak_signals.p1000.n50 <- data.frame((weak_las),(weak_alas),(weak_scad))
      
      RECALL.las.p1000.n50                     <- data.frame(sens_las.p1000.n50,spec_las.p1000.n50,acc_las.p1000.n50,G_las.p1000.n50)
      colnames(RECALL.las.p1000.n50)           <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.alas.p1000.n50                    <- data.frame(sens_alas.p1000.n50,spec_alas.p1000.n50,acc_alas.p1000.n50,G_alas.p1000.n50)
      colnames(RECALL.alas.p1000.n50)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.scad.p1000.n50                    <- data.frame(sens_scad.p1000.n50,spec_scad.p1000.n50,acc_scad.p1000.n50,G_scad.p1000.n50)
      colnames(RECALL.scad.p1000.n50)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL_all.p1000.n50                     <- rbind(RECALL.las.p1000.n50,RECALL.alas.p1000.n50,RECALL.scad.p1000.n50)
      rownames(RECALL_all.p1000.n50)           <- c("Lasso","Adaptive-Lasso","SCAD")
    }
    if (s2==3 & s1==1){
      WR.lasso.fhat.p3000.n50  <-rowMeans(Est.f.LassoWR) 
      WR.alasso.fhat.p3000.n50 <-rowMeans(Est.f.aLassoWR)
      WR.scad.fhat.p3000.n50   <-rowMeans(Est.f.SCADWR)
      
      rmsebeta_las.p3000.n50  <- mean(rmsebeta_las)
      rmsebeta_alas.p3000.n50 <- mean(rmsebeta_alas)
      rmsebeta_scad.p3000.n50 <- mean(rmsebeta_scad)
      
      Rsq_las.p3000.n50  <- mean(Rsq_las)
      Rsq_alas.p3000.n50 <- mean(Rsq_alas)
      Rsq_scad.p3000.n50 <- mean(Rsq_scad)
      
      msef_las.p3000.n50  <- mean(msef_las)
      msef_alas.p3000.n50 <- mean(msef_alas)
      msef_scad.p3000.n50 <- mean(msef_scad)
      
      sens_las.p3000.n50 <- mean(sens_las)
      spec_las.p3000.n50 <- mean(spec_las) 
      acc_las.p3000.n50  <- mean(acc_las)
      G_las.p3000.n50    <- mean(G_las)
      
      sens_alas.p3000.n50 <- mean(sens_alas)
      spec_alas.p3000.n50 <- mean(spec_alas) 
      acc_alas.p3000.n50  <- mean(acc_alas)
      G_alas.p3000.n50    <- mean(G_alas)
      
      sens_scad.p3000.n50 <- mean(sens_scad)
      spec_scad.p3000.n50 <- mean(spec_scad) 
      acc_scad.p3000.n50  <- mean(acc_scad)
      G_scad.p3000.n50    <- mean(G_scad)
      
      RMSE_beta.p3000.n50             <- data.frame(rmsebeta_las.p3000.n50,rmsebeta_alas.p3000.n50,rmsebeta_scad.p3000.n50)
      colnames(RMSE_beta.p3000.n50)   <- c("Lasso","Adaptive-Lasso","SCAD")
      
      RSQ.p3000.n50                   <- data.frame(Rsq_las.p3000.n50,Rsq_alas.p3000.n50,Rsq_scad.p3000.n50)
      colnames(RSQ.p3000.n50)         <- c("Lasso","Adaptive-Lasso","SCAD")
      
      MSEF.p3000.n50                  <- data.frame(msef_las.p3000.n50,msef_alas.p3000.n50,msef_scad.p3000.n50)
      colnames(MSEF.p3000.n50)        <- c("Lasso","Adaptive-Lasso","SCAD")
      
      metrics_all.p3000.n50           <- rbind(RMSE_beta.p3000.n50,RSQ.p3000.n50,MSEF.p3000.n50)
      rownames(metrics_all.p3000.n50) <- c("REMSE","R-square","MSEf") 
      
      stong_signals.p3000.n50 <- data.frame((strong_las),(strong_alas),(strong_scad))
      weak_signals.p3000.n50 <- data.frame((weak_las),(weak_alas),(weak_scad))
      
      RECALL.las.p3000.n50                     <- data.frame(sens_las.p3000.n50,spec_las.p3000.n50,acc_las.p3000.n50,G_las.p3000.n50)
      colnames(RECALL.las.p3000.n50)           <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.alas.p3000.n50                    <- data.frame(sens_alas.p3000.n50,spec_alas.p3000.n50,acc_alas.p3000.n50,G_alas.p3000.n50)
      colnames(RECALL.alas.p3000.n50)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.scad.p3000.n50                    <- data.frame(sens_scad.p3000.n50,spec_scad.p3000.n50,acc_scad.p3000.n50,G_scad.p3000.n50)
      colnames(RECALL.scad.p3000.n50)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL_all.p3000.n50                     <- rbind(RECALL.las.p3000.n50,RECALL.alas.p3000.n50,RECALL.scad.p3000.n50)
      rownames(RECALL_all.p3000.n50)           <- c("Lasso","Adaptive-Lasso","SCAD")
    }
    if (s2==1 & s1==2){
      WR.lasso.fhat.p300.n100  <-rowMeans(Est.f.LassoWR) 
      WR.alasso.fhat.p300.n100 <-rowMeans(Est.f.aLassoWR)
      WR.scad.fhat.p300.n100   <-rowMeans(Est.f.SCADWR)
      
      rmsebeta_las.p300.n100   <- mean(rmsebeta_las)
      rmsebeta_alas.p300.n100 <- mean(rmsebeta_alas)
      rmsebeta_scad.p300.n100 <- mean(rmsebeta_scad)
      
      Rsq_las.p300.n100  <- mean(Rsq_las)
      Rsq_alas.p300.n100 <- mean(Rsq_alas)
      Rsq_scad.p300.n100 <- mean(Rsq_scad)
      
      msef_las.p300.n100  <- mean(msef_las)
      msef_alas.p300.n100 <- mean(msef_alas)
      msef_scad.p300.n100 <- mean(msef_scad)
      
      sens_las.p300.n100 <- mean(sens_las)
      spec_las.p300.n100 <- mean(spec_las) 
      acc_las.p300.n100  <- mean(acc_las)
      G_las.p300.n100    <- mean(G_las)
      
      sens_alas.p300.n100 <- mean(sens_alas)
      spec_alas.p300.n100 <- mean(spec_alas) 
      acc_alas.p300.n100  <- mean(acc_alas)
      G_alas.p300.n100    <- mean(G_alas)
      
      sens_scad.p300.n100 <- mean(sens_scad)
      spec_scad.p300.n100 <- mean(spec_scad) 
      acc_scad.p300.n100  <- mean(acc_scad)
      G_scad.p300.n100    <- mean(G_scad)
      
      RMSE_beta.p300.n100             <- data.frame(rmsebeta_las.p300.n100,rmsebeta_alas.p300.n100,rmsebeta_scad.p300.n100)
      colnames(RMSE_beta.p300.n100)   <- c("Lasso","Adaptive-Lasso","SCAD")
      
      RSQ.p300.n100                   <- data.frame(Rsq_las.p300.n100,Rsq_alas.p300.n100,Rsq_scad.p300.n100)
      colnames(RSQ.p300.n100)         <- c("Lasso","Adaptive-Lasso","SCAD")
      
      MSEF.p300.n100                  <- data.frame(msef_las.p300.n100,msef_alas.p300.n100,msef_scad.p300.n100)
      colnames(MSEF.p300.n100)        <- c("Lasso","Adaptive-Lasso","SCAD")
      
      metrics_all.p300.n100           <- rbind(RMSE_beta.p300.n100,RSQ.p300.n100,MSEF.p300.n100)
      rownames(metrics_all.p300.n100) <- c("REMSE","R-square","MSEf") 
      
      stong_signals.p300.n100 <- data.frame((strong_las),(strong_alas),(strong_scad))
      weak_signals.p300.n100 <- data.frame((weak_las),(weak_alas),(weak_scad))
      
      RECALL.las.p300.n100                     <- data.frame(sens_las.p300.n100,spec_las.p300.n100,acc_las.p300.n100,G_las.p300.n100)
      colnames(RECALL.las.p300.n100)           <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.alas.p300.n100                    <- data.frame(sens_alas.p300.n100,spec_alas.p300.n100,acc_alas.p300.n100,G_alas.p300.n100)
      colnames(RECALL.alas.p300.n100)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.scad.p300.n100                    <- data.frame(sens_scad.p300.n100,spec_scad.p300.n100,acc_scad.p300.n100,G_scad.p300.n100)
      colnames(RECALL.scad.p300.n100)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL_all.p300.n100                     <- rbind(RECALL.las.p300.n100,RECALL.alas.p300.n100,RECALL.scad.p300.n100)
      rownames(RECALL_all.p300.n100)           <- c("Lasso","Adaptive-Lasso","SCAD")
      
    }
    if (s2==2 & s1==2){
      WR.lasso.fhat.p1000.n100  <-rowMeans(Est.f.LassoWR) 
      WR.alasso.fhat.p1000.n100 <-rowMeans(Est.f.aLassoWR)
      WR.scad.fhat.p1000.n100   <-rowMeans(Est.f.SCADWR)
      
      rmsebeta_las.p1000.n100   <- mean(rmsebeta_las)
      rmsebeta_alas.p1000.n100  <- mean(rmsebeta_alas)
      rmsebeta_scad.p1000.n100  <- mean(rmsebeta_scad)
      
      Rsq_las.p1000.n100  <- mean(Rsq_las)
      Rsq_alas.p1000.n100 <- mean(Rsq_alas)
      Rsq_scad.p1000.n100 <- mean(Rsq_scad)
      
      msef_las.p1000.n100  <- mean(msef_las)
      msef_alas.p1000.n100 <- mean(msef_alas)
      msef_scad.p1000.n100 <- mean(msef_scad)
      
      sens_las.p1000.n100 <- mean(sens_las)
      spec_las.p1000.n100 <- mean(spec_las) 
      acc_las.p1000.n100  <- mean(acc_las)
      G_las.p1000.n100    <- mean(G_las)
      
      sens_alas.p1000.n100 <- mean(sens_alas)
      spec_alas.p1000.n100 <- mean(spec_alas) 
      acc_alas.p1000.n100  <- mean(acc_alas)
      G_alas.p1000.n100    <- mean(G_alas)
      
      sens_scad.p1000.n100 <- mean(sens_scad)
      spec_scad.p1000.n100 <- mean(spec_scad) 
      acc_scad.p1000.n100  <- mean(acc_scad)
      G_scad.p1000.n100    <- mean(G_scad)
      
      RMSE_beta.p1000.n100             <- data.frame(rmsebeta_las.p1000.n100,rmsebeta_alas.p1000.n100,rmsebeta_scad.p1000.n100)
      colnames(RMSE_beta.p1000.n100)   <- c("Lasso","Adaptive-Lasso","SCAD")
      
      RSQ.p1000.n100                   <- data.frame(Rsq_las.p1000.n100,Rsq_alas.p1000.n100,Rsq_scad.p1000.n100)
      colnames(RSQ.p1000.n100)         <- c("Lasso","Adaptive-Lasso","SCAD")
      
      MSEF.p1000.n100                  <- data.frame(msef_las.p1000.n100,msef_alas.p1000.n100,msef_scad.p1000.n100)
      colnames(MSEF.p1000.n100)        <- c("Lasso","Adaptive-Lasso","SCAD")
      
      metrics_all.p1000.n100           <- rbind(RMSE_beta.p1000.n100,RSQ.p1000.n100,MSEF.p1000.n100)
      rownames(metrics_all.p1000.n100) <- c("REMSE","R-square","MSEf") 
      
      stong_signals.p1000.n100 <- data.frame((strong_las),(strong_alas),(strong_scad))
      weak_signals.p1000.n100 <- data.frame((weak_las),(weak_alas),(weak_scad))
      
      RECALL.las.p1000.n100                     <- data.frame(sens_las.p1000.n100,spec_las.p1000.n100,acc_las.p1000.n100,G_las.p1000.n100)
      colnames(RECALL.las.p1000.n100)           <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.alas.p1000.n100                    <- data.frame(sens_alas.p1000.n100,spec_alas.p1000.n100,acc_alas.p1000.n100,G_alas.p1000.n100)
      colnames(RECALL.alas.p1000.n100)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.scad.p1000.n100                    <- data.frame(sens_scad.p1000.n100,spec_scad.p1000.n100,acc_scad.p1000.n100,G_scad.p1000.n100)
      colnames(RECALL.scad.p1000.n100)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL_all.p1000.n100                     <- rbind(RECALL.las.p1000.n100,RECALL.alas.p1000.n100,RECALL.scad.p1000.n100)
      rownames(RECALL_all.p1000.n100)           <- c("Lasso","Adaptive-Lasso","SCAD")
      
    }
    if (s2==3 & s1==2){
      WR.lasso.fhat.p3000.n100  <-rowMeans(Est.f.LassoWR) 
      WR.alasso.fhat.p3000.n100 <-rowMeans(Est.f.aLassoWR)
      WR.scad.fhat.p3000.n100   <-rowMeans(Est.f.SCADWR)
      
      rmsebeta_las.p3000.n100   <- mean(rmsebeta_las)
      rmsebeta_alas.p3000.n100  <- mean(rmsebeta_alas)
      rmsebeta_scad.p3000.n100  <- mean(rmsebeta_scad)
      
      Rsq_las.p3000.n100  <- mean(Rsq_las)
      Rsq_alas.p3000.n100 <- mean(Rsq_alas)
      Rsq_scad.p3000.n100 <- mean(Rsq_scad)
      
      msef_las.p3000.n100  <- mean(msef_las)
      msef_alas.p3000.n100 <- mean(msef_alas)
      msef_scad.p3000.n100 <- mean(msef_scad)
      
      sens_las.p3000.n100 <- mean(sens_las)
      spec_las.p3000.n100 <- mean(spec_las) 
      acc_las.p3000.n100  <- mean(acc_las)
      G_las.p3000.n100    <- mean(G_las)
      
      sens_alas.p3000.n100 <- mean(sens_alas)
      spec_alas.p3000.n100 <- mean(spec_alas) 
      acc_alas.p3000.n100  <- mean(acc_alas)
      G_alas.p3000.n100    <- mean(G_alas)
      
      sens_scad.p3000.n100 <- mean(sens_scad)
      spec_scad.p3000.n100 <- mean(spec_scad) 
      acc_scad.p3000.n100  <- mean(acc_scad)
      G_scad.p3000.n100    <- mean(G_scad)
      
      RMSE_beta.p3000.n100             <- data.frame(rmsebeta_las.p3000.n100,rmsebeta_alas.p3000.n100,rmsebeta_scad.p3000.n100)
      colnames(RMSE_beta.p3000.n100)   <- c("Lasso","Adaptive-Lasso","SCAD")
      
      RSQ.p3000.n100                   <- data.frame(Rsq_las.p3000.n100,Rsq_alas.p3000.n100,Rsq_scad.p3000.n100)
      colnames(RSQ.p3000.n100)         <- c("Lasso","Adaptive-Lasso","SCAD")
      
      MSEF.p3000.n100                  <- data.frame(msef_las.p3000.n100,msef_alas.p3000.n100,msef_scad.p3000.n100)
      colnames(MSEF.p3000.n100)        <- c("Lasso","Adaptive-Lasso","SCAD")
      
      metrics_all.p3000.n100           <- rbind(RMSE_beta.p3000.n100,RSQ.p3000.n100,MSEF.p3000.n100)
      rownames(metrics_all.p3000.n100) <- c("REMSE","R-square","MSEf") 
      
      stong_signals.p3000.n100 <- data.frame((strong_las),(strong_alas),(strong_scad))
      weak_signals.p3000.n100 <- data.frame((weak_las),(weak_alas),(weak_scad))
      
      RECALL.las.p3000.n100                     <- data.frame(sens_las.p3000.n100,spec_las.p3000.n100,acc_las.p3000.n100,G_las.p3000.n100)
      colnames(RECALL.las.p3000.n100)           <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.alas.p3000.n100                    <- data.frame(sens_alas.p3000.n100,spec_alas.p3000.n100,acc_alas.p3000.n100,G_alas.p3000.n100)
      colnames(RECALL.alas.p3000.n100)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.scad.p3000.n100                    <- data.frame(sens_scad.p3000.n100,spec_scad.p3000.n100,acc_scad.p3000.n100,G_scad.p3000.n100)
      colnames(RECALL.scad.p3000.n100)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL_all.p3000.n100                     <- rbind(RECALL.las.p3000.n100,RECALL.alas.p3000.n100,RECALL.scad.p3000.n100)
      rownames(RECALL_all.p3000.n100)           <- c("Lasso","Adaptive-Lasso","SCAD")
    }
    if (s2==1 & s1==3){
      WR.lasso.fhat.p300.n200  <-rowMeans(Est.f.LassoWR) 
      WR.alasso.fhat.p300.n200 <-rowMeans(Est.f.aLassoWR)
      WR.scad.fhat.p300.n200   <-rowMeans(Est.f.SCADWR)
      
      rmsebeta_las.p300.n200   <- mean(rmsebeta_las)
      rmsebeta_alas.p300.n200 <- mean(rmsebeta_alas)
      rmsebeta_scad.p300.n200 <- mean(rmsebeta_scad)
      
      Rsq_las.p300.n200  <- mean(Rsq_las)
      Rsq_alas.p300.n200 <- mean(Rsq_alas)
      Rsq_scad.p300.n200 <- mean(Rsq_scad)
      
      msef_las.p300.n200  <- mean(msef_las)
      msef_alas.p300.n200 <- mean(msef_alas)
      msef_scad.p300.n200 <- mean(msef_scad)
      
      sens_las.p300.n200 <- mean(sens_las)
      spec_las.p300.n200 <- mean(spec_las) 
      acc_las.p300.n200  <- mean(acc_las)
      G_las.p300.n200    <- mean(G_las)
      
      sens_alas.p300.n200 <- mean(sens_alas)
      spec_alas.p300.n200 <- mean(spec_alas) 
      acc_alas.p300.n200  <- mean(acc_alas)
      G_alas.p300.n200    <- mean(G_alas)
      
      sens_scad.p300.n200 <- mean(sens_scad)
      spec_scad.p300.n200 <- mean(spec_scad) 
      acc_scad.p300.n200  <- mean(acc_scad)
      G_scad.p300.n200    <- mean(G_scad)
      
      RMSE_beta.p300.n200             <- data.frame(rmsebeta_las.p300.n200,rmsebeta_alas.p300.n200,rmsebeta_scad.p300.n200)
      colnames(RMSE_beta.p300.n200)   <- c("Lasso","Adaptive-Lasso","SCAD")
      
      RSQ.p300.n200                   <- data.frame(Rsq_las.p300.n200,Rsq_alas.p300.n200,Rsq_scad.p300.n200)
      colnames(RSQ.p300.n200)         <- c("Lasso","Adaptive-Lasso","SCAD")
      
      MSEF.p300.n200                  <- data.frame(msef_las.p300.n200,msef_alas.p300.n200,msef_scad.p300.n200)
      colnames(MSEF.p300.n200)        <- c("Lasso","Adaptive-Lasso","SCAD")
      
      metrics_all.p300.n200           <- rbind(RMSE_beta.p300.n200,RSQ.p300.n200,MSEF.p300.n200)
      rownames(metrics_all.p300.n200) <- c("REMSE","R-square","MSEf") 
      
      stong_signals.p300.n200 <- data.frame((strong_las),(strong_alas),(strong_scad))
      weak_signals.p300.n200  <- data.frame((weak_las),(weak_alas),(weak_scad))
      
      RECALL.las.p300.n200                     <- data.frame(sens_las.p300.n200,spec_las.p300.n200,acc_las.p300.n200,G_las.p300.n200)
      colnames(RECALL.las.p300.n200)           <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.alas.p300.n200                    <- data.frame(sens_alas.p300.n200,spec_alas.p300.n200,acc_alas.p300.n200,G_alas.p300.n200)
      colnames(RECALL.alas.p300.n200)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.scad.p300.n200                    <- data.frame(sens_scad.p300.n200,spec_scad.p300.n200,acc_scad.p300.n200,G_scad.p300.n200)
      colnames(RECALL.scad.p300.n200)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL_all.p300.n200                     <- rbind(RECALL.las.p300.n200,RECALL.alas.p300.n200,RECALL.scad.p300.n200)
      rownames(RECALL_all.p300.n200)           <- c("Lasso","Adaptive-Lasso","SCAD")
      
    }
    if (s2==2 & s1==3){
      WR.lasso.fhat.p1000.n200  <-rowMeans(Est.f.LassoWR) 
      WR.alasso.fhat.p1000.n200 <-rowMeans(Est.f.aLassoWR)
      WR.scad.fhat.p1000.n200   <-rowMeans(Est.f.SCADWR)
      
      rmsebeta_las.p1000.n200   <- mean(rmsebeta_las)
      rmsebeta_alas.p1000.n200 <- mean(rmsebeta_alas)
      rmsebeta_scad.p1000.n200 <- mean(rmsebeta_scad)
      
      Rsq_las.p1000.n200  <- mean(Rsq_las)
      Rsq_alas.p1000.n200 <- mean(Rsq_alas)
      Rsq_scad.p1000.n200 <- mean(Rsq_scad)
      
      msef_las.p1000.n200  <- mean(msef_las)
      msef_alas.p1000.n200 <- mean(msef_alas)
      msef_scad.p1000.n200 <- mean(msef_scad)
      
      sens_las.p1000.n200 <- mean(sens_las)
      spec_las.p1000.n200 <- mean(spec_las) 
      acc_las.p1000.n200  <- mean(acc_las)
      G_las.p1000.n200    <- mean(G_las)
      
      sens_alas.p1000.n200 <- mean(sens_alas)
      spec_alas.p1000.n200 <- mean(spec_alas) 
      acc_alas.p1000.n200  <- mean(acc_alas)
      G_alas.p1000.n200    <- mean(G_alas)
      
      sens_scad.p1000.n200 <- mean(sens_scad)
      spec_scad.p1000.n200 <- mean(spec_scad) 
      acc_scad.p1000.n200  <- mean(acc_scad)
      G_scad.p1000.n200    <- mean(G_scad)
      
      stong_signals.p1000.n200 <- data.frame((strong_las),(strong_alas),(strong_scad))
      weak_signals.p1000.n200  <- data.frame((weak_las),(weak_alas),(weak_scad))
      
      RMSE_beta.p1000.n200             <- data.frame(rmsebeta_las.p1000.n200,rmsebeta_alas.p1000.n200,rmsebeta_scad.p1000.n200)
      colnames(RMSE_beta.p1000.n200)   <- c("Lasso","Adaptive-Lasso","SCAD")
      
      RSQ.p1000.n200                   <- data.frame(Rsq_las.p1000.n200,Rsq_alas.p1000.n200,Rsq_scad.p1000.n200)
      colnames(RSQ.p1000.n200)         <- c("Lasso","Adaptive-Lasso","SCAD")
      
      MSEF.p1000.n200                  <- data.frame(msef_las.p1000.n200,msef_alas.p1000.n200,msef_scad.p1000.n200)
      colnames(MSEF.p1000.n200)        <- c("Lasso","Adaptive-Lasso","SCAD")
      
      metrics_all.p1000.n200           <- rbind(RMSE_beta.p1000.n200,RSQ.p1000.n200,MSEF.p1000.n200)
      rownames(metrics_all.p1000.n200) <- c("REMSE","R-square","MSEf") 
      
      RECALL.las.p1000.n200                     <- data.frame(sens_las.p1000.n200,spec_las.p1000.n200,acc_las.p1000.n200,G_las.p1000.n200)
      colnames(RECALL.las.p1000.n200)           <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.alas.p1000.n200                    <- data.frame(sens_alas.p1000.n200,spec_alas.p1000.n200,acc_alas.p1000.n200,G_alas.p1000.n200)
      colnames(RECALL.alas.p1000.n200)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.scad.p1000.n200                    <- data.frame(sens_scad.p1000.n200,spec_scad.p1000.n200,acc_scad.p1000.n200,G_scad.p1000.n200)
      colnames(RECALL.scad.p1000.n200)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL_all.p1000.n200                     <- rbind(RECALL.las.p1000.n200,RECALL.alas.p1000.n200,RECALL.scad.p1000.n200)
      rownames(RECALL_all.p1000.n200)           <- c("Lasso","Adaptive-Lasso","SCAD")
      
    }
    if (s2==3 & s1==3){
      WR.lasso.fhat.p3000.n200  <-rowMeans(Est.f.LassoWR) 
      WR.alasso.fhat.p3000.n200 <-rowMeans(Est.f.aLassoWR)
      WR.scad.fhat.p3000.n200   <-rowMeans(Est.f.SCADWR)
      
      rmsebeta_las.p3000.n200   <- mean(rmsebeta_las)
      rmsebeta_alas.p3000.n200 <- mean(rmsebeta_alas)
      rmsebeta_scad.p3000.n200 <- mean(rmsebeta_scad)
      
      Rsq_las.p3000.n200  <- mean(Rsq_las)
      Rsq_alas.p3000.n200 <- mean(Rsq_alas)
      Rsq_scad.p3000.n200 <- mean(Rsq_scad)
      
      msef_las.p3000.n200  <- mean(msef_las)
      msef_alas.p3000.n200 <- mean(msef_alas)
      msef_scad.p3000.n200 <- mean(msef_scad)
      
      sens_las.p3000.n200 <- mean(sens_las)
      spec_las.p3000.n200 <- mean(spec_las) 
      acc_las.p3000.n200  <- mean(acc_las)
      G_las.p3000.n200    <- mean(G_las)
      
      sens_alas.p3000.n200 <- mean(sens_alas)
      spec_alas.p3000.n200 <- mean(spec_alas) 
      acc_alas.p3000.n200  <- mean(acc_alas)
      G_alas.p3000.n200    <- mean(G_alas)
      
      sens_scad.p3000.n200 <- mean(sens_scad)
      spec_scad.p3000.n200 <- mean(spec_scad) 
      acc_scad.p3000.n200  <- mean(acc_scad)
      G_scad.p3000.n200    <- mean(G_scad)
      
      RMSE_beta.p3000.n200             <- data.frame(rmsebeta_las.p3000.n200,rmsebeta_alas.p3000.n200,rmsebeta_scad.p3000.n200)
      colnames(RMSE_beta.p3000.n200)   <- c("Lasso","Adaptive-Lasso","SCAD")
      
      RSQ.p3000.n200                   <- data.frame(Rsq_las.p3000.n200,Rsq_alas.p3000.n200,Rsq_scad.p3000.n200)
      colnames(RSQ.p3000.n200)         <- c("Lasso","Adaptive-Lasso","SCAD")
      
      MSEF.p3000.n200                  <- data.frame(msef_las.p3000.n200,msef_alas.p3000.n200,msef_scad.p3000.n200)
      colnames(MSEF.p3000.n200)        <- c("Lasso","Adaptive-Lasso","SCAD")
      
      metrics_all.p3000.n200           <- rbind(RMSE_beta.p3000.n200,RSQ.p3000.n200,MSEF.p3000.n200)
      rownames(metrics_all.p3000.n200) <- c("REMSE","R-square","MSEf") 
      
      stong_signals.p3000.n200 <- data.frame((strong_las),(strong_alas),(strong_scad))
      weak_signals.p3000.n200  <- data.frame((weak_las),(weak_alas),(weak_scad))
      
      RECALL.las.p3000.n200                     <- data.frame(sens_las.p3000.n200,spec_las.p3000.n200,acc_las.p3000.n200,G_las.p3000.n200)
      colnames(RECALL.las.p3000.n200)           <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.alas.p3000.n200                    <- data.frame(sens_alas.p3000.n200,spec_alas.p3000.n200,acc_alas.p3000.n200,G_alas.p3000.n200)
      colnames(RECALL.alas.p3000.n200)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.scad.p3000.n200                    <- data.frame(sens_scad.p3000.n200,spec_scad.p3000.n200,acc_scad.p3000.n200,G_scad.p3000.n200)
      colnames(RECALL.scad.p3000.n200)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL_all.p3000.n200                     <- rbind(RECALL.las.p3000.n200,RECALL.alas.p3000.n200,RECALL.scad.p3000.n200)
      rownames(RECALL_all.p3000.n200)           <- c("Lasso","Adaptive-Lasso","SCAD")
      
    }
  }
}
df_list.p300  <- rbind(metrics_all.p300.n50,metrics_all.p300.n100,metrics_all.p300.n200)
df_list.p1000 <- rbind(metrics_all.p1000.n50,metrics_all.p1000.n100,metrics_all.p1000.n200)
df_list.p3000 <- rbind(metrics_all.p3000.n50,metrics_all.p3000.n100,metrics_all.p3000.n200)


df_list.p300
df_list.p1000 
df_list.p3000

strong.signals.p300 <-cbind(stong_signals.p300.n50,stong_signals.p300.n100,stong_signals.p300.n200) 
weak.signals.p300   <-cbind(weak_signals.p300.n50,weak_signals.p300.n100,weak_signals.p300.n200)

strong.signals.p1000 <-cbind(stong_signals.p1000.n50,stong_signals.p1000.n100,stong_signals.p1000.n200) 
weak.signals.p1000   <-cbind(weak_signals.p1000.n50,weak_signals.p1000.n100,weak_signals.p1000.n200)

strong.signals.p3000 <-cbind(stong_signals.p300.n50,stong_signals.p300.n100,stong_signals.p300.n200) 
weak.signals.p3000   <-cbind(weak_signals.p300.n50,weak_signals.p300.n100,weak_signals.p300.n200)

par(mfrow=c(3,1))
boxplot(strong.signals.p300,main="Strong signals for pn=300")
boxplot(strong.signals.p1000,main="Strong signals for pn=1000")
boxplot(strong.signals.p3000,main="Strong signals for pn=3000")

par(mfrow=c(3,1))
boxplot(weak.signals.p300,main="Weak signals for pn=300")
boxplot(weak.signals.p1000,main="Weak signals for pn=1000")
boxplot(weak.signals.p3000,main="Weak signals for pn=3000")

recall.p300 <- data.frame(sens_las.p300.n50,sens_alas.p300.n50,sens_scad.p300.n50,acc_las.p300.n50,acc_alas.p300.n50,acc_scad.p300.n50)
recall.p1000 <- data.frame(sens_las.p1000.n50,sens_alas.p1000.n50,sens_scad.p1000.n50,acc_las.p1000.n50,acc_alas.p1000.n50,acc_scad.p1000.n50)
recall.p3000 <- data.frame(sens_las.p3000.n50,sens_alas.p3000.n50,sens_scad.p3000.n50,acc_las.p3000.n50,acc_alas.p3000.n50,acc_scad.p3000.n50)

recall.p300
recall.p1000
recall.p3000
###################################################################################################
plot(f,type="l",ylim=c(min(f),max(f)),col=1,main="n=200, pn=3000, Scenario 1, rho=0.80",xlab="t",ylab="f(t)")
par(new=TRUE)
plot(WR.lasso.fhat.p3000.n200,type="l",ylim=c(min(f),max(f)),col=2,xlab="t",ylab="f(t)")
par(new=TRUE)
plot(WR.alasso.fhat.p3000.n200+0.1,type="l",ylim=c(min(f),max(f)),col=3,xlab="t",ylab="f(t)")
par(new=TRUE)
plot(WR.scad.fhat.p3000.n200,type="l",ylim=c(min(f),max(f)),col=4,xlab="t",ylab="f(t)")
legend("bottomleft",legend=c("Real f", "Lasso-PSE","aLasso-PSE","SCAD-PSE"),
       col=c(1, 2,3,4),lty=c(1,1,1,1))
grid()

write.csv(df_list.p300,"p300sc1rho2.csv")
write.csv(df_list.p1000,"p800sc1rho2.csv")
write.csv(df_list.p3000,"p1500sc1rho2.csv")

write.csv(recall.p300,"RECALL.p300sc1rho2.csv")
write.csv(recall.p1000,"RECALL.p800sc1rho2.csv")
write.csv(recall.p3000,"RECALL.p1500sc1rho2.csv")


############################################################################################################
#===========================================================================================================
############################################################################################################
#SCENARIO 2 RHO 0.95

#OUTLINE-----------------------------------------------
#Function 1: Data generation for simulation Sudy
#Function 2: Estimators based on Penalty functions
#Function 3: Calculation of evaluation metrics
#------------------------------------------------------
#Requaired packages------------------------------------
#--------------------------------------------------------------------------------------------------------------
ni <- c(50,100,120)
pni <- c(300,350,400)
rhoi <- c(0.80,0.95)              #Number of simulation
for (s1 in 1:3){                  #Sample size index
  for (s2 in 1:3){                #pn index
    if (s1==1){n=50} 
    if (s1==2){n=100}  
    if (s1==3){n=200}
    if (s2==1){pn=300} 
    if (s2==2){pn=600} 
    if (s2==3){pn=1200}
    rho <- 0.95
    sc  <- 2
    #ZERO MATRICES----------------------------------------------------------
    rmsebeta_las    <- 0
    rmsebeta_alas   <- 0
    rmsebeta_scad   <- 0
    
    Rsq_las         <- 0
    Rsq_alas        <- 0
    Rsq_scad        <- 0
    
    msef_las         <- 0
    msef_alas        <- 0
    msef_scad        <- 0
    
    sens_las <- 0
    spec_las <- 0
    acc_las  <- 0
    G_las    <- 0
    sens_alas <- 0
    spec_alas <- 0
    acc_alas  <- 0
    G_alas    <- 0
    sens_scad <- 0
    spec_scad <- 0
    acc_scad  <- 0
    G_scad    <- 0
    
    strong_las  <- 0
    weak_las    <- 0
    strong_alas <- 0
    weak_alas   <- 0 
    strong_scad <- 0 
    weak_scad   <- 0
    
    Est.f.LassoWR  <- matrix(0,n,sim)
    Est.f.aLassoWR <- matrix(0,n,sim)
    Est.f.SCADWR   <- matrix(0,n,sim)
    
    Est.f.Lasso  <- matrix(0,n,sim)
    Est.f.aLasso <- matrix(0,n,sim)
    Est.f.SCAD   <- matrix(0,n,sim)
    
    Est.y.LassoWR  <- matrix(0,n,sim)
    Est.y.aLassoWR <- matrix(0,n,sim)
    Est.y.SCADWR   <- matrix(0,n,sim)
    
    Est.y.Lasso  <- matrix(0,n,sim)
    Est.y.aLasso <- matrix(0,n,sim)
    Est.y.SCAD   <- matrix(0,n,sim)
    #------------------------------------------------------------------------------
    no.seed <- 19345
    aa <- 1
    for (s in 1:sim){
      #Set seed number needs to be changed in every repeat----------------------------
      set.seed(no.seed+aa)
      mydata <- simdata(n,pn,sc,rho)
      set.seed(no.seed+s)
      #TEST AREA (To be Deleted)------------------------------------------------------
      x     <- mydata$X
      y     <- mydata$y
      t     <- mydata$t
      f     <- mydata$f 
      nol   <- 25
      aiclas  <- 0
      aicalas <- 0
      aicscad <- 0
      
      las     <- lasso(x,y)
      alas    <- adalasso(x,y)
      scadm   <- scad(x,y)
      
      lams    <- seq(0.1,1.5,length.out=nol)
      for (j in 1:nol){
        aiclas[j]  <- aicfunc(y,las,lams[j],f)
        aicalas[j]  <- aicfunc(y,alas,lams[j],f)
        aicscad[j]  <- aicfunc(y,scadm,lams[j],f)
      }
      #PLOTS for smoothing parameter selection (AICc)---------------------------------
      dflas  <- data.frame(lams,aiclas)
      dfalas <- data.frame(lams,aicalas)
      dfscad <- data.frame(lams,aicscad)
      #ggplot()+geom_line(data=dflas,aes(x=lams,aiclas,col="Lasso_WR"),lwd=1)+geom_line(data=dfalas,aes(x=lams,aicalas,col="Adapt.-Lasso_WR"),lwd=1)+geom_line(data=dfscad,aes(x=lams,aicscad,col="SCAD_WR"),lwd=1)+xlab("lambda_spline")+ylab("Improved_AIC score")+ggtitle("(e) n=100, pn=3000, rho=0.95")+theme(legend.title = element_blank(),legend.position = "bottom")
      #-------------------------------------------------------------------------------
      for (j2 in 1:nol){
        if (aiclas[j2]==min(aiclas)){
          lam_las <-  lams[j2]
        }
        if (aicalas[j2]==min(aicalas)){
          lam_alas <-  lams[j2]
        }
        if (aicscad[j2]==min(aicscad)){
          lam_scad <-  lams[j2]
        }
      }
      S_las  <- Smatrix(t,2,lam_las)
      S_alas <- Smatrix(t,2,lam_alas)
      S_scad <- Smatrix(t,2,lam_scad)
      
      
      ytil_las  <- (diag(n)-S_las)%*%y
      ytil_alas <- (diag(n)-S_alas)%*%y
      ytil_scad <- (diag(n)-S_scad)%*%y
      
      
      xtil_las  <- (diag(n)-S_las)%*%x
      xtil_alas <- (diag(n)-S_alas)%*%x
      xtil_scad <- (diag(n)-S_scad)%*%x
      
      newX_las  <- (diag(n)-S_las)%*%las$X_S1
      newX_alas <- (diag(n)-S_alas)%*%alas$X_S1
      newX_scad <- (diag(n)-S_scad)%*%scadm$X_S1
      #Estimation based on Penalty functions
      las_beta      <- solve(t(newX_las)%*%newX_las+0.01*diag(las$S1))%*%t(newX_las)%*%ytil_las
      fhat_las      <- S_las%*%(y-las$X_S1%*%las_beta)
      alas_beta     <- solve(t(newX_alas)%*%newX_alas+0.01*diag(alas$S1))%*%t(newX_alas)%*%ytil_alas
      fhat_alas     <- S_alas%*%(y-alas$X_S1%*%alas_beta)
      scad_beta     <- solve(t(newX_scad)%*%newX_scad+0.01*diag(scadm$S1))%*%t(newX_scad)%*%ytil_scad
      fhat_scad     <- S_scad%*%(y-scadm$X_S1%*%scad_beta)
      
      yhat_las  <- las$X_S1%*%las_beta+fhat_las
      yhat_alas <- alas$X_S1%*%alas_beta+fhat_alas
      yhat_scad <- scadm$X_S1%*%scad_beta+fhat_scad
      
      WRobj_las  <- WR(las,S_las,x,y,0.01)
      WRobj_alas <- WR(alas,S_alas,x,y,0.01)
      WRobj_scad <- WR(scadm,S_scad,x,y,0.01)
      
      WR_fhat_las  <- S_las%*%(y-las$X_S1%*%WRobj_las$betaPSE)
      WR_fhat_alas <- S_alas%*%(y-alas$X_S1%*%WRobj_alas$betaPSE)
      WR_fhat_scad <- S_scad%*%(y-scadm$X_S1%*%WRobj_scad$betaPSE)
      
      WR_yhat_las  <- las$X_S1%*%WRobj_las$betaPSE+fhat_las
      WR_yhat_alas <- alas$X_S1%*%WRobj_alas$betaPSE+fhat_alas
      WR_yhat_scad <- scadm$X_S1%*%WRobj_scad$betaPSE
      #---------------------------------------------------------------------
      las_index       <- las$S1
      real_beta_las   <- mydata$beta[las_index]
      alas_index      <- alas$S1
      real_beta_alas  <- mydata$beta[alas_index]
      scad_index      <- scadm$S1
      real_beta_scad  <- mydata$beta[scad_index]
      
      rmsebeta_las[s]  <- rmsebeta(real_beta_las,WRobj_las$betaPSE)
      rmsebeta_alas[s] <- rmsebeta(real_beta_alas,WRobj_alas$betaPSE)
      rmsebeta_scad[s] <- rmsebeta(real_beta_scad,WRobj_scad$betaPSE)
      
      Rsq_las[s]  <- rsq(mydata$y,WR_yhat_las,length(las_index))
      Rsq_alas[s] <- rsq(mydata$y,WR_yhat_alas,length(alas_index))
      Rsq_scad[s] <- rsq(mydata$y,WR_yhat_scad,length(scad_index))
      
      msef_las[s]  <- msef(mydata$f,WR_fhat_las)
      msef_alas[s] <- msef(mydata$f,WR_fhat_alas)
      msef_scad[s] <- msef(mydata$f,WR_fhat_scad)
      
      recall_las  <- recall(real_beta_las,WRobj_las$betaPSE)
      sens_las[s] <- recall_las$sensitivity
      spec_las[s] <- recall_las$specifcity
      acc_las[s]  <- recall_las$accuracy
      G_las[s]    <- recall_las$Gscore
      
      recall_alas  <- recall(real_beta_alas,WRobj_alas$betaPSE)
      sens_alas[s] <- recall_alas$sensitivity
      spec_alas[s] <- recall_alas$specifcity
      acc_alas[s]  <- recall_alas$accuracy
      G_alas[s]    <- recall_alas$Gscore
      
      recall_scad  <- recall(real_beta_scad,WRobj_scad$betaPSE)
      sens_scad[s] <- recall_scad$sensitivity
      spec_scad[s] <- recall_scad$specifcity
      acc_scad[s]  <- recall_scad$accuracy
      G_scad[s]    <- recall_scad$Gscore
      
      strong_las[s]<- length(WRobj_las$betaPSE)
      weak_las[s]  <- WRobj_las$S2
      
      strong_alas[s]<- length(WRobj_alas$betaPSE)
      weak_alas[s]  <- WRobj_alas$S2
      
      strong_scad[s]<- length(WRobj_scad$betaPSE)
      weak_scad[s]  <- WRobj_scad$S2
      #---------------------------------------------------------------------
      Est.f.LassoWR[,s]  <- WR_fhat_las
      Est.f.aLassoWR[,s] <- WR_fhat_alas
      Est.f.SCADWR[,s]   <- WR_fhat_scad
      
      Est.f.Lasso[,s]  <- fhat_las
      Est.f.aLasso[,s] <- fhat_alas
      Est.f.SCAD[,s]   <- fhat_scad
      
      Est.y.Lasso[,s]  <- yhat_las
      Est.y.aLasso[,s] <- yhat_alas
      Est.y.SCAD[,s]   <- yhat_scad
      
      Est.y.LassoWR[,s]  <- WR_yhat_las
      Est.y.aLassoWR[,s] <- WR_yhat_alas
      Est.y.SCADWR[,s]   <- WR_yhat_scad
      
      #plot(f,type="l",ylim=c(min(f),max(f)))
      #par(new=TRUE)
      #plot(WR_fhat_las,type="l",ylim=c(min(f),max(f)),col="red")
      #par(new=TRUE)
      #plot(WR_fhat_alas,type="l",ylim=c(min(f),max(f)),col="blue")
      #par(new=TRUE)
      #plot(WR_fhat_scad,type="l",ylim=c(min(f),max(f)),col="green")
      aa <- aa+1
      message(s, "th Simulation ends for n=",n, " pn=",pn, " and rho=", rho, " in Scenario ",sc)
    }
    if (s2==1 & s1==1){
      WR.lasso.fhat.p300.n50  <-rowMeans(Est.f.LassoWR) 
      WR.alasso.fhat.p300.n50 <-rowMeans(Est.f.aLassoWR)
      WR.scad.fhat.p300.n50   <-rowMeans(Est.f.SCADWR)
      
      rmsebeta_las.p300.n50  <- mean(rmsebeta_las)
      rmsebeta_alas.p300.n50 <- mean(rmsebeta_alas)
      rmsebeta_scad.p300.n50 <- mean(rmsebeta_scad)
      
      Rsq_las.p300.n50  <- mean(Rsq_las)
      Rsq_alas.p300.n50 <- mean(Rsq_alas)
      Rsq_scad.p300.n50 <- mean(Rsq_scad)
      
      msef_las.p300.n50  <- mean(msef_las)
      msef_alas.p300.n50 <- mean(msef_alas)
      msef_scad.p300.n50 <- mean(msef_scad)
      
      sens_las.p300.n50 <- mean(sens_las)
      spec_las.p300.n50 <- mean(spec_las) 
      acc_las.p300.n50  <- mean(acc_las)
      G_las.p300.n50    <- mean(G_las)
      
      sens_alas.p300.n50 <- mean(sens_alas)
      spec_alas.p300.n50 <- mean(spec_alas) 
      acc_alas.p300.n50  <- mean(acc_alas)
      G_alas.p300.n50    <- mean(G_alas)
      
      sens_scad.p300.n50 <- mean(sens_scad)
      spec_scad.p300.n50 <- mean(spec_scad) 
      acc_scad.p300.n50  <- mean(acc_scad)
      G_scad.p300.n50    <- mean(G_scad)
      
      stong_signals.p300.n50 <- data.frame((strong_las),(strong_alas),(strong_scad))
      weak_signals.p300.n50 <- data.frame((weak_las),(weak_alas),(weak_scad))
      
      
      RMSE_beta.p300.n50             <- data.frame(rmsebeta_las.p300.n50,rmsebeta_alas.p300.n50,rmsebeta_scad.p300.n50)
      colnames(RMSE_beta.p300.n50)   <- c("Lasso","Adaptive-Lasso","SCAD")
      
      RSQ.p300.n50                   <- data.frame(Rsq_las.p300.n50,Rsq_alas.p300.n50,Rsq_scad.p300.n50)
      colnames(RSQ.p300.n50)         <- c("Lasso","Adaptive-Lasso","SCAD")
      
      MSEF.p300.n50                  <- data.frame(msef_las.p300.n50,msef_alas.p300.n50,msef_scad.p300.n50)
      colnames(MSEF.p300.n50)        <- c("Lasso","Adaptive-Lasso","SCAD")
      
      metrics_all.p300.n50           <- rbind(RMSE_beta.p300.n50,RSQ.p300.n50,MSEF.p300.n50)
      rownames(metrics_all.p300.n50) <- c("REMSE","R-square","MSEf") 
      
      RECALL.las.p300.n50                     <- data.frame(sens_las.p300.n50,spec_las.p300.n50,acc_las.p300.n50,G_las.p300.n50)
      colnames(RECALL.las.p300.n50)           <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.alas.p300.n50                    <- data.frame(sens_alas.p300.n50,spec_alas.p300.n50,acc_alas.p300.n50,G_alas.p300.n50)
      colnames(RECALL.alas.p300.n50)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.scad.p300.n50                    <- data.frame(sens_scad.p300.n50,spec_scad.p300.n50,acc_scad.p300.n50,G_scad.p300.n50)
      colnames(RECALL.scad.p300.n50)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL_all.p300.n50                     <- rbind(RECALL.las.p300.n50,RECALL.alas.p300.n50,RECALL.scad.p300.n50)
      rownames(RECALL_all.p300.n50)           <- c("Lasso","Adaptive-Lasso","SCAD")
      #-PLOTS of fs---------------------------------------------------------
      plot(f,type="l",ylim=c(min(f),max(f)))
      par(new=TRUE)
      plot(WR.lasso.fhat.p300.n50,type="l",ylim=c(min(f),max(f)),col="red")
      par(new=TRUE)
      plot(WR.alasso.fhat.p300.n50,type="l",ylim=c(min(f),max(f)),col="blue")
      par(new=TRUE)
      plot(WR.scad.fhat.p300.n50,type="l",ylim=c(min(f),max(f)),col="green")
    }
    if (s2==2 & s1==1){
      WR.lasso.fhat.p1000.n50  <-rowMeans(Est.f.LassoWR) 
      WR.alasso.fhat.p1000.n50 <-rowMeans(Est.f.aLassoWR)
      WR.scad.fhat.p1000.n50   <-rowMeans(Est.f.SCADWR)
      
      rmsebeta_las.p1000.n50  <- mean(rmsebeta_las)
      rmsebeta_alas.p1000.n50 <- mean(rmsebeta_alas)
      rmsebeta_scad.p1000.n50 <- mean(rmsebeta_scad)
      
      Rsq_las.p1000.n50  <- mean(Rsq_las)
      Rsq_alas.p1000.n50 <- mean(Rsq_alas)
      Rsq_scad.p1000.n50 <- mean(Rsq_scad)
      
      msef_las.p1000.n50  <- mean(msef_las)
      msef_alas.p1000.n50 <- mean(msef_alas)
      msef_scad.p1000.n50 <- mean(msef_scad)
      
      sens_las.p1000.n50 <- mean(sens_las)
      spec_las.p1000.n50 <- mean(spec_las) 
      acc_las.p1000.n50  <- mean(acc_las)
      G_las.p1000.n50    <- mean(G_las)
      
      sens_alas.p1000.n50 <- mean(sens_alas)
      spec_alas.p1000.n50 <- mean(spec_alas) 
      acc_alas.p1000.n50  <- mean(acc_alas)
      G_alas.p1000.n50    <- mean(G_alas)
      
      sens_scad.p1000.n50 <- mean(sens_scad)
      spec_scad.p1000.n50 <- mean(spec_scad) 
      acc_scad.p1000.n50  <- mean(acc_scad)
      G_scad.p1000.n50    <- mean(G_scad)
      
      RMSE_beta.p1000.n50             <- data.frame(rmsebeta_las.p1000.n50,rmsebeta_alas.p1000.n50,rmsebeta_scad.p1000.n50)
      colnames(RMSE_beta.p1000.n50)   <- c("Lasso","Adaptive-Lasso","SCAD")
      
      RSQ.p1000.n50                   <- data.frame(Rsq_las.p1000.n50,Rsq_alas.p1000.n50,Rsq_scad.p1000.n50)
      colnames(RSQ.p1000.n50)         <- c("Lasso","Adaptive-Lasso","SCAD")
      
      MSEF.p1000.n50                  <- data.frame(msef_las.p1000.n50,msef_alas.p1000.n50,msef_scad.p1000.n50)
      colnames(MSEF.p1000.n50)        <- c("Lasso","Adaptive-Lasso","SCAD")
      
      metrics_all.p1000.n50           <- rbind(RMSE_beta.p1000.n50,RSQ.p1000.n50,MSEF.p1000.n50)
      rownames(metrics_all.p1000.n50) <- c("REMSE","R-square","MSEf") 
      
      stong_signals.p1000.n50 <- data.frame((strong_las),(strong_alas),(strong_scad))
      weak_signals.p1000.n50 <- data.frame((weak_las),(weak_alas),(weak_scad))
      
      RECALL.las.p1000.n50                     <- data.frame(sens_las.p1000.n50,spec_las.p1000.n50,acc_las.p1000.n50,G_las.p1000.n50)
      colnames(RECALL.las.p1000.n50)           <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.alas.p1000.n50                    <- data.frame(sens_alas.p1000.n50,spec_alas.p1000.n50,acc_alas.p1000.n50,G_alas.p1000.n50)
      colnames(RECALL.alas.p1000.n50)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.scad.p1000.n50                    <- data.frame(sens_scad.p1000.n50,spec_scad.p1000.n50,acc_scad.p1000.n50,G_scad.p1000.n50)
      colnames(RECALL.scad.p1000.n50)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL_all.p1000.n50                     <- rbind(RECALL.las.p1000.n50,RECALL.alas.p1000.n50,RECALL.scad.p1000.n50)
      rownames(RECALL_all.p1000.n50)           <- c("Lasso","Adaptive-Lasso","SCAD")
    }
    if (s2==3 & s1==1){
      WR.lasso.fhat.p3000.n50  <-rowMeans(Est.f.LassoWR) 
      WR.alasso.fhat.p3000.n50 <-rowMeans(Est.f.aLassoWR)
      WR.scad.fhat.p3000.n50   <-rowMeans(Est.f.SCADWR)
      
      rmsebeta_las.p3000.n50  <- mean(rmsebeta_las)
      rmsebeta_alas.p3000.n50 <- mean(rmsebeta_alas)
      rmsebeta_scad.p3000.n50 <- mean(rmsebeta_scad)
      
      Rsq_las.p3000.n50  <- mean(Rsq_las)
      Rsq_alas.p3000.n50 <- mean(Rsq_alas)
      Rsq_scad.p3000.n50 <- mean(Rsq_scad)
      
      msef_las.p3000.n50  <- mean(msef_las)
      msef_alas.p3000.n50 <- mean(msef_alas)
      msef_scad.p3000.n50 <- mean(msef_scad)
      
      sens_las.p3000.n50 <- mean(sens_las)
      spec_las.p3000.n50 <- mean(spec_las) 
      acc_las.p3000.n50  <- mean(acc_las)
      G_las.p3000.n50    <- mean(G_las)
      
      sens_alas.p3000.n50 <- mean(sens_alas)
      spec_alas.p3000.n50 <- mean(spec_alas) 
      acc_alas.p3000.n50  <- mean(acc_alas)
      G_alas.p3000.n50    <- mean(G_alas)
      
      sens_scad.p3000.n50 <- mean(sens_scad)
      spec_scad.p3000.n50 <- mean(spec_scad) 
      acc_scad.p3000.n50  <- mean(acc_scad)
      G_scad.p3000.n50    <- mean(G_scad)
      
      RMSE_beta.p3000.n50             <- data.frame(rmsebeta_las.p3000.n50,rmsebeta_alas.p3000.n50,rmsebeta_scad.p3000.n50)
      colnames(RMSE_beta.p3000.n50)   <- c("Lasso","Adaptive-Lasso","SCAD")
      
      RSQ.p3000.n50                   <- data.frame(Rsq_las.p3000.n50,Rsq_alas.p3000.n50,Rsq_scad.p3000.n50)
      colnames(RSQ.p3000.n50)         <- c("Lasso","Adaptive-Lasso","SCAD")
      
      MSEF.p3000.n50                  <- data.frame(msef_las.p3000.n50,msef_alas.p3000.n50,msef_scad.p3000.n50)
      colnames(MSEF.p3000.n50)        <- c("Lasso","Adaptive-Lasso","SCAD")
      
      metrics_all.p3000.n50           <- rbind(RMSE_beta.p3000.n50,RSQ.p3000.n50,MSEF.p3000.n50)
      rownames(metrics_all.p3000.n50) <- c("REMSE","R-square","MSEf") 
      
      stong_signals.p3000.n50 <- data.frame((strong_las),(strong_alas),(strong_scad))
      weak_signals.p3000.n50 <- data.frame((weak_las),(weak_alas),(weak_scad))
      
      RECALL.las.p3000.n50                     <- data.frame(sens_las.p3000.n50,spec_las.p3000.n50,acc_las.p3000.n50,G_las.p3000.n50)
      colnames(RECALL.las.p3000.n50)           <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.alas.p3000.n50                    <- data.frame(sens_alas.p3000.n50,spec_alas.p3000.n50,acc_alas.p3000.n50,G_alas.p3000.n50)
      colnames(RECALL.alas.p3000.n50)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.scad.p3000.n50                    <- data.frame(sens_scad.p3000.n50,spec_scad.p3000.n50,acc_scad.p3000.n50,G_scad.p3000.n50)
      colnames(RECALL.scad.p3000.n50)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL_all.p3000.n50                     <- rbind(RECALL.las.p3000.n50,RECALL.alas.p3000.n50,RECALL.scad.p3000.n50)
      rownames(RECALL_all.p3000.n50)           <- c("Lasso","Adaptive-Lasso","SCAD")
    }
    if (s2==1 & s1==2){
      WR.lasso.fhat.p300.n100  <-rowMeans(Est.f.LassoWR) 
      WR.alasso.fhat.p300.n100 <-rowMeans(Est.f.aLassoWR)
      WR.scad.fhat.p300.n100   <-rowMeans(Est.f.SCADWR)
      
      rmsebeta_las.p300.n100   <- mean(rmsebeta_las)
      rmsebeta_alas.p300.n100 <- mean(rmsebeta_alas)
      rmsebeta_scad.p300.n100 <- mean(rmsebeta_scad)
      
      Rsq_las.p300.n100  <- mean(Rsq_las)
      Rsq_alas.p300.n100 <- mean(Rsq_alas)
      Rsq_scad.p300.n100 <- mean(Rsq_scad)
      
      msef_las.p300.n100  <- mean(msef_las)
      msef_alas.p300.n100 <- mean(msef_alas)
      msef_scad.p300.n100 <- mean(msef_scad)
      
      sens_las.p300.n100 <- mean(sens_las)
      spec_las.p300.n100 <- mean(spec_las) 
      acc_las.p300.n100  <- mean(acc_las)
      G_las.p300.n100    <- mean(G_las)
      
      sens_alas.p300.n100 <- mean(sens_alas)
      spec_alas.p300.n100 <- mean(spec_alas) 
      acc_alas.p300.n100  <- mean(acc_alas)
      G_alas.p300.n100    <- mean(G_alas)
      
      sens_scad.p300.n100 <- mean(sens_scad)
      spec_scad.p300.n100 <- mean(spec_scad) 
      acc_scad.p300.n100  <- mean(acc_scad)
      G_scad.p300.n100    <- mean(G_scad)
      
      RMSE_beta.p300.n100             <- data.frame(rmsebeta_las.p300.n100,rmsebeta_alas.p300.n100,rmsebeta_scad.p300.n100)
      colnames(RMSE_beta.p300.n100)   <- c("Lasso","Adaptive-Lasso","SCAD")
      
      RSQ.p300.n100                   <- data.frame(Rsq_las.p300.n100,Rsq_alas.p300.n100,Rsq_scad.p300.n100)
      colnames(RSQ.p300.n100)         <- c("Lasso","Adaptive-Lasso","SCAD")
      
      MSEF.p300.n100                  <- data.frame(msef_las.p300.n100,msef_alas.p300.n100,msef_scad.p300.n100)
      colnames(MSEF.p300.n100)        <- c("Lasso","Adaptive-Lasso","SCAD")
      
      metrics_all.p300.n100           <- rbind(RMSE_beta.p300.n100,RSQ.p300.n100,MSEF.p300.n100)
      rownames(metrics_all.p300.n100) <- c("REMSE","R-square","MSEf") 
      
      stong_signals.p300.n100 <- data.frame((strong_las),(strong_alas),(strong_scad))
      weak_signals.p300.n100 <- data.frame((weak_las),(weak_alas),(weak_scad))
      
      RECALL.las.p300.n100                     <- data.frame(sens_las.p300.n100,spec_las.p300.n100,acc_las.p300.n100,G_las.p300.n100)
      colnames(RECALL.las.p300.n100)           <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.alas.p300.n100                    <- data.frame(sens_alas.p300.n100,spec_alas.p300.n100,acc_alas.p300.n100,G_alas.p300.n100)
      colnames(RECALL.alas.p300.n100)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.scad.p300.n100                    <- data.frame(sens_scad.p300.n100,spec_scad.p300.n100,acc_scad.p300.n100,G_scad.p300.n100)
      colnames(RECALL.scad.p300.n100)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL_all.p300.n100                     <- rbind(RECALL.las.p300.n100,RECALL.alas.p300.n100,RECALL.scad.p300.n100)
      rownames(RECALL_all.p300.n100)           <- c("Lasso","Adaptive-Lasso","SCAD")
      
    }
    if (s2==2 & s1==2){
      WR.lasso.fhat.p1000.n100  <-rowMeans(Est.f.LassoWR) 
      WR.alasso.fhat.p1000.n100 <-rowMeans(Est.f.aLassoWR)
      WR.scad.fhat.p1000.n100   <-rowMeans(Est.f.SCADWR)
      
      rmsebeta_las.p1000.n100   <- mean(rmsebeta_las)
      rmsebeta_alas.p1000.n100  <- mean(rmsebeta_alas)
      rmsebeta_scad.p1000.n100  <- mean(rmsebeta_scad)
      
      Rsq_las.p1000.n100  <- mean(Rsq_las)
      Rsq_alas.p1000.n100 <- mean(Rsq_alas)
      Rsq_scad.p1000.n100 <- mean(Rsq_scad)
      
      msef_las.p1000.n100  <- mean(msef_las)
      msef_alas.p1000.n100 <- mean(msef_alas)
      msef_scad.p1000.n100 <- mean(msef_scad)
      
      sens_las.p1000.n100 <- mean(sens_las)
      spec_las.p1000.n100 <- mean(spec_las) 
      acc_las.p1000.n100  <- mean(acc_las)
      G_las.p1000.n100    <- mean(G_las)
      
      sens_alas.p1000.n100 <- mean(sens_alas)
      spec_alas.p1000.n100 <- mean(spec_alas) 
      acc_alas.p1000.n100  <- mean(acc_alas)
      G_alas.p1000.n100    <- mean(G_alas)
      
      sens_scad.p1000.n100 <- mean(sens_scad)
      spec_scad.p1000.n100 <- mean(spec_scad) 
      acc_scad.p1000.n100  <- mean(acc_scad)
      G_scad.p1000.n100    <- mean(G_scad)
      
      RMSE_beta.p1000.n100             <- data.frame(rmsebeta_las.p1000.n100,rmsebeta_alas.p1000.n100,rmsebeta_scad.p1000.n100)
      colnames(RMSE_beta.p1000.n100)   <- c("Lasso","Adaptive-Lasso","SCAD")
      
      RSQ.p1000.n100                   <- data.frame(Rsq_las.p1000.n100,Rsq_alas.p1000.n100,Rsq_scad.p1000.n100)
      colnames(RSQ.p1000.n100)         <- c("Lasso","Adaptive-Lasso","SCAD")
      
      MSEF.p1000.n100                  <- data.frame(msef_las.p1000.n100,msef_alas.p1000.n100,msef_scad.p1000.n100)
      colnames(MSEF.p1000.n100)        <- c("Lasso","Adaptive-Lasso","SCAD")
      
      metrics_all.p1000.n100           <- rbind(RMSE_beta.p1000.n100,RSQ.p1000.n100,MSEF.p1000.n100)
      rownames(metrics_all.p1000.n100) <- c("REMSE","R-square","MSEf") 
      
      stong_signals.p1000.n100 <- data.frame((strong_las),(strong_alas),(strong_scad))
      weak_signals.p1000.n100 <- data.frame((weak_las),(weak_alas),(weak_scad))
      
      RECALL.las.p1000.n100                     <- data.frame(sens_las.p1000.n100,spec_las.p1000.n100,acc_las.p1000.n100,G_las.p1000.n100)
      colnames(RECALL.las.p1000.n100)           <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.alas.p1000.n100                    <- data.frame(sens_alas.p1000.n100,spec_alas.p1000.n100,acc_alas.p1000.n100,G_alas.p1000.n100)
      colnames(RECALL.alas.p1000.n100)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.scad.p1000.n100                    <- data.frame(sens_scad.p1000.n100,spec_scad.p1000.n100,acc_scad.p1000.n100,G_scad.p1000.n100)
      colnames(RECALL.scad.p1000.n100)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL_all.p1000.n100                     <- rbind(RECALL.las.p1000.n100,RECALL.alas.p1000.n100,RECALL.scad.p1000.n100)
      rownames(RECALL_all.p1000.n100)           <- c("Lasso","Adaptive-Lasso","SCAD")
      
    }
    if (s2==3 & s1==2){
      WR.lasso.fhat.p3000.n100  <-rowMeans(Est.f.LassoWR) 
      WR.alasso.fhat.p3000.n100 <-rowMeans(Est.f.aLassoWR)
      WR.scad.fhat.p3000.n100   <-rowMeans(Est.f.SCADWR)
      
      rmsebeta_las.p3000.n100   <- mean(rmsebeta_las)
      rmsebeta_alas.p3000.n100  <- mean(rmsebeta_alas)
      rmsebeta_scad.p3000.n100  <- mean(rmsebeta_scad)
      
      Rsq_las.p3000.n100  <- mean(Rsq_las)
      Rsq_alas.p3000.n100 <- mean(Rsq_alas)
      Rsq_scad.p3000.n100 <- mean(Rsq_scad)
      
      msef_las.p3000.n100  <- mean(msef_las)
      msef_alas.p3000.n100 <- mean(msef_alas)
      msef_scad.p3000.n100 <- mean(msef_scad)
      
      sens_las.p3000.n100 <- mean(sens_las)
      spec_las.p3000.n100 <- mean(spec_las) 
      acc_las.p3000.n100  <- mean(acc_las)
      G_las.p3000.n100    <- mean(G_las)
      
      sens_alas.p3000.n100 <- mean(sens_alas)
      spec_alas.p3000.n100 <- mean(spec_alas) 
      acc_alas.p3000.n100  <- mean(acc_alas)
      G_alas.p3000.n100    <- mean(G_alas)
      
      sens_scad.p3000.n100 <- mean(sens_scad)
      spec_scad.p3000.n100 <- mean(spec_scad) 
      acc_scad.p3000.n100  <- mean(acc_scad)
      G_scad.p3000.n100    <- mean(G_scad)
      
      RMSE_beta.p3000.n100             <- data.frame(rmsebeta_las.p3000.n100,rmsebeta_alas.p3000.n100,rmsebeta_scad.p3000.n100)
      colnames(RMSE_beta.p3000.n100)   <- c("Lasso","Adaptive-Lasso","SCAD")
      
      RSQ.p3000.n100                   <- data.frame(Rsq_las.p3000.n100,Rsq_alas.p3000.n100,Rsq_scad.p3000.n100)
      colnames(RSQ.p3000.n100)         <- c("Lasso","Adaptive-Lasso","SCAD")
      
      MSEF.p3000.n100                  <- data.frame(msef_las.p3000.n100,msef_alas.p3000.n100,msef_scad.p3000.n100)
      colnames(MSEF.p3000.n100)        <- c("Lasso","Adaptive-Lasso","SCAD")
      
      metrics_all.p3000.n100           <- rbind(RMSE_beta.p3000.n100,RSQ.p3000.n100,MSEF.p3000.n100)
      rownames(metrics_all.p3000.n100) <- c("REMSE","R-square","MSEf") 
      
      stong_signals.p3000.n100 <- data.frame((strong_las),(strong_alas),(strong_scad))
      weak_signals.p3000.n100 <- data.frame((weak_las),(weak_alas),(weak_scad))
      
      RECALL.las.p3000.n100                     <- data.frame(sens_las.p3000.n100,spec_las.p3000.n100,acc_las.p3000.n100,G_las.p3000.n100)
      colnames(RECALL.las.p3000.n100)           <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.alas.p3000.n100                    <- data.frame(sens_alas.p3000.n100,spec_alas.p3000.n100,acc_alas.p3000.n100,G_alas.p3000.n100)
      colnames(RECALL.alas.p3000.n100)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.scad.p3000.n100                    <- data.frame(sens_scad.p3000.n100,spec_scad.p3000.n100,acc_scad.p3000.n100,G_scad.p3000.n100)
      colnames(RECALL.scad.p3000.n100)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL_all.p3000.n100                     <- rbind(RECALL.las.p3000.n100,RECALL.alas.p3000.n100,RECALL.scad.p3000.n100)
      rownames(RECALL_all.p3000.n100)           <- c("Lasso","Adaptive-Lasso","SCAD")
    }
    if (s2==1 & s1==3){
      WR.lasso.fhat.p300.n200  <-rowMeans(Est.f.LassoWR) 
      WR.alasso.fhat.p300.n200 <-rowMeans(Est.f.aLassoWR)
      WR.scad.fhat.p300.n200   <-rowMeans(Est.f.SCADWR)
      
      rmsebeta_las.p300.n200   <- mean(rmsebeta_las)
      rmsebeta_alas.p300.n200 <- mean(rmsebeta_alas)
      rmsebeta_scad.p300.n200 <- mean(rmsebeta_scad)
      
      Rsq_las.p300.n200  <- mean(Rsq_las)
      Rsq_alas.p300.n200 <- mean(Rsq_alas)
      Rsq_scad.p300.n200 <- mean(Rsq_scad)
      
      msef_las.p300.n200  <- mean(msef_las)
      msef_alas.p300.n200 <- mean(msef_alas)
      msef_scad.p300.n200 <- mean(msef_scad)
      
      sens_las.p300.n200 <- mean(sens_las)
      spec_las.p300.n200 <- mean(spec_las) 
      acc_las.p300.n200  <- mean(acc_las)
      G_las.p300.n200    <- mean(G_las)
      
      sens_alas.p300.n200 <- mean(sens_alas)
      spec_alas.p300.n200 <- mean(spec_alas) 
      acc_alas.p300.n200  <- mean(acc_alas)
      G_alas.p300.n200    <- mean(G_alas)
      
      sens_scad.p300.n200 <- mean(sens_scad)
      spec_scad.p300.n200 <- mean(spec_scad) 
      acc_scad.p300.n200  <- mean(acc_scad)
      G_scad.p300.n200    <- mean(G_scad)
      
      RMSE_beta.p300.n200             <- data.frame(rmsebeta_las.p300.n200,rmsebeta_alas.p300.n200,rmsebeta_scad.p300.n200)
      colnames(RMSE_beta.p300.n200)   <- c("Lasso","Adaptive-Lasso","SCAD")
      
      RSQ.p300.n200                   <- data.frame(Rsq_las.p300.n200,Rsq_alas.p300.n200,Rsq_scad.p300.n200)
      colnames(RSQ.p300.n200)         <- c("Lasso","Adaptive-Lasso","SCAD")
      
      MSEF.p300.n200                  <- data.frame(msef_las.p300.n200,msef_alas.p300.n200,msef_scad.p300.n200)
      colnames(MSEF.p300.n200)        <- c("Lasso","Adaptive-Lasso","SCAD")
      
      metrics_all.p300.n200           <- rbind(RMSE_beta.p300.n200,RSQ.p300.n200,MSEF.p300.n200)
      rownames(metrics_all.p300.n200) <- c("REMSE","R-square","MSEf") 
      
      stong_signals.p300.n200 <- data.frame((strong_las),(strong_alas),(strong_scad))
      weak_signals.p300.n200  <- data.frame((weak_las),(weak_alas),(weak_scad))
      
      RECALL.las.p300.n200                     <- data.frame(sens_las.p300.n200,spec_las.p300.n200,acc_las.p300.n200,G_las.p300.n200)
      colnames(RECALL.las.p300.n200)           <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.alas.p300.n200                    <- data.frame(sens_alas.p300.n200,spec_alas.p300.n200,acc_alas.p300.n200,G_alas.p300.n200)
      colnames(RECALL.alas.p300.n200)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.scad.p300.n200                    <- data.frame(sens_scad.p300.n200,spec_scad.p300.n200,acc_scad.p300.n200,G_scad.p300.n200)
      colnames(RECALL.scad.p300.n200)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL_all.p300.n200                     <- rbind(RECALL.las.p300.n200,RECALL.alas.p300.n200,RECALL.scad.p300.n200)
      rownames(RECALL_all.p300.n200)           <- c("Lasso","Adaptive-Lasso","SCAD")
      
    }
    if (s2==2 & s1==3){
      WR.lasso.fhat.p1000.n200  <-rowMeans(Est.f.LassoWR) 
      WR.alasso.fhat.p1000.n200 <-rowMeans(Est.f.aLassoWR)
      WR.scad.fhat.p1000.n200   <-rowMeans(Est.f.SCADWR)
      
      rmsebeta_las.p1000.n200   <- mean(rmsebeta_las)
      rmsebeta_alas.p1000.n200 <- mean(rmsebeta_alas)
      rmsebeta_scad.p1000.n200 <- mean(rmsebeta_scad)
      
      Rsq_las.p1000.n200  <- mean(Rsq_las)
      Rsq_alas.p1000.n200 <- mean(Rsq_alas)
      Rsq_scad.p1000.n200 <- mean(Rsq_scad)
      
      msef_las.p1000.n200  <- mean(msef_las)
      msef_alas.p1000.n200 <- mean(msef_alas)
      msef_scad.p1000.n200 <- mean(msef_scad)
      
      sens_las.p1000.n200 <- mean(sens_las)
      spec_las.p1000.n200 <- mean(spec_las) 
      acc_las.p1000.n200  <- mean(acc_las)
      G_las.p1000.n200    <- mean(G_las)
      
      sens_alas.p1000.n200 <- mean(sens_alas)
      spec_alas.p1000.n200 <- mean(spec_alas) 
      acc_alas.p1000.n200  <- mean(acc_alas)
      G_alas.p1000.n200    <- mean(G_alas)
      
      sens_scad.p1000.n200 <- mean(sens_scad)
      spec_scad.p1000.n200 <- mean(spec_scad) 
      acc_scad.p1000.n200  <- mean(acc_scad)
      G_scad.p1000.n200    <- mean(G_scad)
      
      stong_signals.p1000.n200 <- data.frame((strong_las),(strong_alas),(strong_scad))
      weak_signals.p1000.n200  <- data.frame((weak_las),(weak_alas),(weak_scad))
      
      RMSE_beta.p1000.n200             <- data.frame(rmsebeta_las.p1000.n200,rmsebeta_alas.p1000.n200,rmsebeta_scad.p1000.n200)
      colnames(RMSE_beta.p1000.n200)   <- c("Lasso","Adaptive-Lasso","SCAD")
      
      RSQ.p1000.n200                   <- data.frame(Rsq_las.p1000.n200,Rsq_alas.p1000.n200,Rsq_scad.p1000.n200)
      colnames(RSQ.p1000.n200)         <- c("Lasso","Adaptive-Lasso","SCAD")
      
      MSEF.p1000.n200                  <- data.frame(msef_las.p1000.n200,msef_alas.p1000.n200,msef_scad.p1000.n200)
      colnames(MSEF.p1000.n200)        <- c("Lasso","Adaptive-Lasso","SCAD")
      
      metrics_all.p1000.n200           <- rbind(RMSE_beta.p1000.n200,RSQ.p1000.n200,MSEF.p1000.n200)
      rownames(metrics_all.p1000.n200) <- c("REMSE","R-square","MSEf") 
      
      RECALL.las.p1000.n200                     <- data.frame(sens_las.p1000.n200,spec_las.p1000.n200,acc_las.p1000.n200,G_las.p1000.n200)
      colnames(RECALL.las.p1000.n200)           <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.alas.p1000.n200                    <- data.frame(sens_alas.p1000.n200,spec_alas.p1000.n200,acc_alas.p1000.n200,G_alas.p1000.n200)
      colnames(RECALL.alas.p1000.n200)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.scad.p1000.n200                    <- data.frame(sens_scad.p1000.n200,spec_scad.p1000.n200,acc_scad.p1000.n200,G_scad.p1000.n200)
      colnames(RECALL.scad.p1000.n200)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL_all.p1000.n200                     <- rbind(RECALL.las.p1000.n200,RECALL.alas.p1000.n200,RECALL.scad.p1000.n200)
      rownames(RECALL_all.p1000.n200)           <- c("Lasso","Adaptive-Lasso","SCAD")
      
    }
    if (s2==3 & s1==3){
      WR.lasso.fhat.p3000.n200  <-rowMeans(Est.f.LassoWR) 
      WR.alasso.fhat.p3000.n200 <-rowMeans(Est.f.aLassoWR)
      WR.scad.fhat.p3000.n200   <-rowMeans(Est.f.SCADWR)
      
      rmsebeta_las.p3000.n200   <- mean(rmsebeta_las)
      rmsebeta_alas.p3000.n200 <- mean(rmsebeta_alas)
      rmsebeta_scad.p3000.n200 <- mean(rmsebeta_scad)
      
      Rsq_las.p3000.n200  <- mean(Rsq_las)
      Rsq_alas.p3000.n200 <- mean(Rsq_alas)
      Rsq_scad.p3000.n200 <- mean(Rsq_scad)
      
      msef_las.p3000.n200  <- mean(msef_las)
      msef_alas.p3000.n200 <- mean(msef_alas)
      msef_scad.p3000.n200 <- mean(msef_scad)
      
      sens_las.p3000.n200 <- mean(sens_las)
      spec_las.p3000.n200 <- mean(spec_las) 
      acc_las.p3000.n200  <- mean(acc_las)
      G_las.p3000.n200    <- mean(G_las)
      
      sens_alas.p3000.n200 <- mean(sens_alas)
      spec_alas.p3000.n200 <- mean(spec_alas) 
      acc_alas.p3000.n200  <- mean(acc_alas)
      G_alas.p3000.n200    <- mean(G_alas)
      
      sens_scad.p3000.n200 <- mean(sens_scad)
      spec_scad.p3000.n200 <- mean(spec_scad) 
      acc_scad.p3000.n200  <- mean(acc_scad)
      G_scad.p3000.n200    <- mean(G_scad)
      
      RMSE_beta.p3000.n200             <- data.frame(rmsebeta_las.p3000.n200,rmsebeta_alas.p3000.n200,rmsebeta_scad.p3000.n200)
      colnames(RMSE_beta.p3000.n200)   <- c("Lasso","Adaptive-Lasso","SCAD")
      
      RSQ.p3000.n200                   <- data.frame(Rsq_las.p3000.n200,Rsq_alas.p3000.n200,Rsq_scad.p3000.n200)
      colnames(RSQ.p3000.n200)         <- c("Lasso","Adaptive-Lasso","SCAD")
      
      MSEF.p3000.n200                  <- data.frame(msef_las.p3000.n200,msef_alas.p3000.n200,msef_scad.p3000.n200)
      colnames(MSEF.p3000.n200)        <- c("Lasso","Adaptive-Lasso","SCAD")
      
      metrics_all.p3000.n200           <- rbind(RMSE_beta.p3000.n200,RSQ.p3000.n200,MSEF.p3000.n200)
      rownames(metrics_all.p3000.n200) <- c("REMSE","R-square","MSEf") 
      
      stong_signals.p3000.n200 <- data.frame((strong_las),(strong_alas),(strong_scad))
      weak_signals.p3000.n200  <- data.frame((weak_las),(weak_alas),(weak_scad))
      
      RECALL.las.p3000.n200                     <- data.frame(sens_las.p3000.n200,spec_las.p3000.n200,acc_las.p3000.n200,G_las.p3000.n200)
      colnames(RECALL.las.p3000.n200)           <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.alas.p3000.n200                    <- data.frame(sens_alas.p3000.n200,spec_alas.p3000.n200,acc_alas.p3000.n200,G_alas.p3000.n200)
      colnames(RECALL.alas.p3000.n200)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL.scad.p3000.n200                    <- data.frame(sens_scad.p3000.n200,spec_scad.p3000.n200,acc_scad.p3000.n200,G_scad.p3000.n200)
      colnames(RECALL.scad.p3000.n200)          <- c("Sensitivity","Specificity","Accuracy","G-score")
      RECALL_all.p3000.n200                     <- rbind(RECALL.las.p3000.n200,RECALL.alas.p3000.n200,RECALL.scad.p3000.n200)
      rownames(RECALL_all.p3000.n200)           <- c("Lasso","Adaptive-Lasso","SCAD")
      
    }
  }
}
df_list.p300  <- rbind(metrics_all.p300.n50,metrics_all.p300.n100,metrics_all.p300.n200)
df_list.p1000 <- rbind(metrics_all.p1000.n50,metrics_all.p1000.n100,metrics_all.p1000.n200)
df_list.p3000 <- rbind(metrics_all.p3000.n50,metrics_all.p3000.n100,metrics_all.p3000.n200)


df_list.p300
df_list.p1000 
df_list.p3000

strong.signals.p300 <-cbind(stong_signals.p300.n50,stong_signals.p300.n100,stong_signals.p300.n200) 
weak.signals.p300   <-cbind(weak_signals.p300.n50,weak_signals.p300.n100,weak_signals.p300.n200)

strong.signals.p1000 <-cbind(stong_signals.p1000.n50,stong_signals.p1000.n100,stong_signals.p1000.n200) 
weak.signals.p1000   <-cbind(weak_signals.p1000.n50,weak_signals.p1000.n100,weak_signals.p1000.n200)

strong.signals.p3000 <-cbind(stong_signals.p300.n50,stong_signals.p300.n100,stong_signals.p300.n200) 
weak.signals.p3000   <-cbind(weak_signals.p300.n50,weak_signals.p300.n100,weak_signals.p300.n200)

par(mfrow=c(3,1))
boxplot(strong.signals.p300,main="Strong signals for pn=300")
boxplot(strong.signals.p1000,main="Strong signals for pn=1000")
boxplot(strong.signals.p3000,main="Strong signals for pn=3000")

par(mfrow=c(3,1))
boxplot(weak.signals.p300,main="Weak signals for pn=300")
boxplot(weak.signals.p1000,main="Weak signals for pn=1000")
boxplot(weak.signals.p3000,main="Weak signals for pn=3000")

recall.p300 <- data.frame(sens_las.p300.n50,sens_alas.p300.n50,sens_scad.p300.n50,acc_las.p300.n50,acc_alas.p300.n50,acc_scad.p300.n50)
recall.p1000 <- data.frame(sens_las.p1000.n50,sens_alas.p1000.n50,sens_scad.p1000.n50,acc_las.p1000.n50,acc_alas.p1000.n50,acc_scad.p1000.n50)
recall.p3000 <- data.frame(sens_las.p3000.n50,sens_alas.p3000.n50,sens_scad.p3000.n50,acc_las.p3000.n50,acc_alas.p3000.n50,acc_scad.p3000.n50)

recall.p300
recall.p1000
recall.p3000
###################################################################################################
plot(f,type="l",ylim=c(min(f),max(f)),col=1,main="n=200, pn=3000, Scenario 1, rho=0.80",xlab="t",ylab="f(t)")
par(new=TRUE)
plot(WR.lasso.fhat.p3000.n200,type="l",ylim=c(min(f),max(f)),col=2,xlab="t",ylab="f(t)")
par(new=TRUE)
plot(WR.alasso.fhat.p3000.n200+0.1,type="l",ylim=c(min(f),max(f)),col=3,xlab="t",ylab="f(t)")
par(new=TRUE)
plot(WR.scad.fhat.p3000.n200,type="l",ylim=c(min(f),max(f)),col=4,xlab="t",ylab="f(t)")
legend("bottomleft",legend=c("Real f", "Lasso-PSE","aLasso-PSE","SCAD-PSE"),
       col=c(1, 2,3,4),lty=c(1,1,1,1))
grid()

write.csv(df_list.p300,"p300sc2rho2.csv")
write.csv(df_list.p1000,"p800sc2rho2.csv")
write.csv(df_list.p3000,"p1500sc2rho2.csv")

write.csv(recall.p300,"RECALL.p300sc2rho2.csv")
write.csv(recall.p1000,"RECALL.p800sc2rho2.csv")
write.csv(recall.p3000,"RECALL.p1500sc2rho2.csv")