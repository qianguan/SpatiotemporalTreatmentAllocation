# logit transformation
logit<-function(x){log(x/(1-x))}
expit<-function(x){exp(x)/(1+exp(x))}

# Compute the values of the Matern correlation function
corfx <- function(d,theta){
  matern(d,exp(theta[1]),exp(theta[2]))
}

# Fit the Bayesian spatiotemporal model and get posterior samples of parameters
ST_Krige<-function(Y,X,trt,A,
                   a=.1,b=.01,
                   sd_beta=100,
                   iters=5000,burn=2000){
  
  library(emulator)
  library(geoR)

  ns       <- nrow(Y)
  nt       <- ncol(Y)
  p        <- dim(X)[2]
  
  m     <- rowSums(A)
  adj   <- apply(A==1,1,which)
  nei1  <- row(A)[A==1]
  nei2  <- col(A)[A==1]
  
  
  # Initial values
  
  pp <- p + 4
  ab <- rep(0,pp)
  beta <- ab[1:p]
  a1 <- ab[p+1]
  b1 <- ab[p+2]
  a2 <- ab[p+3]
  b2 <- ab[p+4]
  
  taue    <- 1/0.02^2
  taus    <- 1
  rho   <- 0.8
  
  theta <- Y
  
  # Pre-compute the determinant for all rho of interest
  
  C        <- diag(1/sqrt(m))%*%A%*%diag(1/sqrt(m))
  lambda   <- eigen(C)$values
  rho.grid <- seq(0.5,0.999,0.001)
  logd     <- rho.grid
  for(j in 1:length(logd)){
    logd[j] <- sum(log(1-rho.grid[j]*lambda))*nt
  }
  
  Q <- taus*(diag(m)-rho*A)
  H <- array(0,c(ns,ns,nt))
  for(t in 1:nt){
    H[,,t] <- (a1+b1*trt[,t])*diag(ns) + (a2+b2*trt[,t])*A/m
  }
  
  # Keep track of stuff
  keep.ab <- matrix(0,iters,pp)
  keepers <- matrix(0,iters,3)
  colnames(keepers) <- c("sig2e","sig2s","rho")
  keep.latent_base <- array(0,c(iters,ns,nt))

  
  for(i in 1:iters){
    
    
    ##############################################:
    #####       THETA (Gibbs)      #######:
    ##############################################:
    for(t in 1:nt){
      if(t==1){
        VVV <- solve(taue*diag(ns) + t(H[,,t+1])%*%Q%*%H[,,t+1])
        MMM <- taue*Y[,t] + t(H[,,t+1])%*%Q%*%(theta[,t+1]-X[,,t+1]%*%beta)
      }
      if(t>1 & t<nt){
        VVV <- solve(taue*diag(ns) + Q + t(H[,,t+1])%*%Q%*%H[,,t+1])
        MMM <- taue*Y[,t] + Q%*%(H[,,t]%*%theta[,t-1]+X[,,t]%*%beta) +
          t(H[,,t+1])%*%Q%*%(theta[,t+1]-X[,,t+1]%*%beta)
      }
      if(t==nt){
        VVV <- solve(taue*diag(ns) + Q)
        MMM <- taue*Y[,t] + Q%*%(H[,,t]%*%theta[,t-1]+X[,,t]%*%beta)
      }
      theta[,t] <- VVV%*%MMM + t(chol(VVV))%*%rnorm(ns)
    }
    
    ##############################################:
    #####       trt PARAMETERS (Gibbs)    #######:
    ##############################################:
    pp <- p + 4
    X_latent <- array(0,c(ns,pp,nt))
    for(t in 2:nt){
      nb_theta <- (A/m)%*%theta[,t-1]
      X_latent[,,t] <- cbind(X[,,t],theta[,t-1],trt[,t]*theta[,t-1],nb_theta,trt[,t]*nb_theta)
    }
    XQX <- 0
    for (t in 2:nt){
      XQX <- XQX + t(X_latent[,,t])%*%Q%*%X_latent[,,t]
    }
    VVV  <- solve(XQX + diag(pp)/sd_beta^2)
    MMM  <- 0
    for(t in 2:nt){
      MMM <- MMM + t(X_latent[,,t])%*%Q%*%theta[,t]
    }
    ab <- VVV%*%MMM + t(chol(VVV))%*%rnorm(pp)
    beta <- ab[1:p]
    a1 <- ab[p+1]
    b1 <- ab[p+2]
    a2 <- ab[p+3]
    b2 <- ab[p+4]
    for(t in 1:nt){
      H[,,t] <- (a1+b1*trt[,t])*diag(ns) + (a2+b2*trt[,t])*A/m
    }
    
    ##############################################:
    #####          VARIANCE (Gibbs)        #######:
    ##############################################:
    epsilon <- matrix(0,ns,t)
    for(t in 2:nt){
      epsilon[,t] <- theta[,t] - H[,,t]%*%theta[,t-1] - X[,,t]%*%beta
    }
    
    TAT   <- sum(epsilon[nei1,]*epsilon[nei2,])
    TMT   <- sum(m*epsilon^2)
    taus  <- rgamma(1,ns*(nt-1)/2+a,(TMT-rho*TAT)/2+b)
    
    taue   <- rgamma(1,ns*nt/2+a,sum((Y-theta)^2)/2+b)
    
    # CAR DEPENDENCE PARAMETER
    
    R    <- 0.5*logd + 0.5*rho.grid*TAT*taus
    rho  <- sample(rho.grid,1,prob=exp(R-max(R)))
    
    Q <- taus*(diag(m)-rho*A)
    
    
    ##############################################:
    #####        KEEP TRACK OF STUFF       #######:
    ##############################################:
    
    keep.ab[i,]    <- ab
    keepers[i,1]   <- 1/sqrt(taue)
    keepers[i,2]   <- 1/sqrt(taus)
    keepers[i,3]   <- rho
    keep.latent_base[i,,] <- theta
    
    if(i%%1000==0){
      print(paste("Done with",i,"of",iters))
    }
  }   
  
  output <- list(ab=keep.ab,keepers=keepers,
                 latent_base=keep.latent_base)
  
  return(output)}


#Highest rate policy
policy_highY <- function(X,Y,nb_Y,alpha){
  X_trt <- rep(0,ns)
  X_trt[sort(Y,decreasing=TRUE, index.return=TRUE)$ix[1:(ns*trt_rate)]] <- 1
  return(X_trt)
}

#Even policy
policy_even <- function(X,Y,nb_Y,alpha){
  X_trt <- rep(trt_rate,ns)
  return(X_trt)
}

#Policy with no resource allocation 
policy_no <- function(X,Y,nb_Y,alpha){
  X_trt <- rep(0,ns)
  return(X_trt)
}

#Policy with resource allocated to everyone
policy_all <- function(X,Y,nb_Y,alpha){
  X_trt <- rep(1,ns)
  return(X_trt)
}

#Random policy
policy_random <- function(X,Y,nb_Y,alpha){
  X_trt <- expit(rnorm(ns))
  X_trt <- X_trt/sum(X_trt)*trt_rate*ns
  return(X_trt)
}

#Linear utility policy
policy_linear <- function(X,Y,nb_Y,alpha){
  alpha_pri <- alpha[1:(ncol(X)+2)]
  alpha01 <- alpha[ncol(X)+3]
  pri <- expit(cbind(X,Y,nb_Y)%*%alpha_pri)
  Q1 <- (diag(m)-Adj_Mat*0.999)
  HHH <- Q1*alpha01 /10
  Amat <- matrix(0,ns,2*ns+1)
  Amat[,1] <- -1/ns
  Amat[1,2:(ns+1)] <- 1
  Amat[1,(ns+2):(2*ns+1)] <- -1
  Aind <- matrix(0,ns+1,2*ns+1)
  Aind[1,] <- c(ns,rep(1,2*ns))
  Aind[2:(ns+1),1] <- 1:ns
  Aind[2,2:(2*ns+1)] <- c(1:ns,1:ns)
  A.optim <- solve.QP.compact(Dmat=2*HHH, dvec=pri, 
                              Amat=Amat,Aind=Aind,
                              bvec=c(-trt_rate,rep(0,ns),rep(-1,ns)))$solution
  return(A.optim)
}


#Quadratic utility policy
policy_quad <- function(X,Y,nb_Y,alpha){
  alpha_pri <- alpha[1:(ncol(X)+2)]
  alpha0 <- alpha[ncol(X)+3]
  pri <- expit(cbind(X,Y,nb_Y)%*%alpha_pri)
  pri[which(pri<0.0001)] <- 0.0001
  Q1 <- diag(as.vector(pri))
  Q2 <- (diag(m)-Adj_Mat*0.999)*alpha0/10
  A.optim <- solve.QP(Dmat=Q1+Q2, dvec=pri, 
                      Amat=cbind(matrix(-1/ns,ns,1),diag(1,ns,ns),diag(-1,ns,ns)),
                      bvec=c(-trt_rate,rep(0,ns),rep(-1,ns)))$solution
  return(A.optim)
}

# Estimated the loss value based on fitted model
get.value.base.pos <- function(alpha, params, X_env, Y_base, nyear, nsamps,policy){
  Ypred <- array(0,c(nsamps,ns,nyear))
  X_env_std <- (X_env-mean(X_env))/sd(X_env)
  samps <- sample(3000:5000,nsamps,replace=FALSE)
  nt <- dim(params$latent_base)[3]
  for(i in 1:nsamps){
    iter <- samps[i]
    curY <- Y_base
    cur_latent <- params$latent_base[iter,,nt]
    cholSigma <- chol(solve(diag(m)-params$keepers[iter,3]*Adj_Mat) * params$keepers[iter,2]^2)
    beta <- params$ab[iter,1:4]
    a1 <- params$ab[iter,5]
    b1 <- params$ab[iter,6]
    a2 <- params$ab[iter,7]
    b2 <- params$ab[iter,8]
    for(year in 1:nyear){
      lagY <- curY
      lag_latent <- cur_latent
      Y_std <- lagY
      nb_std <- (Adj_Mat/m)%*%lagY
      trt <- policy(X_env_std,Y_std,nb_std, alpha)
      cur_latent <- cbind(1,X_env,trt,trt*X_env)%*%beta + 
        (a1+b1*trt)*lag_latent + (a2+b2*trt)*(Adj_Mat/m)%*%lag_latent+
        t(cholSigma)%*%rnorm(ns)
      curY <- cur_latent + rnorm(ns,0,params$keepers[iter,1])
      Ypred[i,,year] <- curY
    }
  }
  return(mean(expit(Ypred)))
}

# Estimate the loss value based on true model
get.value.true <- function(alpha, params, X_env, Y_base, latent_base,nyear, nsamps,policy){
  Ypred <- array(0,c(nsamps,ns,nyear))
  X_env_std <- (X_env-mean(X_env))/sd(X_env)
  for(i in 1:nsamps){
    curY <- Y_base
    cur_latent <- latent_base
    for(year in 1:nyear){
      lagY <- curY
      lag_latent <- cur_latent
      Y_std <- lagY
      nb_std <- (Adj_Mat/m)%*%lagY
      trt <- policy(X_env_std,Y_std,nb_std, alpha)
      cur_latent <- cbind(1,X_env,trt,trt*X_env)%*%params$beta + 
        (params$a1+params$b1*trt)*lag_latent + (params$a2+params$b2*trt)*
        (Adj_Mat/m)%*%lag_latent + t(chol(params$Sigma))%*%rnorm(ns)
      curY <- cur_latent + rnorm(ns,0,params$sige)
      Ypred[i,,year] <- curY
    }
  }
  return(mean(expit(Ypred)))
}


