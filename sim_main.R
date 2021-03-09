setwd("~/malaria")

library(fields)
library(quadprog)
library(geoR)
library(lhs)
library(DiceKriging)
library(DiceOptim)

source("functions.R")
# NUmber of zones
ns <- 100
# Average treatment assignment rate 
trt_rate <- 0.5

# Adjacency matrix
Adj_Mat <- matrix(0,ns,ns)
ngrid <- sqrt(ns)
for (i in 1:ngrid){
  for(j in 1:ngrid){
    pos <- (i-1)*ngrid + j
    if(i==1 & j==1){
      nb_pos <- c(pos+1,pos+ngrid)
    }
    else if (i==1 & j==ngrid){
      nb_pos <- c(pos-1,pos+ngrid)
    }
    else if(i==ngrid & j==1){
      nb_pos <- c(pos-ngrid,pos+1)
    }
    else if(i==ngrid & j==ngrid){
      nb_pos <- c(pos-ngrid,pos-1)
    }
    else if(i==1 & j>1 & j<ngrid){
      nb_pos <- c(pos-1,pos+1,pos+ngrid)
    }
    else if(i==ngrid & j>1 & j<ngrid){
      nb_pos <- c(pos-ngrid,pos-1,pos+1)
    }
    else if(j==1 & i>1 & i<ngrid){
      nb_pos <- c(pos-ngrid,pos+1,pos+ngrid)
    }
    
    else if(j==ngrid & i>1 & i<ngrid){
      nb_pos <- c(pos-ngrid,pos-1,pos+ngrid)
    }
    else{
      nb_pos <- c(pos-ngrid,pos-1,pos+1,pos+ngrid)
    }
    Adj_Mat[pos,nb_pos] <- 1
  }
}
m <- rowSums(Adj_Mat)
S <- as.matrix(expand.grid(x=seq(1,ngrid,length.out=ngrid),y=seq(1,ngrid,length.out=ngrid)))
# Number of years to simulate
nt <- 6

# set true values of model parameters
beta <- c(0.2,0.12,-0.7,-0.1) # int,env,A,A^2,A*env
a1 <- 0.9  # prevY
b1 <- 0.1  # A*prevY

a2 <- 0.1 # nb_prevY
b2 <- -0.1  # A*nb_prevY

rho <- 0.9
sigs <- 0.1
sige <- 0.01
sigs_base <- 0.5
Sigma <- solve(diag(m)-rho*Adj_Mat) * sigs^2

params_or <- list(beta=beta,a1=a1,b1=b1,a2=a2,b2=b2,Sigma=Sigma,sige=sige)

# Simulate one data set

# Simulate environmental covariates
temp_matern <- c(log(2),log(0.5))
Cor_temp <- corfx(rdist(S),temp_matern)
X_env <- rnorm(ns)
X_env <- t(chol(Cor_temp))%*%X_env

# Simulate treatment assignment and disease rate
trt <- matrix(0,ns,nt) 
Y <- matrix(0,ns,nt)
Y[,1] <- t(chol(Sigma))%*%rnorm(ns) /sigs*sigs_base
Y[,1] <- Y[,1]-mean(Y[,1])
for(t in 2:nt){
  X_trt <- 0.1*(t-1)+rnorm(ns,0,0.05)
  X_trt[X_trt<0] <- 0
  X_trt[X_trt>1] <- 1
  trt[,t] <- X_trt
  Y[,t] <- cbind(1,X_env,X_trt,X_trt*X_env)%*%beta +
    (a1+b1*X_trt)*Y[,t-1] + (a2+b2*X_trt)*(Adj_Mat/m)%*%Y[,t-1] + 
    t(chol(Sigma))%*%rnorm(ns)
}
Y_ms <- Y + rnorm(ns*nt)*sige
rate <- expit(Y_ms)

X <-  array(0,c(ns,4,nt))
for(t in 1:nt){
  X[,,t] <- cbind(1,X_env,trt[,t],X_env*trt[,t])
}

# Fit the Bayesian spatiotemporal model
fit <- ST_Krige(Y_ms,X,trt,Adj_Mat)

#Estimate the loss value associated with the policy defined by alpha as the weights
#in the priorty score based on fitted model. 
#Linear uitility function or quadratic utility function can be used
get.value <- function(alpha,nsamps=200,policy=policy_linear){
  value <- get.value.base.pos(alpha, fit, X_env, Y_ms[,nt], nyear=5, nsamps=nsamps,policy=policy)
  value <- value + 0.0001*sum(alpha^2)
  return(value)
}

# Create an initial set of priority score weights from a Latin hypercube design
doe <- optimumLHS(n=100, k=4)
doe[,1:3] <- doe[,1:3]*10-5
doe[,4] <- doe[,4]
n1   <- nrow(doe)
val1 <- rep(0,n1)
# Evaluate the value for each priority score weight vector in the initial set
for(k in 1:n1){
  val1[k] <- get.value(doe[k,])
  if(k%%10==0){
    print(paste("Done with",k,"of",n1))
  }
}

# Fit a Gaussian process regression model
noise.var <- 0.0005
model <- km(y~1, design=data.frame(doe), response=data.frame(y=val1),
            covtype="exp",
            lower=rep(.01,4), upper=rep(1000,4))
# Get optimized priority score weights by sequentially selecting 
# next weight vector that optimizes expected improvement
res <- noisy.optimizer(optim.crit="EQI", optim.param=list(quantile=0.5),
                       model=model, n.ite=200, noise.var=noise.var, funnoise=get.value,
                       lower=c(rep(-10,3),0.001), upper=rep(10,4),
                       control=list(print.level=0),
                       NoiseReEstimate=FALSE, CovReEstimate=TRUE)

# Estimated optimal weights
alpha_opt <- res$best.x
# The loss value assciated with the linear utility policy
val_linear <- get.value.true(res$best.x, params_or, X_env, Y_ms[,nt],Y[,nt], 
                             nyear=5, nsamps=1000,policy=policy_linear)
# # The loss value assciated with the highest rate policy
val_highY <- get.value.true(res$best.x, params_or, X_env, Y_ms[,nt],Y[,nt], 
                             nyear=5, nsamps=1000,policy=policy_highY)
# The loss value assciated with the even policy
val_even <- get.value.true(res$best.x, params_or, X_env, Y_ms[,nt],Y[,nt], 
                           nyear=5, nsamps=1000,policy=policy_even)
# The loss value assciated with the policy with no resource allocation
val_no <- get.value.true(res$best.x, params_or, X_env, Y_ms[,nt],Y[,nt], 
                         nyear=5, nsamps=1000,policy=policy_no)
# The loss value assciated with the policy with resource allocated to everyone
val_all <- get.value.true(res$best.x, params_or, X_env, Y_ms[,nt],Y[,nt], 
                          nyear=5, nsamps=1000,policy=policy_all)
