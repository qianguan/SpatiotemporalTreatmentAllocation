library(rgdal)
library(maptools)
library(spdep)
library(geoR)
library(fields)
library(quadprog)

source("functions.R")

setwd("C:/NCSU/Research/Malaria/real_data_code")

#read Polygon data of health zones in DRC
healthzone <- readOGR(".","HealthZones_InterventionData")

n.cell <- length(healthzone)
healthzone.cell <- SpatialPolygonsDataFrame(healthzone,data.frame(id=1:n.cell),match.ID = F)
for (i in 1:515){
  healthzone.cell@polygons[[i]]@Polygons[[1]]@coords <- 
    round(healthzone.cell@polygons[[i]]@Polygons[[1]]@coords,3)
}
centroids <- getSpPPolygonsLabptSlots(healthzone.cell)
nb <- poly2nb(healthzone.cell,queen=FALSE)

# Create adjacency matrix
Adj_Mat <- matrix(0,n.cell,n.cell)
for (i in 1:n.cell){
  Adj_Mat[i,nb[[i]]] <- 1
}
m <- rowSums(Adj_Mat)
no_nb_list <- which(m==0)

# Manually identify neighbors for some health zones
man_nb <- NULL
man_nb[[1]] <- 56 #57
man_nb[[2]] <- 510 #199
man_nb[[3]] <- 264 #331
man_nb[[4]] <- 291 #333
man_nb[[5]] <- 492 #491
man_nb[[6]] <- 486 #502
man_nb[[7]] <- 486 #503
man_nb[[8]] <- 476 #504
for(i in 1:length(no_nb_list)){
  Adj_Mat[no_nb_list[i],man_nb[[i]]] <- 1
  Adj_Mat[man_nb[[i]],no_nb_list[i]] <- 1
}
m <- rowSums(Adj_Mat)
# remove(healthzone)
# remove(healthzone.cell)

# read PfPR and intervention data
PfPR <- read.csv('PfPR_2007_2015.csv')
ITN <- read.csv('ITN_2007_2015.csv')

# read temperature and precipitation data
tavg <- read.csv('tavg.csv')
prec <- read.csv('prec.csv')

avg_temp <- apply(tavg[,4:15],1,mean)
std_temp <- (avg_temp-mean(avg_temp))/sd(avg_temp)
avg_prec <- apply(prec[,4:15],1,mean)
std_prec <- (avg_prec-mean(avg_prec))/sd(avg_prec)

X_env <- cbind(std_temp,std_prec)

ns <- nrow(PfPR)
nt = 9
X <-  array(0,c(ns,6,nt))
for(t in 1:nt){
  X[,,t] <- cbind(1,X_env,ITN[,t+3],X_env*ITN[,t+3])
}
Y <- logit(PfPR[,4:12])
# fit spatio-temporal model
fit <- ST_Krige(as.matrix(Y),X,as.matrix(ITN[,4:12]),Adj_Mat)

trt_rate <- 0.5

# Posterior draws to estimate posterior distribution of policy parameters
ndraws <- 100
draws <- sample(3000:5000,ndraws,replace=FALSE)

library(quadprog)
library(lhs)
library(DiceOptim, lib.loc="~/R-Library/")
library(DiceKriging)

# Estimated the loss value based on a posterior draw from the fitted model
get.value.base.pos <- function(alpha, params, X_env, Y_base, nyear, nsamps,policy,draw){
  Ypred <- array(0,c(nsamps,ns,nyear))
  for(i in 1:nsamps){
    iter <- draw
    curY <- Y_base
    cur_latent <- params$latent_base[iter,,9]
    cholSigma <- chol(solve(diag(m)-params$keepers[iter,3]*Adj_Mat) * params$keepers[iter,2]^2)
    beta <- params$ab[iter,1:6]
    a1 <- params$ab[iter,7]
    b1 <- params$ab[iter,8]
    a2 <- params$ab[iter,9]
    b2 <- params$ab[iter,10]
    for(year in 1:nyear){
      lagY <- curY
      lag_latent <- cur_latent
      Y_std <- lagY
      nb_std <- (Adj_Mat/m)%*%lagY
      trt <- policy(X_env,Y_std,nb_std, alpha)
      if(sum(is.na(trt))!=0){break}
      cur_latent <- cbind(1,X_env,trt,trt*X_env)%*%beta + 
        (a1+b1*trt)*lag_latent + (a2+b2*trt)*(Adj_Mat/m)%*%lag_latent+
        t(cholSigma)%*%rnorm(ns)
      if(sum(is.na(cur_latent))!=0){break}
      curY <- cur_latent + rnorm(ns,0,params$keepers[iter,1])
      Ypred[i,,year] <- curY
    }
  }
  return(mean(expit(Ypred),na.rm=TRUE))
}

#Store the posterior distribution of the optimal value and optimal policy paramters
val <- rep(0,ndraws)
alpha_opt <- matrix(0,ndraws,5)

for (ii in 1:ndraws){
  draw = draws[ii]
  #Estimate the loss value associated with the policy defined by alpha as the weights
  #in the priorty score based on a posterior draw from the fitted model. 
  #Linear uitility function or quadratic utility function can be used
  get.value <- function(alpha,nsamps=100,policy=policy_quad){
    value <- get.value.base.pos(alpha, fit, X_env, as.matrix(Y[,9]), 
                                nyear=5, nsamps=nsamps,policy=policy,draw=draw)
    value <- value + 0.0001*sum(alpha^2)
    return(value)
  }
  doe <- optimumLHS(n=50, k=5)
  doe[,1:4] <- doe[,1:4]*10-5
  doe[,5] <- doe[,5]
  n1   <- nrow(doe)
  val1 <- rep(0,n1)
  for(k in 1:n1){
    val1[k] <- get.value(doe[k,])
    if(k%%10==0){
      print(paste("Done with",k,"of",n1))
    }
  }
  
  model <- km(y~1, design=data.frame(doe), response=data.frame(y=val1),
              covtype="exp",
              lower=rep(.01,5), upper=rep(1000,5))
  
  res <- noisy.optimizer(optim.crit="EQI", optim.param=list(quantile=0.5),
                         model=model, n.ite=100, noise.var=0.0001, funnoise=get.value,
                         lower=c(rep(-10,4),0.001), upper=rep(10,5),
                         NoiseReEstimate=FALSE, CovReEstimate=TRUE)
  
  val[ii] <- get.value(res$best.x)
  alpha_opt[ii,] <- as.vector(res$best.x)
}