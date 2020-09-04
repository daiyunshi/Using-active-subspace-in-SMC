library(MASS)
library(stats)
library(mvtnorm)
library(mnormt)
library(rgl)
library(ggplot2)
library(plyr)
library(LaplacesDemon)
Q <- 1/2*matrix(c(sqrt(2),-sqrt(2),sqrt(2),sqrt(2)),2,2) ##°´ÁÐÌî³ä
A1 <- matrix(c(1,0,0,0.01),2,2)
A2 <- matrix(c(1,0,0,0.95),2,2)
M1 <- Q%*%A1%*%t(Q)
M2 <- Q%*%A2%*%t(Q)

X <- mvrnorm(1000,as.vector(c(0,0)),diag(2))
X <- t(as.matrix(X))
d <- 0.9                       ##fixed data##

##model
two_para_model <- function(x,i){
  if(i==1){
    1/2*t(x)%*%M1%*%x
  }else{
    1/2*t(x)%*%M2%*%x
  }
}

###################
##misfit function##
###################
f <- function(d=0.9,x,i){
  if(i==1){
    5*as.numeric(d-1/2*t(x)%*%M1%*%x)^2
  }else{
    5*as.numeric(d-1/2*t(x)%*%M2%*%x)^2
  }
}


dif_model <- function(obs_data,x,i){
  if(i==1){
    A <- Q%*%A1%*%t(Q)
    1/0.1*t(t(x)%*%(M1+t(M1)))*as.numeric(obs_data-1/2*t(x)%*%M1%*%x)
  }else{
    1/0.1*t(t(x)%*%(M2+t(M2)))*as.numeric(obs_data-1/2*t(x)%*%M2%*%x)
  }
}

##c1
C_hat_1 <-function(d=0.9,X){
  m <- dif_model(0.9,X[,1],1)%*%t(dif_model(0.9,X[,1],1))
  for(i in 2:1000){
    m <- m + dif_model(0.9,X[,i],1)%*%t(dif_model(0.9,X[,i],1))
  }
  m/1000
} 

C1 <- C_hat_1(0.9,X)
egc1 <- eigen(C1)
eg_value_c1 <- egc1$values
eg_vectors_c1 <- egc1$vectors

###############
##from (2.16)##
###############
ghat<- function(y,M){
  z <- rnorm(M,0,1)
  g <- c()
  for(i in 1:M){
    g[i] <-f(d,as.vector(eg_vectors_c1[,1])*y+as.vector(eg_vectors_c1[,2])*z[i],1)
  }
  return(mean(g))
}
##eq2 unbiased estimate of posterior

ghat_eq2<- function(y,M){
  z <- rnorm(M,0,1)
  g <- c()
  for(i in 1:M){
    g[i] <- exp(-f(d,as.vector(eg_vectors_c1[,1])*y+as.vector(eg_vectors_c1[,2])*z[i],1))
  }
  return(mean(g))
}
#############################
##MCMC with Active Subspace##
#############################
MCMC_Active_Subspace <- function(Runs=10000, sigma=1.0, y0=-0.529,M){
  Y <- c()
  Y[1] <- y0
  U <- runif(Runs)
  for(i in 2:Runs){
    y_prime <- rnorm(1,Y[i-1],sigma)
    mh <- exp(-ghat(y_prime,M)+ghat(Y[i-1],M))*dnorm(y_prime)/dnorm(Y[i-1])
    alpha <- min(1,mh)
    if(alpha>U[i]){
      Y[i] <- y_prime
    }else{
      Y[i] <- Y[i-1]
    }
  }
  return(Y)
}

#############################
##MCMC with Active Subspace##
############################# with eq2
ASMH_eq2 <- function(Runs=10000, sigma=1.0, y0=-0.529,M){
  Y <- c()
  Y[1] <- y0
  U <- runif(Runs)
  for(i in 2:Runs){
    y_prime <- rnorm(1,Y[i-1],sigma)
    mh <- ghat_eq2(y_prime,M)/ghat_eq2(Y[i-1],M)*dnorm(y_prime)/dnorm(Y[i-1])
    alpha <- min(1,mh)
    if(alpha>U[i]){
      Y[i] <- y_prime
    }else{
      Y[i] <- Y[i-1]
    }
  }
  return(Y)
}

eq2_M10 <- ASMH_eq2(100000,0.5,0,10)
eq2_M1 <- ASMH_eq2(100000,0.5,0,1)

mean_eq2_M10 <- c()
for(i in 1:100000){
  mean_eq2_M10[i] <- mean(eq2_M10[1:i])
}

mean_eq2_M1 <- c()
for(i in 1:100000){
  mean_eq2_M1[i] <- mean(eq2_M1[1:i])
}
plot(mean_eq2_M1,type='s',ylim = c(-1,1))
plot(mean_eq2_M10,type='s',ylim = c(-1,1))
##
par(mfcol=c(1,2))
hist(eq2_M10[50000:100000],breaks = 50,col = 'salmon',probability = T,ylim = c(0,0.8))
hist(Y_ASMH_M10_2.38[20000:100000],breaks = 50,col = 'skyblue',probability = T,ylim = c(0,0.8))
Mean_eq2 <- c()
for(i in 1:30){
  y <-  ASMH_eq2(100000,2.38,0,1)
  Mean_eq2[i] <- mean(y)
}
Mean_eq1 <- c()
for(i in 1:30){
  y <-  MCMC_Active_Subspace(100000,2.38,0,1)
  Mean_eq1[i] <- mean(y)
}
t.test(Mean_eq1,Mean_eq2)
####
posterior <- function(X){
  exp(-ghat(as.numeric(as.vector(eg_vectors_c1[,1])%*%X)))*exp(-t(X)%*%solve(diag(2))%*%X)
}
posterior(REX3[,90000+5722])
post_value <- c()
for(i in 90000:100000){
  post_value[i-89999] <- posterior(REX3[,i])
}

plot3d(REX3[1,90000:100000],REX3[2,90000:100000],post_value)

#############
Y_ASMH_M10 <- MCMC_Active_Subspace(Runs = 100000,sigma =0.5 ,y0=0,M=10)
Y_ASMH_M10_2.38 <- MCMC_Active_Subspace(Runs = 100000,sigma =2.38 ,y0=0,M=10)
Y_ASMH_M1 <- MCMC_Active_Subspace(Runs = 100000,sigma =0.5 ,y0=0,M=1)
Y_ASMH_M1_2.38 <- MCMC_Active_Subspace(Runs = 100000,sigma =2.38 ,y0=0,M=1)
Y_ASMH_M100_2.38 <- MCMC_Active_Subspace(Runs = 100000,sigma =2.38 ,y0=0,M=100)
ASMH_M1_M10 <- data.frame(M=factor(rep(c("M=1","M=10"),each=10000)),Y=c(Y_ASMH_M1[90001:100000],Y_ASMH_M10[90001:100000]))
head(ASMH_M1_M10)
#############

#############
hist(eq2_M10[50000:100000],breaks = 100,probability = T)
#############
ggplot(ASMH_M1_M10, aes(x=Y, color=M)) + 
  geom_density(alpha=.5)


hist(Y_ASMH_M10[50000:100000],breaks = 50,probability = T)
hist(Y_ASMH_M1[50000:100000],breaks = 50,probability = T)

plot(Y_ASMH_M10[10000:20000],type='s',col='navyblue')
plot(Y_ASMH_M1[10000:20000],type='s',col='green3')

mean_value_M10 <- c()
for(i in 1:100000){
  mean_value_M10[i]=mean(Y_ASMH_M10[1:i])
}

mean_value_M1 <- c()
for(i in 1:100000){
  mean_value_M1[i]=mean(Y_ASMH_M1[1:i])
}

plot(mean_value_M10[1:50000],type = 's',ylim = c(-1,1),col='yellow green',xlab = 'Runs',ylab = 'Mean Value')
par(new=T)
plot(mean_value_M1[1:50000],type = 's',ylim = c(-1,1),col='blue',xlab = 'Runs',ylab = 'Mean Value')

par(mfcol=c(2,1))
plot(Y_ASMH_M1[10000:20000],type = 's')
plot(Y_ASMH_M10[50000:70000],type='s')

par(mfcol=c(1,1))
hist(Y_ASMH_M10[50000:100000],breaks = 50,probability = T,xlim=c(-2,2),ylim = c(0,0.8))
hist(RYM10[50000:100000],breaks = 50,probability = T,xlim=c(-2,2),ylim = c(0,0.8))


acf1 <- acf(Y_ASMH_M10_2.38[90000:100000],lag.max = 50)
acf2 <- acf(RYM10_2.38[90000:100000],lag.max = 50)

X_space <- matrix(NA,nrow=2,ncol=10000)

Y3_rep <- rep(Y_ASMH_M10_2.38,each=10)
Z2 <- rnorm(1000000)
REX2.38 <- as.vector(eg_vectors_c1[,1])%*%t(as.vector(Y3_rep))+as.vector(eg_vectors_c1[,2])%*%t(as.vector(Z2))
plot(REX2.38[1,c(600000:610000)],type = 's',col="salmon",ylim = c(-4,4),xlab = 'Iteration',ylab = 'Markov chain states')
par(new=T)
plot(REX2.38[2,c(600000:610000)],type = 's',col="skyblue",ylim = c(-4,4),xlab = 'Iteration',ylab = 'Markov chain states')
legend(90000,3,c("x1","x2"),col = c("salmon","skyblue"))

par(mfcol=c(1,2))
acf1 <- acf(REX2.38[1,990000:1000000],lag.max = 150)

acf2 <- acf(EY[1,990000:1000000],lag.max = 150)

