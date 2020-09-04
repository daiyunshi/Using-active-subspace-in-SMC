set.seed(980417)
eASMH <- function(sigma=1.5){
  ##initialization qz~N(0,1)
  x0 <- as.vector(X <- mvrnorm(1,as.vector(c(0,0)),diag(2)))
  y0 <- as.numeric(t(eg_vectors_c1[,1])%*%x0)
  z0j <- as.vector(rnorm(10))
  Z <- matrix(NA,nrow=10,ncol=100000)
  Z[,1] <- z0j
  Y <- c()
  Y[1] <- y0
  rx <- as.vector(eg_vectors_c1[,1])%*%t(as.vector(rep(y0,10)))+as.vector(eg_vectors_c1[,2])%*%t(z0j)
  w0j <- c()
  for (j in 1:10) {
    w0j[j] <- exp(-f(0.9,rx[,j],1))*dmvnorm(rx[,j])/dnorm(z0j[j])
  }
  b <- mean(w0j)
  ##END of initial
  
  for (i in 2:100000) {
    y_prime <- rnorm(1,Y[i-1],sigma)
    z_prime <- as.vector(rnorm(10))
    x_prime <- as.vector(eg_vectors_c1[,1])%*%t(as.vector(rep(y_prime,10)))+as.vector(eg_vectors_c1[,2])%*%t(z_prime)
    w_prime <- c()
    for(m in 1:10){
      w_prime <- exp(-f(0.9,x_prime[,m],1))*dmvnorm(x_prime[,m])/dnorm(z_prime[m])
    }
    b_prime <- mean(w_prime)
    a <- rbinom(1,1,min(1,b_prime/b))
    if(a==1){
      b <- b_prime
      Y[i] <- y_prime
      Z[,i] <- z_prime
    }else{
      Y[i] <- Y[i-1]
      Z[,i] <- Z[,i-1]
    }##end if 
  }##end for
  REX <- as.vector(eg_vectors_c1[,1])%*%t(as.vector(rep(Y,each=10)))+as.vector(eg_vectors_c1[,2])%*%t(as.vector(Z))
  return(REX)
}

eASMH1 <- function(sigma=1.5){
  ##initialization qz~N(0,1)
  x0 <- as.vector(X <- mvrnorm(1,as.vector(c(0,0)),diag(2)))
  y0 <- as.numeric(t(eg_vectors_c1[,1])%*%x0)
  z0j <- as.vector(rnorm(1))
  Z <- matrix(NA,nrow=1,ncol=100000)
  Z[,1] <- z0j
  Y <- c()
  Y[1] <- y0
  rx <- as.vector(eg_vectors_c1[,1])%*%t(as.vector(rep(y0,1)))+as.vector(eg_vectors_c1[,2])%*%t(z0j)
  w0j <- c()
  for (j in 1:1) {
    w0j[j] <- exp(-f(0.9,rx[,j],1))*dmvnorm(rx[,j])/dnorm(z0j[j])
  }
  b <- mean(w0j)
  ##END of initial
  
  for (i in 2:100000) {
    y_prime <- rnorm(1,Y[i-1],sigma)
    z_prime <- as.vector(rnorm(1))
    x_prime <- as.vector(eg_vectors_c1[,1])%*%t(as.vector(rep(y_prime,1)))+as.vector(eg_vectors_c1[,2])%*%t(z_prime)
    w_prime <- c()
    for(m in 1:1){
      w_prime <- exp(-f(0.9,x_prime[,m],1))*dmvnorm(x_prime[,m])/dnorm(z_prime[m])
    }
    b_prime <- mean(w_prime)
    a <- rbinom(1,1,min(1,b_prime/b))
    if(a==1){
      b <- b_prime
      Y[i] <- y_prime
      Z[,i] <- z_prime
    }else{
      Y[i] <- Y[i-1]
      Z[,i] <- Z[,i-1]
    }##end if 
  }##end for
  REX <- as.vector(eg_vectors_c1[,1])%*%t(as.vector(Y))+as.vector(eg_vectors_c1[,2])%*%t(as.vector(Z))
  return(Y)
}

eASMH100 <- function(sigma=1.5){
  ##initialization qz~N(0,1)
  x0 <- as.vector(X <- mvrnorm(1,as.vector(c(0,0)),diag(2)))
  y0 <- as.numeric(t(eg_vectors_c1[,1])%*%x0)
  z0j <- as.vector(rnorm(100))
  Z <- matrix(NA,nrow=100,ncol=100000)
  Z[,1] <- z0j
  Y <- c()
  Y[1] <- y0
  rx <- as.vector(eg_vectors_c1[,1])%*%t(as.vector(rep(y0,100)))+as.vector(eg_vectors_c1[,2])%*%t(z0j)
  w0j <- c()
  for (j in 1:100) {
    w0j[j] <- exp(-f(0.9,rx[,j],1))*dmvnorm(rx[,j])/dnorm(z0j[j])
  }
  b <- mean(w0j)
  ##END of initial
  
  for (i in 2:100000) {
    y_prime <- rnorm(1,Y[i-1],sigma)
    z_prime <- as.vector(rnorm(100))
    x_prime <- as.vector(eg_vectors_c1[,1])%*%t(as.vector(rep(y_prime,100)))+as.vector(eg_vectors_c1[,2])%*%t(z_prime)
    w_prime <- c()
    for(m in 1:100){
      w_prime <- exp(-f(0.9,x_prime[,m],1))*dmvnorm(x_prime[,m])/dnorm(z_prime[m])
    }
    b_prime <- mean(w_prime)
    a <- rbinom(1,1,min(1,b_prime/b))
    if(a==1){
      b <- b_prime
      Y[i] <- y_prime
      Z[,i] <- z_prime
    }else{
      Y[i] <- Y[i-1]
      Z[,i] <- Z[,i-1]
    }##end if 
  }##end for
  return(Y)
}


RYM100 <- eASMH100(sigma = 0.5)
RYM10 <- eASMH(sigma = 0.5)

RYM100_2.38 <- eASMH(sigma = 2.38)
RXM10_2.38 <- eASMH(sigma = 2.38)
RYM1_2.38 <- eASMH1(sigma = 2.38)


RYM1 <- eASMHM1(sigma = 0.5)
REXEASMH <- eASMH(sigma = 0.5)

plot(REXEASMH[1,c(600000:610000)],type = 's',col="salmon",ylim = c(-4,4),xlab = 'Iteration',ylab = 'Markov chain states')
par(new=T)
plot(REXEASMH[2,c(600000:610000)],type = 's',col="skyblue",ylim = c(-4,4),xlab = 'Iteration',ylab = 'Markov chain states')

eASMH_M1_M10 <- data.frame(M=factor(rep(c("M=1","M=10","M=100"),each=50000)),Y=c(RYM1[50001:100000],RYM10[50001:100000],RYM100[50001:100000]))

ggplot(eASMH_M1_M10, aes(x=Y, color=M)) + 
  geom_density(alpha=.5)

par(mfcol=c(3,1))
plot(RYM1[0:10000],type = 's',col='yellowgreen')
plot(RYM10[0:10000],type = 's',col='navyblue')
plot(RYM100[0:10000],type = 's',col='salmon')

mean_M100 <- c()
for(i in 1:100000){
  mean_M100[i]=mean(RYM100[1:i])
}

mean_M10 <- c()
for(i in 1:100000){
  mean_M10[i]=mean(RYM10[1:i])
}

mean_M1 <- c()
for(i in 1:100000){
  mean_M1[i]=mean(RYM1[1:i])
}
par(mfcol=c(1,1))
plot(mean_M100[1:100000],type = 's',ylim = c(-1,1),col='salmon',xlab = 'Runs',ylab = 'Mean Value')
par(new=T)
plot(mean_M10[1:100000],type = 's',ylim = c(-1,1),col='navyblue',xlab = 'Runs',ylab = 'Mean Value')
par(new=T)
plot(mean_M1[1:100000],type = 's',ylim = c(-1,1),col='yellowgreen',xlab = 'Runs',ylab = 'Mean Value')
abline(h=0)


par(mfcol=1,3)
acfM100 <- acf(RYM100_2.38,lag.max = 50)
acfM10 <- acf(RYM10_2.38,lag.max = 50)
acfM1 <- acf(RYM1_2.38,lag.max = 50)