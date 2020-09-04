eASMH_eq1 <- function(sigma=1.5,M){
  ##initialization qz~N(0,1)
  x0 <- as.vector(X <- mvrnorm(1,as.vector(c(0,0)),diag(2)))
  y0 <- as.numeric(t(eg_vectors_c1[,1])%*%x0)
  z0j <- as.vector(rnorm(M))
  Z <- matrix(NA,nrow=M,ncol=100000)
  Z[,1] <- z0j
  Y <- c()
  Y[1] <- y0
  rx <- as.vector(eg_vectors_c1[,1])%*%t(as.vector(rep(y0,M)))+as.vector(eg_vectors_c1[,2])%*%t(z0j)
  w0j <- c()
  for(j in 1:M){
    w0j[j] <-f(d,rx[,j],1)/dnorm(z0j[j])
  }
  b <- exp(-mean(w0j))
  ##END of initial
  
  for (i in 2:100000) {
    y_prime <- rnorm(1,Y[i-1],sigma)
    z_prime <- as.vector(rnorm(M))
    x_prime <- as.vector(eg_vectors_c1[,1])%*%t(as.vector(rep(y_prime,M)))+as.vector(eg_vectors_c1[,2])%*%t(z_prime)
    w_prime <- c()
    for(m in 1:M){
      w_prime <- f(d,x_prime[,m],1)/dnorm(z_prime[m])
    }
    b_prime <- exp(-mean(w_prime))
    a <- rbinom(1,1,min(1,b_prime*dnorm(y_prime)/(b*dnorm(Y[i-1]))))
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

eASMH_eq1_M10_0.5 <- eASMH_eq1(0.5,10)

eASMH_eq1_M10_2.38 <- eASMH_eq1(2.38,10)

eASMH_eq1_M1 <- eASMH_eq1(1.5,1)
par(mfcol=c(1,1))
hist(eASMH_eq1_M10[50000:100000],breaks = 100,probability = T,xlim = c(-2,2))
hist(eASMH_eq1_M1[50000:100000],breaks = 100,probability = T,xlim = c(-2,2))
t.test(eASMH_eq1_M10[90000:100000])
plot(eASMH_eq1_M10_0.5[60001:70000],type = 's')
abtest <- matrix(NA,nrow = 2,ncol=30)
for(i in 1:30){
  y <- eASMH_eq1(0.5,10)
  abtest[1,i] <- mean(y[50000:10000])
  abtest[2,i] <- var(y[50000:10000])
}

t.test(abtest[1,5])

par(mfcol=c(1,2))
hist(RYM10[50000:100000],breaks = 50,xlim = c(-2,2),ylim = c(0,2),xlab = 'y',probability = T,col='salmon')
hist(eASMH_eq1_M10_0.5[50000:100000],breaks = 50,xlim = c(-2,2),ylim = c(0,2),xlab = 'y',probability = T,col='skyblue')

plot(RYM10[60000:70000],type='s',ylim = c(-2,2),ylab = 'Markov Chain States',xlab = 'Runs',col='salmon')
plot(eASMH_eq1_M10_0.5[60000:70000],type='s',ylim = c(-2,2),ylab = 'Markov Chain States',xlab = 'Runs',col='skyblue')
acfeq2 <- acf(RYM10_2.38[50000:100000],lag.max = 50)
acfeq1 <- acf(eASMH_eq1_M10_2.38[50000:100000],lag.max = 50)
