RWMH <- function(sigma,N){
  x <- matrix(NA,nrow = 2,ncol = N)
  x[,1] <- mvrnorm(1,c(0,0),diag(2))
  U <- runif(N)
  for(i in 2:N){
    x_prime <- mvrnorm(1,x[,i-1],sigma*diag(2))
    v <- exp(-f(d,x_prime,1))*dmvnorm(x_prime)/(exp(-f(d,x[,i-1],1))*dmvnorm(x[i-1]))
    a <- min(1,v)
    if(U[i] < a){
      x[,i] <- x_prime
    }else{
      x[,i] <- x[,i-1]
    }
  }
  return(x)
}

vnla <- RWMH(1,1000000)
vnla2.38 <- RWMH(2.38,1000000)
plot(vnla[1,600000:610000],type = 's',ylim = c(-4,4),col='salmon')
par(new=T)
plot(vnla[2,600000:610000],type='s',ylim = c(-4,4),col='skyblue')
IAT(vnla[1,600000:610000])
IATRWMH <- IAT(vnla2.38[1,900000:1000000])
IATREX2.38 <- IAT(REX2.38[1,900000:1000000])
IATEASM <- IAT(RXM10_2.38[1,900000:1000000])
IATRWMH
IATREX2.38
IATEASM