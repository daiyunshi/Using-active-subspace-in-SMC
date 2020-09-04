##Sequencial Importance Sampler
##example with 5-dimensional normal distribution
library(MASS)
library(stats)
library(mvtnorm)
library(mnormt)

N=50 ##nuber of observed data points
set.seed(1966041)
options(digits=5)

##set up prior distribution
pri_mean <- as.matrix(c(200,250,90,180,20))
pri_sigma <- diag(5)
pri_dist <- function(x){
  exp(-1/2*t((t(as.matrix(x))-pri_mean))%*%solve(pri_sigma)%*%(t(as.matrix(x))-pri_mean))
}

##set up posterior distribution
post_mean <- c(196,254,91,170,16)
mat <- matrix(rnorm(25),5,5)
post_sigma <- mat*diag(5)*mat
post_dist <- function(x){
  exp(-1/2*t((t(as.matrix(x))-post_mean))%*%solve(post_sigma)%*%(t(as.matrix(x))-post_mean))
}

##simulate data
mydata <- mvrnorm(N,post_mean,diag(5))
mydata <- as.data.frame(mydata)


##set loglikelihood to avoid likelihood-->0 when N-->+¡Þ
loglikelihood <- function(obs_data,uk_mean,N){
  ##m <- c()
  ##for (i in 1:N) {
  ##  m[i] <- -1/2*t((t(as.matrix(obs_data[i,]))-t(uk_mean)))%*%solve(pri_sigma)%*%(t(as.matrix(obs_data[i,]))-t(uk_mean))
  ##}
  ##return(sum(m))
  mean_matrix <- matrix(t(uk_mean),nrow = N,ncol = length(uk_mean),byrow = T)
  mass_matrix <- -1/2*(as.matrix(obs_data)-mean_matrix)%*%solve(pri_sigma)%*%t(as.matrix(obs_data)-mean_matrix)
  return(sum(diag(mass_matrix)))
}


## the SMC algorithm
SMC_sampler <- function(particle_size){
  ##initialization
  log_weight <- c()
  norm_weight <- c()
  particle_size=10000
  particles <- as.data.frame(mvrnorm(particle_size,pri_mean,pri_sigma))
  for(i in 1:particle_size){
    log_weight[i] <- 0.001*loglikelihood(mydata,particles[i,],N)
  }
  norm_weight <- exp(log_weight)/sum(exp(log_weight))
  resample_index <- sample(1:particle_size,particle_size,replace=T,prob=norm_weight)
  resample_particles <- as.data.frame(matrix(NA,nrow = particle_size,ncol = 5))
  for(i in 1:particle_size){
    resample_particles[i,] <- particles[resample_index[i],]
  }
  ##Random Walk M-H
  MCMC_particles_opt_scale <- as.data.frame(matrix(NA,nrow = particle_size,ncol = 5))
  for (i in 1:particle_size) {
    proposal <- mvrnorm(1,t(resample_particles[i,]),diag(5)*2.38^2/5) 
    U <- runif(1)
    alpha <- min(1,post_dist(proposal)/post_dist(t(resample_particles[i,])))
    if(U < alpha){
      MCMC_particles_opt_scale[i,] <- proposal
    }else{
      MCMC_particles_opt_scale[i,] <- resample_particles[i,]
    }
    MCMC_1 <- MCMC_particles_opt_scale
  }
  #End of Initialization
  
  #Loop, for each iteration power of likelihood + 0.001
  for(j in 2:1000){
    #get weight
    for(i in 1:particle_size){
      log_weight[i] <- 0.001*loglikelihood(mydata,MCMC_particles_opt_scale[i,],N)
    }
    norm_weight <- exp(log_weight)/sum(exp(log_weight))
    ##resample
    resample_index <- sample(1:particle_size,particle_size,replace=T,prob=norm_weight)
    resample_particles <- as.data.frame(matrix(NA,nrow = particle_size,ncol = 5))
    for(i in 1:particle_size){
      resample_particles[i,] <- MCMC_particles_opt_scale[resample_index[i],]
    }
    ##MCMC move with RWM
    MCMC_particles_opt_scale <- as.data.frame(matrix(NA,nrow = particle_size,ncol = 5))
    for (i in 1:particle_size) {
      proposal <- mvrnorm(1,t(resample_particles[i,]),diag(5)*2.38^2/5) 
        U <- runif(1)
      alpha <- min(1,post_dist(proposal)/post_dist(t(resample_particles[i,])))
      if(U < alpha){
        MCMC_particles_opt_scale[i,] <- proposal
      }else{
        MCMC_particles_opt_scale[i,] <- resample_particles[i,]
      }
    }
  }
  resample_particles
}
   
#########
t.test(RYM10_2.38[90000:100000],Y_ASMH_M10_2.38[90000:100000])
t.test(RYM10_2.38[90000:100000])
t.test(Y_ASMH_M10_2.38[90000:100000])
