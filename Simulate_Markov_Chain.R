## Simulate a Markov Chain

source("HMM_functions.R")
Simulate_Markov_Chain <- function()
{
  num_states <- 2
  
  ## specify some transition-matrix:
  tmat <- matrix(0,nrow=num_states,ncol=num_states)
  tmat[1,1] <- 0.999
  tmat[1,2] <- 0.001
  tmat[2,1] <- 0.002
  tmat[2,2] <- 0.998
  delta <- 20000
  
  ## obtain the steady-state-probabilities
  p1 <- tmat[2,1]/(tmat[1,2]+tmat[2,1])
  p2 <- tmat[1,2]/(tmat[1,2]+tmat[2,1])
 
  
  kmat <- matrix(0,nrow=num_states,ncol=num_states)
  ## obtain transition rate matrix:
  for(i in 1:nrow(tmat))
  {
    for(j in 1:ncol(tmat))
    {
      kmat[i,j] <- -log(1.0 -tmat[i,j])*delta
    }
  }
  ## obtain the lifetimes:
  lts <- array(0,dim=c(num_states))
  for(i in 1:num_states)
  {
    lts[i] = 0
    temp <- 0
    for(j in 1:num_states)
    {
      if(j!=i)
      {
        temp <- temp + kmat[i,j]
      }
    }
    lts[i] <- 1.0/temp
  }
  print("expected steady-state probabilities: ")
  print(p1)
  print(p2)
  print("expected lifetimes [s]: ")
  print(lts)
  
  num_pts <- 150000
  
  ## generate initial state with 
  state_sequ <- array(0,dim=c(num_pts))
  obs_sequ <- array(0,dim=c(num_pts))
  times <- seq(from=0,by=1.0/delta,length=num_pts)
  mus <- c(3.8,5.0)#array(0,dim=c(num_states))
  sigmas <- c(0.4,0.4)#array(0,dim=c(num_states))
  
  #set.seed(5)
  set.seed(6)
  
  state_sequ[1] <- sample(x=c(1,2),1,replace=TRUE,prob=c(p1,p2))
  obs_sequ[1] <- rnorm(n=1,mean=mus[state_sequ[1]],sd=sigmas[state_sequ[1]])
  
  print(state_sequ[1])
  
  for(i in 2:num_pts)
  {
    state_sequ[i] <- sample(x=c(1,2),1,replace=TRUE,prob=c(tmat[state_sequ[i-1],]))
    obs_sequ[i] <- rnorm(n=1,mean=mus[state_sequ[i]],sd=sigmas[state_sequ[i]])
  }
  plot(times,obs_sequ,pch=20,cex=0.2,xlab="Zeit [s]",ylab="Kraft [pN]")
  hist(obs_sequ,breaks=100)
  return(obs_sequ)
}

HMM_GaussianMixture_Evaluate_Stan <- function(obs_sequ,num_iters=200)
{
  require("rstan")
  
  stan_data <- list(N=length(obs_sequ),num_states=2,obs_sequ=obs_sequ,
                    alphas=c(1000,1),betas=c(1,1000),
                    mus_mu=c(4.7,4.1),sigmas_mu=c(0.6,0.6),
                    mus_sigma=c(-1,-1),sigmas_sigma=c(1,1))
  
  ## draw initial values from their respective prior distributions:
  trans_mat <- matrix(0,nrow=2,ncol=2)
  
  trans_mat[1,1] <- rbeta(1,1000,1)
  trans_mat[2,1] <- rbeta(1,1,1000)
  trans_mat[1,2] <- (1.0 - trans_mat[1,1])
  trans_mat[2,2] <- (1.0 - trans_mat[2,1])
  
  em_mus <- rnorm(n=2,mean=c(4.7,4.1),sd=c(0.6,0.6))
  em_sigmas <- rlnorm(n=2,meanlog=c(-1,-1),sdlog=c(1,1))
 
  print(em_mus)
  print(em_sigmas)
  print(trans_mat)
  
  stan_init <- list(list(trans_mat=trans_mat,em_mus=em_mus,em_sigmas=em_sigmas))
  
  stan_fit <- stan(file="HMM_normal_probs.stan",data=stan_data,init=stan_init,chain=1,iter=num_iters)
  return(stan_fit)
}


Plot_Prior_and_Posterior <- function(stan_fit)
{
  mat <- as.matrix(stan_fit)
  
  ## the following were the prior parameters:
  alphas <- c(1000,1)
  betas <- c(1,1000)
  mus_mu <- c(4.7,4.1)
  sigmas_mu <- c(0.6,0.6)
  
  mus_sigma <- c(-1,-1)
  sigmas_sigma <- c(1,1)
  
  ## plot priors for the "mus"
  
  xseq <- seq(from=1,to=10,by=0.02)
  yseq_mu1_prior <- dnorm(xseq,mean=mus_mu[1],sd=sigmas_mu[1])
  yseq_mu2_prior <- dnorm(xseq,mean=mus_mu[2],sd=sigmas_mu[2])
  
  plot(xseq,yseq_mu1_prior)
  points(xseq,yseq_mu2_prior)
  
  ## plot priors for the "sigmas"
  xseq_sigmas <- seq(from=0,to=2.5,by=0.01)
  yseq_sigma1_prior <- dlnorm(xseq_sigmas,meanlog=mus_sigma[1],sdlog=sigmas_sigma[1])
  yseq_sigma2_prior <- dlnorm(xseq_sigmas,meanlog=mus_sigma[2],sdlog=sigmas_sigma[2])
  
  
  
  ## plot some of the prior distributions:
  
  
  
  
}

Make_Nice_Plots <- function(stan_fit)
{
  require("invgamma")
  sampling_rate <- 20000
  mat <- as.matrix(stan_fit)
  
  num_samples <- nrow(mat)
  print(num_samples)
  
  tau1 <- array(0,dim=c(num_samples))
  tau2 <- array(0,dim=c(num_samples))
  pi1 <- array(0,dim=c(num_samples))
  pi2 <- array(0,dim=c(num_samples))
  
  for(i in 1:num_samples)
  {
    pi1[i] <- mat[i,2]/(mat[i,2]+mat[i,3])
    pi2[i] <- mat[i,3]/(mat[i,3]+mat[i,2])
    k_12 <- -log(1-mat[i,3])*sampling_rate;
    k_21 <- -log(1-mat[i,2])*sampling_rate;
    tau1[i] <- 1.0/k_12;
    tau2[i] <- 1.0/k_21;
    
    
  }
  
  #trans_mat <- matrix(0,nrow=2,ncol=2)
  
#  trans_mat[1,1] <- rbeta(1,1000,1)
#  trans_mat[2,1] <- rbeta(1,1,1000)
#  trans_mat[1,2] <- (1.0 - trans_mat[1,1])
#  trans_mat[2,2] <- (1.0 - trans_mat[2,1])
  
  ## Draw from prior distributions:
  
  num_samples_prior <- 40000
  
  trans_mat_prior <- array(0,dim=c(num_samples_prior,2,2))
  pi1_prior <- array(0,dim=c(num_samples_prior))
  pi2_prior <- array(0,dim=c(num_samples_prior))
  tau1_prior <- array(0,dim=c(num_samples_prior))
  tau2_prior <- array(0,dim=c(num_samples_prior))
  
  for(i in 1:num_samples_prior)
  {
    trans_mat_prior[i,1,1] <- rbeta(1,1000,1);
    trans_mat_prior[i,2,1] <- rbeta(1,1,1000);
    trans_mat_prior[i,1,2] <- (1.0 - trans_mat_prior[i,1,1])
    trans_mat_prior[i,2,2] <- (1.0 - trans_mat_prior[i,2,1])
    
    pi1_prior[i] <- trans_mat_prior[i,2,1]/(trans_mat_prior[i,1,2]+trans_mat_prior[i,2,1])
    pi2_prior[i] <- trans_mat_prior[i,1,2]/(trans_mat_prior[i,1,2]+trans_mat_prior[i,2,1])
    
    k_12_prior <- -log(1.0-trans_mat_prior[i,1,2])*sampling_rate;
    k_21_prior <- -log(1.0-trans_mat_prior[i,2,1])*sampling_rate;
    
    
    tau1_prior[i] <- 1.0/k_12_prior;
    tau2_prior[i] <- 1.0/k_21_prior;
  }
 
  breaks_pis <- seq(from=0,to=1.00,by=0.02)
  hist1 <- hist(pi1_prior,freq=FALSE,xlim=c(0,1.00),ylim=c(0,18),breaks=breaks_pis)
  hist2 <- hist(pi1,freq=FALSE,xlim=c(0,0.99),ylim=c(0,18),breaks=breaks_pis)
  hist_pi2 <- hist(pi2,freq=FALSE,xlim=c(0,0.99),ylim=c(0,18),breaks=breaks_pis)
  
  plot(hist1,freq=FALSE,add=T,col=rgb(0,0,0,0.25),xlim=c(0,1.00))
  plot(hist2,freq=FALSE,add=T,col=rgb(1,0,0,0.25),xlim=c(0,0.99))
  plot(hist_pi2,freq=FALSE,add=T,col=rgb(0,0,1,0.25),xlim=c(0,0.99))
  df_beta_params <- HMM_Get_Beta_Distribution_Params_By_Moment_Matching(pi1)
  xs <- seq(from=0,to=1,by=1e-4)
 
  points(xs,dbeta(xs,df_beta_params$alpha,df_beta_params$beta),type="l",lty=3,lwd=2,col="red")
  points(xs,dbeta(xs,df_beta_params$beta,df_beta_params$alpha),type="l",lty=3,lwd=2,col="blue")
  ## do some moment matching to fit beta-distribution....
  
  breaks_taus <- seq(from=0,to=1,by=0.005)
  breaks_tau1_prior <- seq(from=0,to=1000,by=0.005)
  
  hist_tau1 <- hist(tau1,freq=FALSE,xlim=c(0,0.2),ylim=c(0,180),breaks=breaks_taus)
  hist_tau2 <- hist(tau2,freq=FALSE,xlim=c(0,0.2),ylim=c(0,180),breaks=breaks_taus)
#  hist_tau1_prior <- hist(tau1_prior[tau1_prior<=1000],freq=FALSE,breaks=breaks_tau1_prior,xlim=c(0,0.2),ylim=c(0,40));
  
  xsequ_tau <- seq(from=0,to=0.2,by=0.0005)
  df_gamma_params <- HMM_Get_Gamma_Distribution_Params_By_Moment_Matching(tau1)
  df_inv_gamma_params_tau1 <- HMM_Get_Inv_Gamma_Distribution_Params_By_Moment_Matching(tau1)
  df_inv_gamma_params_tau2 <- HMM_Get_Inv_Gamma_Distribution_Params_By_Moment_Matching(tau2)
  df_inv_gamma_params_tau1_prior <- HMM_Get_Inv_Gamma_Distribution_Params_By_Moment_Matching(tau1_prior[tau1_prior<0.2])
  
  print(df_inv_gamma_params_tau1)
  plot(hist_tau1,freq=FALSE,add=T,col=rgb(1,0,0,0.25),xlim=c(0,0.2),ylim=c(0,80))
  plot(hist_tau2,freq=FALSE,add=T,col=rgb(0,0,1,0.25),xlim=c(0,0.2),ylim=c(0,80))
 # plot(hist_tau1_prior,freq=FALSE,col=rgb(0,0,0,0.25),add=T)
  #points(xsequ_tau,dgamma(xsequ_tau,shape=df_gamma_params$k_mm,scale=df_gamma_params$theta_mm),type="l",lty=3,lwd=2)
  points(xsequ_tau,dinvgamma(xsequ_tau,shape=df_inv_gamma_params_tau1$alpha_mm,rate=df_inv_gamma_params_tau1$beta_mm),type="l",lty=3,lwd=2,col="red")
  points(xsequ_tau,dinvgamma(xsequ_tau,shape=df_inv_gamma_params_tau2$alpha_mm,rate=df_inv_gamma_params_tau2$beta_mm),type="l",lty=3,lwd=2,col="blue")
  
 # points(xsequ_tau,dinvgamma(xsequ_tau,shape=df_inv_gamma_params_tau1_prior$alpha_mm,rate=df_inv_gamma_params_tau1_prior$beta_mm),type="l",lty=3,lwd=2,col="black")
  
  plot(tau1_prior[tau1_prior<1000])
  
  
  print(df_gamma_params)
  
  
  #hist(pi1_prior)
  #hist(pi1,freq=FALSE)
  

  
  ## first: draw from prior - distribution:
  
  ## then from posterior:
  
  
}
## analysis with 400 steps started at 18:35.... (less than 3 minutes !!!!)
#
#data
#{
#  int<lower=1> N; // number of data points
#  int<lower=2> num_states; // number of states, for this particular example it has to be 2...
#  real obs_sequ[N]; // the observations sequence
  
#  // hyperparameters for the 
#  // rows of the transition matrix
#  real alphas[num_states];
#  real betas[num_states];
  
#  // hyperparameters for the means of
#  // the emission distributions: 
#    real mus_mu[num_states];
#  real sigmas_mu[num_states];
  
#  // hyperparameters for the (log-normally distributed) 
#  // standard-deviations
#  // of the emission distribution
#  real mus_sigma[num_states];
#  real sigmas_sigma[num_states];
#}

#parameters
#{
#  simplex[num_states] trans_mat[num_states];
#  real em_mus[num_states];
#  real em_sigmas[num_states];
#}