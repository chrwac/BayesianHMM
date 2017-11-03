functions
{
  real LogLikObsSequ(real[] obs_sequ,int num_pts,int num_states,vector init_probs,real[] em_mus,real[] em_sigmas,vector[] trans_matrix)
  {
    real log_alphas[num_pts,num_states];
    
    // to speed things up: precalculate the logarithms of the transition matrix:
    real log_trans[num_states,num_states];
    
    real log_term[num_states];
    real log_prob;
    
    for(i in 1:num_states)
    {
      for(j in 1:num_states)
      {
        log_trans[i,j] = log(trans_matrix[i,j]);
      }
    }
    
    // initialization
    for(i in 1:num_states)
    {
      log_alphas[1,i] = log(init_probs[i]) + normal_log(obs_sequ[1],em_mus[i],em_sigmas[i]);
    }
    // induction
    for(t in 1:(num_pts-1))
    {
      for(j in 1:num_states)
      {
        for(i in 1:num_states)
        {
          log_term[i] = log_alphas[t,i]+ log_trans[i,j];
        }
        log_alphas[t+1,j] = normal_log(obs_sequ[t+1],em_mus[j],em_sigmas[j]) + log_sum_exp(log_term);
      }
    }
    // termination (computation of the log-probability)
    log_prob = log_sum_exp(log_alphas[num_pts,]);
    return log_prob;
  }
}

data
{
  int<lower=1> N; // number of data points
  int<lower=2> num_states; // number of states, for this particular example it has to be 2...
  real obs_sequ[N]; // the observations sequence
  
  // hyperparameters for the 
  // rows of the transition matrix
  real alphas[num_states];
  real betas[num_states];
  
  // hyperparameters for the means of
  // the emission distributions: 
  real mus_mu[num_states];
  real sigmas_mu[num_states];
  
  // hyperparameters for the (log-normally distributed) 
  // standard-deviations
  // of the emission distribution
  real mus_sigma[num_states];
  real sigmas_sigma[num_states];
}

parameters
{
  simplex[num_states] trans_mat[num_states];
  real em_mus[num_states];
  real em_sigmas[num_states];
}

model
{
  vector[num_states] pis;
  vector[num_states] t_trans_mat[num_states];
  for(i in 1:num_states)
  {
    //trans_mat[i] ~ beta(alphas[i],betas[i]); // beta-distributed row-elements of the (2x2) transition matrix
    em_mus[i] ~ normal(mus_mu[i],sigmas_mu[i]); // normally distributed mean-values of the emissions
    em_sigmas[i] ~ lognormal(mus_sigma[i],sigmas_sigma[i]); // log-normally distributed standard-deviation of the emissions
  }
  trans_mat[1,1] ~ beta(alphas[1],betas[1]);
  trans_mat[2,1] ~ beta(alphas[2],betas[2]);
  
  t_trans_mat[1,1] = trans_mat[1,1];
  t_trans_mat[1,2] = (1.0 - trans_mat[1,1]);
  t_trans_mat[2,1] = trans_mat[2,1];
  t_trans_mat[2,2] = (1.0 - trans_mat[2,1]);
  
  
  // use detailed balance to express the initial state probabilities as
  // a function of the transition matrix elements:
  pis[1] = trans_mat[2,1]/(trans_mat[1,2]+trans_mat[2,1]);
  pis[2] = trans_mat[1,2]/(trans_mat[1,2]+trans_mat[2,1]);
  
  // increment the logarithm of the posterior probability by the log-likelihood:
  target+=LogLikObsSequ(obs_sequ,N,num_states,pis,em_mus,em_sigmas,t_trans_mat);
}

generated quantities
{
  real p1;
  p1 <- trans_mat[2,1]/(trans_mat[1,2]+trans_mat[2,1]);
}