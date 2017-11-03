functions
{
  real LogProbObsSequ(real[] obs_sequ,int num_pts,int num_states,vector init_probs,real[] em_locations,real[] em_scales,real[] em_shapes,vector[] trans_matrix)
  {
    real log_alphas[num_pts,num_states];
    real log_term[num_states];
    real log_prob;
    
    // initialization
    for(i in 1:num_states)
    {
      log_alphas[1,i] = log(init_probs[i]) + skew_normal_log(obs_sequ[1],em_locations[i],em_scales[i],em_shapes[i]);
    }
    // induction
    for(t in 1:(num_pts-1))
    {
      for(j in 1:num_states)
      {
        for(i in 1:num_states)
        {
          log_term[i] = log_alphas[t,i]+log(trans_matrix[i,j]);
        }
        log_alphas[t+1,j] = skew_normal_log(obs_sequ[t+1],em_locations[j],em_scales[j],em_shapes[j]) + log_sum_exp(log_term);
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
  int<lower=2> num_states; // number of states, has to be 2 at the moment
  real obs_sequ[N];
  vector<lower=0>[num_states] alphas;
  vector<lower=0>[num_states] trans0[num_states];
 
  real location0s[num_states];
  real<lower=0> scale0s[num_states];
  real shape0s[num_states];
}

parameters
{
  simplex[num_states] pis;
  simplex[num_states] trans_mat[num_states];
}

transformed parameters
{
  real locations[num_states];
  real scales[num_states];
  real shapes[num_states];
  for(i in 1:num_states)
  {
    locations[i] = location0s[i];
    scales[i] = scale0s[i];
    shapes[i] = shape0s[i];
  }
}

model
{
  pis ~ dirichlet(alphas);
  for(i in 1:num_states)
  {
    trans_mat[i] ~ dirichlet(trans0[i,]);
  }
  //increment_log_prob(LogProbObsSequ(obs_sequ,N,num_states,pis,locations,scales,shapes,trans_mat));
  target+=LogProbObsSequ(obs_sequ,N,num_states,pis,locations,scales,shapes,trans_mat);
}

generated quantities
{
  real p1;
  p1 <- trans_mat[2,1]/(trans_mat[1,2]+trans_mat[2,1]);
}

