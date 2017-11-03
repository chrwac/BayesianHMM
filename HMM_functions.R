## data loading, EM-learning, Stan-sampling, Viterbi-assignment, max_state_prob - assignment
## for Skewed-Normal Hidden Markov Models

## HMM_Stan_Start_Job_By_ID(filename_js_file,ID)
HMM_Stan_Start_Job_By_ID <- function(filename_js_file,ID)
{
  require("rstan")
  job_specs <- HMM_Read_Job_Specifications(filename_js_file)
  print("ID: ")
  print(ID)
  print("job specs: ")
  print(job_specs)
  print("loading file: ")
  print(toString(job_specs$filename[ID]))
  data_tb_analyzed <- HMM_Load_Data(toString(job_specs$filename[ID]),job_specs$data_thinning_factor[ID])
  
  ## hyperparameters for the initial state probabilities and the state transition matrix elements...
  hps_init_probs <- array(0,dim=c(job_specs$num_states[ID]))
  hps_trans_mat <- matrix(nrow=job_specs$num_states[ID],ncol=job_specs$num_states[ID])
  ## parameters for the 
  location_params <- array(0,dim=c(job_specs$num_states[ID]))
  scale_params <- array(0,dim=c(job_specs$num_states[ID]))
  shape_params <- array(0,dim=c(job_specs$num_states[ID]))

  curr_index <- 6 + 4*job_specs$num_states[ID]+1
  for(i in 1:job_specs$num_states[ID])
  {
    location_params[i] <- job_specs[ID,6+i]
    scale_params[i] <- job_specs[ID,6+1*job_specs$num_states[ID]+i]
    shape_params[i] <- job_specs[ID,6+2*job_specs$num_states[ID]+i]
    hps_init_probs[i] <- job_specs[ID,6+3*job_specs$num_states[ID]+i]
    
    
    for(j in 1:job_specs$num_states[ID])
    {
      hps_trans_mat[i,j] <- job_specs[ID,curr_index]
       curr_index<-curr_index+1;
    }
  }
  print("hyperparameters for the initial state probabilities: ")
  print(hps_init_probs)
  print("hyperparameters for the transition matrix elements: ")
  print(hps_trans_mat)
  
  init_values_stan_hmm <- HMM_Sample_Initial_Values_From_Prior(hyperparams_init_probs=hps_init_probs,hyperparams_trans_matrix=hps_trans_mat,num_states=job_specs$num_states[ID],num_chains=job_specs$num_chains[ID])
  print("initial values for Stan HMM inference: ")
  print(init_values_stan_hmm)
  data_stan_hmm <- list(N=length(data_tb_analyzed),num_states=job_specs$num_states[ID],obs_sequ=data_tb_analyzed,alphas=hps_init_probs,trans0=hps_trans_mat,location0s=location_params,scale0s=scale_params,shape0s=shape_params)
  print("starting stan-sampling: ")
  fit_stan_hmm_sn_emprobs <- stan(file="HMM_SN_emprobs.stan",data=data_stan_hmm,init=init_values_stan_hmm,iter=job_specs$num_total_samples[ID],chains=job_specs$num_chains[ID],thin=job_specs$sample_thinning_factor[ID])
 
  mat_stan_hmm_sn_samples <- as.matrix(fit_stan_hmm_sn_emprobs)
  splitted <- strsplit(x=toString(job_specs$filename[ID]),split=".",fixed=TRUE)
  casename <- splitted[[1]][1]
  filename_hmm_sn_samples <- paste0(casename,"_hmm_sn_samples.csv")
  write.table(x=mat_stan_hmm_sn_samples,file=filename_hmm_sn_samples,sep=";")
  HMM_SN_Stan_Evaluate_Samples(mat_stan_hmm_sn_samples=mat_stan_hmm_sn_samples,obs_seq=data_tb_analyzed,num_states=job_specs$num_states[ID],casename=casename)
  return(fit_stan_hmm_sn_emprobs)
}

## HMM_SN_Stan_Evaluate_Samples(mat_stan_hmm_sn_samples,obs_seq,num_states=2,casename="test_case")
HMM_SN_Stan_Evaluate_Samples <- function(mat_stan_hmm_sn_samples,obs_seq,num_states=2,casename="test_case")
{
  ## calculate viterbi sequence assignment and save it:
  print("assigning viterbi sequence: ")
  vss <- HMM_SN_Assign_Viterbi_Sequence_From_Mean_Posterior(obs_seq,mat_stan_hmm_sn_samples,num_states)
  df_viterbi <- data.frame(x=obs_seq,vss=vss)
  filename_hmm_sn_vss_plot <-  paste0(casename,"_hmm_sn_vss.pdf")
  filename_hmm_sn_vss <- paste0(casename,"_hmm_sn_vss.csv")
  write.table(x=df_viterbi,file=filename_hmm_sn_vss)
  pdf(file=filename_hmm_sn_vss_plot,width=14,height=11,useDingbats=FALSE)
  HMM_Plot_Trace_With_State_Assignment(obs_seq=obs_seq,state_seq=vss,dt=1,lab_x="index",lab_y="deflection [nm]",case_str = casename)
  dev.off()
  ## calculate the most probable state assignment (point-wise!, not Viterbi) and save it
  print("assigning most probable state (point-wise) sequence: ")
  mss <- HMM_SN_Assign_Most_Probable_States_From_Mean_Posterior(obs_seq,mat_stan_hmm_sn_samples,num_states)
  df_mss <- data.frame(x=obs_seq,mss=mss) 
  filename_hmm_sn_mss_plot <- paste0(casename,"_hmm_sn_mss.pdf")
  filename_hmm_sn_mss <- paste0(casename,"_hmm_sn_mss.csv")
  write.table(x=df_mss,file=filename_hmm_sn_mss)
  pdf(file=filename_hmm_sn_mss_plot,width=14,height=11,useDingbats=FALSE)
  HMM_Plot_Trace_With_State_Assignment(obs_seq=obs_seq,state_seq=mss,dt=1,lab_x="index",lab_y="deflection [nm]",case_str = casename)
  dev.off()
  filename_state1_ss_prob_plot <- paste0(casename,"_hmm_sn_s1_ssp.pdf")
  filename_state1_ss_prob <- paste0(casename,"_hmm_sn_s1_ssp.csv")
  pdf(file=filename_state1_ss_prob_plot,width=14,height=11,useDingbats=FALSE)
  df_descr_state1 <- HMM_SN_Plot_Steady_State_Probability(mat_stan_hmm_sn_samples,num_states,1)
  dev.off()
  ## write table with information(description of the beta-distribution)
  write.table(x=df_descr_state1,file=filename_state1_ss_prob)
}  

## Stan_Sample_HMM_Load_Data(filename,thinning_factor)
HMM_Load_Data <- function(filename,thinning_factor=1)
{
  data_hmm <- read.table(file=filename,header=TRUE)
  data_hmm_thinned <- data_hmm[seq(from=1,to=length(data_hmm[,1]),by=thinning_factor),1]
  return(data_hmm_thinned)
}

## HMM_Set_Initial_Values(filename,thinning_factor)
HMM_Set_Initial_Values <- function(init_state_probs,init_trans_matrix)
{
  init_ll <- list(list(pis=init_state_probs,trans_mat=trans_matrix))
  return(init_ll)
}

## HMM_Read_Job_Specifications(filename_case_file)
HMM_Read_Job_Specifications <- function(filename_js_file)
{
  #require("gdata")
  job_specifications <- read.csv(file=filename_js_file,header=TRUE,sep=";")
  return(job_specifications)
}

## HMM_Sample_Initial_Values_From_Prior(hyperparams_init_probs,hyperparams_trans_matrix,num_states,num_chains)
## samples Stan-input for the initial values from the specified prior distributions (initial state probabilities and a transition matrix)
## the initial state probability - vector is sampled from a Dirichlet distribution with hyperparameters "hyperparams_init_probs"
## each row of the (state) transition matrix is sampled from a Dirichtlet distribution. The hyperparameters for a particular
## row of the (state) transition matrix are the entries in the corresponding row of the "hyperparams_trans_matrix" 
## a list of "num_chains" lists is created and returned
HMM_Sample_Initial_Values_From_Prior <- function(hyperparams_init_probs,hyperparams_trans_matrix,num_states,num_chains)
{
  require("gtools")
  init_list_hmm_sn_emprobs <- list()
  trans_mat_sampled <- matrix(nrow=num_states,ncol=num_states)
  for(i in 1:num_chains)
  {
    pis_sampled <- as.vector(rdirichlet(n=1,hyperparams_init_probs))
    for(j in 1:num_states)
    {
      trans_mat_sampled[j,] <- as.vector(rdirichlet(n=1,hyperparams_trans_matrix[j,]))
    }
    init_list_hmm_sn_emprobs[[i]] <- list(pis=pis_sampled,trans_mat=trans_mat_sampled)
  }
  return(init_list_hmm_sn_emprobs)
}

## AssignViterbiSequenceFromMeanPosterior(obs_seq,stan_hmm_fit)
## input: an observation sequence "obs_seq" and a stan-fit-object
## output: the Viterbi sequence for the sequence computed
## with the mean posterior values
HMM_SN_Assign_Viterbi_Sequence_From_Mean_Posterior <- function(obs_seq,mat_hmm_sn_stan_samples,num_states)
{
  pi_means <- matrix(nrow=num_states,ncol=1);
  trans_mat_means <- matrix(nrow=num_states,ncol=num_states);
  location_means <- matrix(nrow=num_states,ncol=1);
  scale_means <- matrix(nrow=num_states,ncol=1);
  shape_means <- matrix(nrow=num_states,ncol=1);
  
  for(i in 1:num_states)
  {
    pi_means[i] <- mean(mat_hmm_sn_stan_samples[,i]);
  }
  curr_index <- num_states+1;
  for(i in 1:num_states)
  {
    for(j in 1:num_states)
    {
      trans_mat_means[j,i] <- mean(mat_hmm_sn_stan_samples[,curr_index]);
      curr_index <- curr_index+1;
    }
  }
  for(i in 1:num_states)
  {
    location_means[i] <- mean(mat_hmm_sn_stan_samples[,num_states+num_states^2 + i])
    scale_means[i] <- mean(mat_hmm_sn_stan_samples[,2*num_states+num_states^2 + i])
    shape_means[i] <-  mean(mat_hmm_sn_stan_samples[,3*num_states+num_states^2 + i])
  }
  
  print(pi_means);
  print(trans_mat_means);
  print(location_means);
  print(scale_means);
  print(shape_means);
  viterbi_seq <- HMM_SN_Assign_Viterbi_Sequence(obs_seq,num_states,pi_means,location_means,scale_means,shape_means,trans_mat_means);
  return(viterbi_seq);
}

## AssignViterbiSequence(obs_seq, num_states,init_probs,em_locations,em_sds,transition_matrix)
## input: an observation sequence "obs_seq", the number of states "num_states", initial probabilities
## "init_probs", the number of state "num_states", the initial probabilities "init_probs", and
## parameters describing the emission probability distribution, in this case "em_locations" und "em_sds"
#AssignViterbiSequence <- function(obs_seq,num_states,init_probs,em_locations,em_sds,em_shapes,transition_matrix)
HMM_SN_Assign_Viterbi_Sequence <- function(obs_seq,num_states,init_probs,em_locations,em_scales,em_shapes,transition_matrix)
{
  npnts <- length(obs_seq);
  log_delta <- matrix(nrow=npnts,ncol=num_states);
  psi <- matrix(nrow=npnts,ncol=num_states);
  
  ## initialization:
  for(i in 1:num_states)
  {
    log_delta[1,i] <- log(init_probs[i]) + HMM_Return_Log_Em_Prob(x=obs_seq[1],c(em_locations[i],em_scales[i],em_shapes[i]));
    psi[1,i] <- 0.0;
  }
  for(t in 2:npnts)
  {
    for(j in 1:num_states)
    {
      log_max <- log_delta[t-1,1] + log(transition_matrix[1,j]);
      arg_max <- 1;
      for(i in 2:num_states)
      {
        log_temp <- log_delta[t-1,i] + log(transition_matrix[i,j]);
        if(log_temp > log_max)
        {
          log_max <- log_temp;
          arg_max <- i;
        }
      }
      log_delta[t,j] <- log_max + HMM_Return_Log_Em_Prob(x=obs_seq[t],c(em_locations[j],em_scales[j],em_shapes[j]));
      psi[t,j] <- arg_max;
    }
  }
  
  ## backtracking
  q <- matrix(nrow=npnts,ncol=1);
  argmax <- 1;
  log_p_max <- log_delta[npnts,1];
  for(i in 2:num_states)
  {
    if(log_delta[npnts,i]>log_p_max)
    {
      log_p_max <- log_delta[npnts,i];
      argmax <- i;
    }
  }
  q[npnts] <- argmax;
  for(t in (npnts-1):1)
  {
    q[t] <- psi[t+1,q[t+1]];
  }
  return(q);
}

## for now, quick'n dirty stuff...
HMM_Norm_Assign_Viterbi_Sequence_From_MP <- function(stan_fit,obs_sequ)
{
  mat <- as.matrix(stan_fit)
  num_states <- 2 ## for now...
  trans_mat <- matrix(0,nrow=num_states,ncol=num_states)
  trans_mat[1,1] <- mean(mat[,1])
  trans_mat[2,1] <- mean(mat[,2])
  trans_mat[1,2] <- mean(mat[,3])
  trans_mat[2,2] <- mean(mat[,4])
  
  pi1 <- trans_mat[2,1]/(trans_mat[1,2]+trans_mat[2,1])
  pi2 <- trans_mat[1,2]/(trans_mat[1,2]+trans_mat[2,1])
  init_probs <- c(pi1,pi2)
  mus <- c(mean(mat[,5]),mean(mat[,6]))
  sigmas <- c(mean(mat[,7]),mean(mat[,8]))
  ssequ <- HMM_Norm_Assign_Viterbi_Sequence(obs_sequ,2,init_probs,mus,sigmas,trans_mat)
  return(ssequ)
}

## HMM_SN_Assign_Viterbi_Sequence_Norm(obs_seq, num_states,init_probs,em_locations,em_sds,transition_matrix)
## input: an observation sequence "obs_seq", the number of states "num_states", initial probabilities
## "init_probs", the number of state "num_states", the initial probabilities "init_probs", and
## parameters describing the emission probability distribution, in this case "em_locations" und "em_sds"
#AssignViterbiSequence <- function(obs_seq,num_states,init_probs,em_locations,em_sds,em_shapes,transition_matrix)
## difference to the version above is, that the emission probabilities are assumed to be distributed normally....
HMM_Norm_Assign_Viterbi_Sequence <- function(obs_seq,num_states,init_probs,mus,sigmas,transition_matrix)
{
  npnts <- length(obs_seq);
  log_delta <- matrix(nrow=npnts,ncol=num_states);
  psi <- matrix(nrow=npnts,ncol=num_states);
  
  ## initialization:
  for(i in 1:num_states)
  {
    log_delta[1,i] <- log(init_probs[i]) + HMM_Return_Log_Em_Prob_Norm(x=obs_seq[1],c(mus[i],sigmas[i]));
    #print(i)
    #print(log_delta[1,i])
    psi[1,i] <- 0.0;
  }
  for(t in 2:npnts)
  {
    for(j in 1:num_states)
    {
      log_max <- log_delta[t-1,1] + log(transition_matrix[1,j]);
      arg_max <- 1;
      
     # print(log_max)
      i<- 2
      
      #while(i<=num_states)
      
      for(i in 2:num_states)
      {
       # print(t)
        #print(i)
        #print(log_delta[t-1,i])
        log_temp <- log_delta[t-1,i] + log(transition_matrix[i,j]);
        #print(log_temp)
        if(log_temp > log_max)
        {
          log_max <- log_temp;
          arg_max <- i;
        }
        i <- i+1
      }
      log_delta[t,j] <- log_max + HMM_Return_Log_Em_Prob_Norm(x=obs_seq[t],c(mus[j],sigmas[j]));
      psi[t,j] <- arg_max;
    }
  }
  
  ## backtracking
  q <- matrix(nrow=npnts,ncol=1);
  argmax <- 1;
  log_p_max <- log_delta[npnts,1];
  for(i in 2:num_states)
  {
    if(log_delta[npnts,i]>log_p_max)
    {
      log_p_max <- log_delta[npnts,i];
      argmax <- i;
    }
  }
  q[npnts] <- argmax;
  for(t in (npnts-1):1)
  {
    q[t] <- psi[t+1,q[t+1]];
  }
  return(q);
}

## PlotTraceWithStateAssignment(obs_seq,state_seq,dt)
## input: an observation-sequence "obs_seq", and a state-sequence
## "state_seq" (which could e.g. be inferred via the Viterbi-algorithm)
## output: creates a plot with the observation sequence colored according
## to the state assignment
HMM_Plot_Trace_With_State_Assignment <- function(obs_seq,state_seq,dt,lab_x="time [s]",lab_y="deflection [nm]",case_str)
{
  time_seq <- seq(from=0,to=(length(obs_seq)-1)*dt,by=dt);
  states = unique(state_seq);
  minval <- min(obs_seq);
  maxval <- max(obs_seq);
  
  amp <- (maxval-minval)
  minval <- minval - 0.1*amp
  maxval <- maxval + 0.1*amp
  color_table <- c("blue","red","green","yellow","black");
  plot(time_seq[state_seq==states[1]],obs_seq[state_seq==states[1]],col=color_table[1],ylim=c(minval,maxval),xlab=lab_x,ylab=lab_y,pch=".")
  for(i in 2:length(states))
  {
    points(time_seq[state_seq==states[i]],obs_seq[state_seq==states[i]],col=color_table[i],pch=".")
  }
  #filename_graph <- paste(case_string,"_",assignment_type,"_StateAssignment.pdf")
  #dev.copy2pdf(file=filename_graph,width=11,height=7,useDingbats=FALSE);
}

## HMM_Return_Log_Em_Prob(x,params)
## input: a position "x" and parameters describing
## the emission probability distribution "params"
HMM_Return_Log_Em_Prob <- function(x,params)
{
  require("sn");
  ##return(dnorm(x,mean=params[1],sd=params[2],log=TRUE));
  return(dsn(x, xi=params[1], omega=params[2], alpha=params[3], log=TRUE));
}

## HMM_Return_Log_Em_Prob_Norm(x,params)
## input: a position "x" and parameters describing
## the emission probability distribution "params"
HMM_Return_Log_Em_Prob_Norm <- function(x,params)
{
  return(dnorm(x,mean=params[1],sd=params[2],log=TRUE));
}

## HMM_SN_Assign_Most_Probable_States_From_Mean_Posterior(obs_seq,mat_hmm_sn_stan_samples,num_states)
## Calculate gamma-values and
## assign the states according to
## the maximum gamma value
## for each individual time-point...
HMM_SN_Assign_Most_Probable_States_From_Mean_Posterior <- function(obs_seq,mat_hmm_sn_stan_samples,num_states)
{
  pi_means <- matrix(nrow=num_states,ncol=1);
  trans_mat_means <- matrix(nrow=num_states,ncol=num_states);
  location_means <- matrix(nrow=num_states,ncol=1);
  scale_means <- matrix(nrow=num_states,ncol=1);
  shape_means <- matrix(nrow=num_states,ncol=1);
  
  for(i in 1:num_states)
  {
    pi_means[i] <- mean(mat_hmm_sn_stan_samples[,i]);
  }
  curr_index <- num_states+1;
  for(i in 1:num_states)
  {
    for(j in 1:num_states)
    {
      trans_mat_means[j,i] <- mean(mat_hmm_sn_stan_samples[,curr_index]);
      curr_index <- curr_index+1;
    }
  }
  for(i in 1:num_states)
  {
    location_means[i] <- mean(mat_hmm_sn_stan_samples[,num_states+num_states^2 + i])
    scale_means[i] <- mean(mat_hmm_sn_stan_samples[,2*num_states+num_states^2 + i])
    shape_means[i] <-  mean(mat_hmm_sn_stan_samples[,3*num_states+num_states^2 + i])
  }
  
  print(pi_means);
  print(trans_mat_means);
  print(location_means);
  print(scale_means);
  print(shape_means);
  
  mps <- HMM_SN_Assign_Most_Probable_States(obs_seq,num_states,pi_means,location_means,scale_means,shape_means,trans_mat_means);
  return(mps);
}

## HMM_SN_Assign_Most_Probable_States(obs_seq,num_states,init_probs,em_locations,em_scales,em_shapes,trans_matrix)
## calculates a state sequence by assigning at each data point
## the most probable state......(not Viterbi!)
HMM_SN_Assign_Most_Probable_States <- function(obs_seq,num_states,init_probs,em_locations,em_scales,em_shapes,trans_matrix)
{
  require("matrixStats");
  npnts <- length(obs_seq);
  mps <- matrix(nrow=npnts,ncol=1);
  log_alphas <- HMM_SN_Calculate_Log_Alphas(obs_seq,npnts,num_states,init_probs,em_locations,em_scales,em_shapes,trans_matrix);
  log_betas <- HMM_SN_Calculate_Log_Betas(obs_seq,npnts,num_states,init_probs,em_locations,em_scales,em_shapes,trans_matrix);
  
  ## loop over all data points
  for(t in (1:npnts))
  {
    argmax <- 1;
    log_max <- log_alphas[t,1] + log_betas[t,1];
    for(i in (2:num_states))
    {
      log_curr <- log_alphas[t,i] + log_betas[t,i];
      if(log_curr>log_max)
      {
        log_max <- log_curr;
        argmax <- i;
      }
    }
    mps[t] <- argmax;
  }
  return(mps);
}

## CalculateLogAlphas(obs_sequ,num_ptscons,num_states,init_probs,em_means,em_sds,trans_matrix)
## input: observation sequence "obs_sequ", the number of points to be considered "num_ptscons",
## the number of states "num_states", the initial state probabilities "init_probs", the 
## mean values of the (assumed to be normal) emission probabilities "em_means" and the standard
## deviations of the emission probabilities "em_sds", as well as the transition matrix
## "trans_matrix"
HMM_SN_Calculate_Log_Alphas <- function(obs_sequ,num_ptscons,num_states,init_probs,em_locations,em_scales,em_shapes,trans_matrix)
{
  require("matrixStats");
  log_alphas <- matrix(nrow=num_ptscons,ncol=num_states)
  for(i in 1:num_states)
  {
    log_alphas[1,i] <- log(init_probs[i])+HMM_Return_Log_Em_Prob(x=obs_sequ[1],c(em_locations[i],em_scales[i],em_shapes[i]));
  }
  temp_log_sums <- matrix(nrow=num_states,ncol=1);
  for(t in 1:(num_ptscons-1))
  {
    for(j in 1:num_states)
    {
      for(i in 1:num_states)
      {
        temp_log_sums[i] <- log_alphas[t,i] + log(trans_matrix[i,j]);
      }
      log_alphas[t+1,j] <- HMM_Return_Log_Em_Prob(x=obs_sequ[t+1],c(em_locations[i],em_scales[i],em_shapes[i])) + logSumExp(temp_log_sums);
    }
  }
  return(log_alphas)
}

## HMM_SN_Calculate_Log_Betas(obs_seq,num_ptscons,num_states,init_probs,em_locations,em_scales,em_shapes,trans_matrix)
## input: the observation sequence "obs_sequ", the number of points that shall be considered "num_ptscons"
## the initial state probabilities "init_probs", the mean values of the gaussian emission probabilities
## "em_means", the corresponding standard deviations "em_sds" and the transition matrix "trans_matrix"
HMM_SN_Calculate_Log_Betas <- function(obs_seq,num_ptscons,num_states,init_probs,em_locations,em_scales,em_shapes,trans_matrix)
{
  require("matrixStats");
  log_betas <- matrix(nrow=num_ptscons,ncol=num_states)
  for(i in 1:num_states)
  {
    log_betas[num_ptscons,i] <- 0;
  }
  temp_log_sums <- matrix(nrow=num_states,ncol=1);
  for(t in (num_ptscons-1):1)
  {
    for(i in 1:num_states)
    {
      for(j in 1:num_states)
      {
        temp_log_sums[j] <- log(trans_matrix[i,j]) + HMM_Return_Log_Em_Prob(x=obs_seq[t+1],c(em_locations[j],em_scales[j],em_shapes[j])) + log_betas[t+1,j];
      }
      log_betas[t,i] <- logSumExp(temp_log_sums);
    }
  }
  return(log_betas);
}

## HMM_Calculate_Steady_State_Probs_N2_Simple(trans_matrix)
## Calculates the Steady State Probabilities for a given transition matrix
## "trans_matrix", for a 2-state system
HMM_Calculate_Steady_State_Probs_N2_Simple <- function(trans_matrix)
{
  if((nrow(trans_matrix)!=ncol(trans_matrix)) || (nrow(trans_matrix)!=2))
  {
    return;
  }
  else
  {
    x <- matrix(0,nrow=2,ncol=1)
    x[1] <- trans_matrix[2,1]/(trans_matrix[1,2]+trans_matrix[2,1])
    x[2] <- 1.0 - x[1]
    return(x)
  }
}

## HMM_Calculate_Steady_State_Probs_N2(trans_matrix)
## calculates the solution of detailed balance equation
## for a 2-state transition matrix
## so, basically solves the two equations (trivial:
## P_1 w_12 = P_2 w_21
## P_1 + P_2 = 1
HMM_Calculate_Steady_State_Probs_N2 <- function(trans_matrix)
{
  if((nrow(trans_matrix)!=ncol(trans_matrix)) || (nrow(trans_matrix)!=2))
  {
    return;
  }
  else
  {
    A <- matrix(1,nrow=2,ncol=2)
    b <- matrix(0,nrow=2,ncol=1)
    A[1,1] = trans_matrix[1,2]
    A[1,2] = -trans_matrix[2,1]
    b[2] <- 1;
    x <- solve(A,b);
    return(x)
  }
}

## HMM_Calculate_Steady_State_Probs(trans_matrix)
## solves the detailed balance equations (APPROXIMATELX) for the
## given transition matrix "trans_matrix"
## use the function HMM_Calculate_Steady_State_Probs_N2(trans_matrix)
## in the case of num_states=2
## example for num_states=3: 
## (I) P_1 w_12 = P_2 w_21
## (II) P_1 w_13 = P_3 w_31
## (III) P_2 w_23 = P_3 w_32
## (IV) P_1 + P_2 + P_3 = 1
## can be expressed in the form: Ax = b and solved (approximatelx !, it is overdetermined)
## by using the "Moore-Penrose-Inverse / Pseudo-inverse"
HMM_Calculate_Steady_State_Probs <- function(trans_matrix)
{
  require("MASS")
  num_states <- nrow(trans_matrix)
  num_equations <- ((num_states-1)*num_states)/2 + 1
 
  A <- matrix(0,nrow=num_equations,ncol=num_states)
  b <- matrix(0,nrow=num_equations,ncol=1)
  A[num_equations,] <- matrix(1,nrow=1,ncol=num_states)
  b[num_equations,1] <- 1

  ## initialize values of 
  curr_index <- 1
  for(i in 1:num_states)
  {
    for(j in (i+1):num_states)
    {
      if(j<=num_states && (j>i))
      {
        A[curr_index,i] = trans_matrix[i,j]
        A[curr_index,j] = -trans_matrix[j,i]
        curr_index <- curr_index +1
      }
    }
  }
  Ag <- ginv(A);
  sol <- Ag%*%b
  exp_b <- A%*%sol
  list_sol <- list(sol=sol,exp_b=exp_b)
  return(list_sol)
}

## HMM_SN_Make_Steady_State_Probability_Plot(mat_hmm_sn_stan_samples, index_state)
##
HMM_SN_Plot_Steady_State_Probability <- function(mat_hmm_sn_stan_samples,num_states,index_state)
{
  ## obtain vector of steady-state probabilites from the samples:
  steady_state_prob_vectors <- HMM_SN_Calculate_Steady_State_Probs_From_Samples(mat_hmm_sn_stan_samples,num_states)
  ## assign a beta distribution and create a histogram for each state:
  ## Obtain parameters of beta distribution by "moment matching"
  bd_params <- HMM_Get_Beta_Distribution_Params_By_Moment_Matching(steady_state_prob_vectors[,index_state])
  hist(steady_state_prob_vectors[,index_state],breaks=6,freq=FALSE,col="blue",xlab="probability of state occupation",ylab="probability density",main="occup. probability")
    
  x_minimum = min(steady_state_prob_vectors[,index_state])
  x_maximum = max(steady_state_prob_vectors[,index_state])
  dx = (x_maximum-x_minimum)
  ## Plot beta-distribution and normal distribution:  
  beta_plot <- Plot_Beta_Distribution(bd_params,num_points=500)
  norm_plot <- Plot_Normal_Distribution(mu=mean(steady_state_prob_vectors[,index_state]),sigma=sd(steady_state_prob_vectors[,index_state]),x_min=x_minimum-0.2*dx,x_max=x_maximum+0.2*dx,num_points=200)
  points(beta_plot$x_beta_dist,beta_plot$y_beta_dist,type="l",lwd=2,col="red")
  points(norm_plot$x_norm_dist,norm_plot$y_norm_dist,type="l",lwd=2,lty=2)
    
  ## plot 10pc percentile, 90pc percentile and mean value of the beta distribution:
  ten_pc_quantile <- qbeta(0.1,bd_params$alpha,bd_params$beta)
  ninety_pc_quantile <- qbeta(0.9,bd_params$alpha,bd_params$beta)
    
  xpts_tpc <- c(ten_pc_quantile,ten_pc_quantile)
  ypts_tpc <- c(0,dbeta(ten_pc_quantile,bd_params$alpha,bd_params$beta))
  xpts_npc <- c(ninety_pc_quantile,ninety_pc_quantile)
  ypts_npc <- c(0,dbeta(ninety_pc_quantile,bd_params$alpha,bd_params$beta))
    
  mw <- bd_params$alpha/(bd_params$alpha+bd_params$beta)
  xpts_mean <- c(mw,mw)
  ypts_mean <- c(0,dbeta(mw,bd_params$alpha,bd_params$beta))
    
  points(xpts_tpc,ypts_tpc,type="l",col="red",lwd=2)
  points(xpts_npc,ypts_npc,type="l",col="red",lwd=2)
  points(xpts_mean,ypts_mean,type="l",col="red",lwd=2)
  df_descr_steady_state_probs <- data.frame(state_index=index_state,alpha=bd_params$alpha,beta=bd_params$beta,tpc_qtl=ten_pc_quantile,mean=mw,npc_qtl=ninety_pc_quantile)
  return(df_descr_steady_state_probs)
}

## HMM_SN_Calculate_Steady_State_Probs_From_Samples(mat_hmm_sn_stan_samples,num_states)
HMM_SN_Calculate_Steady_State_Probs_From_Samples <- function(mat_hmm_sn_stan_samples,num_states)
{
  ## extract the matrix information:
  ## loop over the samples:
  num_samples <- nrow(mat_hmm_sn_stan_samples)
  steady_state_prob_vectors <- matrix(0,nrow=num_samples,ncol=num_states)
  
  print("number of samples: ")
  print(num_samples)
  for(i in 1:num_samples)
  {
    trans_matrix_loc <- matrix(as.numeric(mat_hmm_sn_stan_samples[i,(num_states+1):(num_states+num_states^2)]),nrow=num_states,ncol=num_states)
    print(trans_matrix_loc)
    if(num_states==2)
    {
      steady_state_prob_vectors[i,] <- HMM_Calculate_Steady_State_Probs_N2_Simple(trans_matrix_loc)##_Simple(trans_matrix_loc)
      
    }
    else
    {
      steady_state_prob_vectors[i,] <- HMM_Calculate_Steady_State_Probs(trans_matrix_loc)
    }
  }
  return(steady_state_prob_vectors)
}

HMM_Get_Beta_Distribution_Params_By_Moment_Matching <- function(samples)
{
  mu <- mean(samples)
  sigma <- sd(samples)
  
  alpha_plus_beta <- (mu*(1-mu)/(sigma^2)) -1.0;
  alpha <- alpha_plus_beta*mu
  beta <- alpha_plus_beta - alpha
  df_beta_params <- data.frame(alpha=alpha,beta=beta)
  return(df_beta_params)
}

HMM_Get_Inv_Gamma_Distribution_Params_By_Moment_Matching <- function(samples)
{
  samp_mean <- mean(samples)
  samp_var <- var(samples)
  
  alpha_mm <- (samp_mean^2)/samp_var + 2
  beta_mm <- samp_mean*(((samp_mean^2)/samp_var) + 1)
  
  df_inv_gamma_params <- data.frame(alpha_mm=alpha_mm,beta_mm=beta_mm)
  return(df_inv_gamma_params)
}

HMM_Get_Gamma_Distribution_Params_By_Moment_Matching <- function(samples)
{
  samp_mean <- mean(samples)
  samp_var <- var(samples)
  
  ## method of moments:
  
  theta_mm <- samp_var/samp_mean
  k_mm <- (samp_mean^2)/samp_var
  df_gamma_params <- data.frame(theta_mm=theta_mm,k_mm=k_mm)
}

Plot_Beta_Distribution <- function(df_beta_params,num_points)
{
  x_values <- seq(from=0,to=1,length=num_points)
  y_values <- matrix(nrow=length(x_values),ncol=1)
  for(i in 1:length(x_values))
  {
    y_values[i] <- dbeta(x_values[i],df_beta_params$alpha,df_beta_params$beta)
  }
  df_beta_plot <- data.frame(x_beta_dist=x_values,y_beta_dist=y_values)
  return(df_beta_plot)
}

## Plot_Normal_Distribution(mu,sigma,x_min,x_max,num_points)
Plot_Normal_Distribution <- function(mu,sigma,x_min,x_max,num_points)
{
  x_values <- seq(from=x_min,to=x_max,length=num_points)
  y_values <- matrix(nrow=num_points,ncol=1)
  for(i in 1:num_points)
  {
    y_values[i] <- dnorm(x_values[i],mean=mu,sd=sigma)
  }
  df_norm_plot <- data.frame(x_norm_dist=x_values,y_norm_dist=y_values)
  return(df_norm_plot)
}

## HMM_Calculate_Lifetime_Likelihood(t_m,tau,n)
## calculates the likelihood of observing a mean value (t_m) for
## the lifetime of n exponentially distributed random variables,
## when the true underlying lifetime is tau 
## sum of n exponentially distributed RVs is Erlang-distributed,
## tiny parameter transformation leads to tbe implemented Likelihood function
HMM_Calculate_Lifetime_Likelihood <- function(t_m,tau,n)
{
  factor1 <- (1.0/tau)^n
  factor2 <- (n*t_m)^(n-1)
  factor3 <- exp(-(n*t_m)/tau)
  factor4 <- gamma(n)
  lik <- n*factor1*factor2*factor3/factor4;
  return(lik)
}

## HMM_Calculate_Lifetime_LogLikelihood(t_m,tau,n)
## calculates the log likelihood for the mean value (t_m)
## out of n exponentially distributed RVs with underlying 
## lifetime tau, see description above
HMM_Calculate_Lifetime_LogLikelihood <- function(t_m,tau,n)
{
  log_lik <- log(n) - n*log(tau)+(n-1)*log(n*t_m)-n*t_m/tau - lgamma(n)
  return(log_lik)
}

## some functions for local (on Mac) evaluation of HMM-data:
HMM_Evaluate_Multiple_Files <- function()
{
  require("gtools")
  require("rstan")
  filename_case_list <- "/Users/christianwachauf/Documents/Daten/data_stack_faki/AuswertungDeflectionTraces/list_of_cases_to_evaluate.csv"
  locte <- read.csv(file=filename_case_list,sep=";")
  print(locte)
  ##print(nrow(locte))
  num_rows <- nrow(locte)
  
  ## loop over input files, case indices,...
  #for(i in 1:num_rows)#length(nrow(locte)))
  for(i in 3:num_rows)
  {
    print(toString(locte$folder_name[i]))
    complete_path_curr_data <- paste0("/Users/christianwachauf/Documents/Daten/data_stack_faki/DeflectionTraces/",toString(locte$folder_name[i]),"/")
    case_indices <- toString(locte$case_indices[i])
    
    min_index <- as.numeric(locte$index_min[i])
    max_index <- as.numeric(locte$index_max[i])
    
    i#ndices <- c(0,2,3,4,5,6,7,8,9)
    
    #j<-min_index;
    for(j in 2:max_index)
    {
      filename_scaffold <- paste0(toString(locte$folder_name[i]) ,"_",case_indices,"extext",toString(j));
      filename_scaffold2 <- paste0(toString(locte$folder_name[i]) ,"_",case_indices,"ext",toString(j));
      filename_post_mean <- paste0(filename_scaffold,"_post_mean.csv")
      filename_data <- paste0(filename_scaffold2,".txt")
      
      path_post_mean <- paste0(complete_path_curr_data,filename_post_mean)
      path_data <- paste0(complete_path_curr_data,filename_data)
      
      print(path_post_mean)
      print(path_data)
      
      df_post_mean <- read.csv(file=path_post_mean,header=TRUE,sep=";")
      print(df_post_mean)
      
      
      ## data to pass to stan
      alphas0 <- c(df_post_mean$phi1_mean*20,df_post_mean$phi2_mean*20)
      location_params <- c(df_post_mean$loc1_mean,df_post_mean$loc2_mean)
      scale_params <- c(df_post_mean$scales1_mean,df_post_mean$scales2_mean)
      shape_params <- c(df_post_mean$shape1_mean,df_post_mean$shape2_mean)
      
      
      init_init_probs <- as.vector(rdirichlet(1,alphas0))
      init_trans_matrix <- matrix(nrow=2,ncol=2)
      init_trans_matrix[1,1] <- 0.999
      init_trans_matrix[1,2] <- 0.001
      init_trans_matrix[2,1] <- 0.001
      init_trans_matrix[2,2] <- 0.999
      
      ## nochmal überarbeiten: Sinnvolle Auswahl für die Hyperparameter
      ## der Transition-Matrix
      hps_trans_mat <- matrix(nrow=2,ncol=2)
      hps_trans_mat[1,1] <- 99
      hps_trans_mat[1,2] <- 1
      hps_trans_mat[2,1] <- 1
      hps_trans_mat[2,2] <- 99
      init_values_stan_hmm <- list(list(pis=init_init_probs,trans_mat=init_trans_matrix))
      
      print("initial values: ")
      print(init_values_stan_hmm)
      
      df_data <- read.csv(file=path_data,header=TRUE,sep=";")
      curr_data_subset <- df_data$bspext[seq(from=1,to=length(df_data$bspext),by=10)]
      
      num_iters <- 1000;
      num_chains <- 1;
      thinning_factor<-1;
    
      data_stan_hmm <- list(N=length(curr_data_subset),num_states=2,obs_sequ=curr_data_subset,alphas=alphas0,trans0=hps_trans_mat,location0s=location_params,scale0s=scale_params,shape0s=shape_params)
      fit_stan_hmm_sn_emprobs <- stan(file="HMM_SN_emprobs.stan",data=data_stan_hmm,init=init_values_stan_hmm,iter=num_iters,chains=num_chains,thin=thinning_factor)
      summary(fit_stan_hmm_sn_emprobs)
      print(summary(fit_stan_hmm_sn_emprobs))
      ## dump samples and dump mean-posterior values:
      filename_samples_hmm <- paste0(filename_scaffold2,"_HMM_samples.csv");
      path_samples_hmm <- paste0(complete_path_curr_data,filename_samples_hmm)
      mat_samples <- as.matrix(fit_stan_hmm_sn_emprobs)
      print(head(mat_samples))
      print("writing samples to: ")
      print(path_samples_hmm)
      write.table(x=mat_samples,file=path_samples_hmm,sep=";")
      
      #mean posterior values:
      df_mp_hmm_values <- HMM_Mean_Posteriors_From_Samples(fit_stan_hmm_sn_emprobs)
      filename_mpvs_hmm <- paste0(filename_scaffold2,"_HMM_mpvs.csv");
      path_mpvs_hmm <- paste0(complete_path_curr_data,filename_mpvs_hmm)
  
      write.table(x=df_mp_hmm_values,file=path_mpvs_hmm,sep=";")
      
      ## calculate Viterbi-sequence, save it to a file and
      ## plot Viterbi-sequence to ".pdf" 
      ## calculate viterbi sequence assignment and save it:
      ##print("assigning viterbi sequence: ")
      ##vss <- HMM_SN_Assign_Viterbi_Sequence_From_Mean_Posterior(obs_seq,mat_stan_hmm_sn_samples,num_states)
      print("calculating Viterbi sequence: ")
      vss <- HMM_SN_Assign_Viterbi_Sequence_From_Mean_Posterior(curr_data_subset,mat_samples,2)
      df_viterbi <- data.frame(x=curr_data_subset,vss=vss)
      filename_hmm_sn_vss <- paste0(filename_scaffold2,"_hmm_sn_vss.csv")
      path_hmm_sn_vss <- paste0(complete_path_curr_data,filename_hmm_sn_vss)
      write.table(x=df_viterbi,file=path_hmm_sn_vss,sep=";")
      
      ## Create pdf-file with deflection-plot (Viterbi state sequence)
      
      filename_hmm_sn_vss_plot <-  paste0(filename_scaffold2,"_hmm_sn_vss.pdf")
      path_hmm_sn_vss_plot <- paste0(complete_path_curr_data,filename_hmm_sn_vss_plot)
      pdf(file=path_hmm_sn_vss_plot,width=14,height=11,useDingbats=FALSE)
      casename <- filename_scaffold2
      HMM_Plot_Trace_With_State_Assignment(obs_seq=curr_data_subset,state_seq=vss,dt=1,lab_x="index",lab_y="deflection [nm]",case_str = casename)
      dev.off()
      
      
      ## plot steady-state probability (samples of it) of state1
      filename_hmm_state1_prob_plot <- paste0(filename_scaffold2,"_hmm_sn_s1_ssp.pdf")
      filename_hmm_state1_prob <- paste0(filename_scaffold2,"_hmm_sn_s1_ssp.csv")
      path_hmm_state1_prob_plot <- paste0(complete_path_curr_data,filename_hmm_state1_prob_plot)
      path_hmm_state1_prob <- paste0(complete_path_curr_data,filename_hmm_state1_prob)
      
      
      
      
      pdf(file=path_hmm_state1_prob_plot,width=14,height=11,useDingbats=FALSE)
      df_descr_state1 <- HMM_SN_Plot_Steady_State_Probability(mat_samples,2,1)
      dev.off()
      ## write table with information(description of the beta-distribution)
      write.table(x=df_descr_state1,file=path_hmm_state1_prob,sep=";")
      
      #write.table(x=df_viterbi,file=filename_hmm_sn_vss)
      
      #pdf(file=filename_hmm_sn_vss_plot,width=14,height=11,useDingbats=FALSE)
      #HMM_Plot_Trace_With_State_Assignment(obs_seq=curr_data_subset,state_seq=vss,dt=1,lab_x="index",lab_y="deflection [nm]",case_str = casename)
      #dev.off()
      
    }
  }
}

HMM_Mean_Posteriors_From_Samples <- function(fit_stan_hmm_sn_emprobs)
{
    mat <- as.matrix(fit_stan_hmm_sn_emprobs)
    phi1_mean <- mean(mat[,1])
    phi2_mean <- mean(mat[,2])
    
    mat11_mean <- mean(mat[,3])
    mat12_mean <- mean(mat[,4])
    mat21_mean <- mean(mat[,5])
    mat22_mean <- mean(mat[,6])
    
    df_post_mean_hmm <- data.frame(phi1_mean,phi2_mean,mat11_mean,mat12_mean,mat21_mean,mat22_mean)
    return(df_post_mean_hmm)
}