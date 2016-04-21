## Load necessary pacakges
require(mvtnorm)

## Function 1: sampler_ESrepl
sampler_ESrepl <- function (data,burnin=10000,iterations=20000,prior,type,Htype=1,assessment_values="auto",d_diff=0){
  phi <- d_diff 
  # Sample values
  m_sample  <- sum(data[,2]*data[,1])/sum(data[,2]) # mean group 1
  sd_sample <- aggregate(data[,1],list(data[,2]),sd) # sd per group in sample
  sd_pool   <- sqrt(sum((table(data[,2])-1)*(sd_sample[,2])^2)/sum(table(data[,2])-1)) # pooled sd
  d_sample  <- (m_sample-sum(data[,3]*data[,1])/sum(data[,3]))/sd_pool # effect size (Cohen's d) in sample
  
  # Proposal distribution MH algorithm fixed scale and shape
  shape_post_prop <- (nrow(data)/2)+prior[3,1]
  scale_post_prop <- prior[3,2]+(1/2)*sum((data[,1]-m_sample+data[,3]*d_sample*sd_pool)^2)
  proposal_parameters <- c(shape_post_prop,scale_post_prop)
  
  # Start values for all parameters
  m_start   <- m_sample
  d_start   <- d_sample
  var_start <- 1/rgamma(1,shape=proposal_parameters[1],rate=proposal_parameters[2])
  
  # Create matrix for iteration values
  iter_val <- matrix(NA,nrow=Htype*iterations+1,ncol=3)
  colnames(iter_val) <- paste(c("m","d","var"))
  
  # Place start values at first row
  iter_val[1,] <- c(m_start,d_start,var_start)
  
  # Create matrix for MH-algorithm values
  iter_MH <- matrix(NA,nrow=Htype*iterations,ncol=5)
  colnames(iter_MH) <- paste(c("cand_it","a_it","a_m","aq_m","a_j"))
  
  if(Htype==1){
    prior[2,2] <- sqrt(prior[2,2]^2+d_diff^2)
    
  }
  
  else {
    priord <- matrix(c(prior[2,1]-d_diff,prior[2,2],prior[2,1]+d_diff,prior[2,2]),nrow=2,ncol=2,byrow=T)
    prior[2,] <- priord[1,]}
  
  for (i in 1:(Htype*iterations)){
    if(Htype==2) {if (i >= iterations) {prior[2,] <- priord[2,]}}
    #### SAMPLE M and D from two univariate normal distributions (Gibbs sampler)
    m_mean_post <- (((1/iter_val[i,3])*(sum(data[,1]+data[,3]*iter_val[i,2]*sqrt(iter_val[i,3]))))+(prior[1,1]/(prior[1,2]^2)))/((nrow(data)/iter_val[i,3])+(1/prior[1,2]^2))
    m_sd_post <- sqrt(1/((nrow(data)/iter_val[i,3])+(1/prior[1,2]^2))    )
    m_it <- rnorm(1,m_mean_post,m_sd_post)
    iter_val[i+1,1] <- m_it
    
    d_mean_post <- ((prior[2,1]/prior[2,2]^2)-((sum(data[,3]*(data[,1]-iter_val[i+1,1]))/sqrt(iter_val[i,3]))))/(sum(data[,3])+(1/prior[2,2]^2))
    d_sd_post <- sqrt(1/(sum(data[,3])+(1/prior[2,2]^2)))
    d_it <- rnorm(1,d_mean_post,d_sd_post)
    iter_val[i+1,2] <- d_it 
    
    #### SAMPLE variance from an unknown posterior distribution (Metropolis-Hastings-algorithm)
    var_it_cand <- 1/rgamma(1,shape=proposal_parameters[1],rate=proposal_parameters[2])
    iter_MH[i,1] <- var_it_cand
    
    density_proposal_cand <- dgamma(1/var_it_cand,shape=proposal_parameters[1],rate=proposal_parameters[2])
    density_proposal_prev <- dgamma(1/iter_val[i,3],shape=proposal_parameters[1],rate=proposal_parameters[2])
    
    density_posterior_cand <- (var_it_cand^(((-1)*nrow(data)/2)-prior[3,1]-1)*exp((-1/(var_it_cand))*(prior[3,2]+(1/2)*sum((data[,1]-iter_val[i+1,1]+data[,3]*iter_val[i+1,2]*sqrt(var_it_cand))^2))))
    density_posterior_prev <- ((iter_val[i,3])^(((-1)*nrow(data)/2)-prior[3,1]-1)*exp((-1/(iter_val[i,3]))*(prior[3,2]+(1/2)*sum((data[,1]-iter_val[i+1,1]+data[,3]*iter_val[i+1,2]*sqrt(iter_val[i,3]))^2))))
    
    r_it <- (density_proposal_prev/density_proposal_cand)*(density_posterior_cand/density_posterior_prev)
    if(is.na(r_it)==T) {r_it<-0}
    
    a_it <- min(1,r_it)
    iter_MH[i,2] <- a_it
    pick_cand_it <- rbinom(1,1,prob=a_it)
    
    if(pick_cand_it==0){iter_val[i+1,3]<-iter_val[i,3]} else {iter_val[i+1,3]<-var_it_cand}   
  }
  
  
  ## posterior means
  if(Htype==1){
    m_m_post <- mean(iter_val[(burnin+2):(Htype*iterations+1),1])
    m_sd_post <- sd(iter_val[(burnin+2):(Htype*iterations+1),1])
    d_m_post <- mean(iter_val[(burnin+2):(Htype*iterations+1),2])
    d_sd_post <- sd(iter_val[(burnin+2):(Htype*iterations+1),2])
    var_m_post <- mean(iter_val[(burnin+2):(Htype*iterations+1),3])
    var_sd_post <- sd(iter_val[(burnin+2):(Htype*iterations+1),3])
  }
  else
  {
    m_m_post <- mean(c(iter_val[(burnin+2):(iterations+1),1],iter_val[((Htype-1)*burnin+iterations+1):(Htype*iterations+1),1]))
    m_sd_post <-  sd(c(iter_val[(burnin+2):(iterations+1),1],iter_val[((Htype-1)):(Htype*iterations+1),1]))
    d_m_post <- mean(c(iter_val[(burnin+2):(iterations+1),2],iter_val[((Htype-1)*burnin+iterations+1):(Htype*iterations+1),2]))
    d_sd_post <-  sd(c(iter_val[(burnin+2):(iterations+1),2],iter_val[((Htype-1)*burnin+iterations+1):(Htype*iterations+1),2]))
    var_m_post <-mean(c(iter_val[(burnin+2):(iterations+1),3],iter_val[((Htype-1)*burnin+iterations+1):(Htype*iterations+1),3]))
    var_sd_post <-sd(c(iter_val[(burnin+2):(iterations+1),3],iter_val[((Htype-1)*burnin+iterations+1):(Htype*iterations+1),3]))
    
    ## posterior means mixture parts
    
    ###### IDEA: AGGREGATE function ###########
    m_m_post_mix1 <-  mean(iter_val[(burnin+2):(iterations+1),1])
    m_sd_post_mix1 <- sd(iter_val[(burnin+2):(iterations+1),1])
    d_m_post_mix1 <-  mean(iter_val[(burnin+2):(iterations+1),2])
    d_sd_post_mix1 <- sd(iter_val[(burnin+2):(iterations+1),2])
    var_m_post_mix1 <-mean(iter_val[(burnin+2):(iterations+1),3])
    var_sd_post_mix1 <- sd(iter_val[(burnin+2):(iterations+1),3])
    
    m_m_post_mix2 <- mean(iter_val[((Htype-1)*burnin+iterations+1):(Htype*iterations+1),1])
    m_sd_post_mix2 <- sd(iter_val[((Htype-1)*burnin+iterations+1):(Htype*iterations+1),1])
    d_m_post_mix2 <- mean(iter_val[((Htype-1)*burnin+iterations+1):(Htype*iterations+1),2])
    d_sd_post_mix2 <- sd(iter_val[((Htype-1)*burnin+iterations+1):(Htype*iterations+1),2])
    var_m_post_mix2 <- mean(iter_val[((Htype-1)*burnin+iterations+1):(Htype*iterations+1),3])
    var_sd_post_mix2 <- sd(iter_val[((Htype-1)*burnin+iterations+1):(Htype*iterations+1),3])
  }
  
  ## Save posterior means and sds
  posterior_summary <- matrix(c(m_m_post,m_sd_post,d_m_post,d_sd_post,var_m_post,var_sd_post),nrow=3,ncol=2,byrow=T)
  colnames(posterior_summary) <- paste(c("mean","sd"))
  rownames(posterior_summary) <- paste(c("m","d","var"))
  
  if(Htype==2){
    posterior_summary_mix1 <- matrix(c(m_m_post_mix1,m_sd_post_mix1,d_m_post_mix1,d_sd_post_mix1,var_m_post_mix1,var_sd_post_mix1),nrow=3,ncol=2,byrow=T)
    posterior_summary_mix2 <- matrix(c(m_m_post_mix2,m_sd_post_mix2,d_m_post_mix2,d_sd_post_mix2,var_m_post_mix2,var_sd_post_mix2),nrow=3,ncol=2,byrow=T)
  }
  
  ## Output
  results <- list("Gibbs"=iter_val,"MH"=iter_MH[,(1:2)],"posterior_summary"=posterior_summary,"prior"=prior)
  
  ### If type="orig" return results and exit function, if type="repl" calculate
  ### marginal likelihood of the data using Chib and Jeliazkov (2001)
  if (type=="orig"){return(results)} else {
    
    if (assessment_values=="auto"){assessment_values <- posterior_summary[,1]}
    
    # log f (y|m*,d*,var*)
    Chib_data_ln <- sum(dnorm(data[,1],mean=(assessment_values[1]-(data[,3]*assessment_values[2]*sqrt(assessment_values[3]))),sd=sqrt(assessment_values[3]),log=T))
    # log pi (m*,d*,var*)
    if(Htype==1){
      Chib_prior_ln <- dnorm(assessment_values[1],mean=prior[1,1],sd=prior[1,2],log=T)+dnorm(assessment_values[2],mean=prior[2,1],sd=prior[2,2],log=T)+dgamma(1/assessment_values[3],shape=prior[3,1],rate=prior[3,2],log=T)
    }
    else{
      Chib_prior_ln <- dnorm(assessment_values[1],mean=prior[1,1],sd=prior[1,2],log=T)+log((1/2)*dnorm(assessment_values[2],mean=priord[1,1],sd=priord[1,2],log=F)+(1/2)*dnorm(assessment_values[2],mean=priord[2,1],sd=priord[2,2],log=F))+dgamma(1/assessment_values[3],shape=prior[3,1],rate=prior[3,2],log=T)
    }
    # log p (m*,d*|y,var*) from bivariate normal distribution
    if (Htype==1){
      mu_matrix <- rep(NA,2)
      A <- ( (posterior_summary[3,1])^(-1)*sum(data[,1]) + prior[1,1]*prior[1,2]^(-2) + (posterior_summary[3,1])^(-1/2)*prior[2,1]*prior[2,2]^(-2)*sum(data[,3]*(sum(data[,3]+prior[2,2]^(-2)))^(-1)) ) * (nrow(data)*posterior_summary[3,1]^(-1) + prior[1,2]^(-2))^(-1)
      B <- posterior_summary[3,1]^(-1) * (nrow(data)*posterior_summary[3,1]^(-1) + prior[1,2]^(-2))^(-1) * (sum(data[,3])+prior[2,2]^(-2))^(-1)
      mu_matrix[1] <- (1-B*sum(data[,3])^2)^(-1) * (A-B*sum(data[,3])*sum(data[,3]*data[,1]))
      mu_matrix[2] <- ( (prior[2,1]*prior[2,2]^(-2)) - (posterior_summary[3,1]^(-1/2) * sum(data[,3]*(data[,1]-((1-B*sum(data[,3])^2)^(-1) * (A-B*sum(data[,3])*sum(data[,3]*data[,1])))))  ) ) * (sum(data[,3])+prior[2,2]^(-2))^(-1)
      
      information_matrix <- matrix(c((nrow(data)/assessment_values[3])+(1/prior[1,2]^2),-sum(data[,3])/sqrt(assessment_values[3]),-sum(data[,3])/sqrt(assessment_values[3]),(sum(data[,3])+(1/prior[2,2]^2))),nrow=2,ncol=2)
      sigma_matrix <- solve(information_matrix)
      Chib_conditional_ln <- dmvnorm(x=c(assessment_values[1],assessment_values[2]),mean=mu_matrix,sigma=sigma_matrix,log=T)                   
    }
    else{
      mu_matrix_mix1 <- rep(NA,2)
      A_mix1 <- ( (posterior_summary_mix1[3,1])^(-1)*sum(data[,1]) + prior[1,1]*prior[1,2]^(-2) + (posterior_summary_mix1[3,1])^(-1/2)*priord[1,1]*prior[2,2]^(-2)*sum(data[,3]*(sum(data[,3]+prior[2,2]^(-2)))^(-1)) ) * (nrow(data)*posterior_summary_mix1[3,1]^(-1) + prior[1,2]^(-2))^(-1)
      B_mix1 <- posterior_summary_mix1[3,1]^(-1) * (nrow(data)*posterior_summary_mix1[3,1]^(-1) + prior[1,2]^(-2))^(-1) * (sum(data[,3])+prior[2,2]^(-2))^(-1)
      mu_matrix_mix1[1] <- (1-B_mix1*sum(data[,3])^2)^(-1) * (A_mix1-B_mix1*sum(data[,3])*sum(data[,3]*data[,1]))
      mu_matrix_mix1[2] <- ( (priord[1,1]*prior[2,2]^(-2)) - (posterior_summary_mix1[3,1]^(-1/2) * sum(data[,3]*(data[,1]-((1-B_mix1*sum(data[,3])^2)^(-1) * (A_mix1-B_mix1*sum(data[,3])*sum(data[,3]*data[,1])))))  ) ) * (sum(data[,3])+prior[2,2]^(-2))^(-1)
      
      mu_matrix_mix2 <- rep(NA,2)
      A_mix2 <- ( (posterior_summary_mix2[3,1])^(-1)*sum(data[,1]) + prior[1,1]*prior[1,2]^(-2) + (posterior_summary_mix2[3,1])^(-1/2)*priord[2,1]*prior[2,2]^(-2)*sum(data[,3]*(sum(data[,3]+prior[2,2]^(-2)))^(-1)) ) * (nrow(data)*posterior_summary_mix2[3,1]^(-1) + prior[1,2]^(-2))^(-1)
      B_mix2 <- posterior_summary_mix2[3,1]^(-1) * (nrow(data)*posterior_summary_mix2[3,1]^(-1) + prior[1,2]^(-2))^(-1) * (sum(data[,3])+prior[2,2]^(-2))^(-1)
      mu_matrix_mix2[1] <- (1-B_mix2*sum(data[,3])^2)^(-1) * (A_mix2-B_mix2*sum(data[,3])*sum(data[,3]*data[,1]))
      mu_matrix_mix2[2] <- ( (priord[2,1]*prior[2,2]^(-2)) - (posterior_summary_mix2[3,1]^(-1/2) * sum(data[,3]*(data[,1]-((1-B_mix2*sum(data[,3])^2)^(-1) * (A_mix2-B_mix2*sum(data[,3])*sum(data[,3]*data[,1])))))  ) ) * (sum(data[,3])+prior[2,2]^(-2))^(-1)
      
      information_matrix <- matrix(c((nrow(data)/assessment_values[3])+(1/prior[1,2]^2),-sum(data[,3])/sqrt(assessment_values[3]),-sum(data[,3])/sqrt(assessment_values[3]),(sum(data[,3])+(1/prior[2,2]^2))),nrow=2,ncol=2)
      sigma_matrix <- solve(information_matrix)
      
      Chib_conditional_ln <- log((1/2)*dmvnorm(x=c(assessment_values[1],assessment_values[2]),mean=mu_matrix_mix1,sigma=sigma_matrix,log=F)+(1/2)*dmvnorm(x=c(assessment_values[1],assessment_values[2]),mean=mu_matrix_mix2,sigma=sigma_matrix,log=F))
    }
    
    # log p (var*|y)
    ### M runs ###
    for(m in (burnin+1):(Htype*iterations)){
      density_proposal_star <- dgamma(1/assessment_values[3],shape=proposal_parameters[1],rate=proposal_parameters[2])
      density_proposal_it <- dgamma(1/iter_MH[m,1],shape=proposal_parameters[1],rate=proposal_parameters[2])
      
      density_posterior_star <- (assessment_values[3]^(((-1)*nrow(data)/2)-prior[3,1]-1)*exp((-1/(assessment_values[3]))*(prior[3,2]+(1/2)*sum((data[,1]-iter_val[m+1,1]+data[,3]*iter_val[m+1,2]*sqrt(assessment_values[3]))^2))))
      density_posterior_it <- ((iter_MH[m,1])^(((-1)*nrow(data)/2)-prior[3,1]-1)*exp((-1/(iter_MH[m,1]))*(prior[3,2]+(1/2)*sum((data[,1]-iter_val[m+1,1]+data[,3]*iter_val[m+1,2]*sqrt(iter_MH[m,1]))^2))))
      
      r_m <- (density_proposal_it/density_proposal_star)*(density_posterior_star/density_posterior_it)
      
      a_m <- min(1,r_m)
      iter_MH[m,3] <- a_m
      
      aq_m <- a_m*density_proposal_star
      iter_MH[m,4] <- aq_m
    }
    
    ### J runs ###
    for(j in (burnin+1):(Htype*iterations)){
      density_proposal_it <- dgamma(1/iter_MH[j,1],shape=proposal_parameters[1],rate=proposal_parameters[2])
      density_proposal_star <- dgamma(1/assessment_values[3],shape=proposal_parameters[1],rate=proposal_parameters[2])
      
      density_posterior_it <- ((iter_MH[j,1])^(((-1)*nrow(data)/2)-prior[3,1]-1)*exp((-1/(iter_MH[j,1]))*(prior[3,2]+(1/2)*sum((data[,1]-iter_val[j+1,1]+data[,3]*iter_val[j+1,2]*sqrt(iter_MH[j,1]))^2))))
      density_posterior_star <- (assessment_values[3]^(((-1)*nrow(data)/2)-prior[3,1]-1)*exp((-1/(assessment_values[3]))*(prior[3,2]+(1/2)*sum((data[,1]-iter_val[j+1,1]+data[,3]*iter_val[j+1,2]*sqrt(assessment_values[3]))^2))))
      
      r_t <- (density_proposal_star/density_proposal_it)*(density_posterior_it/density_posterior_star)
      
      a_t <- min(1,r_t)
      iter_MH[j,5] <- a_t
    }
    
    
    
    if(Htype==1){
      Chib_marginal_ln <- log(mean(iter_MH[(burnin+1):(iterations),4])/mean(iter_MH[(burnin+1):(iterations),5]))
    }
    else{
      Chib_marginal_ln <- log((1/2)*(mean(iter_MH[(burnin+1):(iterations),4])/mean(iter_MH[(burnin+1):(iterations),5]))+(1/2)*(mean(iter_MH[((iterations+1):(Htype*iterations)),4])/mean(iter_MH[((iterations+1):(Htype*iterations)),5])))
    }
    
    Chib_ln_marginal_result <- Chib_data_ln+Chib_prior_ln-Chib_conditional_ln-Chib_marginal_ln
    
    Chib_ln <- c(Chib_data_ln,Chib_prior_ln,Chib_conditional_ln,Chib_marginal_ln)
    
    results_repl <-   list("Gibbs"=iter_val,"MH"=iter_MH,"posterior_summary"=posterior_summary,"Chib"=Chib_ln,"marginal"=Chib_ln_marginal_result,"phi"=phi,"prior"=prior,"iterations"=iterations,"burnin"=burnin)
    
  }
  return(results_repl)
}


## Function 2: BF_repl_es
BF_repl_es <- function(data_orig,data_repl,ddiff,iter) {
      prior_orig <- matrix(c(0,100,0,100,.001,.001),nrow=3,ncol=2,byrow=T)
  orig <- sampler_ESrepl(data=data_orig,prior=prior_orig,type="orig",iterations=iter,d_diff=0)
      prior_repl <- matrix(c(0,100,orig$posterior_summary[2,1],orig$posterior_summary[2,2],.001,.001),nrow=3,ncol=2,byrow=T)
  repl_H1 <- sampler_ESrepl(data=data_repl,prior=prior_repl,Htype=1,type="repl",d_diff=0,iterations=iter)
  repl_H2 <- sampler_ESrepl(data=data_repl,prior=prior_repl,Htype=2,type="repl",d_diff=ddiff,iterations=iter)
  repl_H3 <- sampler_ESrepl(data=data_repl,prior=prior_repl,Htype=1,type="repl",d_diff=ddiff,iterations=iter)

  BF12 <- exp(repl_H1$marginal-repl_H2$marginal)
  BF13 <- exp(repl_H1$marginal-repl_H3$marginal)
  BF23 <- exp(repl_H2$marginal-repl_H3$marginal)
  
  BF_results <- list("BF12"=BF12,"BF13"=BF13,"BF23"=BF23,"orig"=orig,"repl_H1"=repl_H1,"repl_H2"=repl_H2,"repl_H3"=repl_H3)
  return(BF_results)
}

## Function 3: Prior distributions plot
plot_prior <- function(output){
  # ddiff <- output$repl_H1$phi
  phi <- output$repl_H2$phi
  curve(dnorm(x,mean=output$repl_H1$prior[2,1],sd=output$repl_H1$prior[2,2]),xlim=c((output$repl_H1$prior[2,1]-3*output$repl_H1$prior[2,2]),(output$repl_H1$prior[2,1]+3*output$repl_H1$prior[2,2])),xlab=expression(delta),ylab="density")
  curve((1/2)*dnorm(x,mean=output$repl_H1$prior[2,1]-phi,sd=output$repl_H1$prior[2,2])+(1/2)*dnorm(x,mean=output$repl_H1$prior[2,1]+phi,sd=output$repl_H1$prior[2,2]),lty=2,add=T)
  curve(dnorm(x,mean=output$repl_H3$prior[2,1],sd=output$repl_H3$prior[2,2]),lty=3,add=T)
  legend((output$repl_H1$prior[2,1]+1.5*output$repl_H1$prior[2,2]),dnorm(output$repl_H1$prior[2,1],mean=output$repl_H1$prior[2,1],sd=output$repl_H1$prior[2,2]),c("H_1","H_2","H_3"),lty=c(1,2,3))
}