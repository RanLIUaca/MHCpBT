#' threshold_gibbs
#' 
#' this is the main function of MHCpBT, input sequences and corresponding binding scores
#' @param data1 an array whose element is a string, e.g, ('ATGC','TTGG')
#' @param Y the binding scores
#' @param motif_len the length of the motif
#' @param dict named vector for letters
#' @param burn_in the burn_in setting for Gibbs
#' @param end_times the total iterations
#' @param result_path where to store the results
#' @param only_get_func whether only loading the functions
#' @param input_prior whether inputting priors by users
#' @param prior different prior settings
#' @import stringr
#' @import ggplot2
#' @import ggseqlogo
#' @return a dataframe with three columns: pep, real_score, est_score
#' @export
threshold_gibbs<-function(data1, Y, motif_len, dict, burn_in=2000, end_times=3000, result_path,
                          only_get_func = F,input_prior = F,prior = 0){
  
  order_index = order(as.numeric(Y),decreasing = T)
  data1 = data1[order_index]
  Y = Y[order_index]
  
  digit_dict = as.character(1:length(dict))
  names(digit_dict) = dict
  
  # print(digit_dict)
  interval_sel = 1
  if(dir.exists(result_path)==0){dir.create(result_path)}
  
  dir_names = result_path
  pdf(paste(dir_names,'result.pdf',sep=''))
  
  # prior 
  J = motif_len
  n = length(data1)
  beta_0 = rep(0,J+1)
  beta = beta_0
  
  
  alpha_0 = 1
  T_0 = 1
  Sigma_0 = diag(rep(10,J+1))
  # c('0.5_1','1_0.5','5i','15i')
  if(input_prior==T){
    print(prior)
    if(prior == '5i'){
      Sigma_0 = diag(rep(5,J+1))
    }
    if(prior == '15i'){
      Sigma_0 = diag(rep(15,J+1))
    }
    if(prior == '1_0.5'){
      alpha_0 = 1
      T_0 = 0.5
    }
    if(prior == '0.5_1'){
      alpha_0 = 0.5
      T_0 = 1
    }
  }
  
  
  lb=c(-Inf,rep(0,nrow(Sigma_0)-1))
  
  eta_0 = matrix(rep(1,length(dict)*J),nrow = length(dict))
  gamma_0 = rep(1,length(dict))
  
  
  
  {
    H <<- function(sel_sequence){
      
      H = rep(0,length(dict))
      tmp_index = 1:str_length(sel_sequence)
      
      sel_sequence = str_sub(sel_sequence,tmp_index,tmp_index)
      sel_sequence = str_replace_all(sel_sequence, digit_dict)
      for (i in sel_sequence) {
        tmp_i = as.integer(i)
        H[tmp_i] = H[tmp_i] + 1
      }
      
      return(H)
    }
    
    # calculate the theta
    cal_theta<<- function(sel_element){
      theta = matrix(nrow = length(dict), ncol = J)
      for (i in 1:J) {
        hakj = simplify2array(lapply(sel_element,function(x) str_sub(x,i,i)))
        temp = H(paste(hakj,collapse = ''))
        theta[,i] = as.character(temp/sum(temp))
      }
      row.names(theta) = dict
      
      return(theta)
    }
    
    
    # calculate the X
    cal_X<<-function(sel_element, theta){
      X = rep(1,length(sel_element))
      for (i in 1:ncol(theta)) {
        hakj = simplify2array(lapply(sel_element,function(x) str_sub(x,i,i)))
        
        ttmp = theta[,i]
        names(ttmp) = dict
        mode(ttmp) = "character"
        temp = str_replace_all(hakj, ttmp)
        X = cbind(X,temp)
      }
      mode(X) = 'double'
      return(X)
    }
    
    sel<<-function(data1,sel_index){
      
      sel_element = c()
      if(mode(nrow(sel_index)) == 'NULL' ){
        
        # select one col sel(data1,A[,j])
        if(length(sel_index)==length(data1)){
          tmp_col = c()
          for (j in 1:length(sel_index)) {
            tmp_index = sel_index[j]
            tmp = str_sub(data1[j],tmp_index,tmp_index)
            tmp_col = paste(tmp_col, tmp,sep='')
          }
          return(tmp_col)
        }
        
        
      }else{
        # select all A sel(data1,A)
        for (i in 1:nrow(sel_index)) {
          tmp_row = c()
          for (j in 1:length(sel_index[i,])) {
            tmp_index = sel_index[i,j]
            tmp = str_sub(data1[i],tmp_index,tmp_index)
            tmp_row = paste(tmp_row, tmp,sep='')
          }
          sel_element = c(sel_element,tmp_row)
        }
        return(sel_element)
      }
      
      
    }
    
    
    V<<-function(X,beta,sigma2){
      temp0 = t(Y-X%*%beta)%*%(Y-X%*%beta)
      V = exp(-temp0/(2*sigma2))
      return(V)
    }
    
    log_joint_post<<-function(A,theta,theta_0,beta,sigma,d,print_part = F){
      binding_index = 1:d
      sel_element = sel(data1,A)
      X = cal_X(sel_element, theta)
      
      temp0 = t(Y-X%*%beta)%*%(Y-X%*%beta)
      joint_post = -n*log(sigma) - temp0/(2*sigma**2)
      
      if(print_part == T){
        print(paste('linear lik part=',joint_post))
      }
      
      ### the closest to zero
      temp = 0
      for (j in 1:ncol(theta)) {
        haj = H(paste(sel(data1[binding_index],A[binding_index,j]),collapse = ''))
        temp = temp + sum((haj+eta_0[j]-1)*log(theta[,j]))
      }
      joint_post = joint_post + temp
      if(print_part == T){
        print(paste('theta part=',temp))
      }
      
      hall = H(paste(data1, collapse = ''))
      ha = H(paste(sel(data1[binding_index],A[binding_index,]),collapse = ''))
      
      hac = hall - ha
      
      # update theta_0
      temp = sum((hac+gamma_0-1)*log(theta_0))
      joint_post = joint_post + temp
      
      if(print_part == T){
        print(paste('theta0 part=',temp))
      }
      
      temp = sum(TruncatedNormal::dtmvnorm(beta,beta_0,Sigma_0,lb=lb,log=T))
      joint_post = joint_post + temp
      if(print_part == T){
        print(paste('beta part=',temp))
      }
      
      
      temp = sum(dgamma((1/sigma)**2,alpha_0/2,rate = T_0/2,log=T))
      joint_post = joint_post+temp
      if(print_part == T){
        print(paste('sigma2=',temp))
      }
      
      
      return(joint_post)
    }
    
    
    
  }# useful function
  if(only_get_func == T){
    return(0)
  }
  

  d_index = 5:(n-1)
  
  
  
  {### initial value
    a = c()
    for (i in 1:length(data1)) {
      a = c(a,sample(1:(str_length(data1[i])-J+1),1,replace=T))
    }
    
    
    A = matrix(nrow = length(data1), ncol = J)
    for (i in 1:length(a)) {
      A[i,] = seq(a[i],length.out = J, by=1)
    }
    
    # d = floor(n/2)
    d = max(floor(n/5),5)
    binding_index = 1:d
    
    data_binding = data1[1:d]
    Ad = A[1:d,]
    theta = matrix(rep(0,J*length(dict)),nrow = length(dict))
    for (j in 1:J) {
      haj = H(paste(sel(data_binding,Ad[,j]),collapse = ''))
      theta[,j] = MCMCpack::rdirichlet(1,haj+eta_0[,j])
    }
    row.names(theta) = dict
    
    
  }# initial value
  
  
  
  
  iter = 1
  
  
  
  # all_ar = c()
  
  result_theta = 0
  result_beta = 0
  result_score = 0
  result_theta_0 = 0
  result_sigma2 = 0
  result_a = rep(0,length(data1))
  result_d_trace = c()
  
  result_coe = c()
  result_error = c()
  result_MSE = c()
  result_log_joint = c()
  result_mean_joint = c()
  
  
  repeat{
    binding_index = 1:d
    sel_element = sel(data1,A)
    X = cal_X(sel_element, theta)
    
    # update sigma2
    rr = t(Y-X%*%beta)%*%(Y-X%*%beta)
    Isigma2 = rgamma(1,(alpha_0+n)/2,rate = (T_0+rr)/2)
    sigma2 = 1/Isigma2
    
    
    # non-corr prior
    temp_XX = t(X)%*%X
    tttmp = MASS::ginv(MASS::ginv(Sigma_0)+temp_XX/sigma2)
    mu = tttmp%*%(MASS::ginv(Sigma_0)%*%beta_0+t(X)%*%Y/sigma2)
    Sigma = tttmp
    
    
    
    # non-negative
    lb=c(-Inf,rep(0,nrow(Sigma)-1))
    beta = TruncatedNormal::rtmvnorm(1,mu,Sigma,lb=lb)
    
    
    data_binding = data1[binding_index]
    Ad = A[binding_index,]
    hall = H(paste(data1, collapse = ''))
    ha = H(paste(sel(data_binding,Ad),collapse = ''))
    hac = hall - ha
    
    # update theta_0
    theta_0 = MCMCpack::rdirichlet(1,hac+gamma_0)
    
    # update theta
    theta_star = theta
    print(iter)
    for (j in 1:J) {
      haj = H(paste(sel(data_binding,Ad[,j]),collapse = ''))
      theta_star[,j] = MCMCpack::rdirichlet(1,haj+eta_0[,j])
      X_star = cal_X(sel_element, theta_star)
      
      # update beta
      # non-corr prior
      new_XX = t(X_star)%*%X_star
      new_tttmp = MASS::ginv(MASS::ginv(Sigma_0)+new_XX/sigma2)
      new_mu = new_tttmp%*%(MASS::ginv(Sigma_0)%*%beta_0+t(X_star)%*%Y/sigma2)
      new_Sigma = new_tttmp
      
      # non-negative
      lb=c(-Inf,rep(0,nrow(Sigma)-1))
      beta_star = TruncatedNormal::rtmvnorm(1,new_mu,new_Sigma,lb=lb)
      
      pa = log_joint_post(A,theta_star,theta_0,beta_star,sqrt(sigma2),d)
      pb = log_joint_post(A,theta,theta_0,beta,sqrt(sigma2),d)
      
      temp_ratio = -sum(log(MCMCpack::ddirichlet(theta_star[,j],haj+eta_0[,j]))) + 
        sum(log(MCMCpack::ddirichlet(theta[,j],haj+eta_0[,j])))
      
      temp_ratio = temp_ratio - sum(TruncatedNormal::dtmvnorm(beta_star, new_mu, new_Sigma, lb,log=T)) +
        sum(TruncatedNormal::dtmvnorm(beta, mu, Sigma, lb,log=T))
      
      temp_ar = pa - pb + temp_ratio
      ar = min(0,temp_ar) 
      # print(paste(delta,'first_shift_rate=',exp(ar)))
      u = runif(1,0,1)
      u = log(u)
      if(u<=ar){
        # print(paste('j=',j,'shift success'))
        theta = theta_star
        beta = beta_star
        X = X_star
        mu = new_mu
        Sigma = new_Sigma
      }
      else{
        theta_star[,j] = theta[,j]
      }
    }
    
    
    
    # update A
    nb_aprob = 0
    for (k in 1:n) {
      if(k<=d){
        # the length of possible positions
        lpp = str_length(data1[k]) - J + 1
        # the posterior probability of possible positions
        ppp = c()
        for (i in 1:lpp) {
          temp_Ak = matrix(i:(i+J-1),nrow=1)
          temp_sel = sel(data1[k],temp_Ak)
          temp_Xk = cal_X(temp_sel, theta)
          
          ratio = theta/theta_0[1,]
          row.names(ratio) = dict
          hakj = cal_X(temp_sel,ratio)
          hakj = prod(hakj)
          
          temp0 = (Y[k]-temp_Xk%*%beta)**2
          ttpp = exp(-temp0/(2*sigma2))
          
          result = ttpp*hakj
          ppp = c(ppp,result)
          # print(ppp)
        }
        # print(ppp)
        ppp = ppp/sum(ppp)
        #print(ppp)
        new_ak = sample(1:lpp,1, prob = ppp)
        A[k,] = seq(new_ak,length.out = J, by=1)
      }
      else{
        # the length of possible positions
        lpp = str_length(data1[k]) - J + 1
        # the posterior probability of possible positions
        ppp = c()
        for (i in 1:lpp) {
          
          temp_Ak = matrix(i:(i+J-1),nrow=1)
          temp_sel = sel(data1[k],temp_Ak)
          temp_Xk = cal_X(temp_sel, theta)
          
          temp0 = (Y[k]-temp_Xk%*%beta)**2
          ttpp = exp(-temp0/(2*sigma2))
          
          result = ttpp
          ppp = c(ppp,result)
        }
        ppp = ppp/sum(ppp)
        new_ak = sample(1:lpp,1, prob = ppp)
        nb_aprob = nb_aprob + log(ppp[new_ak])
        A[k,] = seq(new_ak,length.out = J, by=1)
      }
    }
    
    Ad = A[binding_index,]
    Xd = cal_X(sel(data1[binding_index],A[binding_index,]), theta)
    X = cal_X(sel(data1,A), theta)
    Yd = Y[binding_index]
    
    
    # non-binding shift
    if(iter%%5==0&iter<=end_times&iter>0){
      # update non-binding A, theta, beta
      theta_star = theta
      # print(iter)
      temp_ratio = 0
      for (j in 1:J) {
        haj = H(paste(sel(data_binding,Ad[,j]),collapse = ''))
        theta_star[,j] = MCMCpack::rdirichlet(1,haj+eta_0[,j])
        theta_star[,j] = theta_star[,j]/sum(theta_star[,j])
        temp_ratio = temp_ratio-sum(log(MCMCpack::ddirichlet(theta_star[,j],haj+eta_0[,j]))) + 
          sum(log(MCMCpack::ddirichlet(theta[,j],haj+eta_0[,j])))
      }
      sel_element_binding = sel(data_binding,Ad)
      X_star = cal_X(sel_element_binding, theta_star)
      
      # update beta
      # non-corr prior
      new_XX = t(X_star)%*%X_star
      new_tttmp = MASS::ginv(MASS::ginv(Sigma_0)+new_XX/sigma2)
      new_mu = new_tttmp%*%(MASS::ginv(Sigma_0)%*%beta_0+t(X_star)%*%Yd/sigma2)
      new_Sigma = new_tttmp
      
      # non-negative
      lb=c(-Inf,rep(0,nrow(Sigma)-1))
      beta_star = TruncatedNormal::rtmvnorm(1,new_mu,new_Sigma,lb=lb)
      
      old_XX = t(Xd)%*%Xd
      old_tttmp = MASS::ginv(MASS::ginv(Sigma_0)+old_XX/sigma2)
      old_mu = old_tttmp%*%(MASS::ginv(Sigma_0)%*%beta_0+t(Xd)%*%Yd/sigma2)
      old_Sigma = old_tttmp
      
      # print(temp_ratio)
      temp_ratio = temp_ratio - sum(TruncatedNormal::dtmvnorm(beta_star, new_mu, new_Sigma, lb,log=T)) +
        sum(TruncatedNormal::dtmvnorm(beta, old_mu, old_Sigma, lb,log=T))
      # print(temp_ratio)
      new_A = A
      new_nb_aprob = 0
      for (k in (d+1):n) {
        # the length of possible positions
        lpp = str_length(data1[k]) - J + 1
        # the posterior probability of possible positions
        ppp = c()
        for (i in 1:lpp) {
          
          temp_Ak = matrix(i:(i+J-1),nrow=1)
          temp_sel = sel(data1[k],temp_Ak)
          temp_Xk = cal_X(temp_sel, theta_star)
          
          temp0 = (Y[k]-temp_Xk%*%beta_star)**2
          ttpp = exp(-temp0/(2*sigma2))
          
          result = ttpp
          ppp = c(ppp,result)
        }
        ppp = ppp/sum(ppp)
        new_ak = sample(1:lpp,1, prob = ppp)
        new_nb_aprob = new_nb_aprob + log(ppp[new_ak])
        new_A[k,] = seq(new_ak,length.out = J, by=1)
      }
      
      temp_ratio = temp_ratio - new_nb_aprob + nb_aprob
      pa = log_joint_post(new_A,theta_star,theta_0,beta_star,sqrt(sigma2),d)
      pb = log_joint_post(A,theta,theta_0,beta,sqrt(sigma2),d)
      temp_ar = pa - pb + temp_ratio
      ar = min(0,temp_ar) 
      # print(paste(delta,'first_shift_rate=',exp(ar)))
      u = runif(1,0,1)
      u = log(u)
      if(u<=ar){
        print(paste('non-binding shift success,prob=',exp(ar)))
        theta = theta_star
        A = new_A
        beta = beta_star
        sel_element = sel(data1,A)
        X = cal_X(sel_element, theta)
        nb_aprob = new_nb_aprob
        old_XX = new_XX
        old_mu = new_mu
        old_Sigma = new_Sigma
      }
    }
    
    
    # Mix step
    if(iter%%5==0&iter<=end_times&iter>0){
      delta = sample(c(-1,1),size = 1,prob=rep(1/2,2))
      print(paste('delta =',delta))
      new_Ad = matrix(nrow = d, ncol = motif_len)
      tmp_Ad = matrix(nrow = d, ncol = motif_len)
      
      p_tmp_Ad = matrix(nrow = d, ncol = motif_len)
      tmp_A = A[1:d,ncol(A)]+1
      max_length = simplify2array(lapply(data1[1:d],function(x) str_length(x)))
      tmp_A[tmp_A>max_length] = max_length[tmp_A>max_length]
      for (i in 1:d) {
        p_tmp_Ad[i,] = seq(tmp_A[i]-motif_len+1,length.out = motif_len, by=1)
      }
      
      n_tmp_Ad = matrix(nrow = d, ncol = motif_len)
      tmp_A = A[1:d,1]-1
      tmp_A[tmp_A<1] = 1
      for (i in 1:d) {
        n_tmp_Ad[i,] = seq(tmp_A[i],length.out = motif_len, by=1)
      }
      
      if(delta>0){
        tmp_Ad = p_tmp_Ad
        other_Ad = n_tmp_Ad
      }
      if(delta<0){
        tmp_Ad = n_tmp_Ad
        other_Ad = p_tmp_Ad
      }
      
      
      theta_star = theta
      temp_ratio = 0
      ratio_theta = 0
      for (j in 1:J) {
        haj = H(paste(sel(data_binding,tmp_Ad[,j]),collapse = ''))
        other_haj = H(paste(sel(data_binding,other_Ad[,j]),collapse = ''))
        theta_star[,j] = MCMCpack::rdirichlet(1,haj+eta_0[,j])
        theta_star[,j] = theta_star[,j]/sum(theta_star[,j])
        ratio_theta = ratio_theta-sum(log(MCMCpack::ddirichlet(theta_star[,j],haj+eta_0[,j]) + 
                                            MCMCpack::ddirichlet(theta_star[,j],other_haj+eta_0[,j])))
      }
      
      sel_element_binding = sel(data_binding,tmp_Ad)
      
      hall = H(paste(data1, collapse = ''))
      ha = H(paste(sel_element_binding,collapse = ''))
      hac = hall - ha
      
      # update theta_0
      theta_0_star = MCMCpack::rdirichlet(1,hac+gamma_0)
      
      ratio_theta0 = - sum(log(MCMCpack::ddirichlet(theta_0_star,hac+gamma_0))) 
      
      # update binding A
      new_b_aprob = 0
      for (k in 1:d) {
        # the length of possible positions
        lpp = str_length(data1[k]) - J + 1
        # the posterior probability of possible positions
        ppp = c()
        for (i in 1:lpp) {
          temp_Ak = matrix(i:(i+J-1),nrow=1)
          temp_sel = sel(data1[k],temp_Ak)
          temp_Xk = cal_X(temp_sel, theta_star)
          
          ratio = theta_star/theta_0_star[1,]
          row.names(ratio) = dict
          hakj = cal_X(temp_sel,ratio)
          hakj = prod(hakj)
          result = hakj
          ppp = c(ppp,result)
        }
        ppp = ppp/sum(ppp)
        new_ak = sample(1:lpp,1, prob = ppp)
        new_b_aprob = new_b_aprob +log(ppp[new_ak])
        new_Ad[k,] = seq(new_ak,length.out = J, by=1)
      }
      
      b_aprob = 0
      for (k in 1:d) {
        # the length of possible positions
        lpp = str_length(data1[k]) - J + 1
        # the posterior probability of possible positions
        ppp = c()
        for (i in 1:lpp) {
          temp_Ak = matrix(i:(i+J-1),nrow=1)
          temp_sel = sel(data1[k],temp_Ak)
          temp_Xk = cal_X(temp_sel, theta)
          
          ratio = theta/theta_0[1,]
          row.names(ratio) = dict
          hakj = cal_X(temp_sel,ratio)
          hakj = prod(hakj)
          result = hakj
          ppp = c(ppp,result)
        }
        ppp = ppp/sum(ppp)
        b_aprob = b_aprob +log(ppp[A[k,1]])
      }
      ratio_Ad = - new_b_aprob + b_aprob
      temp_ratio = temp_ratio + ratio_Ad
      
      sel_element_binding = sel(data_binding, new_Ad)
      X_star = cal_X(sel_element_binding, theta_star)
      
      
      # reverse Ad
      rtmp_Ad = matrix(nrow = d, ncol = motif_len)
      
      rp_Ad = matrix(nrow = d, ncol = motif_len)
      tmp_A = new_Ad[1:d,ncol(A)]+1
      max_length = simplify2array(lapply(data1[1:d],function(x) str_length(x)))
      tmp_A[tmp_A>max_length] = max_length[tmp_A>max_length]
      for (i in 1:d) {
        rp_Ad[i,] = seq(tmp_A[i]-motif_len+1,length.out = motif_len, by=1)
      }
      
      rn_Ad = matrix(nrow = d, ncol = motif_len)
      tmp_A = new_Ad[1:d,1]-1
      tmp_A[tmp_A<1] = 1
      for (i in 1:d) {
        rn_Ad[i,] = seq(tmp_A[i],length.out = motif_len, by=1)
      }
      
      if(-delta>0){
        rtmp_Ad = rp_Ad
        rother_Ad = rn_Ad
      }
      if(-delta<0){
        rtmp_Ad = rn_Ad
        rother_Ad = rp_Ad
      }
      
      for (j in 1:J) {
        old_haj = H(paste(sel(data_binding,rtmp_Ad[,j]),collapse = ''))
        other_haj = H(paste(sel(data_binding,rother_Ad[,j]),collapse = ''))
        ratio_theta = ratio_theta+sum(log(MCMCpack::ddirichlet(theta[,j],old_haj+eta_0[,j])+
                                            MCMCpack::ddirichlet(theta[,j],other_haj+eta_0[,j])))
      }
      
      temp_ratio = ratio_theta
      
      old_ha = H(paste(sel(data1[binding_index],rtmp_Ad),collapse = ''))
      old_hac = hall - old_ha
      ratio_theta0 = ratio_theta0 + sum(log(MCMCpack::ddirichlet(theta_0,old_hac+gamma_0)))
      temp_ratio = temp_ratio + ratio_theta0
      
      # update beta
      new_XX = t(X_star)%*%X_star
      new_tttmp = MASS::ginv(MASS::ginv(Sigma_0)+new_XX/sigma2)
      new_mu = new_tttmp%*%(MASS::ginv(Sigma_0)%*%beta_0+t(X_star)%*%Yd/sigma2)
      new_Sigma = new_tttmp
      lb=c(-Inf,rep(0,nrow(Sigma)-1))
      beta_star = TruncatedNormal::rtmvnorm(1,new_mu,new_Sigma,lb=lb)
      
      
      # update sigma2
      new_sel_element = sel(data1[1:d],new_Ad)
      tmp_X = cal_X(new_sel_element, theta_star)
      rr = t(Yd-tmp_X%*%beta_star)%*%(Yd-tmp_X%*%beta_star)
      old_rr = t(Yd-Xd%*%beta)%*%(Yd-Xd%*%beta)
      Isigma2 = rgamma(1,(alpha_0+n)/2,rate = (T_0+rr)/2)
      sigma2_star = 1/Isigma2
      ratio_sigma2 = - sum(dgamma(1/sigma2_star,(alpha_0+n)/2,rate = (T_0+rr)/2,log=T)) +
        sum(dgamma(1/sigma2, (alpha_0+n)/2,rate = (T_0+old_rr)/2,log=T))
      temp_ratio = temp_ratio + ratio_sigma2
      
      old_XX = t(Xd)%*%Xd
      old_tttmp = MASS::ginv(MASS::ginv(Sigma_0)+old_XX/sigma2_star)
      old_mu = old_tttmp%*%(MASS::ginv(Sigma_0)%*%beta_0+t(Xd)%*%Yd/sigma2_star)
      old_Sigma = old_tttmp
      
      ratio_beta = - sum(TruncatedNormal::dtmvnorm(beta_star, new_mu, new_Sigma, lb,log=T)) +
        sum(TruncatedNormal::dtmvnorm(beta, old_mu, old_Sigma, lb,log=T))
      temp_ratio = temp_ratio + ratio_beta
      
      
      # update non-binding locations
      new_A = A
      new_A[1:d,] = new_Ad
      new_nb_aprob = 0
      for (k in (d+1):n) {
        # the length of possible positions
        lpp = str_length(data1[k]) - J + 1
        # the posterior probability of possible positions
        ppp = c()
        for (i in 1:lpp) {
          
          temp_Ak = matrix(i:(i+J-1),nrow=1)
          temp_sel = sel(data1[k],temp_Ak)
          temp_Xk = cal_X(temp_sel, theta_star)
          
          temp0 = (Y[k]-temp_Xk%*%beta_star)**2
          ttpp = exp(-temp0/(2*sigma2_star))
          result = ttpp
          ppp = c(ppp,result)
        }
        ppp = ppp/sum(ppp)
        new_ak = sample(1:lpp,1, prob = ppp)
        new_nb_aprob = new_nb_aprob + log(ppp[new_ak])
        new_A[k,] = seq(new_ak,length.out = J, by=1)
      }
      ratio_And = - new_nb_aprob + nb_aprob
      temp_ratio = temp_ratio + ratio_And
      
      
      
      # target distribution
      pa = log_joint_post(new_A,theta_star,theta_0_star,beta_star,sqrt(sigma2_star),d,print_part = F)
      pb = log_joint_post(A,theta,theta_0,beta,sqrt(sigma2),d,print_part = F)
      
      temp_ar = pa - pb + temp_ratio
      ar = min(0,temp_ar) 
      u = runif(1,0,1)
      u = log(u)
      if(u<=ar){
        print(paste('Mix step',delta ,'shift success,prob=',exp(ar)))
        old_mu = new_mu
        old_Sigma = new_Sigma
        sigma2 = sigma2_star
        nb_aprob = new_nb_aprob
        theta = theta_star
        theta_0 = theta_0_star
        beta = beta_star
        A = new_A
        sel_element = sel(data1,A)
        X = cal_X(sel_element, theta)
      }
    }
    
    # update d
    # possible d = 5,...,n
    pp = c()
    for (k in d_index) {
      temp = 0
      temp_binding = 1:k
      for (j in 1:motif_len) {
        temp = temp + sum(H(paste(sel(data1[temp_binding],A[temp_binding,j]),
                                  collapse = ''))*log(theta[,j]))
      }
      if(k<n){
        temp = temp + sum(H(paste(sel(data1[-temp_binding],A[-temp_binding,]),
                                  collapse = ''))*log(theta_0))
      }
      
      pp = c(pp, temp)
    }
    pp = pp - max(pp)
    pp = exp(pp)
    pp = pp/sum(pp)
    d = sample(d_index,size=1,replace=T,prob=pp)
    binding_index = 1:d
    
    
    
    
    # final distribution
    
    if(iter%%interval_sel==0){
      ttemp = X%*%beta
      temp_score = ttemp
      temp_coe = cor(x = as.numeric(temp_score), y = as.numeric(Y), method = c("spearman"))
      temp_MSE = sum((temp_score-Y)**2)/length(Y)
      
      
      temp_joint = log_joint_post(A,theta,theta_0,beta,sqrt(sigma2),d,print_part = F)
      
      result_log_joint = c(result_log_joint, temp_joint)
      result_MSE = c(result_MSE, temp_MSE)
      
      result_coe = c(result_coe, temp_coe)
      write.csv(X,paste(dir_names,'temp_X.csv',sep=''))
      write.csv(beta,paste(dir_names,'temp_beta.csv',sep=''))
      print(paste('threshold=',d))
      
      if(iter>burn_in&iter%%interval_sel==0){
        result_a = result_a + A[,1]
        result_d_trace = c(result_d_trace, d)
        result_beta = result_beta + beta
        result_theta = result_theta + theta
        result_theta_0 = result_theta_0 + theta_0 
        result_score = result_score + temp_score
        result_sigma2 = result_sigma2 + sigma2
      }
    }
    
    iter = iter + 1
    if(iter>end_times){
      break
    }
  }
  
  result_beta = result_beta/((end_times-burn_in)/interval_sel)
  result_theta = result_theta/((end_times-burn_in)/interval_sel)
  result_theta_0 = result_theta_0/((end_times-burn_in)/interval_sel)
  
  g_result_beta <<- result_beta
  g_result_theta <<- result_theta
  g_result_theta_0 <<- result_theta_0
  
  result_score = result_score/((end_times-burn_in)/interval_sel)
  result_sigma2 = result_sigma2/((end_times-burn_in)/interval_sel)
  result_a = result_a/((end_times-burn_in)/interval_sel)
  
  for (l in 1:length(result_a)) {
    if(result_a[l] - trunc(result_a[l])>0.5){
      result_a[l] = trunc(result_a[l])+1
    }else{
      result_a[l] = trunc(result_a[l])
    }
  }
  
  result_A = matrix(nrow = length(data1), ncol = J)
  for (i in 1:length(a)) {
    result_A[i,] = seq(result_a[i],length.out = J, by=1)
  }
  
  
  result_d = mean(result_d_trace)
  if(result_d - trunc(result_d)>0.5){
    result_d = trunc(result_d)+1
  }else{
    result_d = trunc(result_d)
  }
  result_binding_index = 1:result_d
  result_threshold = Y[result_d]
  est_label = rep(0,n)
  est_label[1:result_d] = 1
  
  
  write.csv(result_threshold,paste(dir_names,'result_threshold.csv',sep=''))
  
  # the prob for each seq to be binding
  dd = rep(0,n)
  whole_index = 1:n
  for (ii in result_d_trace) {
    temp_index = whole_index<=ii
    dd = dd + rep(1,n)*temp_index
  }
  dd = dd/length(result_d_trace)
  plot(dd,type='l',main = 'Binding_prob',
       xlab = 'seq_index', ylab = 'Prob')
  {
    sel_element = sel(data1[result_binding_index],result_A[result_binding_index,])
    temp_fig = ggseqlogo(sel_element)
    temp_fig = temp_fig + labs(title = "final_logo")+theme(plot.title = element_text(hjust = 0.5))
    print(temp_fig)
    
    
    # 
    write.csv(result_d_trace,paste(dir_names,'d_trace.csv',sep=''))
    write.csv(result_log_joint,paste(dir_names,'log_joint_post_path.csv',sep=''))
    write.csv(result_coe,paste(dir_names,'coe_path.csv',sep=''))
    write.csv(result_MSE,paste(dir_names,'MSE_path.csv',sep=''))
    
    plot(result_log_joint,type='l',main = 'Joint_path',
         xlab = 'Iteration', ylab = 'joint_post')
    
    
    plot(result_coe,type='l',main = 'Coefficient_path',
         xlab = 'Iteration', ylab = 'Coefficient')
    
    plot(result_MSE,type='l',main = 'MSE_path',
         xlab = 'Iteration', ylab = 'MSE')
    
    result_X = sel(data1,result_A)
    result_X = cal_X(result_X, result_theta)
    
    write.csv(result_X,paste(dir_names,'result_X.csv',sep=''))
    write.csv(result_A,paste(dir_names,'result_A.csv',sep=''))
    write.csv(result_beta,paste(dir_names,'result_beta.csv',sep=''))
    write.csv(result_theta,paste(dir_names,'result_theta.csv',sep=''))
    write.csv(result_sigma2,paste(dir_names,'result_sigma2.csv',sep=''))
    write.csv(result_theta_0,paste(dir_names,'result_theta_0.csv',sep=''))
    write.csv(result_score,paste(dir_names,'result_score.csv',sep=''))
    
    
    
    est_label = rep(0,n)
    est_label[1:result_d] = 1
    final_result = cbind(data1,Y,result_score,dd)
    colnames(final_result) = c('pep','score','est','binding_prob')
    write.csv(final_result,paste(dir_names,'final_result.csv',sep=''))
    
    # plot
    final_result = as.data.frame(final_result,stringsAsFactors=F)
    final_result$score = as.numeric(final_result$score)
    final_result$est = as.numeric(final_result$est)
    
    final_result = final_result[order(final_result$score,decreasing = F),]
    p<-ggplot(data=final_result,aes(x=score,y=est))+geom_point()+ geom_smooth(method="lm")+
      geom_hline(aes(yintercept=result_threshold), colour="red", linetype="dashed")+
      geom_vline(aes(xintercept=result_threshold), colour="red", linetype="dashed")+
      labs(title = "real_est")+theme(plot.title = element_text(hjust = 0.5))
    print(p)
    
    
    pearson = cor(x = final_result$score, y = final_result$est,
                  method = c("pearson"))
    
    kendall = cor(x = final_result$score, y = final_result$est,
                  method = c("kendall"))
    
    spearman = cor(x = final_result$score, y = final_result$est,
                   method = c("spearman"))
    
    result_coefficient = c(pearson =pearson, kendall = kendall,spearman = spearman)
    write.csv(result_coefficient,paste(dir_names,'result_coefficient.csv',sep=''))
  } #analysis output save
  dev.off()
  return(final_result)
}








