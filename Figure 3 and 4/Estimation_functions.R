##y=log(AAR)
##x=Temperature
##z=Age

##Function for least square estimation  
LSE=function(y,x,z)
{
  log_g_hat=pava(y,decreasing =  T)
  model=nls(y-log_g_hat~log(1+v*x),start = list(v=sum((exp(y-log_g_hat)-1)*x)/sum(x*x)) )
  gamma_hat=as.numeric(coef(model))
  Q=mean((y-log(1+gamma_hat*x)-log_g_hat)^2)
  for( i in 2:10000)
  {
    log_g_hat=pava(y-log(1+gamma_hat*x),decreasing =  T)
    model=nls(y-log_g_hat~log(1+v*x),start = list(v=gamma_hat))
    gamma_hat=as.numeric(coef(model))  
    Q[i]=mean((y-log(1+gamma_hat*x)-log_g_hat)^2)
    if(i>3 & abs(Q[i]-Q[i-1])<10^(-10) )
    {
      break
    }
  }
  Estimates=list(gamma_hat=gamma_hat,
                 sig2_hat=Q[i],
                 log_g_hat=log_g_hat)
  return(Estimates)
  
}

## Function for smooth estimation using LSE
Estimates=function(h,y,x,z)
{
  
  A=LSE(y,x,z)
  n=length(y)
  tilde_log_g=kreg(z, A$log_g_hat, kernel = c("epanechnikov"), bandwidth = h,grid=z)$y  
  model=nls(y-tilde_log_g~log(1+v*x),start = list(v=A$gamma_hat) )
  tilde_gamma=as.numeric(coef(model))
  tilde_res=y-log(1+tilde_gamma*x)-tilde_log_g
  tilde_phi=sum(tilde_res[2:n]*tilde_res[1:(n-1)])/sum(tilde_res[2:n]^2)
  tilde_sig2=sum((tilde_res[2:n]-tilde_phi*tilde_res[1:(n-1)])^2)/n  
  tilde_innov=tilde_res[2:n]-tilde_phi*tilde_res[1:(n-1)]
  Estimates=list(tilde_gamma=tilde_gamma,
                 tilde_log_g=tilde_log_g,
                 tilde_res=tilde_res,
                 tilde_phi=tilde_phi,
                 tilde_sig2=tilde_sig2,
                 tilde_innov=tilde_innov)
  return(Estimates)
  
}  

##Function for bootstrap samples and uncertainity computation
Estimates_Bootstrap=function(x,z,tilde_gamma,tilde_log_g,tilde_innov,tilde_phi)
{
  n=length(z)
  M=1000
  tilde_gamma_star=matrix(nrow=M,ncol=1)
  tilde_log_g_star=matrix(nrow=M,ncol=n)
  tilde_phi_star=matrix(nrow=M,ncol=1)
  tilde_sig2_star=matrix(nrow=M,ncol=1)
  predicted_Adjusted_AAR_star=matrix(nrow=M,ncol=n)
  predicted_log_AAR_star=matrix(nrow=M,ncol=n)
  j=0
  while (j <M)
  {
  tilde_innov_star=sample(tilde_innov,n,replace = T)
  if (var(tilde_innov_star)>var(tilde_innov))
  {
    eps_star=arima.sim(list(order=c(1,0,0),ar=tilde_phi),
                     n,innov =  tilde_innov_star,
                     n.start=1000,
                     start.innov=sample(tilde_innov,1000,replace = T)
                     )
    y_star=log(1+tilde_gamma*x)+tilde_log_g+eps_star
    E_star=Estimates(h,y_star,x,z)
    tilde_gamma_star[(j+1)]=E_star$tilde_gamma
    tilde_log_g_star[(j+1),]=E_star$tilde_log_g
    tilde_phi_star[(j+1)]=E_star$tilde_phi
    tilde_sig2_star[(j+1)]=E_star$tilde_sig2
    resid_star=y_star-log(1+E_star$tilde_gamma*x)-E_star$tilde_log_g
    predicted_resid_star=E_star$tilde_phi*c(0,resid_star[1:(n-1)])
    predicted_Adjusted_AAR_star[(j+1),]=(1+E_star$tilde_gamma*x)*exp(predicted_resid_star)
    predicted_log_AAR_star[(j+1),]=log((1+E_star$tilde_gamma*x))+E_star$tilde_log_g+predicted_resid_star
    j=j+1
  }
  }
  
  gamma_quantile_025_star= quantile(tilde_gamma_star,0.025)
  gamma_quantile_975_star= quantile(tilde_gamma_star,0.975)
  gamma_var_star=var(tilde_gamma_star)
  
  log_g_quantile_025_star=colQuantiles(tilde_log_g_star,probs=0.025)
  log_g_quantile_975_star=colQuantiles(tilde_log_g_star,probs=0.975)
  log_g_var_star=colVars(tilde_log_g_star)
  
  phi_quantile_025_star= quantile(tilde_phi_star,0.025)
  phi_quantile_975_star= quantile(tilde_phi_star,0.975)
  phi_var_star=var(tilde_phi_star)
  
  sig2_quantile_025_star= quantile(tilde_sig2_star,0.025)
  sig2_quantile_975_star= quantile(tilde_sig2_star,0.975)
  sig2_var_star=var(tilde_sig2_star)
  
  predicted_Adjusted_AAR_quantile_025_star=colQuantiles(predicted_Adjusted_AAR_star,probs=0.025)
  predicted_Adjusted_AAR_quantile_975_star=colQuantiles(predicted_Adjusted_AAR_star,probs=0.975)
  
  predicted_log_AAR_quantile_025_star=colQuantiles(predicted_log_AAR_star,probs=0.025)
  predicted_log_AAR_quantile_975_star=colQuantiles(predicted_log_AAR_star,probs=0.975)
  
  
  Bootstrap_estimates=list(   gamma_quantile_025_star=gamma_quantile_025_star,
                              gamma_quantile_975_star=gamma_quantile_975_star,
                              gamma_var_star=gamma_var_star,
                              log_g_quantile_025_star=log_g_quantile_025_star,
                              log_g_quantile_975_star=log_g_quantile_975_star,
                              log_g_var_star=log_g_var_star,
                              phi_quantile_025_star=phi_quantile_025_star,
                              phi_quantile_975_star=phi_quantile_975_star,
                              phi_var_star=phi_var_star,
                              sig2_quantile_025_star=sig2_quantile_025_star,
                              sig2_quantile_975_star=sig2_quantile_975_star,
                              sig2_var_star=sig2_var_star,
                              predicted_Adjusted_AAR_quantile_025_star=predicted_Adjusted_AAR_quantile_025_star,
                              predicted_Adjusted_AAR_quantile_975_star=predicted_Adjusted_AAR_quantile_975_star,
                              predicted_log_AAR_quantile_025_star=predicted_log_AAR_quantile_025_star,
                              predicted_log_AAR_quantile_975_star=predicted_log_AAR_quantile_975_star
                            )
  
  return(Bootstrap_estimates)
  
}
