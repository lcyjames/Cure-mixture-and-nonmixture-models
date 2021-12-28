#wmcmEM is the estimation procedure of the Weibull Mixture Cure Model based on the EM algorithm
#wnmcmEM is the estimation procedure of the Weibull Non-Mixture Cure Model based on the EM algorithm
#Yi is a vector of the right censoring times, with size n 
#cen is the corresponding censoring indicator with 1 being cases and 0 being censored, with size n
#X is a covariate matrix for the incidence component with size n times the number of covariates
#Z is a covariate matrix for the latency component with size n times the number of covariates
#X and Z should not contain a column of 1 (i.e. no intercept required) but they only contain the covariates from the data
#trace=FALSE by default. For tracking the converging path of the parameter estimation, set trace=TRUE 
#tolerance is the converging criteria typically assigned to be 0.0001

wmcmEM<-function(Yi,cen,X,Z,trace=FALSE,tolerance=10^{-4}){
  n         <-length(Yi)
  X         <-cbind(rep(1,n),X)
  Z         <-cbind(rep(1,n),Z)
  p         <-ncol(X)
  q         <-ncol(Z)
  b.d       <-1
  alpha.d   <-rep(0,p)         
  beta.d    <-rep(0,q) 
  times<-0
  
  repeat{
    alphaX    <-apply(matrix(rep(alpha.d,n),nrow=n,byrow=TRUE)*X,MARGIN = 1,sum)
    betaZ     <-apply(matrix(rep(beta.d,n) ,nrow=n,byrow=TRUE)*Z,MARGIN = 1,sum)
    
    #E-step
    GSI <-exp(alphaX)/(1+exp(alphaX))
    Su  <-exp(-(Yi)^b.d*exp(betaZ))
    Eu  <-cen+(1-cen)*(GSI*Su)/(1-GSI+GSI*Su)
    
    b.dummy<-b.d; alpha.dummy<-alpha.d; beta.dummy <-beta.d
    
    alphaX      <-apply(matrix(rep(alpha.d,n),nrow=n,byrow=TRUE)*X,MARGIN = 1,sum)
    GSI         <-exp(alphaX)/(1+exp(alphaX))
    GSId1       <-exp(alphaX)/(1+exp(alphaX))^2
    
    dQ1a   <-sapply(c(1:p),function(o){sum(Eu*X[,o]-GSI*X[,o])})
    dQ1a2  <-sapply(c(1:p),function(o1){sapply(c(1:p),function(o2){-sum(X[,o1]*X[,o2]*GSId1)})})
    alpha.d<- c(alpha.d)-solve(dQ1a2)%*%dQ1a
    
    dLp_beta<-c(sum(cen*(1/b.d+log(Yi))-Eu*Yi^{b.d}*log(Yi)*exp(betaZ)),
                sapply(c(1:q),function(o){sum(cen*Z[,o]-(Eu*Yi^b.d*exp(betaZ)*Z[,o]))}))
    
    dLp2<-matrix(0,nrow=q+1, ncol=q+1)
    dLp2[1,1]<-sum((-cen/b.d^2)-Eu*Yi^{b.d}*(log(Yi))^2*exp(betaZ))
    dLp2[1,c(2:(q+1))]<-dLp2[c(2:(q+1)),1] <-sapply(c(1:q),function(o){-sum(Z[,o]*Eu*(Yi)^{b.d}*log(Yi)*exp(betaZ))})
    dLp2[c(2:(q+1)),c(2:(q+1))]<-sapply(c(1:q),function(o1){sapply(c(1:q),function(o2){-sum(Eu*Yi^{b.d}*exp(betaZ)*Z[,o1]*Z[,o2])})})
    
    update <-as.numeric(c(b.d,beta.d)-solve(dLp2)%*%dLp_beta)
    b.d    <-update[1]; beta.d<-update[2:(q+1)]
    
    dist<-max(abs(c(b.d-b.dummy, alpha.d-alpha.dummy, beta.d-beta.dummy)))
    
    times<-times+1
    
    if(trace==TRUE){print(round(c(times, b.d, alpha.d, beta.d, dist),4))}
    if(dist<tolerance){
      est     <-c(b.d, alpha.d, beta.d)
      GSI     <-exp(alphaX)/(1+exp(alphaX))
      GSId1   <-exp(alphaX)/(1+exp(alphaX))^2
      GSId2   <- -(exp(alphaX)*(exp(alphaX)-1))/(1+exp(alphaX))^3
      den     <-(1-GSI+GSI*exp(-Yi^b.d*exp(betaZ)))
      
      #Fisher information calculated based on the score function
      score.matrix<-cbind( (cen*(b.d^{-1}+log(Yi)-Yi^{b.d}*log(Yi)*exp(betaZ))+
                              (1-cen)*den^{-1}*GSI*exp(-Yi^b.d*exp(betaZ))*-Yi^b.d*log(Yi)*exp(betaZ)),
                           sapply(c(1:p),function(o){cen*GSI^{-1}*GSId1*X[,o]+(1-cen)*den^{-1}*(-GSId1+GSId1*exp(-Yi^b.d*exp(betaZ)))*X[,o]}),
                           sapply(c(1:q),function(o){cen*(1-Yi^{b.d}*exp(betaZ))*Z[,o]+(1-cen)*den^{-1}*(GSI*exp(-Yi^b.d*exp(betaZ))*(-Yi^b.d*exp(betaZ)))*Z[,o]}))
      
      #apply(score.matrix,2,sum)
      hess<-matrix(0,nrow=(p+q+1), ncol=(p+q+1))
      for(J in c(1:n)){
        hess<-hess+score.matrix[J,]%o%score.matrix[J,]
      }
      se          <-sqrt(diag(solve(hess)))
      
      est_L_95    <-est-1.96*se
      est_R_95    <-est+1.96*se
      exp_est     <-NA.vec<-as.numeric(rep(NA,3))
      if((p==1)&(q==1)){coef<-c("shape_Weibull","intercept_logit","intercept_surv")}else if(p==1){
        exp_est<-append(exp_est,exp(beta.d[2:q]),p+2);NA.vec<-append(NA.vec,rep(1,q-1),p+2);
        coef<-c("shape_Weibull","intercept_logit","intercept_surv",paste0("beta",c(1:(q-1))))}else if(q==1){
          exp_est<-append(exp_est,exp(alpha.d[2:p]),2);NA.vec<-append(NA.vec,rep(1,p-1),2);
          coef<-c("shape_Weibull","intercept_logit",paste0("alpha",c(1:(p-1))),"intercept_surv")}else{
            exp_est<-c(NA,NA,exp(alpha.d[2:p]),NA,exp(beta.d[2:q]));NA.vec<-c(NA,NA,rep(1,p-1),NA,rep(1,q-1));
            coef<-c("shape_Weibull","intercept_logit",paste0("alpha",c(1:(p-1))),"intercept_surv",paste0("beta",c(1:(q-1))))}
      
      exp_est_L_95<-exp(est_L_95)*NA.vec
      exp_est_R_95<-exp(est_R_95)*NA.vec
      pv          <-pchisq(q = (est/se)^2,df = 1,lower.tail = FALSE)*NA.vec
      
      results<-cbind(coef,round(data.frame(est=est, se=se, est_L_95=est_L_95, est_R_95=est_R_95, 
                                           exp_est=exp_est, exp_est_L_95=exp_est_L_95,exp_est_R_95=exp_est_R_95,pv=pv),4))
      
      return(results)
      break}
  }
}

wnmcmEM<-function(Yi,cen,X,Z,trace=FALSE,tolerance=10^{-4}){
  n         <-length(Yi)
  X         <-cbind(rep(1,n),X)
  Z         <-cbind(rep(1,n),Z)
  p         <-ncol(X)
  q         <-ncol(Z)
  b.d       <-1
  alpha.d   <-rep(0,p)         
  beta.d    <-rep(0,q) 
  times<-0
  
  repeat{
    alphaX    <-apply(matrix(rep(alpha.d,n),nrow=n,byrow=TRUE)*X,MARGIN = 1,sum)
    betaZ     <-apply(matrix(rep(beta.d,n) ,nrow=n,byrow=TRUE)*Z,MARGIN = 1,sum)
    
    b.dummy<-b.d; alpha.dummy<-alpha.d; beta.dummy <-beta.d
    
    #E-step
    Eyi<- exp(alphaX)*exp(-Yi^b.d*exp(betaZ))+cen
    
    #round(grad(func=comp.like, x = c(b.d,alpha.d,beta.d),Yi=Yi,cen=cen,X=X,Z=Z,Eyi=Eyi),4)
    #round(hessian(func=comp.like, x = c(b.d,alpha.d,beta.d),Yi=Yi,cen=cen,X=X,Z=Z,Eyi=Eyi),4)
    
    #M-step
    dl_db     <-sum(cen*(-Yi^b.d*log(Yi)*exp(betaZ)+b.d^{-1}+log(Yi))+(Eyi-cen)*(-Yi^b.d*log(Yi)*exp(betaZ)))
    dl_dalpha <-sapply(c(1:p),function(o){sum(X[,o]*(Eyi-exp(alphaX)))})
    dl_dbeta  <-sapply(c(1:q),function(o){sum(Z[,o]*(cen*(1-Yi^b.d*exp(betaZ))+(Eyi-cen)*(-Yi^b.d*exp(betaZ))))})
    
    dl_alpha2 <-sapply(c(1:p),function(o1){sapply(c(1:p),function(o2){sum(-exp(alphaX)*X[,o1]*X[,o2])})})
    alpha.d   <-alpha.d-solve(dl_alpha2)%*%dl_dalpha
    
    dL2<-matrix(0,nrow=q+1,ncol=q+1)
    dL2[1,1]                    <-sum(cen*(-Yi^b.d*(log(Yi))^2*exp(betaZ)-b.d^{-2})+(Eyi-cen)*(-Yi^b.d*(log(Yi))^2*exp(betaZ)))
    dL2[c(2:(q+1)),c(2:(q+1))]  <-sapply(c(1:q),function(o1){sapply(c(1:q),function(o2){
      sum(Z[,o1]*Z[,o2]*(cen*(-Yi^b.d*exp(betaZ))+(Eyi-cen)*(-Yi^b.d*exp(betaZ)))) })})
    dL2[1,c(2:(q+1))]<-
      dL2[c(2:(q+1)),1]<-sapply(c(1:q),function(o){sum(Z[,o]*(cen*(-Yi^b.d*log(Yi)*exp(betaZ))+(Eyi-cen)*(-Yi^b.d*log(Yi)*exp(betaZ))))})      
    
    update<-c(b.d,beta.d)-solve(dL2)%*%c(dl_db,dl_dbeta)
    b.d   <-update[1]
    beta.d<-update[2:(q+1)]
    
    dist<-max(abs(c(b.d-b.dummy, alpha.d-alpha.dummy, beta.d-beta.dummy)))
    
    times<-times+1
    
    if(trace==TRUE){print(round(c(times, b.d, alpha.d, beta.d, dist),4))}
    if(dist<tolerance){
      est       <-c(b.d, alpha.d, beta.d)
      alphaX    <-apply(matrix(rep(alpha.d,n),nrow=n,byrow=TRUE)*X,MARGIN = 1,sum)
      betaZ     <-apply(matrix(rep(beta.d,n) ,nrow=n,byrow=TRUE)*Z,MARGIN = 1,sum)
      
      #grad(func = obs.like.N,x = est, p=p, X=X, Z=Z, cen=cen, Yi=Yi)
      
      dSt_db    <- -exp(alphaX+betaZ)*log(Yi)*exp(-Yi^b.d*exp(betaZ))*(Yi^b.d)
      dSt_dalpha<- -exp(alphaX)*(1-exp(-Yi^b.d*exp(betaZ)))
      dSt_dbeta <- -exp(alphaX+betaZ)*exp(-Yi^b.d*exp(betaZ))*(Yi^b.d)
      
      score.matrix<-cbind((cen*(-Yi^b.d*log(Yi)*exp(betaZ)+b.d^{-1}+log(Yi))+dSt_db),
                          sapply(c(1:p),function(o){X[,o]*(cen+dSt_dalpha)}),
                          sapply(c(1:q),function(o){Z[,o]*(cen*(1-Yi^b.d*exp(betaZ))+dSt_dbeta)}))
      
      hess<-matrix(0,nrow=(p+q+1), ncol=(p+q+1))
      for(J in c(1:n)){
        hess<-hess+score.matrix[J,]%o%score.matrix[J,]
      }
      se          <-sqrt(diag(solve(hess)))
      
      est_L_95    <-est-1.96*se
      est_R_95    <-est+1.96*se
      exp_est     <-NA.vec<-as.numeric(rep(NA,3))
      if((p==1)&(q==1)){coef<-c("shape_Weibull","intercept_logit","intercept_surv")}else if(p==1){
        exp_est<-append(exp_est,exp(beta.d[2:q]),p+2);NA.vec<-append(NA.vec,rep(1,q-1),p+2);
        coef<-c("shape_Weibull","intercept_logit","intercept_surv",paste0("beta",c(1:(q-1))))}else if(q==1){
          exp_est<-append(exp_est,exp(alpha.d[2:p]),2);NA.vec<-append(NA.vec,rep(1,p-1),2);
          coef<-c("shape_Weibull","intercept_logit",paste0("alpha",c(1:(p-1))),"intercept_surv")}else{
            exp_est<-c(NA,NA,exp(alpha.d[2:p]),NA,exp(beta.d[2:q]));NA.vec<-c(NA,NA,rep(1,p-1),NA,rep(1,q-1));
            coef<-c("shape_Weibull","intercept_logit",paste0("alpha",c(1:(p-1))),"intercept_surv",paste0("beta",c(1:(q-1))))}
      
      exp_est_L_95<-exp(est_L_95)*NA.vec
      exp_est_R_95<-exp(est_R_95)*NA.vec
      pv          <-pchisq(q = (est/se)^2,df = 1,lower.tail = FALSE)*NA.vec
      
      results<-cbind(coef,round(data.frame(est=est, se=se, est_L_95=est_L_95, est_R_95=est_R_95, 
                                           exp_est=exp_est, exp_est_L_95=exp_est_L_95,exp_est_R_95=exp_est_R_95,pv=pv),4))
      
      return(results)
      break}
  }
}
