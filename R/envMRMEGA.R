#'@title The weighted environment-adjusted linear regression.
#'@description Fit the linear regression with weights using GWAS data across cohorts.
#'The contribution of each GWAS is weighted by the estimated inverse variance of the reference allele effect at the corresponding variant.
#'@param Y BETA values of the genetic variants across all cohorts.
#'@param X PCs plus environment covariates.
#'@param W The vector of weights, each component referring to the estimated inverse variance of the allele effect at the corresponding variants. The length of the vector is the number of cohorts.
#'@param n_pc The number of PCs.
#'@return Output the data frame containing information of fitting the regression model: the estimated coefficients, standard errors, the deviance of the different constrained  models and the corresponding degrees of freedom.
#'@author Siru Wang

lr_env_w<-function(Y,X,W,n_pc){
  #Y referes to cond.BETA
  #X refers to env+pcs
  #W refers to invSE2
  #n_pc: the number of PCs
  X=as.matrix(X)
  if((length(Y)-dim(X)[2])<1){
    df_lr_w=NA
  }else{
    m0<-lm(Y~1,weights=W)
    m1<-lm(Y~X,weights=W)
    #note this is temporary revised 250805
    m_env<-lm(Y~X[,-(((dim(X)[2])-n_pc+1):(dim(X)[2]))],weights=W)
    m_pc<-lm(Y~X[,(((dim(X)[2])-n_pc+1):(dim(X)[2]))],weights=W)

    RSS=sum(W*(m1$residuals^2))
    TSS0=sum(W*Y^2)
    TSS=sum(W*m0$residuals^2)
    #note this is temporary revised 250805
    RSS_env=sum(W*(m_env$residuals^2))
    RSS_pc=sum(W*(m_pc$residuals^2))

    #sum_m1<-summary(m1)
    #coef_se<-c(t(sum_m1$coefficients[,1:2]))
    #names(coef_se)=sprintf(rep(c("beta%s","se%s"),NROW(sum_m1$coefficients)),rep(0:(NROW(sum_m1$coefficients)-1),each=2))
    if(all(!is.na(m1$coefficients))){
      sum_m1<-summary(m1)
      coef_se<-c(t(sum_m1$coefficients[,1:2]))
    }else{
      coef_se<-rep(NA,2*(dim(X)[2]+1))
      sum_m1<-summary(m1)
      coef_se[seq(1,2*(dim(X)[2]+1),2)]<-m1$coefficients
      coef_se[seq(2,2*(dim(X)[2]+1),2)][!is.na(m1$coefficients)]<-sum_m1$coefficients[,2]
    }

    names(coef_se)=sprintf(rep(c("beta%s","se%s"),dim(X)[2]+1),rep(0:dim(X)[2],each=2))
    coef.se=as.data.frame(as.list(coef_se),col.names=names(coef_se))
    #revised 250805
    df_lr_w=data.frame(coef.se,RSS=RSS,df_RSS=m1$df.residual,TSS0=TSS0,TSS=TSS,
                       RSS_env=RSS_env,RSS_pc=RSS_pc)
  }
  return(df_lr_w)
}


#'@title The weighted linear regression.
#'@description Fit the linear regression with weights using GWAS data across cohorts.
#'The contribution of each GWAS is weighted by the estimated inverse variance of the reference allele effect at the corresponding variant.
#'@param Y BETA values of the genetic variants across all cohorts.
#'@param X PCs.
#'@param W The vector of weights, each component referring to the estimated inverse variance of the allele effect at the corresponding variants. The length of the vector is the number of cohorts.
#'@return Output the data frame containing information of fitting the regression model: the estimated coefficients, standard errors, the deviance of the different constrained  models and the corresponding degrees of freedom.
#'@author Siru Wang

lr_w<-function(Y,X,W){
  #Y referes to BETA
  #X refers to pcs
  #W refers to SE

  X=as.matrix(X)
  if((length(Y)-dim(X)[2])<1){
    df_lr_w=NULL
  }else{
    m0<-lm(Y~1,weights=W)
    m1<-lm(Y~X,weights=W)

    RSS=sum(W*(m1$residuals^2))
    TSS0=sum(W*Y^2)
    TSS=sum(W*m0$residuals^2)

    if(all(!is.na(m1$coefficients))){
      sum_m1<-summary(m1)
      coef_se<-c(t(sum_m1$coefficients[,1:2]))
    }else{
      coef_se<-rep(NA,2*(dim(X)[2]+1))
      sum_m1<-summary(m1)
      coef_se[seq(1,2*(dim(X)[2]+1),2)]<-m1$coefficients
      coef_se[seq(2,2*(dim(X)[2]+1),2)][!is.na(m1$coefficients)]<-sum_m1$coefficients[,2]
    }

    names(coef_se)=sprintf(rep(c("beta%s","se%s"),dim(X)[2]+1),rep(0:dim(X)[2],each=2))

    #sum_m1<-summary(m1)
    #coef_se<-c(t(sum_m1$coefficients[,1:2]))
    #names(coef_se)=sprintf(rep(c("beta%s","se%s"),NROW(sum_m1$coefficients)),rep(0:(NROW(sum_m1$coefficients)-1),each=2))
    coef.se=as.data.frame(as.list(coef_se),col.names=names(coef_se))
    df_lr_w=data.frame(coef.se,RSS=RSS,df_RSS=m1$df.residual,TSS0=TSS0,TSS=TSS)
  }
  return(df_lr_w)
}

#'@title environment-adjusted MR-MEGA approach or MR-MEGA approach.
#'@description Fit the (environment-adjusted) meta-regression model using the normalized statistic data across all cohorts.
#'@param beta_pop The estimated allelic effects of the target SNP conditioned on the SNP set across all the cohorts.
#'@param invse2_pop The inverse of the standard errors of the target SNP conditioned on the SNP set across all the cohorts.
#'@param cohort_count_filt The number of presence for the target SNP across all the cohorts.
#'@param env The study-level environment factors across all the cohorts. Each row refers to one population and each column refers to one environment covariate.
#'For MR-MEGA approach, env=NULL.
#'@param PCs The axes of genetic variation, which can also be called the principal components (PCs). Each row refers to one population.
#'Note: Each \code{env} row and \code{PCs} row should correspond to same population.
#'@param ncores The number of cores which would be used for running in parallel.
#'@returns Output a file containing names of genetic variants, estimated coefficients, standard errors, chisq value of the association, the number of degrees of freedom of the association, p-value of the association, chisq value of the heterogeneity due to different ancestry, ndf of the heterogeneity due to different ancestry, p-value of the heterogeneity due to different ancestry, chisq value of the residual heterogeneity, ndf of the residual heterogeneity,  p-value of the residual heterogeneity.
#'
#'@author Siru Wang
#'@import parallel
#'@import doParallel
#'@import foreach
#'@importFrom stats lm pchisq
#'@export
MR_mega_run<-function(beta_pop,invse2_pop,cohort_count_filt,env,PCs,ncores){

  if(!is.null(env)&&!is.matrix(env)){env=as.matrix(env)}

  if(!is.matrix(PCs)){PCs=as.matrix(PCs)}

  if(!is.null(env)){
    covs=cbind(env,PCs)
  }else{
    covs=PCs
  }

  pcCount=dim(PCs)[2]
  #Run linear regression in parallel and return results of hyporthese testings
  ncores <- min(c(ncores, parallel::detectCores(logical = TRUE)))
  cat("Linear regression model for each genetic variant, and the number of cores are ",ncores,"\n")
  if(ncores>1){
    cl<-makeCluster(ncores,type="FORK")#shared memory
    registerDoParallel(cl)
    lr_out<-foreach(i_snp=1:dim(beta_pop)[1],.combine=rbind)%dopar%{
      check_w=!is.na(invse2_pop[i_snp,-1])

      if((length(unlist(beta_pop[i_snp,-1])[check_w])-(pcCount+ifelse(is.null(env),0,dim(env)[2])))<=2){
        #revised 250805
        if(!is.null(env)){
          exp_lr=paste0("lr_out=data.frame(",paste0(rep(c("beta","se"),dim(covs)[2]+1),rep(0:dim(covs)[2],each=2),"=NA",collapse=","),",RSS=NA,df_RSS=NA,TSS0=NA,TSS=NA,RSS_env=NA,RSS_pc=NA,ndf_association=NA,ndf_chisq_anceheter=NA)")
          eval(parse(text=exp_lr))
        }else{
          exp_lr=paste0("lr_out=data.frame(",paste0(rep(c("beta","se"),dim(covs)[2]+1),rep(0:dim(covs)[2],each=2),"=NA",collapse=","),",RSS=NA,df_RSS=NA,TSS0=NA,TSS=NA,ndf_association=NA,ndf_chisq_anceheter=NA)")
          eval(parse(text=exp_lr))
        }
      }else if(!is.null(env)){
        lr_out<-lr_env_w(unlist(beta_pop[i_snp,-1])[check_w],as.matrix(covs)[check_w,],unlist(invse2_pop[i_snp,-1])[check_w],pcCount)
        lr_out<-data.frame(lr_out,ndf_association=dim(covs)[2]+1,ndf_chisq_anceheter=dim(covs)[2])
      }else{
        lr_out<-lr_w(unlist(beta_pop[i_snp,-1])[check_w],as.matrix(covs)[check_w,],unlist(invse2_pop[i_snp,-1])[check_w])
        lr_out<-data.frame(lr_out,ndf_association=dim(covs)[2]+1,ndf_chisq_anceheter=dim(covs)[2])
      }
      return(lr_out)
    }
    stopCluster(cl)
  }else{
    lr_out<-foreach(i_snp=1:dim(beta_pop)[1],.combine=rbind)%do%{
      #cat("i_snp=",i_snp,"\n")
      check_w=!is.na(invse2_pop[i_snp,-1])
      #cat("the length of Y is ",length(unlist(beta_pop[i_snp,-1])[check_w])," and the number of cov is ",dim(as.matrix(covs)[check_w,])[2],"\n")

      if((length(unlist(beta_pop[i_snp,-1])[check_w])-(pcCount+ifelse(is.null(env),0,dim(env)[2])))<=2){
        if(!is.null(env)){
          exp_lr=paste0("lr_out=data.frame(",paste0(rep(c("beta","se"),dim(covs)[2]+1),rep(0:dim(covs)[2],each=2),"=NA",collapse=","),",RSS=NA,df_RSS=NA,TSS0=NA,TSS=NA,RSS_env=NA,RSS_pc=NA,ndf_association=NA,ndf_chisq_anceheter=NA)")
          eval(parse(text=exp_lr))
        }else{
          exp_lr=paste0("lr_out=data.frame(",paste0(rep(c("beta","se"),dim(covs)[2]+1),rep(0:dim(covs)[2],each=2),"=NA",collapse=","),",RSS=NA,df_RSS=NA,TSS0=NA,TSS=NA,ndf_association=NA,ndf_chisq_anceheter=NA)")
          eval(parse(text=exp_lr))
        }
      }else if(!is.null(env)){
        lr_out<-lr_env_w(unlist(beta_pop[i_snp,-1])[check_w],as.matrix(covs)[check_w,],unlist(invse2_pop[i_snp,-1])[check_w],pcCount)
        lr_out<-data.frame(lr_out,ndf_association=dim(covs)[2]+1,ndf_chisq_anceheter=dim(covs)[2])
      }else{
        lr_out<-lr_w(unlist(beta_pop[i_snp,-1])[check_w],as.matrix(covs)[check_w,],unlist(invse2_pop[i_snp,-1])[check_w])
        lr_out<-data.frame(lr_out,ndf_association=dim(covs)[2]+1,ndf_chisq_anceheter=dim(covs)[2])
      }
      return(lr_out)
    }
  }

  TSS_RSS=abs(lr_out$TSS-lr_out$RSS)
  RSS=abs(lr_out$RSS)
  TSS0_RSS=abs(lr_out$TSS0-lr_out$RSS)
  #Note this is temporary 250805
  #RSS_env_RSS=abs(lr_out$RSS_env-lr_out$RSS)
  #RSS_pc_RSS=abs(lr_out$RSS_pc-lr_out$RSS)
  #revised 250805
  if(!is.null(env)){
    RSS_env_RSS=abs(lr_out$RSS_env-lr_out$RSS)
    RSS_pc_RSS=abs(lr_out$RSS_pc-lr_out$RSS)
    pResEnv=pchisq(RSS_env_RSS,pcCount,lower.tail = FALSE)
    pResPC=pchisq(RSS_pc_RSS,dim(env)[2],lower.tail = FALSE)
  }


  pModelHet=pchisq(TSS_RSS,dim(covs)[2],lower.tail = FALSE)
  #df_pModelHet=dim(covs)[2]
  pResidHet=pchisq(RSS,lr_out$df_RSS,lower.tail = FALSE)
  #df_pResidHet=lr_out$df_RSS
  pModelTest=pchisq(TSS0_RSS,dim(covs)[2]+1,lower.tail = FALSE)
  ##df_pModelTest=dim(covs)[2]+1

  logBF=rep(NA,length(cohort_count_filt$count_pop))
  logBF[(cohort_count_filt$count_pop-2)>dim(covs)[2]]=0.5*(TSS0_RSS[(cohort_count_filt$count_pop-2)>dim(covs)[2]]-(dim(covs)[2]+1)*log(cohort_count_filt$count_pop[(cohort_count_filt$count_pop-2)>dim(covs)[2]]))

  #revised 250805
  if(!is.null(env)){
    logout=data.frame(MARKERNAME=beta_pop$MARKERNAME,
                      chisq_association=TSS0_RSS,
                      ndf_association=lr_out$ndf_association,#dim(covs)[2]+1,
                      pvalue_association=pModelTest,
                      chisq_anc_env_het=TSS_RSS,
                      ndf_chisq_anc_env_het=lr_out$ndf_chisq_anceheter,#dim(covs)[2],
                      pvalue_anc_env_het=pModelHet,
                      chisq_residual=RSS,
                      ndf_residual=lr_out$df_RSS,
                      pvalue_residual=pResidHet,
                      logBF=logBF,
                      #                #Note thisis temporary
                      chisq_env_het=RSS_pc_RSS,
                      ndf_env_het=dim(env)[2],
                      pvalue_env_het=pResPC,
                      chisq_anc_het=RSS_env_RSS,
                      ndf_anc_het=pcCount,
                      pvalue_anc_het=pResEnv)
  }else{

  logout=data.frame(MARKERNAME=beta_pop$MARKERNAME,
                    chisq_association=TSS0_RSS,
                    ndf_association=lr_out$ndf_association,#dim(covs)[2]+1,
                    pvalue_association=pModelTest,
                    chisq_anc_het=TSS_RSS,
                    ndf_chisq_anc_het=lr_out$ndf_chisq_anceheter,#dim(covs)[2],
                    pvalue_anc_het=pModelHet,
                    chisq_residual=RSS,
                    ndf_residual=lr_out$df_RSS,
                    pvalue_residual=pResidHet,
                    logBF=logBF
                    #                #Note thisis temporary
                    #                chisq_env=RSS_pc_RSS,
                    #                ndf_env=dim(envs)[2],
                    #                pvalue_env=pResPC,
                    #                chisq_PC=RSS_env_RSS,
                    #                ndf_PC=pcCount,
                    #                pvalue_PC=pResEnv

  )
}

  return(logout)
}
