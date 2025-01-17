#'@title (env-)MR-MEGA fine-mapping approach (env-MR-MEGAfm)
#'@description Fine-mapping genetic associations across multiple cohorts, allowing for multiple causal variants.
#'@param gwas.list A list of length K_g which contains the pre-processed GWAS files. Each component contains one GWAS file comprising these required columns:
#'"MARKERNAME","CHROMOSOME","POSITION","EA","NEA","EAF","N","BETA","SE".
#'@param ld.list A list of length K_ld which contains the pre-processed LD structures. Each component contains one LD structure.
#'@param which.ld A character vector of length K_g. Each component of the vector corresponds to one LD structure. The length of which.ld should equal to the number of gwas files.
#'@param meta.file A data frame obtained from the output of MR-MEGA method or env-MR-MEGA method.
#'@param PCs The axes of genetic variation, which can also be called the principal components (PCs). Each row refers to one population and each column refers to one axes of genetic variation.
#'The order of matrix row in \code{PCs} should be consistent with the order of gwas names.
#'@param env The matrix which contains the study-level environmental factors. Each row refers to on population and each column refers to one environment covariate.
#'The order of matrix row in \code{env} should be consistent with the order of gwas names. If MR-MEGAfm is utilized, env=NULL.
#'@param out_loc Path to save the output of (env-)MR-MEGAfm. By default, out_loc=NULL.
#'@param ncores The the number of cores which would be used for running in parallel.
#'@param collinear A threshold to filter out the target SNP in high LD with the SNP set. If the squared multiple correlation between the target SNP exceeds the threshold,
#'such as 0.9, the target SNP is ignored.
#'@param pvalue_cutoff A cutoff p-value to identify the potential associated SNP with the smallest p-value exceeding the cutoff p-value. By default, pvalue_cutoff=5e-8.
#'@param cred.thr Credible threshold for the credible set for each selected potential SNP. By default, cred.thr=0.99 refers to 99% credible sets.
#'@param actual.geno An indicator to specify whether the true cohort-level LD structure is applied. If actual.geno is TRUE, the inputted LD structures are derived from the true individual-level genotype data.
#'@return A list of length 2. One component contains the selected potential associated SNP. Another component contains their credible sets with the specified thresholds. By default, the credible threshold is set to 0.99.
#'@export
#'@import parallel
#'@import doParallel
#'@import foreach
#'@importFrom dplyr left_join
#'@author Siru Wang
env_MR_MEGA_fm<-function(gwas.list,ld.list,which.ld,meta.file,PCs,env,out_loc=NULL,ncores=1,collinear=0.9,pvalue_cutoff=5e-8,cred.thr=0.99,actual.geno=TRUE){

  if(!is.list(gwas.list)){
    stop("The gwas files should be inputted as list!")
  }
  n_cohort=length(gwas.list)

  if(!is.list(ld.list)){
    stop("The LD.list should be inputted as list comprising the genotype across populations.")
  }
  cat("(Collinearity cutoff = ",collinear,")\n")

  cat("The number of cores are ",ncores,", which are used for running in parallel\n")
  cat("The threshold for credible set is set to ", cred.thr,".\n")
  #####Initialization of gwas files and the corresponding genotype matrix############

  sel.set=NULL
  remn.set=NULL
  meta.temp<-meta.file
  trace.cond<-0
  #################Stepwise conditioning selection approach############################

  while(any(meta.temp$pvalue_association[!is.na(meta.temp$pvalue_association)]<pvalue_cutoff)&&trace.cond<8){

    idx.na=is.na(meta.temp$pvalue_association)
    minp.assoc=min(meta.temp$pvalue_association[!idx.na])
    idx.minp=(meta.temp$pvalue_association[!idx.na]==minp.assoc)

    if(sum(as.numeric(idx.minp))>1){
      meta.minp=meta.temp[!idx.na,][idx.minp,]
      sel.set=c(sel.set,meta.minp$MARKERNAME[which.max(meta.minp$chisq_association)])
    }else{
      sel.set=c(sel.set,meta.temp$MARKERNAME[which.min(meta.temp$pvalue_association)])
    }

    print("The set of the selected SNPs is")
    print(sel.set)

    remn.set=meta.temp$MARKERNAME[!(meta.temp$MARKERNAME%in%sel.set)]
    trace.cond<-trace.cond+1
    cat("trace.cond is ",trace.cond,"\n")

    if(trace.cond>1){
      selgwas.list=lapply(gwas.list,function(ml,ss){
        loc=match(ss,ml$MARKERNAME)
        return(ml[loc[!is.na(loc)],])
      },sel.set)

      selld.list=vector("list",n_cohort)
      for(i.pop in 1:n_cohort){
        loc=match(sel.set,colnames(ld.list[[which.ld[i.pop]]]))
        selld.list[[i.pop]]=ld.list[[which.ld[i.pop]]][loc[!is.na(loc)],loc[!is.na(loc)]]
      }#revised 240901


      checkld<-lapply(ld.list,function(ll,rs){
        ldcoh.sub=ll[(colnames(ll)%in%rs),(colnames(ll)%in%rs)]
        ldcoh.sub=round(ldcoh.sub,3)
        if(is.matrix(ldcoh.sub)){
          ldcoh.chk=(ldcoh.sub>=collinear)
          ldcoh.chk[lower.tri(ldcoh.chk,diag=TRUE)]=0
        }
        return(any(colSums(ldcoh.chk)>0))
      },remn.set)

      checkld=unlist(checkld)
      if(any(checkld)){break}

      #selgeno.list=lapply(geno.list,function(gl,ss){
      #  loc=match(ss,colnames(gl));
      #  return(gl[,loc[!is.na(loc)]])
      #},sel.set)#revised 240227

      #selgeno.list=vector("list",n_cohort)
      #for(i.pop in 1:n_cohort){
      #  loc=match(sel.set,colnames(geno.list[[uniLD[i.pop]]]))
      #  selgeno.list[[i.pop]]=geno.list[[uniLD[i.pop]]][,loc[!is.na(loc)]]
      #}#revised 240228


      condgwas.list<-Calc_condBF(selgwas.list,selld.list,meta.snp=NULL,sel.set,n_cohort,env,PCs,out_loc=NULL,ncores,cred.thr=NULL,actual.geno,collinear)
      condgwas.list=as.data.frame(condgwas.list)

      if(any(condgwas.list$pvalue_association[!is.na(condgwas.list$pvalue_association)]>pvalue_cutoff)){
        cat("Whether selected SNPs would be excluded in the joint analysis",condgwas.list$pvalue_association[!is.na(condgwas.list$pvalue_association)]>pvalue_cutoff,"\n")
        sel.set<-sel.set[!(sel.set%in%condgwas.list$MARKERNAME[which.max(condgwas.list$pvalue_association)])]

      }

    }

    ncores <- min(c(ncores,parallel::detectCores(logical = TRUE)))
    cl<-makeCluster(ncores,type="FORK")#shared memory
    registerDoParallel(cl)
    newgwas.list<-foreach(i.pop=1:n_cohort,.packages=c("foreach","dplyr"),.export=c("meta_sel_condest","meta_sel_joinest","Get_GenInverse"))%dopar%{

      cat("i.pop is ",i.pop,"\n")

      gwas.file=gwas.list[[i.pop]]
      #gwas.ld=ld.list[[i.pop]] #revised 240227
      gwas.ld=ld.list[[which.ld[i.pop]]] #revised 240901
      #gwas.geno=geno.list[[uniLD[i.pop]]]#revised 240228

      cond.upd=foreach(i.remn=remn.set,.combine=rbind,.export=c("meta_sel_condest","meta_sel_joinest","Get_GenInverse"))%do%{
        #cat("i.remn is ",i.remn,"\n")
        if(all(!gwas.file$MARKERNAME%in%i.remn)|all(!colnames(gwas.ld)%in%i.remn)){
          beta.se=data.frame(MARKERNAME=i.remn,cond.BETA=NA,cond.invSE2=NA)
        }else{

          loc.gwas=match(sel.set,gwas.file$MARKERNAME)
          loc.ld=match(sel.set,colnames(gwas.ld))
          sel.gwas=gwas.file[loc.gwas[!is.na(loc.gwas)],]
          add.gwas=gwas.file[gwas.file$MARKERNAME%in%i.remn,]

          if((length(loc.ld[!is.na(loc.ld)])>0)&(length(loc.gwas[!is.na(loc.gwas)])>0)){

            addsel.ld=matrix(gwas.ld[,colnames(gwas.ld)%in%i.remn],ncol=1,dimnames=list(rownames(gwas.ld),i.remn))
            addsel.ld=matrix(addsel.ld[loc.ld[!is.na(loc.ld)]],ncol=1,dimnames=list(rownames(addsel.ld)[loc.ld[!is.na(loc.ld)]],i.remn))

            if(length(loc.ld[!is.na(loc.ld)])>1){
              sel.ld=gwas.ld[loc.ld[!is.na(loc.ld)],loc.ld[!is.na(loc.ld)]]
            }else{
              sel.ld=matrix(gwas.ld[loc.ld[!is.na(loc.ld)],loc.ld[!is.na(loc.ld)]],dimnames=list(rownames(gwas.ld)[loc.ld[!is.na(loc.ld)]],colnames(gwas.ld)[loc.ld[!is.na(loc.ld)]]))
            }

            if(dim(sel.ld)[2]!=dim(sel.gwas)[1]){
              if(dim(sel.ld)[2]<dim(sel.gwas)[1]){
                loc=match(colnames(sel.ld),sel.gwas$MARKERNAME)
                sel.gwas=sel.gwas[loc[!is.na(loc)],]
              }else{
                loc=match(sel.gwas$MARKERNAME,colnames(sel.ld))
                sel.ld=sel.ld[loc[!is.na(loc)],loc[!is.na(loc)]]
              }
            }
            beta.se=meta_sel_condest(add.gwas,sel.gwas,addsel.ld,sel.ld,collinear,actual.geno)
          }else{
            loc=match(i.remn,colnames(gwas.ld))

            if(!is.na(gwas.ld[loc[!is.na(loc)],loc[!is.na(loc)]])){#241004
              if(length(loc[!is.na(loc)])>1){
                add.ld=gwas.ld[loc[!is.na(loc)],loc[!is.na(loc)]]
              }else{
                add.ld=matrix(gwas.ld[loc[!is.na(loc)],loc[!is.na(loc)]],dimnames=list(rownames(gwas.ld[[i.pop]])[loc[!is.na(loc)]],colnames(gwas.ld[[i.pop]])[loc[!is.na(loc)]]))
              }

              if(dim(add.ld)[2]!=dim(add.gwas)[1]){
                if(dim(add.ld)[2]<dim(add.gwas)[1]){
                  loc=match(colnames(add.ld),add.gwas$MARKERNAME)
                  add.gwas=add.gwas[loc[!is.na(loc)],]
                }else{
                  loc=match(add.gwas$MARKERNAME,colnames(add.ld))
                  add.ld=add.ld[loc[!is.na(loc)],loc[!is.na(loc)]]
                }
              }

              jbeta.se=meta_sel_joinest(sel.gwas=add.gwas,sel.ld=add.ld,actual.geno)
              beta.se=data.frame(MARKERNAME=jbeta.se$MARKERNAME,cond.BETA=jbeta.se$Join.BETA,cond.invSE2=jbeta.se$Join.invSE2)

            }else{
              beta.se=data.frame(MARKERNAME=i.remn,cond.BETA=NA,cond.invSE2=NA)
            }


          }
        }

        return(beta.se)
      }
      gwas.file.cond<-subset(gwas.file,select=c(MARKERNAME,POSITION,CHROMOSOME))
      gwas.file.cond<-left_join(gwas.file.cond,cond.upd,by="MARKERNAME")

      return(gwas.file.cond)
    }
    stopCluster(cl)
    names(newgwas.list)=names(gwas.list)

    exp=paste("beta_pop=data.frame(MARKERNAME=remn.set,",paste("Pop",1:n_cohort,"=NA",sep ="",collapse=","),")",sep="")
    eval(parse(text=exp))
    exp2=paste("invse2_pop=data.frame(MARKERNAME=remn.set,",paste("Pop",1:n_cohort,"=NA",sep ="",collapse=","),")",sep="")
    eval(parse(text=exp2))
    cohort_count_filt=data.frame(MARKERNAME=remn.set,count_pop=0)
    for(i.pop in 1:n_cohort){
      cat("i.pop is ",i.pop,"\n")
      cover=beta_pop$MARKERNAME%in%newgwas.list[[i.pop]]$MARKERNAME
      loc=match(beta_pop$MARKERNAME,newgwas.list[[i.pop]]$MARKERNAME)
      beta_pop[,i.pop+1][cover]=newgwas.list[[i.pop]]$cond.BETA[loc[complete.cases(loc)]]
      invse2_pop[,i.pop+1][cover]=newgwas.list[[i.pop]]$cond.invSE2[loc[complete.cases(loc)]]
      cohort_count_filt$count_pop[cover][which(!is.na(beta_pop[,i.pop+1][cover]))]=cohort_count_filt$count_pop[cover][which(!is.na(beta_pop[,i.pop+1][cover]))]+1

    }
    meta.temp=MR_mega_run(beta_pop,invse2_pop,cohort_count_filt,env,PCs,ncores)

  }
  cat("There are ",length(sel.set)," potential associated SNPs!","\n")
  condBF_tgSNP<-Calc_condBF(gwas.list,ld.list,which.ld,meta.snp=meta.file$MARKERNAME,sel.set,n_cohort,env,PCs,out_loc,ncores,cred.thr,actual.geno,collinear=1)

  return(list(sel.set=sel.set,tgCS.thr=condBF_tgSNP))
}

#'@title Conditional Bayesian factors
#'@description Calculate the Bayesian factor for each target SNP conditioned on the SNP set
#'@param gwas.list  A list of length K_g which contains the pre-processed GWAS files. Each component contains one GWAS file comprising these required columns:
#'"MARKERNAME","CHROMOSOME","POSITION","EA","NEA","EAF","N","BETA","SE".
#'@param ld.list A list of length K_ld which contains the pre-processed LD structures. Each component contains one LD structure.
#'@param which.ld A character vector of length K_g. Each component of the vector corresponds to one LD structure. The length of which.ld should equal to the number of gwas files.
#'@param meta.snp A character vector containing all SNP names present in the output of MR-MEGA method or env-MR-MEGA method.
#'@param sel.set A set contains all potential associated SNPs.
#'@param n_cohort The number of cohorts.
#'@param env The study-level environment factors across all the cohorts. Each row refers to one population and each column refers to one environment covariate.
#'For MR-MEGA approach, env=NULL.
#'@param PCs The axes of genetic variation, which can also be called the principal components (PCs). Each row refers to one population.
#'Note: Each \code{env} row and \code{PCs} row should correspond to same population.
#'@param out_loc Path to save pre-processed GWAS files and LD structures. By default, out_loc=NULL.
#'@param ncores The the number of cores which would be used for running in parallel.
#'@param cred.thr Credible threshold for the credible set for each selected potential SNP. By default, cred.thr=0.99 refers to 99% credible sets.
#'@param actual.geno An indicator to specify whether the true cohort-level LD structure is applied.
#'@param collinear A threshold to filter out the target SNP in high LD with the SNP set. If the squared multiple correlation between the target SNP exceeds the threshold,
#'such as 0.9, the target SNP is ignored.
#'@return Output the results of (env-)MR-MEGAfm conditioned on the subset of the selected potential SNP set.
#'@import doParallel
#'@import parallel
#'@import foreach
#'@importFrom dplyr left_join
#'@importFrom data.table setorder fwrite
#'@author Siru Wang
#'@export
Calc_condBF<-function(gwas.list,ld.list,which.ld,meta.snp,sel.set,n_cohort,env,PCs,out_loc,ncores,cred.thr,actual.geno,collinear){

  if(!is.null(meta.snp)){
    condmeta.list=vector(mode="list",length=length(sel.set))
    names(condmeta.list)=sel.set
  }else{
    condmeta.list=data.frame()
  }

  for(tg.sel in sel.set){
    cat("The target selected genetic variant is", tg.sel,"\n")

    if(length(sel.set)==1){
      cond.sel=sel.set
    }else{
      cond.sel<-sel.set[!(sel.set%in%tg.sel)]
    }

    if(!is.null(meta.snp)){
      remn.set=meta.snp[!(meta.snp%in%cond.sel)]
    }else{
      remn.set=tg.sel
    }

    #############################Generate new GWAS files condition on the set of selected SNPs#####################################
    ncores <- min(c(ncores,parallel::detectCores(logical = TRUE)))
    cl<-makeCluster(ncores,type="FORK")#shared memory
    registerDoParallel(cl)
    newgwas.list<-foreach(i.pop=1:n_cohort,.packages="foreach",.export=c("meta_sel_condest","meta_sel_joinest","Get_GenInverse"))%dopar%{

      cat("i.pop is ",i.pop,"\n")
      gwas.file=gwas.list[[i.pop]]
      #gwas.ld=ld.list[[i.pop]]#revised 240327
      gwas.ld=ld.list[[which.ld[i.pop]]]#revised 240901


      cond.upd=foreach(i.remn=remn.set,.combine=rbind,.export=c("meta_sel_condest","meta_sel_joinest","Get_GenInverse"))%do%{
        #    cat("i.remn is ",i.remn,"\n")
        if(all(!gwas.file$MARKERNAME%in%i.remn)|all(!colnames(gwas.ld)%in%i.remn)){
          beta.se=data.frame(MARKERNAME=i.remn,cond.BETA=NA,cond.invSE2=NA)
        }else{
          #sel.meta=meta.file[meta.file$MARKERNAME%in%cond.sel,]
          loc.gwas=match(cond.sel,gwas.file$MARKERNAME)
          loc.ld=match(cond.sel,colnames(gwas.ld))
          sel.gwas=gwas.file[loc.gwas[!is.na(loc.gwas)],]
          add.gwas=gwas.file[gwas.file$MARKERNAME%in%i.remn,]
          if((length(loc.ld[!is.na(loc.ld)])>0)&(length(loc.gwas[!is.na(loc.gwas)])>0)){

            addsel.ld=matrix(gwas.ld[,colnames(gwas.ld)%in%i.remn],ncol=1,dimnames=list(rownames(gwas.ld),i.remn))
            addsel.ld=matrix(addsel.ld[loc.ld[!is.na(loc.ld)]],ncol=1,dimnames=list(rownames(addsel.ld)[loc.ld[!is.na(loc.ld)]],i.remn))

            if(length(loc.ld[!is.na(loc.ld)])>1){
              sel.ld=gwas.ld[loc.ld[!is.na(loc.ld)],loc.ld[!is.na(loc.ld)]]
            }else{
              sel.ld=matrix(gwas.ld[loc.ld[!is.na(loc.ld)],loc.ld[!is.na(loc.ld)]],dimnames=list(rownames(gwas.ld)[loc.ld[!is.na(loc.ld)]],colnames(gwas.ld)[loc.ld[!is.na(loc.ld)]]))
            }

            if(dim(sel.ld)[2]!=dim(sel.gwas)[1]){
              if(dim(sel.ld)[2]<dim(sel.gwas)[1]){
                loc=match(colnames(sel.ld),sel.gwas$MARKERNAME)
                sel.gwas=sel.gwas[loc[!is.na(loc)],]
              }else{
                loc=match(sel.gwas$MARKERNAME,colnames(sel.ld))
                sel.ld=sel.ld[loc[!is.na(loc)],loc[!is.na(loc)]]
              }
            }
            beta.se=meta_sel_condest(add.gwas,sel.gwas,addsel.ld,sel.ld,collinear,actual.geno)
          }else{
            loc=match(i.remn,colnames(gwas.ld))

            if(!is.na(gwas.ld[loc[!is.na(loc)],loc[!is.na(loc)]])){#241004
              if(length(loc[!is.na(loc)])>1){
                add.ld=gwas.ld[loc[!is.na(loc)],loc[!is.na(loc)]]
              }else{
                add.ld=matrix(gwas.ld[loc[!is.na(loc)],loc[!is.na(loc)]],dimnames=list(rownames(gwas.ld[[i.pop]])[loc[!is.na(loc)]],colnames(gwas.ld[[i.pop]])[loc[!is.na(loc)]]))
              }

              if(dim(add.ld)[2]!=dim(add.gwas)[1]){
                if(dim(add.ld)[2]<dim(add.gwas)[1]){
                  loc=match(colnames(add.ld),add.gwas$MARKERNAME)
                  add.gwas=add.gwas[loc[!is.na(loc)],]
                }else{
                  loc=match(add.gwas$MARKERNAME,colnames(add.ld))
                  add.ld=add.ld[loc[!is.na(loc)],loc[!is.na(loc)]]
                }
              }

              jbeta.se=meta_sel_joinest(sel.gwas=add.gwas,sel.ld=add.ld,actual.geno)
              beta.se=data.frame(MARKERNAME=jbeta.se$MARKERNAME,cond.BETA=jbeta.se$Join.BETA,cond.invSE2=jbeta.se$Join.invSE2)
            }else{
              beta.se=data.frame(MARKERNAME=i.remn,cond.BETA=NA,cond.invSE2=NA)
            }

          }
        }
        return(beta.se)
      }

      gwas.file.cond<-subset(gwas.file,select=c(MARKERNAME,POSITION,CHROMOSOME))
      gwas.file.cond<-left_join(gwas.file.cond,cond.upd,by="MARKERNAME")

      return(gwas.file.cond)
    }
    stopCluster(cl)
    ################################Perform environment-adjusted MR-MEGA############################
    names(newgwas.list)=names(gwas.list)
    exp=paste("beta_pop=data.frame(MARKERNAME=remn.set,",paste("Pop",1:n_cohort,"=NA",sep ="",collapse=","),")",sep="")
    eval(parse(text=exp))
    exp2=paste("invse2_pop=data.frame(MARKERNAME=remn.set,",paste("Pop",1:n_cohort,"=NA",sep ="",collapse=","),")",sep="")
    eval(parse(text=exp2))
    cohort_count_filt=data.frame(MARKERNAME=remn.set,count_pop=0)
    for(i.pop in 1:n_cohort){
      cat("i.pop is ",i.pop,"\n")
      cover=beta_pop$MARKERNAME%in%newgwas.list[[i.pop]]$MARKERNAME
      loc=match(beta_pop$MARKERNAME,newgwas.list[[i.pop]]$MARKERNAME)
      beta_pop[,i.pop+1][cover]=newgwas.list[[i.pop]]$cond.BETA[loc[complete.cases(loc)]]

      invse2_pop[,i.pop+1][cover]=newgwas.list[[i.pop]]$cond.invSE2[loc[complete.cases(loc)]]
      cohort_count_filt$count_pop[cover][which(!is.na(beta_pop[,i.pop+1][cover]))]=cohort_count_filt$count_pop[cover][which(!is.na(beta_pop[,i.pop+1][cover]))]+1
      #fwrite(newmeta.list[[i.pop]],file=paste0(out_loc,"/",names(newmeta.list)[i.pop],"_cond",trace.cond,".txt"))
    }

    meta.temp=MR_mega_run(beta_pop,invse2_pop,cohort_count_filt,env,PCs,ncores)
    #PPcv=vector(3)

    logsum <- function(x) {
      my.max <- max(x)
      my.res <- my.max + log(sum(exp(x - my.max )))
      return(my.res)
    }

    #meta.snp=NULL
    if(!is.null(meta.snp)){
      PPs=data.frame(MARKERNAME=remn.set,PPs=NA)
      na.BF=(!is.na(meta.temp$logBF))
      PPs$PPs[na.BF]=exp(meta.temp$logBF[na.BF]-logsum(meta.temp$logBF[na.BF]))#exp(meta.temp$logBF[na.BF])/sum(exp(meta.temp$logBF[na.BF]))
      #PPcv
      meta.temp=left_join(meta.temp,PPs,by="MARKERNAME")
      meta.desc=setorder(meta.temp,-PPs,na.last=TRUE)

      naidx=is.na(meta.desc$PPs)
      meta.desc=meta.desc[!naidx,]

      if(!is.null(cred.thr)){
        csm=cumsum(meta.desc$PPs)
        condmeta.list[[tg.sel]]=meta.desc[1:(which(csm>cred.thr)[1]),]
      }else{
        condmeta.list[[tg.sel]]=meta.desc
      }
    }else{
      meta.desc=meta.temp
      condmeta.list=rbind(condmeta.list,meta.desc)
    }

    if(!is.null(out_loc)){
      #condmeta.list=rbind(condmeta.list,meta.desc)
      fwrite(meta.desc,file=paste0(out_loc,"/envmeta_",tg.sel,"_cvPPs.txt"),sep="\t")
    }

  }
  return(condmeta.list)
}
