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
#'@param actual.geno An indicator to specify whether the true cohort-level LD structure is applied. If actual.geno is TRUE, the inputted LD structures are derived from the true individual-level genotype data. By default, actual.geno=FALSE is recommended in GCTA-COJO software.
#'@return A list of length 2. One component contains the selected potential associated SNP. Another component contains their credible sets with the specified thresholds. By default, the credible threshold is set to 0.99.
#'@export
#'@import parallel
#'@import doParallel
#'@import foreach
#'@importFrom dplyr left_join
#'@importFrom stats setNames
#'@author Siru Wang
env_MR_MEGA_fm<-function(gwas.list,ld.list,which.ld,meta.file,PCs,env,out_loc=NULL,ncores=1,collinear=0.9,pvalue_cutoff=5e-8,cred.thr=0.99,actual.geno=FALSE){

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

  #revised 250415
  ncores <- min(c(ncores,parallel::detectCores(logical = TRUE)))
  cl<-makeCluster(ncores,type="FORK")#shared memory
  registerDoParallel(cl)
  gwas.list<-foreach(i.pop=1:n_cohort, .final = function(x) {setNames(x,names(gwas.list))})%do%{
    cat("i.pop is ",i.pop,"\n")
    gwas.file=gwas.list[[i.pop]]

    #calculate Vp revised 250415
    h.buf=2*gwas.file$EAF*(1-gwas.file$EAF)
    Vp.init=gwas.file$N*h.buf*gwas.file$SE*gwas.file$SE+h.buf*gwas.file$BETA*gwas.file$BETA*gwas.file$N/(gwas.file$N-1)
    Vp.med=median(Vp.init)
    gwas.file=data.frame(gwas.file,Vp.med=Vp.med)

    if(!isTRUE(actual.geno)){
      cat("Whether the inputted LD is cohort-matched LD",isTRUE(actual.geno),"\n")
      Neff=(Vp.med-h.buf*gwas.file$BETA*gwas.file$BETA)/(h.buf*gwas.file$SE*gwas.file$SE)+1
      gwas.file=data.frame(gwas.file,Neff=Neff)
    }

    return(gwas.file)
  }
  stopCluster(cl)

  sel.set=NULL
  remn.set=NULL
  meta.temp<-meta.file
  trace.cond<-0
  #################Stepwise conditioning selection approach############################

  while(any(meta.temp$pvalue_association[!is.na(meta.temp$pvalue_association)]<pvalue_cutoff)&&trace.cond<5){

    idx.na=is.na(meta.temp$pvalue_association)
    minp.assoc=min(meta.temp$pvalue_association[!idx.na])
    idx.minp=(meta.temp$pvalue_association[!idx.na]==minp.assoc)

    if(sum(as.numeric(idx.minp))>1){
      meta.minp=meta.temp[!idx.na,][idx.minp,]
      new.sel=meta.minp$MARKERNAME[which.max(meta.minp$chisq_association)]
      #revised 250805
      #new.sel.p=data.frame(sel=new.sel,cond.p=meta.minp$pvalue_association[which.max(meta.minp$chisq_association)])
    }else{
      new.sel=meta.temp$MARKERNAME[which.min(meta.temp$pvalue_association)]
      #250805
      #new.sel.p=data.frame(sel=new.sel,cond.p=meta.temp$pvalue_association[which.min(meta.temp$pvalue_association)])
    }

    sel.set=c(sel.set,new.sel)
    #revised 250805
    #condp.sel.set=rbind(condp.sel.set,new.sel.p)

    remn.set=meta.temp$MARKERNAME[!(meta.temp$MARKERNAME%in%sel.set)]
    trace.cond<-trace.cond+1
    cat("At the ",trace.cond,"iteration, the set of the selected SNPs include \n")
    print(sel.set)

    if(trace.cond>1){
      selgwas.list=lapply(gwas.list,function(ml,ss){
        loc=match(ss,ml$MARKERNAME)
        return(ml[loc[!is.na(loc)],])
      },sel.set)

      selld.list=vector("list",length(ld.list))
      names(selld.list)=names(ld.list)
      for(i.pop in 1:length(ld.list)){
        loc=match(sel.set,colnames(ld.list[[i.pop]]))
        selld.list[[i.pop]]=ld.list[[i.pop]][loc[!is.na(loc)],loc[!is.na(loc)]]
      }#revised 250404

      #selld.list=vector("list",n_cohort)
      #for(i.pop in 1:n_cohort){
      #  loc=match(sel.set,colnames(ld.list[[which.ld[i.pop]]]))
      #  selld.list[[i.pop]]=ld.list[[which.ld[i.pop]]][loc[!is.na(loc)],loc[!is.na(loc)]]
      #}#revised 240901


      condgwas.list<-Calc_condBF(selgwas.list,selld.list,which.ld,meta.snp=NULL,sel.set,n_cohort,env,PCs,out_loc=NULL,ncores,cred.thr=NULL,actual.geno,collinear)
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

      #cat("i.pop is ",i.pop,"\n")

      gwas.file=gwas.list[[i.pop]]
      #gwas.ld=ld.list[[i.pop]] #revised 240227
      gwas.ld=ld.list[[which.ld[i.pop]]] #revised 240901
      #gwas.geno=geno.list[[uniLD[i.pop]]]#revised 240228

      cond.upd=foreach(i.remn=remn.set,.combine=rbind,.export=c("meta_sel_condest","Get_GenInverse"))%do%{

        if(all(!gwas.file$MARKERNAME%in%i.remn)|all(!colnames(gwas.ld)%in%i.remn)|any(!(sel.set%in%gwas.file$MARKERNAME))|any(!(sel.set%in%colnames(gwas.ld)))){# revised 250310
          #revised 250317
          if(any(sel.set%in%rownames(gwas.ld))&(i.remn%in%colnames(gwas.ld))){
            addsel.maxld=max(abs(gwas.ld[rownames(gwas.ld)%in%sel.set,colnames(gwas.ld)%in%i.remn]),na.rm=TRUE)
          }else{
            addsel.maxldi=NA
          }
          beta.se=data.frame(MARKERNAME=i.remn,cond.BETA=NA,cond.SE=NA,cond.invSE2=NA,mul.cor=NA,ld.coh=addsel.maxld)
        }else{

          #revised 250317
          if(any(sel.set%in%rownames(gwas.ld))&(i.remn%in%colnames(gwas.ld))){
            addsel.maxld=max(abs(gwas.ld[rownames(gwas.ld)%in%sel.set,colnames(gwas.ld)%in%i.remn]),na.rm=TRUE)
          }else{
            addsel.maxld=NA
          }

          loc.gwas=match(sel.set,gwas.file$MARKERNAME)
          loc.ld=match(sel.set,colnames(gwas.ld))
          sel.gwas=gwas.file[loc.gwas[!is.na(loc.gwas)],]
          add.gwas=gwas.file[gwas.file$MARKERNAME%in%i.remn,]

          addsel.ld=matrix(gwas.ld[,colnames(gwas.ld)%in%i.remn],ncol=1,dimnames=list(rownames(gwas.ld),i.remn))
          addsel.ld=matrix(addsel.ld[loc.ld[!is.na(loc.ld)]],ncol=1,dimnames=list(rownames(addsel.ld)[loc.ld[!is.na(loc.ld)]],i.remn))

          if(length(loc.ld[!is.na(loc.ld)])>1){
            sel.ld=gwas.ld[loc.ld[!is.na(loc.ld)],loc.ld[!is.na(loc.ld)]]
          }else{
            sel.ld=matrix(gwas.ld[loc.ld[!is.na(loc.ld)],loc.ld[!is.na(loc.ld)]],dimnames=list(rownames(gwas.ld)[loc.ld[!is.na(loc.ld)]],colnames(gwas.ld)[loc.ld[!is.na(loc.ld)]]))
          }

          if(any(is.na(sel.ld))|any(is.na(addsel.ld))|any(is.na(sel.gwas$BETA))|any(is.na(sel.gwas$SE))){
            beta.se=data.frame(MARKERNAME=i.remn,cond.BETA=NA,cond.SE=NA,cond.invSE2=NA,mul.cor=NA,ld.coh=addsel.maxld)
          }else{
            #revised 250415
            Vp.med=add.gwas$Vp.med
            beta.se=meta_sel_condest(add.gwas,sel.gwas,addsel.ld,sel.ld,Vp.med,collinear,actual.geno)
            beta.se=data.frame(beta.se,ld.coh=addsel.maxld)
          }

        }

        return(beta.se)
      }
      gwas.file.cond<-subset(gwas.file,select=c(MARKERNAME,POSITION,CHROMOSOME,EA,NEA))
      gwas.file.cond<-left_join(gwas.file.cond,cond.upd,by="MARKERNAME")

      return(gwas.file.cond)
    }
    stopCluster(cl)
    names(newgwas.list)=names(gwas.list)

    #Assigned to the same reference panel
    print("All GWAS are assigned to the same reference allele")
    print("At the same, the correponding LD structures should be adjusted accordingly.")
    #print("Set the ancestry 1 as the reference allele......")

    for(i_pop in 2:length(newgwas.list)){
      #print(paste("Check ancestry",i_pop,sep=""))
      for(j_pop in (i_pop-1):1){
        cover=newgwas.list[[i_pop-j_pop]]$MARKERNAME%in%newgwas.list[[i_pop]]$MARKERNAME
        loc=match(newgwas.list[[i_pop-j_pop]]$MARKERNAME,newgwas.list[[i_pop]]$MARKERNAME)
        same_ea=(newgwas.list[[i_pop]]$EA[loc[complete.cases(loc)]]!=newgwas.list[[i_pop-j_pop]]$EA[cover])

        if(any(same_ea)){
          newgwas.list[[i_pop]]$EA[loc[complete.cases(loc)]][same_ea]=newgwas.list[[i_pop-j_pop]]$EA[cover][same_ea]
          newgwas.list[[i_pop]]$NEA[loc[complete.cases(loc)]][same_ea]=newgwas.list[[i_pop-j_pop]]$NEA[cover][same_ea]
          #newgwas.list[[i_pop]]$EAF[loc[complete.cases(loc)]][same_ea]=1-newgwas.list[[i_pop]]$EAF[loc[complete.cases(loc)]][same_ea]
          newgwas.list[[i_pop]]$cond.BETA[loc[complete.cases(loc)]][same_ea]=-newgwas.list[[i_pop]]$cond.BETA[loc[complete.cases(loc)]][same_ea]
        }
      }
    }


    exp=paste("beta_pop=data.frame(MARKERNAME=remn.set,",paste("Pop",1:n_cohort,"=NA",sep ="",collapse=","),")",sep="")
    eval(parse(text=exp))
    exp2=paste("invse2_pop=data.frame(MARKERNAME=remn.set,",paste("Pop",1:n_cohort,"=NA",sep ="",collapse=","),")",sep="")
    eval(parse(text=exp2))

    #250317
    exp3=paste("mcor_pop=data.frame(MARKERNAME=remn.set,",paste("Pop",1:n_cohort,"=NA",sep ="",collapse=","),")",sep="")
    eval(parse(text=exp3))
    exp4=paste("ld_pop=data.frame(MARKERNAME=remn.set,",paste("Pop",1:n_cohort,"=NA",sep ="",collapse=","),")",sep="")
    eval(parse(text=exp4))

    cohort_count_filt=data.frame(MARKERNAME=remn.set,count_pop=0)
    for(i.pop in 1:n_cohort){
      #cat("i.pop is ",i.pop,"\n")
      cover=beta_pop$MARKERNAME%in%newgwas.list[[i.pop]]$MARKERNAME
      loc=match(beta_pop$MARKERNAME,newgwas.list[[i.pop]]$MARKERNAME)
      beta_pop[,i.pop+1][cover]=newgwas.list[[i.pop]]$cond.BETA[loc[complete.cases(loc)]]
      #revised 250317
      mcor_pop[,i.pop+1][cover]=newgwas.list[[i.pop]]$mul.cor[loc[complete.cases(loc)]]
      ld_pop[,i.pop+1][cover]=newgwas.list[[i.pop]]$ld.coh[loc[complete.cases(loc)]]

      invse2_pop[,i.pop+1][cover]=newgwas.list[[i.pop]]$cond.invSE2[loc[complete.cases(loc)]]
      cohort_count_filt$count_pop[cover][which(!is.na(beta_pop[,i.pop+1][cover]))]=cohort_count_filt$count_pop[cover][which(!is.na(beta_pop[,i.pop+1][cover]))]+1

    }

    #revised 250317
    mcor.thr=apply(mcor_pop[,-1],1,function(mp){
      if(all(is.na(mp))){
        thr=TRUE
      }else{
        thr=(max(mp,na.rm=TRUE)<collinear)
      }
      return(thr)})
    tgcor.thr=apply(ld_pop[,-1],1,function(lp){
      if(all(is.na(lp))){
        thr=TRUE
      }else{
        thr=(max(lp,na.rm=TRUE)<collinear)
      }
      return(thr)
    })

    meta.temp=MR_mega_run(beta_pop[mcor.thr&tgcor.thr,],invse2_pop[mcor.thr&tgcor.thr,],cohort_count_filt[mcor.thr&tgcor.thr,],env,PCs,ncores)

  }
  cat("There are ",length(sel.set)," potential associated SNPs!","\n")
  condBF_tgSNP<-Calc_condBF(gwas.list,ld.list,which.ld,meta.snp=meta.file$MARKERNAME,sel.set,n_cohort,env,PCs,out_loc,ncores,cred.thr,actual.geno,collinear=2)

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

  logsum <- function(x) {
    my.max <- max(x)
    my.res <- my.max + log(sum(exp(x - my.max )))
    return(my.res)
  }

  if(!is.null(meta.snp)){
    condmeta.list=vector(mode="list",length=length(sel.set))
    names(condmeta.list)=sel.set
  }else{
    condmeta.list=data.frame()
  }

  for(tg.sel in sel.set){
    cat("The target selected genetic variant is", tg.sel,"\n")

    if(!is.null(meta.snp)){
      cond.sel<-sel.set[!(sel.set%in%tg.sel)]

      if(length(cond.sel)!=0){
        remn.set=meta.snp[!(meta.snp%in%cond.sel)]
      }else{
        na.BF=(!is.na(meta.file$logBF))
        meta.file$PPs[na.BF]=exp(meta.file$logBF[na.BF]-logsum(meta.file$logBF[na.BF]))
        meta.desc=setorder(meta.file,-PPs,na.last=TRUE)

        naidx=is.na(meta.desc$PPs)
        meta.desc=meta.desc[!naidx,]

        #revised 250319
        exp=paste("ld_pop=data.frame(MARKERNAME=meta.desc$MARKERNAME,",paste("Pop",1:n_cohort,"=NA",sep ="",collapse=","),")",sep="")
        eval(parse(text=exp))

        for(pop.ld in 1:n_cohort){
          tg.ld=ld.list[[which.ld[pop.ld]]]
          if((tg.sel%in%colnames(tg.ld))&any(meta.desc$MARKERNAME%in%rownames(tg.ld))){
            cover=meta.desc$MARKERNAME%in%rownames(tg.ld)
            loc.ld=match(meta.desc$MARKERNAME,rownames(tg.ld))
            ld_pop[,pop.ld+1][cover]=tg.ld[rownames(tg.ld)[loc.ld[complete.cases(loc.ld)]],colnames(tg.ld)%in%tg.sel]
          }
        }

        ld_max=data.frame(MARKERNAME=ld_pop$MARKERNAME,maxld=apply(ld_pop[,-1],1,function(lp){
          if(!all(is.na(lp))){
            ld.m=lp[which.max(abs(lp))]
          }else{
            ld.m=NA
          }
          return(ld.m)}))

        meta.desc=left_join(meta.desc,ld_max,by="MARKERNAME")
        if(!is.null(cred.thr)){
          csm=cumsum(meta.desc$PPs)
          #condmeta.list[[tg.sel]]=meta.desc[1:(which(csm>cred.thr)[1]),]
          #revised 250314
          meta.desc=data.frame(meta.desc,cs.inc=rep(c(1,0),c(which(csm>cred.thr)[1],dim(meta.desc)[1]-which(csm>cred.thr)[1])))
        }else{
          #revised 250314
          meta.desc=data.frame(meta.desc,cs.inc=rep(c(1,0),c(dim(meta.desc)[1],0)))
          condmeta.list[[tg.sel]]=meta.desc
        }


        if(!is.null(out_loc)){
          #condmeta.list=rbind(condmeta.list,meta.desc)
          fwrite(meta.desc,file=paste0(out_loc,"/envmeta_",tg.sel,"_Neff_cvPPs.txt"),sep="\t")
        }
        next
      }
    }else{
      cond.sel<-sel.set[!(sel.set%in%tg.sel)]# revised 250220
      remn.set=tg.sel
    }

    #############################Generate new GWAS files condition on the set of selected SNPs#####################################
    ncores <- min(c(ncores,parallel::detectCores(logical = TRUE)))
    cl<-makeCluster(ncores,type="FORK")#shared memory
    registerDoParallel(cl)
    newgwas.list<-foreach(i.pop=1:n_cohort,.packages="foreach",.export=c("meta_sel_condest","Get_GenInverse"))%dopar%{

      #cat("i.pop is ",i.pop,"\n")
      gwas.file=gwas.list[[i.pop]]
      #gwas.ld=ld.list[[i.pop]]#revised 240327
      gwas.ld=ld.list[[which.ld[i.pop]]]#revised 240901


      cond.upd=foreach(i.remn=remn.set,.combine=rbind,.export=c("meta_sel_condest","Get_GenInverse"))%do%{
        #    cat("i.remn is ",i.remn,"\n")
        if(all(!gwas.file$MARKERNAME%in%i.remn)|all(!colnames(gwas.ld)%in%i.remn)|any(!(cond.sel%in%gwas.file$MARKERNAME))|any(!(cond.sel%in%colnames(gwas.ld)))){
          #revised 250314
          if((tg.sel%in%rownames(gwas.ld))&(i.remn%in%colnames(gwas.ld))){
            addtg.ld=gwas.ld[rownames(gwas.ld)%in%tg.sel,colnames(gwas.ld)%in%i.remn]
          }else{
            addtg.ld=NA
          }
          beta.se=data.frame(MARKERNAME=i.remn,cond.BETA=NA,cond.SE=NA,cond.invSE2=NA,mul.cor=NA,ld.coh=addtg.ld)
        }else{

          #revised 250314
          if((tg.sel%in%rownames(gwas.ld))&(i.remn%in%colnames(gwas.ld))){
            addtg.ld=gwas.ld[rownames(gwas.ld)%in%tg.sel,colnames(gwas.ld)%in%i.remn]
          }else{
            addtg.ld=NA
          }

          #sel.meta=meta.file[meta.file$MARKERNAME%in%cond.sel,]
          loc.gwas=match(cond.sel,gwas.file$MARKERNAME)
          loc.ld=match(cond.sel,colnames(gwas.ld))
          sel.gwas=gwas.file[loc.gwas[!is.na(loc.gwas)],]
          add.gwas=gwas.file[gwas.file$MARKERNAME%in%i.remn,]

          addsel.ld=matrix(gwas.ld[,colnames(gwas.ld)%in%i.remn],ncol=1,dimnames=list(rownames(gwas.ld),i.remn))
          addsel.ld=matrix(addsel.ld[loc.ld[!is.na(loc.ld)]],ncol=1,dimnames=list(rownames(addsel.ld)[loc.ld[!is.na(loc.ld)]],i.remn))

          if(length(loc.ld[!is.na(loc.ld)])>1){
            sel.ld=gwas.ld[loc.ld[!is.na(loc.ld)],loc.ld[!is.na(loc.ld)]]
          }else{
            sel.ld=matrix(gwas.ld[loc.ld[!is.na(loc.ld)],loc.ld[!is.na(loc.ld)]],dimnames=list(rownames(gwas.ld)[loc.ld[!is.na(loc.ld)]],colnames(gwas.ld)[loc.ld[!is.na(loc.ld)]]))
          }

          if(any(is.na(sel.ld))|any(is.na(addsel.ld))|any(is.na(sel.gwas$BETA))|any(is.na(sel.gwas$SE))){
            beta.se=data.frame(MARKERNAME=i.remn,cond.BETA=NA,cond.SE=NA,cond.invSE2=NA,mul.cor=NA,ld.coh=addtg.ld)
          }else{
            #revised 250415
            Vp.med=add.gwas$Vp.med
            beta.se=meta_sel_condest(add.gwas,sel.gwas,addsel.ld,sel.ld,Vp.med,collinear,actual.geno)
            #revised 250317
            beta.se=data.frame(beta.se,ld.coh=addtg.ld)
          }
        }
        return(beta.se)
      }

      gwas.file.cond<-subset(gwas.file,select=c(MARKERNAME,POSITION,CHROMOSOME,EA,NEA))
      gwas.file.cond<-left_join(gwas.file.cond,cond.upd,by="MARKERNAME")

      return(gwas.file.cond)
    }
    stopCluster(cl)
    ################################Perform environment-adjusted MR-MEGA############################
    names(newgwas.list)=names(gwas.list)

    #Assigned to the same reference panel
    print("All GWAS are assigned to the same reference allele")
    print("At the same, the correponding LD structures should be adjusted accordingly.")
    #print("Set the ancestry 1 as the reference allele......")

    for(i_pop in 2:length(newgwas.list)){
      #print(paste("Check ancestry",i_pop,sep=""))
      for(j_pop in (i_pop-1):1){
        cover=newgwas.list[[i_pop-j_pop]]$MARKERNAME%in%newgwas.list[[i_pop]]$MARKERNAME
        loc=match(newgwas.list[[i_pop-j_pop]]$MARKERNAME,newgwas.list[[i_pop]]$MARKERNAME)
        same_ea=(newgwas.list[[i_pop]]$EA[loc[complete.cases(loc)]]!=newgwas.list[[i_pop-j_pop]]$EA[cover])

        if(any(same_ea)){
          newgwas.list[[i_pop]]$EA[loc[complete.cases(loc)]][same_ea]=newgwas.list[[i_pop-j_pop]]$EA[cover][same_ea]
          newgwas.list[[i_pop]]$NEA[loc[complete.cases(loc)]][same_ea]=newgwas.list[[i_pop-j_pop]]$NEA[cover][same_ea]
          #newgwas.list[[i_pop]]$EAF[loc[complete.cases(loc)]][same_ea]=1-newgwas.list[[i_pop]]$EAF[loc[complete.cases(loc)]][same_ea]
          newgwas.list[[i_pop]]$cond.BETA[loc[complete.cases(loc)]][same_ea]=-newgwas.list[[i_pop]]$cond.BETA[loc[complete.cases(loc)]][same_ea]
        }
      }
    }

    exp=paste("beta_pop=data.frame(MARKERNAME=remn.set,",paste("Pop",1:n_cohort,"=NA",sep ="",collapse=","),")",sep="")
    eval(parse(text=exp))
    exp2=paste("invse2_pop=data.frame(MARKERNAME=remn.set,",paste("Pop",1:n_cohort,"=NA",sep ="",collapse=","),")",sep="")
    eval(parse(text=exp2))

    #revised 250314
    exp3=paste("ld_pop=data.frame(MARKERNAME=remn.set,",paste("Pop",1:n_cohort,"=NA",sep ="",collapse=","),")",sep="")
    eval(parse(text=exp3))
    #revised 250317
    exp4=paste("mcor_pop=data.frame(MARKERNAME=remn.set,",paste("Pop",1:n_cohort,"=NA",sep ="",collapse=","),")",sep="")
    eval(parse(text=exp4))

    cohort_count_filt=data.frame(MARKERNAME=remn.set,count_pop=0)
    for(i.pop in 1:n_cohort){
      #cat("i.pop is ",i.pop,"\n")
      cover=beta_pop$MARKERNAME%in%newgwas.list[[i.pop]]$MARKERNAME
      loc=match(beta_pop$MARKERNAME,newgwas.list[[i.pop]]$MARKERNAME)
      beta_pop[,i.pop+1][cover]=newgwas.list[[i.pop]]$cond.BETA[loc[complete.cases(loc)]]

      #revised 250314
      ld_pop[,i.pop+1][cover]=newgwas.list[[i.pop]]$ld.coh[loc[complete.cases(loc)]]
      #revised 250317
      mcor_pop[,i.pop+1][cover]=newgwas.list[[i.pop]]$mul.cor[loc[complete.cases(loc)]]

      invse2_pop[,i.pop+1][cover]=newgwas.list[[i.pop]]$cond.invSE2[loc[complete.cases(loc)]]
      cohort_count_filt$count_pop[cover][which(!is.na(beta_pop[,i.pop+1][cover]))]=cohort_count_filt$count_pop[cover][which(!is.na(beta_pop[,i.pop+1][cover]))]+1
      #fwrite(newmeta.list[[i.pop]],file=paste0(out_loc,"/",names(newmeta.list)[i.pop],"_cond",trace.cond,".txt"))
    }

    meta.temp=MR_mega_run(beta_pop,invse2_pop,cohort_count_filt,env,PCs,ncores)
    if(!is.null(meta.snp)){
      #revised 250319
      meta.temp=MR_mega_run(beta_pop,invse2_pop,cohort_count_filt,env,PCs,ncores)

      ld_max=data.frame(MARKERNAME=ld_pop$MARKERNAME,maxld=apply(ld_pop[,-1],1,function(lp){
        if(!all(is.na(lp))){
          ld.m=lp[which.max(abs(lp))]
        }else{
          ld.m=NA
        }
        return(ld.m)}))

      PPs=data.frame(MARKERNAME=remn.set,PPs=NA)
      na.BF=(!is.na(meta.temp$logBF))
      PPs$PPs[na.BF]=exp(meta.temp$logBF[na.BF]-logsum(meta.temp$logBF[na.BF]))#exp(meta.temp$logBF[na.BF])/sum(exp(meta.temp$logBF[na.BF]))
      #PPcv
      meta.temp=left_join(meta.temp,PPs,by="MARKERNAME")
      meta.temp=left_join(meta.temp,ld_max,by="MARKERNAME")

      meta.desc=setorder(meta.temp,-PPs,na.last=TRUE)
      naidx=is.na(meta.desc$PPs)
      meta.desc=meta.desc[!naidx,]

      if(!is.null(cred.thr)){
        csm=cumsum(meta.desc$PPs)
        #condmeta.list[[tg.sel]]=meta.desc[1:(which(csm>cred.thr)[1]),]
        #revised 250314
        meta.desc=data.frame(meta.desc,cs.inc=rep(c(1,0),c(which(csm>cred.thr)[1],dim(meta.desc)[1]-which(csm>cred.thr)[1])))
      }else{
        #condmeta.list[[tg.sel]]=meta.desc
        #revised 250314
        meta.desc=data.frame(meta.desc,cs.inc=rep(c(1,0),c(dim(meta.desc)[1],0)))
      }
      condmeta.list[[tg.sel]]=meta.desc[meta.desc$cs.inc==1,]
    }else{
      #revised 250319
      mcor.thr=apply(mcor_pop[,-1],1,function(mp){
        if(all(is.na(mp))){
          thr=TRUE
        }else{
          thr=(max(mp,na.rm=TRUE)<collinear)
        }
        return(thr)})

      meta.temp=MR_mega_run(beta_pop[mcor.thr,],invse2_pop[mcor.thr,],cohort_count_filt[mcor.thr,],env,PCs,ncores)

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
