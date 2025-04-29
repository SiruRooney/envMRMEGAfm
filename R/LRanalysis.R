#'@title The inverse of the matrix
#'@description Calculate the inverse of the matrix.
#'@param V.item The matrix.
#'@return Output the inverse of the matrix.
#'
Get_GenInverse<-function(V.item){
  temp <- eigen(V.item)
  idx <- which(abs(temp$values) > 10^-5)
  #cat("idx is ",idx,"\n")
  if(length(idx)>1){
    V.item.inv <- temp$vectors[, idx] %*% diag(1/temp$values[idx]) %*%
      t(temp$vectors[, idx])
  }else{
    V.item.inv <- temp$vectors[, idx] %*% as.matrix(1/temp$values[idx]) %*%
      t(temp$vectors[, idx])
  }

  return(V.item.inv)
}


#'@title Conditional analysis of the linear regression model
#'@description Estimate the allelic effect size and the standard errors of the target SNP conditional on the SNP set.
#'@param add.gwas  A dataframe which contains the GWAS information for the target SNP. GWAS information comprise "MARKERNAME","CHROMOSOME","POSITION","EA","NEA","EAF","N","BETA","SE".
#'@param sel.gwas  A dataframe which contains the GWAS information for SNP set.
#'@param addsel.ld A LD matrix containing the correlation matrix between the target SNP and the SNP set.
#'@param sel.ld A LD matrix containing the correlation between the SNP set.
#'@param Vp.med The median of variance of phenotype. Refer to GCTA-COJO.
#'@param collinear A threshold to filter out the target SNP in high LD with the SNP set. If the squared multiple correlation between the target SNP exceeds the threshold,
#'such as 0.9, the target SNP is ignored.
#'@param actual.geno An indicator to specify whether the true cohort-level LD structure is applied. If actual.geno is TRUE, the inputted LD structures are derived from the true individual-level genotype data.
#' If actual.geno is FALSE, the the inputted LD structures are derived from the approximation LD strucure.
#'@return Output a data.frame containing the estimated allelic effect and  the inverse of squared standard error for the target SNP.
#'@author Siru Wang
#'@importFrom stats median

meta_sel_condest<-function(add.gwas,sel.gwas,addsel.ld,sel.ld,Vp.med,collinear,actual.geno){

  if(!is.matrix(addsel.ld)){addsel.ld=as.matrix(addsel.ld)}
  if(!is.matrix(sel.ld)){sel.ld=as.matrix(sel.ld)}


  MSX.sel=2*sel.gwas$EAF*(1-sel.gwas$EAF)
  if(!is.matrix(MSX.sel)){MSX.sel<-as.matrix(MSX.sel)}

  MSX.add=2*add.gwas$EAF*(1-add.gwas$EAF)
  if(!is.matrix(MSX.add)){MSX.add<-as.matrix(MSX.add)}

  if(isTRUE(actual.geno)){
    Nd.sel=sel.gwas$N
    Nd.add=add.gwas$N
    if(!is.matrix(Nd.sel)){Nd.sel<-as.matrix(Nd.sel)}
    if(!is.matrix(Nd.add)){Nd.add<-as.matrix(Nd.add)}

    if(length(Nd.sel)!=1){
      Ndmin=apply(Nd.sel,1,function(x,y){apply(y,1,function(yy,xx){min(yy,xx)},x)},Nd.sel)
      Ndmin.insc=apply(Nd.sel,1,function(x,y){min(x,y)},Nd.add)
    }else{
      Ndmin=Nd.sel
      Ndmin.insc=min(Nd.sel,Nd.add)
    }
    #MSX.selN=matrixcalc::hadamard.prod(Ndmin*sqrt(MSX.sel%*%t(MSX.sel)),sel.ld)
    MSX.selN=Ndmin*sqrt(MSX.sel%*%t(MSX.sel))*sel.ld
    #add.selN=matrixcalc::hadamard.prod(Ndmin.insc*sqrt(MSX.add%*%t(MSX.sel)),t(addsel.ld))
    add.selN=Ndmin.insc*sqrt(MSX.add%*%t(MSX.sel))*t(addsel.ld)

    MSX.addN=MSX.add*Nd.add
  }else{
    #revised 250415
    Nd.sel=sel.gwas$Neff
    Nd.add=add.gwas$Neff
    if(!is.matrix(Nd.sel)){Nd.sel<-as.matrix(Nd.sel)}
    if(!is.matrix(Nd.add)){Nd.add<-as.matrix(Nd.add)}


    if(length(Nd.sel)!=1){
      Ndmin=apply(Nd.sel,1,function(x,y){apply(y,1,function(yy,xx){min(yy,xx)},x)},Nd.sel)
      Ndmin.insc=apply(Nd.sel,1,function(x,y){min(x,y)},Nd.add)
    }else{
      Ndmin=Nd.sel
      Ndmin.insc=min(Nd.sel,Nd.add)
    }

    #MSX.selN=matrixcalc::hadamard.prod(Ndmin*sqrt(MSX.sel%*%t(MSX.sel)),sel.ld)
    MSX.selN=Ndmin*sqrt(MSX.sel%*%t(MSX.sel))*sel.ld
    #add.selN=matrixcalc::hadamard.prod(Ndmin.insc*sqrt(MSX.add%*%t(MSX.sel)),t(addsel.ld))
    add.selN=Ndmin.insc*sqrt(MSX.add%*%t(MSX.sel))*t(addsel.ld)

    MSX.addN=MSX.add*Nd.add
  }


  Nd.medsel=median(Nd.sel)
  if(dim(MSX.selN)[2]!=1){
    inv.selN=Get_GenInverse(MSX.selN)
  }else{
    inv.selN=1/MSX.selN
  }
  inv.addN=1/MSX.addN

  if(dim(sel.ld)[2]>1){
    inv.sel.ld=Get_GenInverse(sel.ld)
  }else{
    inv.sel.ld=1/sel.ld
  }

  sq.mul.cor=as.vector(round(t(addsel.ld)%*%inv.sel.ld%*%addsel.ld,6))
  #cat("The squared multiple correlation is ",round(t(addsel.ld)%*%inv.sel.ld%*%addsel.ld,6),"\n")

  if(!is.na(sq.mul.cor)){
    if(sq.mul.cor<collinear){
      D1BETA1=MSX.sel*Nd.sel*sel.gwas$BETA
      cond.BETA=add.gwas$BETA-inv.addN%*%add.selN%*%inv.selN%*%D1BETA1
      cond.BETA=round(cond.BETA,6)
      #print("Calculating the residual variance in the conditional analysis....")
      bJ=inv.selN%*%D1BETA1
      if((Nd.medsel-length(bJ))<1){
        stop("no degree of freedom is left for the residues. The model is over-fitted. Please specify a more stringent p-value cut-off.")
      }

      if(isTRUE(actual.geno)){
        Ve=sum(MSX.sel*Nd.sel*bJ*sel.gwas$BETA)
        jma_Ve=((Nd.medsel-1)*Vp.med-Ve)/(Nd.medsel-dim(sel.gwas)[1])
        if(jma_Ve<0){
          stop("Residual variance is out of boundary. The model is over-fitted. Please specify a more stringent p-value cut-off.")
        }
        #cond.SE2=(jma_Ve-(MSX.add*Nd.add*cond.BETA*add.gwas$BETA)/(Nd.add-dim(sel.gwas)[1]-1))*((MSX.add*Nd.add)-add.selN%*%inv.selN%*%t(add.selN))/(MSX.add*Nd.add)^2
        cond.SE2=(jma_Ve-(MSX.add*Nd.add*cond.BETA*add.gwas$BETA)/(Nd.add-dim(sel.gwas)[1]-1))*(1/(MSX.add*Nd.add))
      }else{
        jma_Ve=Vp.med
        #cond.SE2=jma_Ve*((MSX.add*Nd.add)-add.selN%*%inv.selN%*%t(add.selN))/(MSX.add*Nd.add)^2
        cond.SE2=jma_Ve*1/(MSX.add*Nd.add)
      }

      cond.invSE2=round(cond.SE2^{-1})
      beta.se=data.frame(MARKERNAME=add.gwas$MARKERNAME,cond.BETA=cond.BETA,cond.SE=sqrt(cond.SE2),cond.invSE2=cond.invSE2,mul.cor=sq.mul.cor)
    }else{
      beta.se=data.frame(MARKERNAME=add.gwas$MARKERNAME,cond.BETA=NA,cond.SE=NA,cond.invSE2=NA,mul.cor=sq.mul.cor)
    }
  }else{
    beta.se=data.frame(MARKERNAME=add.gwas$MARKERNAME,cond.BETA=NA,cond.SE=NA,cond.invSE2=NA,mul.cor=sq.mul.cor)
  }

  return(beta.se)
}

#'@title Joint analysis of a linear regression model
#'@description Estimate the allelic effect size and the standard errors of the target SNP in the single SNP linear model.
#'@param sel.gwas A dataframe which contains the GWAS information for SNP set. GWAS information comprise "MARKERNAME","CHROMOSOME","POSITION","EA","NEA","EAF","N","BETA","SE".
#'@param sel.ld A LD matrix containing the correlation between the SNP set.
#'@param actual.geno An indicator to specify whether the true cohort-level LD structure is applied. If actual.geno is TRUE, the inputted LD structures are derived from the true individual-level genotype data.
#'@return Output a data.frame containing the estimated allelic effect, standard error, the inverse of squared standard error and the p-value of association test for the target SNP.
#'@importFrom stats median
#'@author Siru Wang
meta_sel_joinest<-function(sel.gwas,sel.ld,actual.geno){
  h.bufsel=2*sel.gwas$EAF*(1-sel.gwas$EAF)
  Vp.sel=sel.gwas$N*h.bufsel*sel.gwas$SE*sel.gwas$SE+h.bufsel*sel.gwas$BETA*sel.gwas$BETA*sel.gwas$N/(sel.gwas$N-1)
  Vp.medsel=median(Vp.sel)

  MSX.sel=2*sel.gwas$EAF*(1-sel.gwas$EAF)
  if(!is.matrix(MSX.sel)){MSX.sel=as.matrix(MSX.sel)}

  if(isTRUE(actual.geno)){
    Nd.sel=sel.gwas$N
    if(!is.matrix(Nd.sel)){Nd.sel=as.matrix(Nd.sel)}

    if(length(Nd.sel)!=1){
      Ndmin=apply(Nd.sel,1,function(x,y){apply(y,1,function(yy,xx){min(yy,xx)},x)},Nd.sel)
    }else{
      Ndmin=Nd.sel
    }

  }else{


    Nd.sel=(Vp.medsel-MSX.sel*sel.gwas$BETA*sel.gwas$BETA)/(MSX.sel*sel.gwas$SE*sel.gwas$SE)+1
    if(!is.matrix(Nd.sel)){Nd.sel=as.matrix(Nd.sel)}

    if(length(Nd.sel)!=1){
      Ndmin=apply(Nd.sel,1,function(x,y){apply(y,1,function(yy,xx){min(yy,xx)},x)},Nd.sel)
    }else{
      Ndmin=Nd.sel
    }
  }

  #MSX.selN=matrixcalc::hadamard.prod(Ndmin*sqrt(MSX.sel%*%t(MSX.sel)),sel.ld)
  MSX.selN=Ndmin*sqrt(MSX.sel%*%t(MSX.sel))*sel.ld
  Nd.medsel=median(Nd.sel)
  D1BETA1=MSX.sel*Nd.sel*sel.gwas$BETA

  if(dim(MSX.selN)[2]!=1){
    inv.selN=Get_GenInverse(MSX.selN)
  }else{
    inv.selN=1/MSX.selN
  }
  bJ_pop=inv.selN%*%D1BETA1
  bJ_se2_pop=diag(inv.selN)

  if(isTRUE(actual.geno)){
    Ve=sum(MSX.sel*Nd.sel*bJ_pop*sel.gwas$BETA)
    jma_Ve=((Nd.medsel-1)*Vp.medsel-Ve)/(Nd.medsel-dim(sel.gwas)[1])
    if(jma_Ve<0){
      stop("Residual variance is out of boundary. The model is over-fitted. Please specify a more stringent p-value cut-off.")
    }
  }else{
    jma_Ve=Vp.medsel
  }

  bJ_se2_pop=bJ_se2_pop*jma_Ve
  bJ_pop[(bJ_se2_pop<=1e-30)]=0
  bJ_se2_pop[(bJ_se2_pop<=1e-30)]=0

  bJ_se_pop=sqrt(bJ_se2_pop)
  chisq_bJ=bJ_pop[(bJ_se_pop>1e-30)]/bJ_se_pop[(bJ_se_pop>1e-30)]

  bJ_invse2_pop=rep(NA,length(bJ_pop))
  bJ_invse2_pop[(bJ_se2_pop>1e-30)]=round(bJ_se2_pop[(bJ_se2_pop>1e-30)]^{-1})

  p_assoc_jsel=rep(NA,length(bJ_pop))
  p_assoc_jsel[(bJ_se_pop>1e-30)]=pchisq(chisq_bJ*chisq_bJ,1,lower.tail=FALSE)
  ############################################################################################################
  beta.se=data.frame(MARKERNAME=sel.gwas$MARKERNAME,Join.BETA=bJ_pop,Join.SE=bJ_se_pop,Join.invSE2=bJ_invse2_pop,p_assoc_join=p_assoc_jsel)
  return(beta.se)

}
