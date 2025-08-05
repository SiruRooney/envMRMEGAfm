#'@title Read-in GWAS files and LD structures
#'@description Read-in GWAS files and LD structures, and assign the first-read-in GWAS file as the reference panel. The following GWAS files and ld structures need to be adjusted based on
#'EA and NEA of the first GWAS file.
#'@param gwas.loc A file contains a vector of length K_g. Each component of the vector refers to the location of GWAS file.
#'@param ld.loc A file contains a vector of length K_ld. Each component of the vector refers to the location of LD strucrue file.
#'@param which.ld A character vector of length K_g. Each component of the vector corresponds to one LD structure. The length of which.ld should equal to the number of gwas files.
#'@param cohort_name A character vector of length K_g. Each component of the vector corresponds to one GWAS file.
#'@param ld_name A character vector of length K_ld. Each component of the vector corresponds to one LD structure.
#'@param out_loc Path to save pre-processed GWAS files and LD structures. By default, out_loc=NULL.
#'@return The output is a list with two main components. One component contains the pre-processed GWAS file list, and another component contains the pre-processed LD structures.
#'@author Siru Wang
#'@import data.table
#'@importFrom stats complete.cases
#'@export
read_gwasld<-function(gwas.loc,ld.loc,which.ld,cohort_name,ld_name,out_loc=NULL){
  if(length(cohort_name)!=length(gwas.loc)){
    stop("The length of cohort_name should be consistent with the number of gwas files!")
  }

  if(length(ld_name)!=length(ld.loc)){
    stop("The length of ld.loc should be consistent with the number of LD structure files!")
  }

  gwas.list=vector(mode="list",length=length(gwas.loc))
  if(!is.null(cohort_name)){
    names(gwas.list)=cohort_name
  }

  ld.list=vector(mode="list",length=length(ld.loc))
  if(!is.null(ld_name)){
    names(ld.list)=ld_name
  }
  sum_stat_nd=c("MARKERNAME","CHROMOSOME","POSITION","EA","NEA","EAF","N","BETA","SE")

  if (!is.null(gwas.loc)){
    for (i in 1:length(gwas.loc)){
      print(paste("i in filesize ",i,sep=""))

      print(paste("Reading gwas file ",gwas.loc[i],sep=""))

      gwas.list[[i]]=fread(gwas.loc[i],select=sum_stat_nd)

      smallsp<-(gwas.list[[i]]$N<=10)
      if(any(smallsp)){
        gwas.list[[i]]=gwas.list[[i]][!smallsp,]
        cat("There are ",sum(smallsp)," SNPs of which the corresponding sample sizes are invalid (<=10)\n")
      }
      cat("There are ",dim(gwas.list[[i]])[1]," SNPs in this GWAS file\n")
    }
  }else{
    stop("The option gwas.loc should contain the locations of several summary statistics files")
  }

  if(!is.null(ld.loc)){

    for(i in 1:length(ld.loc)){

      print(paste("i in filesize ",i,sep=""))
      print(paste("Reading LD file ",ld.loc[i],sep=""))

      ld.list[[i]]=fread(ld.loc[i])
      if(!is.matrix(ld.list[[i]])){
        ld.list[[i]]=as.matrix(ld.list[[i]]);
        print(dim(ld.list[[i]]))
      }
      rownames(ld.list[[i]])=colnames(ld.list[[i]])
    }
  }else{
    stop("The option ld.loc should contain the locations of LD structure files")
  }

  # Assign each GWAS file to the corresponding LD structure
  print("Assign each GWAS file to the corresponding LD structure......")

  if(length(gwas.list)!=length(which.ld)){
    stop("The length of which.ld should be consistent with the length of gwas.list!")
  }

  if(any(!which.ld%in%names(ld.list))){
    stop("All components of which.ld come from the names of ld.list!")
  }

  for(ld.nm in names(ld.list)){

    ind.ref=which(which.ld==ld.nm)
    snp_pop<-lapply(gwas.list[ind.ref],function(ssp){
      return(unique(ssp$MARKERNAME))
    })

    snp_pop_tot=unique(unlist(snp_pop))
    print(class(ld.list[[ld.nm]]))# expect matrix
    col.bad=!(colnames(ld.list[[ld.nm]])%in%snp_pop_tot)
    row.bad=!(rownames(ld.list[[ld.nm]])%in%snp_pop_tot)


    ld.list[[ld.nm]]=ld.list[[ld.nm]][!row.bad,][,!col.bad]

    print(dim(ld.list[[ld.nm]]))
  }

  ##Output the GWAS file after flipping the EA and NEA.

  if(!is.null(out_loc)){
    if(!is.null(cohort_name)){
      for(i in 1:length(gwas.loc)){
        print(paste0("Writing the updated GWAS file into ",out_loc,"/",cohort_name[i],".txt"))
        fwrite(gwas.list[[i]],paste0(out_loc,"/Ref_",cohort_name[i],".txt"),sep="\t")
      }
    }else{
      stop("If gwas files need to outputted, the cohort_name should not be null!")
    }

    if(!is.null(ld_name)){
      for(i in 1:length(ld.loc)){
        print(paste0("Writing the updated GWAS file into ",out_loc,"/",ld_name[i],".txt"))
        fwrite(ld.list[[i]],paste0(out_loc,"/Ref_",ld_name[i],".txt"),sep="\t")
      }
    }else{
      stop("If LD files need to be outputted, the ld_name should not be null!")
    }
  }

  return(list(gwas.list=gwas.list,ld.list=ld.list))
}
