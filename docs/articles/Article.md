# The script for running fine-mapping methods on Africa data

## Introduction

Here, we provide step-by-step R and slurm scripts that describe how to
run the proposed fine-mapping methods, MR-MEGAfm/env-MR-MEGAfm and
SuSiEx on the Africa data. For the three fine-mapping methods, they
require summary-level data and cohort-level LD matrices. Additionally,
(env-)MR-MEGAfm also requires MR-MEGA/env-MR-MEGA meta-analysis outputs
as inputs. In the following, we present the scripts for running the
three fine-mapping methods, which are shown in the slurm file.

## The R script for running (env-)MR-MEGAfm

### Input preparation

Before submitting slurm files, we write all locations of inputs as txt
files and create a new directory for storing the fine-mapping outputs.
Specifically, all required inputs for (env-)MR-MEGAfm include several
cohort-level summary stats data and LD matrices. The files
“AFR_gwas_loc.txt” and “AFR_ld_loc.txt” list all GWAS summary data
locations and the corresponding LD matrices locations.

1.  Make a directory containing all inputs of env-MR-MEGAfm or
    MR-MEGAfm.

- PCs (the previous PCs used in env-MR-MEGA/MR-MEGA meta-analysis).

- The document which contains GWAS file locations.

- The document which contains LD structure locations.

- Environment covariates (for env-MR-MEGA).

- The output files obtained from env-MR-MEGA meta-analysis or MR-MEGA
  meta-analysis.

2.  Make a directory containing outputs of env-MR-MEGAfm or MR-MEGAfm.

3.  Make sure that each GWAS file has mandatory column headers:
    MARKERNAME, CHROMOSOME, POSITION, EA, NEA, EAF, N, BETA, SE.

``` r
library(data.table)

# out_loc refers to an output directory
out_loc=as.character(Sys.getenv("out_loc"))
cat("The output location is ",out_loc,"\n")

# meta_loc refers to the input directory containing output of env-MR-MEGA meta-analysis, PCs, env, gwasfile.loc, ld.loc
meta_loc=as.character(Sys.getenv("meta_loc"))
cat("The meta-analysis location is ",meta_loc,"\n")

# gwasfile.loc contains GWAS file locations
gwasfile.loc=read.table(paste0(meta_loc,"/AFR_gwas_loc.txt"))
gwasfile.loc=as.matrix(gwasfile.loc)

# ld.loc contains the location of the cohort-level LD matrices.
ld.loc=read.table(paste0(meta_loc,"/AFR_ld_loc.txt"))
ld.loc=as.matrix(ld.loc)

shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(doParallel))
shhh(library(foreach))
```

### Read in inputs of (env-)MR-MEGAfm methods

In Africa data application of env-MR-MEGA, we used the twelve
sex-stratified Africa cohorts from AWI-Gen east Africa, AWI-Gen west
Africa and AWI-Gen south Africa, Uganda, Zulu case, Zulu control. The
corresponding five cohort-level LD matrices consist of AWI-Gen east
Africa LD, AWI-Gen west Africa LD, AWI-Gen south Africa LD, Uganda LD
and Zulu LD. It is noted that, no matter which LD patterns are used for
fine-mapping, the option `actual.geno` is set to FALSE, which is
recommended in GCTA-cojo software.

``` r
library(envMRMEGAfm)
cohort_name=c("east_f_chr19", "east_m_chr19", "south_f_chr19", "south_m_chr19", "uganda_f_chr19", "uganda_m_chr19", "west_f_chr19", "west_m_chr19", "zulu_f_case_chr19","zulu_f_control_chr19", "zulu_m_case_chr19", "zulu_m_control_chr19")

#Name LD structures  as ld_name
ld_name=c("east_chr19","south_chr19","uganda_chr19","west_chr19","zulu_chr19","zulu_chr19")


#If run env-MR-MEGAfm, read in the output of env-MR-MEGA analysis
envmeta.file=fread(paste0(meta_loc,"AFR_envmrmega_chr19.txt"))
#If run MR-MEGAfm, read in the output of MR-MEGA analysis
meta.file=fread(paste0(meta_loc,"AFR_mrmega_chr19.txt"))

# which.ld aims at linking each LD to the correct gwas file.
#For example, gwas files are listed in order as "esn_female","esn_male","gwd_female","gwd_male","lwk_female","lwk_male","msl_female","msl_male","yri_female","yri_male" (the order is same as cohort_name)
#If we use the true cohort-level LD, the names of LD (ld_name) are "esn","gwd","lwk","msl","yri".
#We make both esn_female and esn_male correspond to esn LD, gwd_female and gwd_male correspond to gwd LD.....
#So we set which.ld =rep(c("esn","gwd","lwk","msl","yri"),each=2)
gwasld.list=read_gwasld(gwasfile.loc,ld.loc,which.ld=rep(c("east_chr19","south_chr19","uganda_chr19","west_chr19","zulu_chr19", "zulu_chr19"),each=2),cohort_name,ld_name=ld_name,out_loc=NULL)

gwas.list=gwasld.list$gwas.list
ld_ref.list=gwasld.list$ld.list
rm(gwasld.list)

#Read in PC file stored in the inputted directory "meta_loc"
PCs=read.table(paste0(meta_loc,"/pcs.txt"))
```

### Run MR-MEGAfm and env-MR-MEGAfm analysis

``` r
library(envMRMEGAfm)

#If run env-MR-MEGAfm, read in env file stored in the inputted directory "meta_loc"
env=read.table(paste0(meta_loc,"/envs.txt"))
envMRMEGAfm.out=env_MR_MEGA_fm(gwas.list=gwas.list,ld.list=ld_ref.list,which.ld=rep(c("east_chr19","south_chr19","uganda_chr19","west_chr19","zulu_chr19", "zulu_chr19"),each=2),meta.file=envmeta.file,PCs,env=env,out_loc=out_loc,ncores=4,collinear=0.9,pvalue_cutoff=5e-8,cred.thr=0.99)


#If run MR-MEGAfm, just set env=NULL
env=NULL
MRMEGAfm.out=env_MR_MEGA_fm(gwas.list=gwas.list,ld.list=ld_ref.list,which.ld=rep(c("east_chr19","south_chr19","uganda_chr19","west_chr19","zulu_chr19", "zulu_chr19"),each=2),meta.file=envmeta.file,PCs,env=env,out_loc=out_loc,ncores=4,collinear=0.9,pvalue_cutoff=5e-8,cred.thr=0.99)
```

## The slurm script for calling the R script of running (env-)MR-MEGAfm

The R script named `env-MR-MEGAfm-running.R` contains above R commands
for running (env-)MR-MEGAfm, and then is called by
`Rscript env-MR-MEGAfm-running.R` in the corresponding slurm script,
which is submitted to HPC.

## The slurm script for running SuSiex

More details about using SuSiEx can be found in
[SuSiEx](https://github.com/getian107/SuSiEx) github ([Yuan et al.
2024](#ref-yuan2024fine)).

``` r
/home/sw2077/rds/hpc-work/Susiex/SuSiEx/bin/SuSiEx \
        --sst_file=/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_susiex/chr19_gwas_in/chr19_gwas_extracted/east_f_assoc2_chr19.txt,/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_susiex/chr19_gwas_in/chr19_gwas_extracted/east_m_assoc2_chr19.txt,/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_susiex/chr19_gwas_in/chr19_gwas_extracted/south_f_assoc2_chr19.txt,/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_susiex/chr19_gwas_in/chr19_gwas_extracted/south_m_assoc2_chr19.txt,/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_susiex/chr19_gwas_in/chr19_gwas_extracted/west_f_assoc2_chr19.txt,/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_susiex/chr19_gwas_in/chr19_gwas_extracted/west_m_assoc2_chr19.txt,/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_susiex/chr19_gwas_in/chr19_gwas_extracted/uganda_f_assoc2_chr19.txt,/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_susiex/chr19_gwas_in/chr19_gwas_extracted/uganda_m_assoc2_chr19.txt,/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_susiex/chr19_gwas_in/chr19_gwas_extracted/zulu_f_case_assoc2_chr19.txt,/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_susiex/chr19_gwas_in/chr19_gwas_extracted/zulu_m_case_assoc2_chr19.txt,/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_susiex/chr19_gwas_in/chr19_gwas_extracted/zulu_f_control_assoc2_chr19.txt,/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_susiex/chr19_gwas_in/chr19_gwas_extracted/zulu_m_control_assoc2_chr19.txt \
        --n_gwas=961,805,3040,2225,1984,1883,3501,2618,1224,375,675,298 \
        --ref_file=/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_plk/east_chr19,/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_plk/east_chr19,/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_plk/south_chr19,/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_plk/south_chr19,/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_plk/west_chr19,/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_plk/west_chr19,/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_plk/ugr_chr19_region,/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_plk/ugr_chr19_region,/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_plk/zulu_chr19_region,/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_plk/zulu_chr19_region,/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_plk/zulu_chr19_region,/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_plk/zulu_chr19_region \
        --ld_file=/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_plk/east_chr19,/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_plk/east_chr19,/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_plk/south_chr19,/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_plk/south_chr19,/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_plk/west_chr19,/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_plk/west_chr19,/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_plk/ugr_chr19_region,/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_plk/ugr_chr19_region,/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_plk/zulu_chr19_region,/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_plk/zulu_chr19_region,/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_plk/zulu_chr19_region,/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_plk/zulu_chr19_region \
        --out_dir=/home/sw2077/rds/hpc-work/StepwsCond/Africadta_app/rs741_susiex/chr19_gwas_out \
        --out_name=SuSiEx.chr19.cs99 \
        --level=0.99 \
        --chr=19 \
        --bp=44912079,45912079 \
        --snp_col=2,2,2,2,2,2,2,2,2,2,2,2 \
        --chr_col=1,1,1,1,1,1,1,1,1,1,1,1 \
        --bp_col=3,3,3,3,3,3,3,3,3,3,3,3 \
        --a1_col=4,4,4,4,4,4,4,4,4,4,4,4 \
        --a2_col=5,5,5,5,5,5,5,5,5,5,5,5 \
        --eff_col=8,8,8,8,8,8,8,8,8,8,8,8 \
        --se_col=9,9,9,9,9,9,9,9,9,9,9,9 \
        --pval_col=10,10,10,10,10,10,10,10,10,10,10,10 \
        --plink=../utilities/plink\
        --mult-step=True \
        --min_purity=0.0001 \
        --keep-ambig=True \
        --threads=16
```

## Reference

Yuan, Kai, Ryan J Longchamps, Antonio F Pardiñas, Mingrui Yu, Tzu-Ting
Chen, Shu-Chin Lin, Yu Chen, et al. 2024. “Fine-Mapping Across Diverse
Ancestries Drives the Discovery of Putative Causal Variants Underlying
Human Complex Traits and Diseases.” *Nature Genetics* 56 (9): 1841–50.
