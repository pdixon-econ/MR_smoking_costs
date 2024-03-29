

#Smoking analysis 

rm(list=ls(all=TRUE))

###################################

#Variant level analysis

###################################

#New version - check everything is ok with this one

library(TwoSampleMR)
library(dplyr)
library(ggplot2)


setwd("//…Individual genotypes/")

##First, read in my file name with the non-standard headings

#INITIATION 1



exp_dat<-read_exposure_data("ss1_smoking_initiation_mrbase.csv",
                            sep = ",",
                            snp_col = "snp",
                            beta_col = "betaexposure",
                            se_col = "seexposure",
                            effect_allele_col = "effect_alleleexposure",
                            other_allele_col = "other_alleleexposure",
                            eaf_col = "eafexposure",
                            pval_col = "pvalexposure",
                            #units_col = "unitsexposure",
                            #gene_col = "Gene",
                            #samplesize_col = "samplesizeexposure"
)

##Clump this file - 

exp_dat <- clump_data(exp_dat)

#read in cost data

cost_out_dat<-read_outcome_data(
  filename="ss1_smoking_initiation_mrbase.csv",
  snps = exp_dat$SNP,
  sep = ",",
  snp_col = "snp",
  beta_col = "beta__outcome",
  se_col = "se__outcome",
  effect_allele_col = "effect_alleleexposure",
  other_allele_col = "other_alleleexposure",
  eaf_col = "eafexposure",
  pval_col = "p_outcome",
  #units_col = "unit_output",
  #gene_col = "Gene",
  #samplesize_col = "samplesize_out"
)
     




##Harmonise data 

dat <- harmonise_data(
  exposure_dat = exp_dat, 
  outcome_dat = cost_out_dat
)

##Run MR Mix run_mrmix(dat)



initiation.res.summary <- mr(dat, method_list=c("mr_ivw", "mr_egger_regression","mr_penalised_weighted_median","mr_weighted_mode"))
initiation.resegger <- mr(dat, method_list=c("mr_ivw", "mr_egger_regression"))
initiation.resmedian <- mr(dat, method_list=c("mr_ivw","mr_simple_median", "mr_weighted_median", "mr_penalised_weighted_median"))	
initiation.resmode <- mr(dat, method_list=c("mr_ivw", "mr_simple_mode", "mr_weighted_mode"))			
initiation.resivw <- mr(dat, method_list = "mr_ivw")


##Steiger test for directionality - 


samplesize.exposure<-rep(149995,10)
samplesize.outcome<-rep(149995,10)
dat$samplesize.outcome<-samplesize.outcome
dat$samplesize.exposure<-samplesize.exposure

steiger1 <- directionality_test(dat)
steiger1

# Heterogeneity and Egger intercept
initiation.mr_het <- mr_heterogeneity(dat)
initiation.mr_het
initiation.mr_egger_int <- mr_pleiotropy_test(dat)
initiation.mr_egger_int




# single SNP analyses
initiation.res_single <- mr_singlesnp(dat, all_method=c("mr_ivw", "mr_egger_regression", "mr_penalised_weighted_median", "mr_weighted_mode"))
initiation.res_single

# leave one out analyses
initiation.res_loo <- mr_leaveoneout(dat)
initiation.res_loo

####List results

initiation.res.summary 


#######GRAPHICS 

##Plain scatter 


#Scatter plot

initiation.main.scatter<- mr_scatter_plot(initiation.res.summary, dat)
initiation.main.scatter
ggsave(initiation.main.scatter[[1]], file="main_scatter.png", width=7, height=7)

initiation.ivw.scatter.only<-mr_scatter_plot(initiation.resivw,dat)
initiation.ivw.scatter.only
ggsave(initiation.ivw.scatter.only[[1]], file="ivw_scatter.png", width=7, height=7)

#FOREST PLOT 


initiation.res_single <- mr_singlesnp(dat)
initiation.forest.plot <- mr_forest_plot(initiation.res_single)
initiation.forest.plot[[1]]
ggsave(initiation.forest.plot[[1]], file="initiationforest.png", width=7, height=7)


#############################################################################
##Smoking_initiation split 2


##First, read in my file name with the non-standard headings

exp_dat<-read_exposure_data("ss2_smoking_initiation_mrbase.csv",
                            sep = ",",
                            snp_col = "snp",
                            beta_col = "betaexposure",
                            se_col = "seexposure",
                            effect_allele_col = "effect_alleleexposure",
                            other_allele_col = "other_alleleexposure",
                            eaf_col = "eafexposure",
                            pval_col = "pvalexposure",
                            #units_col = "unitsexposure",
                            #gene_col = "Gene",
                            #samplesize_col = "samplesizeexposure"
)

##Clump this file - 

exp_dat <- clump_data(exp_dat)

#read in cost data

cost_out_dat<-read_outcome_data(
  filename="ss2_smoking_initiation_mrbase.csv",
  snps = exp_dat$SNP,
  sep = ",",
  snp_col = "snp",
  beta_col = "beta__outcome",
  se_col = "se__outcome",
  effect_allele_col = "effect_alleleexposure",
  other_allele_col = "other_alleleexposure",
  eaf_col = "eafexposure",
  pval_col = "p_outcome",
  #units_col = "unit_output",
  #gene_col = "Gene",
  #samplesize_col = "samplesize_out"
)



##Harmonise data 

dat <- harmonise_data(
  exposure_dat = exp_dat, 
  outcome_dat = cost_out_dat
)


smoking_initiation.ss2.res.summary <- mr(dat, method_list=c("mr_ivw", "mr_egger_regression","mr_penalised_weighted_median","mr_weighted_mode"))


smoking_initiation.ss2.res.summary

write.csv(smoking_initiation.ss2.res.summary,'smoking_initiationss2ressummary.csv')


# Heterogeneity and Egger intercept
smoking_initiation.ss2.mr_het <- mr_heterogeneity(dat)
smoking_initiation.ss2.mr_het
smoking_initiation.ss2.mr_egger_int <- mr_pleiotropy_test(dat)
smoking_initiation.ss2.mr_egger_int



##Steiger test for directionality - 

samplesize.exposure<-rep(150050,11)
samplesize.outcome<-rep(150050,11)
dat$samplesize.outcome<-samplesize.outcome
dat$samplesize.exposure<-samplesize.exposure

steiger1 <- directionality_test(dat)


##############################


###LIFETIME SMOKING - 1



##First, read in my file name with the non-standard headings

exp_dat<-read_exposure_data("ss1_csi_mrbase.csv",
                            sep = ",",
                            snp_col = "snp",
                            beta_col = "betaexposure",
                            se_col = "seexposure",
                            effect_allele_col = "effect_alleleexposure",
                            other_allele_col = "other_alleleexposure",
                            eaf_col = "eafexposure",
                            pval_col = "pvalexposure",
                            #units_col = "unitsexposure",
                            #gene_col = "Gene",
                            #samplesize_col = "samplesizeexposure"
)

##Clump this file - 

exp_dat <- clump_data(exp_dat)

#read in cost data

cost_out_dat<-read_outcome_data(
  filename="ss1_csi_mrbase.csv",
  snps = exp_dat$SNP,
  sep = ",",
  snp_col = "snp",
  beta_col = "beta__outcome",
  se_col = "se__outcome",
  effect_allele_col = "effect_alleleexposure",
  other_allele_col = "other_alleleexposure",
  eaf_col = "eafexposure",
  pval_col = "p_outcome",
  #units_col = "unit_output",
  #gene_col = "Gene",
  #samplesize_col = "samplesize_out"
)



##Harmonise data 

dat <- harmonise_data(
  exposure_dat = exp_dat, 
  outcome_dat = cost_out_dat
)


lifetime_smoking.ss1.res.summary <- mr(dat, method_list=c("mr_ivw", "mr_egger_regression","mr_penalised_weighted_median","mr_weighted_mode"))


lifetime_smoking.ss1.res.summary

write.csv(lifetime_smoking.ss1.res.summary,'csiss1ressummary.csv')



# Heterogeneity and Egger intercept
lifetime_smoking.ss1.mr_het <- mr_heterogeneity(dat)
lifetime_smoking.ss1.mr_het
lifetime_smoking.ss1.mr_egger_int <- mr_pleiotropy_test(dat)
lifetime_smoking.ss1.mr_egger_int


##Steiger test for directionality - 


samplesize.exposure<-rep(149995,17)
samplesize.outcome<-rep(149995,17)
dat$samplesize.outcome<-samplesize.outcome
dat$samplesize.exposure<-samplesize.exposure

steiger1 <- directionality_test(dat)
steiger1


#############################################################################
##Lifetime_smoking split 2


##First, read in my file name with the non-standard headings

exp_dat<-read_exposure_data("ss2_csi_mrbase.csv",
                            sep = ",",
                            snp_col = "snp",
                            beta_col = "betaexposure",
                            se_col = "seexposure",
                            effect_allele_col = "effect_alleleexposure",
                            other_allele_col = "other_alleleexposure",
                            eaf_col = "eafexposure",
                            pval_col = "pvalexposure",
                            #units_col = "unitsexposure",
                            #gene_col = "Gene",
                            #samplesize_col = "samplesizeexposure"
)

##Clump this file - 

exp_dat <- clump_data(exp_dat)

#read in cost data

cost_out_dat<-read_outcome_data(
  filename="ss2_csi_mrbase.csv",
  snps = exp_dat$SNP,
  sep = ",",
  snp_col = "snp",
  beta_col = "beta__outcome",
  se_col = "se__outcome",
  effect_allele_col = "effect_alleleexposure",
  other_allele_col = "other_alleleexposure",
  eaf_col = "eafexposure",
  pval_col = "p_outcome",
  #units_col = "unit_output",
  #gene_col = "Gene",
  #samplesize_col = "samplesize_out"
)



##Harmonise data 

dat <- harmonise_data(
  exposure_dat = exp_dat, 
  outcome_dat = cost_out_dat
)

lifetime_smoking.ss2.res.summary <- mr(dat, method_list=c("mr_ivw", "mr_egger_regression","mr_penalised_weighted_median","mr_weighted_mode"))


lifetime_smoking.ss2.res.summary

write.csv(lifetime_smoking.ss2.res.summary,'lifetime_smokingss2ressummary.csv')


# Heterogeneity and Egger intercept
lifetime_smoking.ss2.mr_het <- mr_heterogeneity(dat)
lifetime_smoking.ss2.mr_het
lifetime_smoking.ss2.mr_egger_int <- mr_pleiotropy_test(dat)
lifetime_smoking.ss2.mr_egger_int


##Steiger test for directionality - 


samplesize.exposure<-rep(150050,15)
sampelsize.outcome<-rep(150050,15)
dat$samplesize.outcome<-samplesize.outcome
dat$samplesize.exposure<-samplesize.exposure

steiger1 <- directionality_test(dat)
steiger1

#########

#RISK

#This is the MRB_tot phenotype, not MRB_score

########


###RISK - 1



##First, read in my file name with the non-standard headings

exp_dat<-read_exposure_data("ss1_mrb_mrbase.csv",
                            sep = ",",
                            snp_col = "snp",
                            beta_col = "betaexposure",
                            se_col = "seexposure",
                            effect_allele_col = "effect_alleleexposure",
                            other_allele_col = "other_alleleexposure",
                            eaf_col = "eafexposure",
                            pval_col = "pvalexposure",
                            #units_col = "unitsexposure",
                            #gene_col = "Gene",
                            #samplesize_col = "samplesizeexposure"
)

##Clump this file - 

exp_dat <- clump_data(exp_dat)

#read in cost data

cost_out_dat<-read_outcome_data(
  filename="ss1_mrb_mrbase.csv",
  snps = exp_dat$SNP,
  sep = ",",
  snp_col = "snp",
  beta_col = "beta__outcome",
  se_col = "se__outcome",
  effect_allele_col = "effect_alleleexposure",
  other_allele_col = "other_alleleexposure",
  eaf_col = "eafexposure",
  pval_col = "p_outcome",
  #units_col = "unit_output",
  #gene_col = "Gene",
  #samplesize_col = "samplesize_out"
)



##Harmonise data 

dat <- harmonise_data(
  exposure_dat = exp_dat, 
  outcome_dat = cost_out_dat
)


mrb.ss1.res.summary <- mr(dat, method_list=c("mr_ivw", "mr_egger_regression","mr_penalised_weighted_median","mr_weighted_mode"))


mrb.ss1.res.summary

write.csv(mrb.ss1.res.summary,'mrbss1ressummary.csv')



# Heterogeneity and Egger intercept
mrb.ss1.mr_het <- mr_heterogeneity(dat)
mrb.ss1.mr_het
mrb.ss1.mr_egger_int <- mr_pleiotropy_test(dat)
mrb.ss1.mr_egger_int



##Steiger test for directionality - 

samplesize.exposure<-rep(150000,10)
sampelsize.outcome<-rep(150000,10)
dat$samplesize.outcome<-samplesize.outcome
dat$samplesize.exposure<-samplesize.exposure

steiger1 <- directionality_test(dat)


#############################################################################
##Mrb split 2


##First, read in my file name with the non-standard headings

exp_dat<-read_exposure_data("ss2_mrb_mrbase.csv",
                            sep = ",",
                            snp_col = "snp",
                            beta_col = "betaexposure",
                            se_col = "seexposure",
                            effect_allele_col = "effect_alleleexposure",
                            other_allele_col = "other_alleleexposure",
                            eaf_col = "eafexposure",
                            pval_col = "pvalexposure",
                            #units_col = "unitsexposure",
                            #gene_col = "Gene",
                            #samplesize_col = "samplesizeexposure"
)

##Clump this file - 

exp_dat <- clump_data(exp_dat)

#read in cost data

cost_out_dat<-read_outcome_data(
  filename="ss2_mrb_mrbase.csv",
  snps = exp_dat$SNP,
  sep = ",",
  snp_col = "snp",
  beta_col = "beta__outcome",
  se_col = "se__outcome",
  effect_allele_col = "effect_alleleexposure",
  other_allele_col = "other_alleleexposure",
  eaf_col = "eafexposure",
  pval_col = "p_outcome",
  #units_col = "unit_output",
  #gene_col = "Gene",
  #samplesize_col = "samplesize_out"
)



##Harmonise data 

dat <- harmonise_data(
  exposure_dat = exp_dat, 
  outcome_dat = cost_out_dat
)

mrb.ss2.res.summary <- mr(dat, method_list=c("mr_ivw", "mr_egger_regression","mr_penalised_weighted_median","mr_weighted_mode"))


mrb.ss2.res.summary

write.csv(mrb.ss2.res.summary,'mrbss2ressummary.csv')


# Heterogeneity and Egger intercept
mrb.ss2.mr_het <- mr_heterogeneity(dat)
mrb.ss2.mr_het
mrb.ss2.mr_egger_int <- mr_pleiotropy_test(dat)
mrb.ss2.mr_egger_int



##Steiger test for directionality - 

samplesize.exposure<-rep(150000,10)
sampelsize.outcome<-rep(150000,10)
dat$samplesize.outcome<-samplesize.outcome
dat$samplesize.exposure<-samplesize.exposure

steiger1 <- directionality_test(dat)

###################################################################################################

#MULTIVARIABLE MR

###################################################################################################


#install.packages("remotes")

library(remotes)
install_github("WSpiller/MVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)

library(MVMR)

##getting data



##Sample 1


dat_mvmr <- read.csv(file = "…mvmr_analysis_sample1.csv")

head(dat_mvmr)



##format these data

F.data.csi <- format_mvmr(BXGs = dat_mvmr[,c(4,6)],
                      BYG = dat_mvmr[,2],
                      seBXGs = dat_mvmr[,c(5,7)],
                      seBYG = dat_mvmr[,3],
                      RSID = dat_mvmr[,1])


head(F.data.csi)

F.data.ini <- format_mvmr(BXGs = dat_mvmr[,c(4,8)],
                          BYG = dat_mvmr[,2],
                          seBXGs = dat_mvmr[,c(5,9)],
                          seBYG = dat_mvmr[,3],
                          RSID = dat_mvmr[,1])


head(F.data.ini)

##upload the gencov file (experimental, 26 march 2021)

#gencov<-snpcov_mvmr()

###Testing for weak instruments

sres.csi <- strength_mvmr(r_input = F.data.csi, gencov = 0)
sres.ini <- strength_mvmr(r_input = F.data.ini, gencov = 0)

#Comment - not strong enough for MVMR...

##Test for horizontal pleiotropy
pres.csi <- pleiotropy_mvmr(r_input = F.data.csi, gencov = 0)
pres.ini <- pleiotropy_mvmr(r_input = F.data.ini, gencov = 0)

##Comment - need to add the actual covariance matrices in here to check...

##Estimate causal effects

res.csi <- ivw_mvmr(r_input = F.data.csi)
res.ini <- ivw_mvmr(r_input = F.data.ini)

#########################
##Sample 2


dat_mvmr <- read.csv(file = "…mvmr_analysis_sample2.csv")

head(dat_mvmr)


##format these data

F.data.csi.2 <- format_mvmr(BXGs = dat_mvmr[,c(4,6)],
                            BYG = dat_mvmr[,2],
                            seBXGs = dat_mvmr[,c(5,7)],
                            seBYG = dat_mvmr[,3],
                            RSID = dat_mvmr[,1])


head(F.data.csi.2)

F.data.ini.2 <- format_mvmr(BXGs = dat_mvmr[,c(4,8)],
                            BYG = dat_mvmr[,2],
                            seBXGs = dat_mvmr[,c(5,9)],
                            seBYG = dat_mvmr[,3],
                            RSID = dat_mvmr[,1])


head(F.data.ini.2)

###Testing for weak instruments

sres.csi.2 <- strength_mvmr(r_input = F.data.csi.2, gencov = 0)
sres.ini.2 <- strength_mvmr(r_input = F.data.ini.2, gencov = 0)

#Comment - not strong enough for MVMR...

##Test for horizontal pleiotropy
pres.csi.2 <- pleiotropy_mvmr(r_input = F.data.csi.2, gencov = 0)
pres.in.2i <- pleiotropy_mvmr(r_input = F.data.ini.2, gencov = 0)

##Comment - need to add the actual covariance matrices in here to check...

##Estimate causal effects

res.csi.2 <- ivw_mvmr(r_input = F.data.csi.2)
res.ini.2 <- ivw_mvmr(r_input = F.data.ini.2)

