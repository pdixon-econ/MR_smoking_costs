///STATA OBSERVATIONAL

cd "…GRS"


use ss_analysis_smoking, clear


cap replace MRB_tot="" if MRB_tot=="NA"

cap destring MRB_tot, replace


//Some descriptive analysis on the dataset 

su age_at_recruitment sex
su age_at_recruitment sex if csi>0

codebook sex
tab sex smk_init


//Simple observational regressions

reg cost_person i.sex age_at_recruitment  i.assessment_centre eco_tdi smk_init


reg cost_person i.sex age_at_recruitment  i.assessment_centre eco_tdi csi 

reg cost_person i.sex age_at_recruitment  i.assessment_centre eco_tdi csi if csi>0 

*note csi=0 is when smk_init=1

reg cost_person i.sex age_at_recruitment  i.assessment_centre eco_tdi MRB_score

reg cost_person i.sex age_at_recruitment  i.assessment_centre eco_tdi MRB_tot


////////////////////////////////// 

*percentage of variance explained by GRS in the phenotype


foreach num of numlist 1/2 {
    
	logistic  smk_init grs_init_`num' 
	
}

foreach num of numlist 1/2 {
    
	reg csi grs_csi_`num' 
	
}


foreach num of numlist 1/2 {
    
	reg MRB_score grs_mrb_`num' 
	
}


foreach num of numlist 1/2 {
    
	reg MRB_tot grs_mrb_tot_`num' 
	
}


//////////////////////////////////////
STATA CAUSAL

SMOKING INITIATION

*/

cd "…GRS"

//Import text file of PRS data for each of the three phenotypes

import delimited "smk_init_PRS_sample1.csv ",  clear
rename score1_sum grs_init_1
save smk_init_PRS_sample1.dta, replace

import delimited "smk_init_PRS_sample2.csv ",  clear
rename score1_sum grs_init_2
save smk_init_PRS_sample2.dta, replace

append using smk_init_PRS_sample1.dta



rename iid ieu_id
sort ieu_id

drop fid

save smk_init_PRS_sample12, replace


//csi
cd "…GRS"

import delimited "csi_PRS_sample1.csv ",  clear
rename score1_sum grs_csi_1
save csi_PRS_sample1.dta, replace

import delimited "csi_PRS_sample2.csv ",  clear
rename score1_sum grs_csi_2
save csi_PRS_sample2.dta, replace

append using csi_PRS_sample1.dta



rename iid ieu_id
sort ieu_id

drop fid

save csi_PRS_sample12, replace

//mrb

import delimited "mrb_score_PRS_sample1.csv ",  clear
rename score1_sum grs_mrb_1
save mrb_PRS_sample1.dta, replace

import delimited "mrb_score_PRS_sample2.csv ",  clear
rename score1_sum grs_mrb_2
save mrb_PRS_sample2.dta, replace

append using mrb_PRS_sample1.dta



rename iid ieu_id
sort ieu_id
drop fid

save mrb_PRS_sample12, replace


//mrb_tot

import delimited "mrb_tot_PRS_sample1.csv ",  clear
rename score1_sum grs_mrb_tot_1
save mrb_tot_PRS_sample1.dta, replace

import delimited "mrb_tot_PRS_sample2.csv ",  clear
rename score1_sum grs_mrb_tot_2
save mrb_tot_PRS_sample2.dta, replace

append using mrb_tot_PRS_sample1.dta



rename iid ieu_id
sort ieu_id
drop fid

save mrb_tot_PRS_sample12, replace


//

merge 1:1 ieu_id using mrb_PRS_sample12

drop _merge

merge 1:1 ieu_id using csi_PRS_sample12

drop _merge

merge 1:1 ieu_id using smk_init_PRS_sample12

drop _merge

save grs_smoking_ss.dta, replace

//Merge with cost data

cd "…Many traits analysis"


use  neid_pca_excl_pheno_grs, clear

cd "…GRS"


drop _merge


//Merge with the cost-pheno data (ipd)

merge 1:1 ieu_id using "…grs_smoking_ss.dta"


keep if _merge==3

rename _merge split_sample_merge_smoking

*Drop the old grs scores for avoidance of doubt and to make the code tractable 


drop grs_alcohol_intake-grs_type_2_diabetes




save ss_analysis_pd_merged_smoking, replace


use ss_analysis_pd_merged_smoking, clear

drop pca*

merge 1:1 ieu_id using "…split_pheno_covars", nogen force

drop sample


//Keep the analysis set

keep if cost_person_year!=.

//Put new id at start of the data window

sort ieu_id

//Add in the sample designation

merge 1:1 ieu_id using "…_sample_list.dta"

//Need to rename the sample coding here


tab _merge

keep if _merge==3

drop _merge

sort ieu_id

save ss_analysis_smoking, replace


//////////////////////////2SLS GRS regressions here

use  "….dta", clear
rename linkedid ieu_id
cap drop smk_init
cap gen smk_init=0 if csi!=.
replace smk_init=1 if csi>1 & csi<.
sort ieu_id

cd "…GRS"

merge ieu_id using ss_analysis_smoking


save ss_analysis_smoking, replace

use ss_analysis_smoking,clear



foreach num of numlist 1/2 {
parmby "ivreg2 cost_person (smk_init = grs_init_`num') age i.sex i.assessment_centre pc1-pc40 , first robust endog(smk_init) partial(i.assessment_centre)", saving(parms_grs_init_`num'.dta, replace)
}

foreach num of numlist 1/2 {
	
use parms_grs_init_`num', clear

keep in 1

save, replace

}



use parms_grs_init_1, clear

append using parms_grs_init_2.dta

drop parmseq

save parms_grs_init_meta.dta, replace

use parms_grs_init_meta.dta

metan estimate stderr, nograph


scalar effct=r(ES)
gen effect=effct

scalar ste=r(seES)
gen stderror=ste

scalar pv=r(p_z)
gen pvaluw=pv

scalar uci=r(ci_upp)
gen uci=uci

scalar lci=r(ci_low)
gen lci=lci


 keep in 1
 
 drop parm-_WT
 
 gen exposure="initiation"
 

 
save meta_split_grs_init, replace


//Continuous smoking



use ss_analysis_smoking,clear


foreach num of numlist 1/2 {
parmby "ivreg2 cost_person (csi = grs_csi_`num') age i.sex i.assessment_centre pc1-pc40 if sample==`num' ,first  robust endog(csi) partial(i.assessment_centre)", saving(parms_grs_csi_`num'.dta, replace)
}

foreach num of numlist 1/2 {
	
use parms_grs_csi_`num', clear

keep in 1

save, replace

}



use parms_grs_csi_1, clear

append using parms_grs_csi_2.dta

drop parmseq

save parms_grs_csi_meta.dta, replace

use parms_grs_csi_meta.dta

metan estimate stderr, nograph


scalar effct=r(ES)
gen effect=effct

scalar ste=r(seES)
gen stderror=ste

scalar pv=r(p_z)
gen pvaluw=pv

scalar uci=r(ci_upp)
gen uci=uci

scalar lci=r(ci_low)
gen lci=lci


 keep in 1
 
 drop parm-_WT
 
 gen exposure="csi"
 

 
save meta_split_grs_csi, replace


//mrb_score - risk




use ss_analysis_smoking,clear


foreach num of numlist 1/2 {
parmby "ivreg2 cost_person (MRB_score = grs_mrb_`num') age i.sex i.assessment_centre pc1-pc40 if sample==`num', first robust endog(MRB_score) partial(i.assessment_centre)", saving(parms_grs_mrb_`num'.dta, replace)
}

foreach num of numlist 1/2 {
	
use parms_grs_mrb_`num', clear

keep in 1

save, replace

}



use parms_grs_mrb_1, clear

append using parms_grs_mrb_2.dta

drop parmseq

save parms_grs_mrb_meta.dta, replace

use parms_grs_mrb_meta.dta

metan estimate stderr, nograph


scalar effct=r(ES)
gen effect=effct

scalar ste=r(seES)
gen stderror=ste

scalar pv=r(p_z)
gen pvaluw=pv

scalar uci=r(ci_upp)
gen uci=uci

scalar lci=r(ci_low)
gen lci=lci


 keep in 1
 
 drop parm-_WT
 
 gen exposure="mrb"
 

 
save meta_split_grs_mrb, replace


///


//mrb_tot - risk




use ss_analysis_smoking,clear


foreach num of numlist 1/2 {
parmby "ivreg2 cost_person (MRB_tot = grs_mrb_tot_`num') age i.sex i.assessment_centre pc1-pc40 if sample==`num', first robust endog(MRB_tot) partial(i.assessment_centre)", saving(parms_grs_mrb_tot_`num'.dta, replace)
}

foreach num of numlist 1/2 {
	
use parms_grs_mrb_tot_`num', clear

keep in 1

save, replace

}



use parms_grs_mrb_tot_1, clear

append using parms_grs_mrb_tot_2.dta

drop parmseq

save parms_grs_mrb_tot_meta.dta, replace

use parms_grs_mrb_tot_meta.dta

metan estimate stderr, nograph


scalar effct=r(ES)
gen effect=effct

scalar ste=r(seES)
gen stderror=ste

scalar pv=r(p_z)
gen pvaluw=pv

scalar uci=r(ci_upp)
gen uci=uci

scalar lci=r(ci_low)
gen lci=lci


 keep in 1
 
 drop parm-_WT
 
 gen exposure="mrb_tot"
 

 
save meta_split_grs_mrb_tot, replace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//1. Preparing data for SUMMARY MR to be ANALYZED IN MR BASE 

////INDIVIDUAL SNP ANALYSIS 


cd "…Individual genotypes"


import delimited ind_genotypes_all_snps_new.csv, clear

local bases "g c t a"
foreach i of local bases { 
rename (*_`i') (*)

}
drop fid
rename iid ieu_id

save all_ind_genotypes_new.dta, replace

use all_ind_genotypes_new, clear

//check the new file if it is ss_analysis new or something else  

merge 1:1 ieu_id using "…_analysis_new.dta", force

keep if cost_person!=.

keep if _merge==3

drop _merge

/////

save all_ind_geno_analysis_new.dta, replace


//////////

*INITIATION

//////////


import delim "…Smk_init_summ_stats_MRBase_sample1.csv", clear
cap rename ïsnp snp
drop in 11/17
list snp

////////

use all_ind_geno_analysis_new, clear



keep if sample == 2 //snps from sample 1, cost data from sample 2

local snps_1 "rs12996225 rs4856463 rs7665036 rs986391 rs240957 rs10093628 rs10767646 rs7105462 rs7173514 rs62068572"




global tflist ""
 global modseq=0
foreach X of local snps_1 {
global modseq=$modseq+1
tempfile tfcur
 parmby "regr cost_person `X' pc1-pc40 i.sex age i.assessment_centre, robust", command format(estimate min95 max95 %8.2f p %8.1e) idn($modseq) saving(`"`tfcur'"',replace) flist(tflist)
 }
 
 drop _all
append using $tflist
sort idnum command  parmseq
describe
by idnum command:list parm  estimate min95 max95 p,noobs


keep if parm!="_cons"
keep parm estimate stderr p

rename estimate beta_`num'_outcome
rename stderr se_`num'_outcome
rename parm snp

*Now need to keep only the coefficients from the SNPs

keep if strpos(snp,"rs")

//save files 

save cost_ini_ss1_betas_new.dta, replace


///Now sample 2


*Now, lookup genome wide significant SNPs so that they can be included in the local macros below. 

import delim "…_init_summ_stats_MRBase_sample2.csv", clear
cap rename ïsnp snp
drop in 12/17
list snp

///////////////



use all_ind_geno_analysis_new, clear


keep if sample == 1 //snps from sample 2, cost data from sample 1

local snps_2 "rs7567570 rs6745444	rs222446	rs458806	rs2240294	rs1043450	rs7113596	rs35534970	rs28681284	rs9747332	rs45577732"




global tflist ""
 global modseq=0
foreach X of local snps_2 {
global modseq=$modseq+1
tempfile tfcur
 parmby "regr cost_person `X' pc1-pc40 i.sex age i.assessment_centre, robust", command format(estimate min95 max95 %8.2f p %8.1e) idn($modseq) saving(`"`tfcur'"',replace) flist(tflist)
 }
 
 drop _all
append using $tflist
sort idnum command  parmseq
describe
by idnum command:list parm  estimate min95 max95 p,noobs


keep if parm!="_cons"
keep parm estimate stderr p

rename estimate beta_`num'_outcome
rename stderr se_`num'_outcome
rename parm snp

*Now need to keep only the coefficients from the SNPs

keep if strpos(snp,"rs")

save cost_ini_ss2_betas_new.dta, replace



///////////////////////////

*CSI

//////////////////////////

import delim "…_summ_stats_MRBase_sample1.csv", clear
cap rename ïsnp snp
list snp

////////

use all_ind_geno_analysis_new, clear


keep if sample == 2 //snps from sample 1, cost data from sample 2

local snps_1 "rs499257	rs2890772	rs4856591	rs326341	rs6852351	rs17159727	rs986391	rs16879271	rs10226228	rs10233018	rs10093628	rs113382419	rs11030088	rs9919670	rs6590701	rs4763463	rs7173514"




global tflist ""
 global modseq=0
foreach X of local snps_1 {
global modseq=$modseq+1
tempfile tfcur
 parmby "regr cost_person `X' pc1-pc40 i.sex age i.assessment_centre, robust", command format(estimate min95 max95 %8.2f p %8.1e) idn($modseq) saving(`"`tfcur'"',replace) flist(tflist)
 }
 
 drop _all
append using $tflist
sort idnum command  parmseq
describe
by idnum command:list parm  estimate min95 max95 p,noobs


keep if parm!="_cons"
keep parm estimate stderr p

rename estimate beta_`num'_outcome
rename stderr se_`num'_outcome
rename parm snp

*Now need to keep only the coefficients from the SNPs

keep if strpos(snp,"rs")

save cost_csi_ss1_betas_new.dta, replace



///Now sample 2


import delim "…CSI_summ_stats_MRBase_sample2.csv", clear
cap rename ïsnp snp
drop in 16/17
list snp


///////////////



use all_ind_geno_analysis_new, clear


keep if sample == 1 //snps from sample 2, cost data from sample 1

local snps_2 "rs10922907	rs7559547	rs1863161	rs16824949	rs263771	rs3866330	rs17657924	rs12553882	rs56116178	rs7948789	rs12897150	rs12901436	rs28669908	rs159058	rs151176846"




global tflist ""
 global modseq=0
foreach X of local snps_2 {
global modseq=$modseq+1
tempfile tfcur
 parmby "regr cost_person `X' pc1-pc40 i.sex age i.assessment_centre, robust", command format(estimate min95 max95 %8.2f p %8.1e) idn($modseq) saving(`"`tfcur'"',replace) flist(tflist)
 }
 
 drop _all
append using $tflist
sort idnum command  parmseq
describe
by idnum command:list parm  estimate min95 max95 p,noobs


keep if parm!="_cons"
keep parm estimate stderr p

rename estimate beta_`num'_outcome
rename stderr se_`num'_outcome
rename parm snp

*Now need to keep only the coefficients from the SNPs

keep if strpos(snp,"rs")

save cost_csi_ss2_betas_new.dta, replace





/////////////////////////////////////////

//MRB - comment - this is only for the MRB_tot as the MRB_score doesn't have sufficient SNPs 

////////////////////////////////////////


import delim "…MRB_tot_summ_stats_MRBase_sample1.csv", clear
cap rename ïsnp snp
list snp

////////

use all_ind_geno_analysis_new, clear


keep if sample == 2 //snps from sample 1, cost data from sample 2

local snps_1 " rs113642272 rs145410819 rs70600" 





global tflist ""
 global modseq=0
foreach X of local snps_1 {
global modseq=$modseq+1
tempfile tfcur
 parmby "regr cost_person `X' pc1-pc40 i.sex age i.assessment_centre, robust", command format(estimate min95 max95 %8.2f p %8.1e) idn($modseq) saving(`"`tfcur'"',replace) flist(tflist)
 }
 
 drop _all
append using $tflist
sort idnum command  parmseq
describe
by idnum command:list parm  estimate min95 max95 p,noobs


keep if parm!="_cons"
keep parm estimate stderr p

rename estimate beta_`num'_outcome
rename stderr se_`num'_outcome
rename parm snp

*Now need to keep only the coefficients from the SNPs

keep if strpos(snp,"rs")

save cost_mrb_tot_ss1_betas_new.dta, replace



///Now sample 2


*Now, lookup genome wide significant SNPs so that they can be included in the local macros below. 

import delim "…MRB_tot_summ_stats_MRBase_sample2.csv", clear
cap rename ïsnp snp
drop in 3
list snp


///////////////



use all_ind_geno_analysis_new, clear


keep if sample == 1 //snps from sample 2, cost data from sample 1

local snps_2 "rs9787523 rs12901069"




global tflist ""
 global modseq=0
foreach X of local snps_2 {
global modseq=$modseq+1
tempfile tfcur
 parmby "regr cost_person `X' pc1-pc40 i.sex age i.assessment_centre, robust", command format(estimate min95 max95 %8.2f p %8.1e) idn($modseq) saving(`"`tfcur'"',replace) flist(tflist)
 }
 
 drop _all
append using $tflist
sort idnum command  parmseq
describe
by idnum command:list parm  estimate min95 max95 p,noobs


keep if parm!="_cons"
keep parm estimate stderr p

rename estimate beta_`num'_outcome
rename stderr se_`num'_outcome
rename parm snp

*Now need to keep only the coefficients from the SNPs

keep if strpos(snp,"rs")

save cost_mrb_tot_ss2_betas_new.dta, replace





//////GET DATA READY FOR MRBASE  



//INITIATION 1
import delim "…_init_snps_clumped_r2_0.001.csv", clear


gen outcome="Cost per person-year"
replace exposure="Smoking initiation"

merge 1:1 snp using "…cost_ini_ss1_betas_new.dta" , nogen


drop v1 

rename  p p_outcome 

save "…ss1_smoking_initiation_mrbase", replace

//Allele flipping

// Swap values of effect_alleleexposure and other_alleleexposure
gen temp_effect = effect_alleleexposure
gen temp_other = other_alleleexposure
replace effect_alleleexposure = temp_other if betaexposure < 0
replace other_alleleexposure = temp_effect if betaexposure < 0
drop temp_effect temp_other

// Subtract eafexposure from 1 whenever betaexposure is negative
replace eafexposure = 1 - eafexposure if betaexposure < 0

// Change the value of betaexposure to be positive whenever it is negative
replace betaexposure = -betaexposure if betaexposure < 0



//

export delimited "…ss1_smoking_initiation_mrbase.csv", replace


//INITIATION 2 

import delim "…_init_snps_clumped_r2_0.001_sample2.csv", clear

gen outcome="Cost per person-year"
replace exposure="Smoking initiation"

merge 1:1 snp using "…cost_ini_ss2_betas_new.dta" , nogen


drop v1

rename  p p_outcome 

save "…ss2_smoking_initiation_mrbase", replace

//Allele flipping

// Swap values of effect_alleleexposure and other_alleleexposure
gen temp_effect = effect_alleleexposure
gen temp_other = other_alleleexposure
replace effect_alleleexposure = temp_other if betaexposure < 0
replace other_alleleexposure = temp_effect if betaexposure < 0
drop temp_effect temp_other

// Subtract eafexposure from 1 whenever betaexposure is negative
replace eafexposure = 1 - eafexposure if betaexposure < 0

// Change the value of betaexposure to be positive whenever it is negative
replace betaexposure = -betaexposure if betaexposure < 0


export delimited "…ss2_smoking_initiation_mrbase.csv", replace

//////////////////////////////////////////////////

//LIFETIME SMOKING 1

import delim "…CSI_snps_clumped_r2_0.001.csv", clear


gen outcome="Cost per person-year"
replace exposure="Lifetime smoking"

merge 1:1 snp using "…cost_csi_ss1_betas_new.dta" , nogen


drop v1 

rename  p p_outcome 

save "…ss1_csi_mrbase", replace

//Allele flipping

// Swap values of effect_alleleexposure and other_alleleexposure
gen temp_effect = effect_alleleexposure
gen temp_other = other_alleleexposure
replace effect_alleleexposure = temp_other if betaexposure < 0
replace other_alleleexposure = temp_effect if betaexposure < 0
drop temp_effect temp_other

// Subtract eafexposure from 1 whenever betaexposure is negative
replace eafexposure = 1 - eafexposure if betaexposure < 0

// Change the value of betaexposure to be positive whenever it is negative
replace betaexposure = -betaexposure if betaexposure < 0


export delimited "Y…\ss1_csi_mrbase.csv", replace


//LIFETIME SMOKING 2 

import delim .Y:\projects\ieu2\p1\013\working\data\Smoking\GWAS hits\CSI_snps_clumped_r2_0.001_sample2.csv", clear

gen outcome="Cost per person-year"
replace exposure="Lifetime smoking"


merge 1:1 snp using "…cost_csi_ss2_betas_new.dta" , nogen


drop v1

rename  p p_outcome 

save "…ss2_csi_mrbase", replace

//Allele flipping

// Swap values of effect_alleleexposure and other_alleleexposure
gen temp_effect = effect_alleleexposure
gen temp_other = other_alleleexposure
replace effect_alleleexposure = temp_other if betaexposure < 0
replace other_alleleexposure = temp_effect if betaexposure < 0
drop temp_effect temp_other

// Subtract eafexposure from 1 whenever betaexposure is negative
replace eafexposure = 1 - eafexposure if betaexposure < 0

// Change the value of betaexposure to be positive whenever it is negative
replace betaexposure = -betaexposure if betaexposure < 0


export delimited "…_csi_mrbase.csv", replace

//////////////////////////////////////////////////

//RISK - mrb tot

//RISK 1

import delim "…mrb_snps_tot_clumped_r2_0.001.csv", clear


gen outcome="Cost per person-year"
replace exposure="Risk"

merge 1:1 snp using "…\cost_mrb_tot_ss1_betas_new.dta" , nogen


drop v1 

rename  p p_outcome 

save "…ss1_mrb_mrbase", replace

//Allele flipping

// Swap values of effect_alleleexposure and other_alleleexposure
gen temp_effect = effect_alleleexposure
gen temp_other = other_alleleexposure
replace effect_alleleexposure = temp_other if betaexposure < 0
replace other_alleleexposure = temp_effect if betaexposure < 0
drop temp_effect temp_other

// Subtract eafexposure from 1 whenever betaexposure is negative
replace eafexposure = 1 - eafexposure if betaexposure < 0

// Change the value of betaexposure to be positive whenever it is negative
replace betaexposure = -betaexposure if betaexposure < 0


export delimited "…ss1_mrb_mrbase.csv", replace


//RISK 2 

import delim "…_snps_tot_clumped_r2_0.001_sample2.csv", clear

gen outcome="Cost per person-year"

replace exposure="Risk"

merge 1:1 snp using "…cost_mrb_tot_ss2_betas_new.dta" , nogen


drop v1

rename  p p_outcome 

save "…ss2_mrb_mrbase", replace

//Allele flipping

// Swap values of effect_alleleexposure and other_alleleexposure
gen temp_effect = effect_alleleexposure
gen temp_other = other_alleleexposure
replace effect_alleleexposure = temp_other if betaexposure < 0
replace other_alleleexposure = temp_effect if betaexposure < 0
drop temp_effect temp_other

// Subtract eafexposure from 1 whenever betaexposure is negative
replace eafexposure = 1 - eafexposure if betaexposure < 0

// Change the value of betaexposure to be positive whenever it is negative
replace betaexposure = -betaexposure if betaexposure < 0


export delimited "…ss2_mrb_mrbase.csv", replace

///////////////////////////////////////////////////////////////////////////////////////

///MVMR 


cd “…All SNPs for MVMR"




//Initiation

import delimited "Smk_init_all_snps_MRBase_sample1.csv", clear
rename beta init_beta
rename se init_se
cap rename ïsnp  snp
drop chr effect_allele other_allele eaf pval
drop in 28/29

save smk_init_all_snps_MRBase_sample1,replace

//Lifetime smoking

import delimited "CSI_all_snps_MRBase_sample1.csv", clear
rename beta csi_beta
rename se csi_se
cap rename ïsnp  snp
drop chr effect_allele other_allele eaf pval
drop in 28/29

save CSI_all_snps_MRBase_sample1,replace

//Risk tot

import delimited "MRB_tot_all_snps_MRBase_sample1.csv", clear
rename beta MRB_tot_beta
rename se MRB_tot_se
cap rename ïsnp  snp
drop chr effect_allele other_allele eaf pval
drop in 28/29

save MRB_tot_all_snps_MRBase_sample1,replace

//Merge these files
 merge 1:1 snp using CSI_all_snps_MRBase_sample1, nogen
 merge 1:1 snp using smk_init_all_snps_MRBase_sample1, nogen
 
 save mvmr_sample1, replace

 
 ///Sample 2
 
 
 

//Initiation

import delimited "Smk_init_all_snps_MRBase_sample2.csv", clear
rename beta init_beta
rename se init_se
cap rename ïsnp  snp
drop chr effect_allele other_allele eaf pval

save smk_init_all_snps_MRBase_sample2,replace

//Lifetime smoking

import delimited "CSI_all_snps_MRBase_sample2.csv", clear
rename beta csi_beta
rename se csi_se
cap rename ïsnp  snp
drop chr effect_allele other_allele eaf pval

save CSI_all_snps_MRBase_sample2,replace

//Risk tot

import delimited "MRB_tot_all_snps_MRBase_sample2.csv", clear
rename beta MRB_tot_beta
rename se MRB_tot_se
cap rename ïsnp  snp
drop chr effect_allele other_allele eaf pval

save MRB_tot_all_snps_MRBase_sample2,replace

//Merge these files
 merge 1:1 snp using CSI_all_snps_MRBase_sample2, nogen
 merge 1:1 snp using smk_init_all_snps_MRBase_sample2, nogen
 
 save mvmr_sample2, replace

 
 
 
 
 
 //SAMPLE 1 SNPS, COST FROM SAMPLE 2
 
//Getting outcome data - SNPs from sample 1, cost data from sample 2
cd "…Individual genotypes"

use all_ind_geno_analysis_new, clear

cd "…All SNPs for MVMR"


keep if sample == 2 //snps from sample 1, cost data from sample 2

local snps_1 "rs499257	rs12996225	rs2890772	rs4856463	rs113642272	rs4856591	rs326341	rs6852351	rs7665036	rs17159727	rs986391	rs16879271	rs240957	rs10226228	rs10233018	rs10093628	rs113382419	rs10767646	rs11030088	rs9919670	rs7105462	rs6590701	rs4763463	rs7173514	rs145410819	rs70600	rs62068572"



global tflist ""
 global modseq=0
foreach X of local snps_1 {
global modseq=$modseq+1
tempfile tfcur
 parmby "regr cost_person `X' pc1-pc40 i.sex age i.assessment_centre, robust", command format(estimate min95 max95 %8.2f p %8.1e) idn($modseq) saving(`"`tfcur'"',replace) flist(tflist)
 }
 
 drop _all
append using $tflist
sort idnum command  parmseq
describe
by idnum command:list parm  estimate min95 max95 p,noobs


keep if parm!="_cons"
keep parm estimate stderr p

rename estimate beta_`num'_outcome
rename stderr se_`num'_outcome
rename parm snp

*Now need to keep only the coefficients from the SNPs

keep if strpos(snp,"rs")


save cost_mvmr_ss1_betas.dta, replace

use "…cost_mvmr_ss1_betas.dta", clear

merge 1:1 snp using "…mvmr_sample1", nogen

drop p

 

export delimited "…mvmr_analysis_sample1.csv", replace





 
 //SAMPLE 2 SNPS, COST FROM SAMPLE 1
 
//Getting outcome data - SNPs from sample 1, cost data from sample 2
cd "…Individual genotypes"

use all_ind_geno_analysis_new, clear

cd "…All SNPs for MVMR"


keep if sample == 1 //snps from sample 2, cost data from sample 1

local snps_1 "rs1043450 rs10922907 rs12553882 rs12897150 rs12901069 rs12901436 rs138534356 rs151176846 rs159058 rs16824949 rs17657924 rs1863161 rs222446 rs2240294 rs263771 rs28669908 rs28681284 rs35534970 rs3866330 rs45577732 rs458806 rs56116178 rs6745444 rs7113596 rs7559547 rs7567570 rs7948789 rs9747332 rs9787523"


global tflist ""
 global modseq=0
foreach X of local snps_1 {
global modseq=$modseq+1
tempfile tfcur
 parmby "regr cost_person `X' pc1-pc40 i.sex age i.assessment_centre, robust", command format(estimate min95 max95 %8.2f p %8.1e) idn($modseq) saving(`"`tfcur'"',replace) flist(tflist)
 }
 
 drop _all
append using $tflist
sort idnum command  parmseq
describe
by idnum command:list parm  estimate min95 max95 p,noobs


keep if parm!="_cons"
keep parm estimate stderr p

rename estimate beta_`num'_outcome
rename stderr se_`num'_outcome
rename parm snp

*Now need to keep only the coefficients from the SNPs

keep if strpos(snp,"rs")


save cost_mvmr_ss2_betas.dta, replace

use cost_mvmr_ss2_betas, clear


merge 1:1 snp using "…mvmr_sample2", nogen

drop p


//flip alleles

foreach var in MRB_tot_beta  csi_beta  init_beta  {
    replace `var' = -`var' if `var' < 0
}


export delimited "...mvmr_analysis_sample2.csv", replace
