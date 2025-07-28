# age and smoking EWAS
# QC phase 1
# checking for EWAS signal associated with age and smoking
# QUESTION - how are we coding smoking for child cohorts?
  # this will be maternal smoking (will flag different sites - add to check specific ones in output)
  # we can add an if statement looking for the maternal smoking var
# QUESTION - are there any cohorts that may be newborns with little maternal smoking? if so we may need a different approach
  # this could maybe be done by checking epigenetic age looks reasonable instead of these EWAS
  # cohorts to check: check senegal and india (mumbai) studies. and ?malawi. they probs have smoking data?
  # if they don't have smoking data then maybe we can do an if statement in the bash script and run a different R script

# Args
arguments <- commandArgs(T)

beta_file <- arguments[1]
cov_file <- arguments[2] # [winsorized_covariates_file from check_phenotypes]
meth_array <- arguments[3]
cellcount_file <- arguments[4] 
study_name <- arguments[5]
study_specific_vars <- arguments[6]
smoking_ewas_stats <- arguments[7]
smoking_ewas_report <- arguments[8]
age_ewas_stats <- arguments[9]
age_ewas_report <- arguments[10]

suppressPackageStartupMessages(library(meffil))

message("Reading in data and matching up samples across files")#######################################
load(cov_file)
load(beta_file)
load(cellcount_file)
participants <- as.character(intersect(colnames(norm.beta),covar$IID))
covar <- covar[covar$IID%in%participants,]
norm.beta <- norm.beta[,participants]

celltypes <- colnames(cellcount_file) # replace this - what is the name of the cell counts df??
# also is IID a column in the cellcounts file?

message("Setting up EWAS")#######################################

# 1. smoking ewas

# in case some cohorts only have smoking quantity,
# recode quantity into a factor if they don't have categories (0 = never, 1+= yes)
  # QUESTION - check that numerical smoking vars will indeed be called Smoking_numeric
  # rather than pack-years or something

if ("Smoking_factor" %in% colnames(covar)) {
  message("Smoking_factor variable already exists.")
} else if ("Smoking_numeric" %in% colnames(covar)) {
  # Create Smoking_factor based on Smoking_numeric
  covar$Smoking_factor <- ifelse(covar$Smoking_numeric == 0, "No", "Yes")
  message("Smoking_factor variable created based on Smoking_numeric.")
} else {
  message("Smoking_numeric column not found. Cannot create Smoking_factor. Please check your phenotype data")
}

# QUESTION - is there a config file option that says what array the methylation data is on?
featureset <- meffil:::guess.featureset(rownames(norm.beta))
ewas_threshold <- ifelse(meth_array == '450k', 2.4e-7,
                ifelse(meth_array %in% c('epic', 'epic2'), 9e-8, NA))
if(is.na(ewas_threshold)){
  message("ERROR: No EWAS threshold has been selected as the methylation array type does
          not match the expected values of 450k or EPIC. Please check your config file")
  # QUESTION - now we can either stop the script, or we can do a bonferroni correction based on the n of probes
}

ewas.parameters <- meffil.ewas.parameters(sig.threshold=ewas_threshold,  ## EWAS p-value threshold
                                          max.plots=10, ## plot at most 10 CpG sites
                                          qq.inflation.method="median",  ## measure inflation using median
                                          model="all") ## select default EWAS model; 

ewas_covars_smoking <- c("Age_numeric","Sex_factor",celltypes,study_specific_vars) # need to update these var names?
# QUESTION : might we need to stratify smoking EWAS by sex in some cohorts? - 
  # answer:No not at the moment

message("Starting smoking EWAS")#######################################

if ("Smoking_factor" %in% colnames(df)) {
  # make sure meth and covs are in same order
  participants <- as.character(covs$IID)
  meth.temp <- norm.beta[,participants]
  ewas.smoking <- meffil.ewas(meth.temp, variable=pheno$Smoking.factor, covariates=pheno.temp[,ewas_covars_smoking], sva=T, isva=F, random.seed=23) 
  # save out the ewas summary stats:
  ewas.out <- ewas.smoking$analyses
  save(ewas.out, file=paste0(smoking_ewas_stats,"_",study_name,".Robj"))
  # generate and save html report of EWAS:
  ewas.summary<-meffil.ewas.summary(ewas.smoking,meth.temp,parameters=ewas.parameters)                              
  meffil.ewas.report(ewas.summary, output.file=paste0(smoking_ewas_report,"_",study_name,".html"))
  # let's have a quick glance at how many hits we get:
  hits <- ewas.smoking$analyses$none$table
  message("There were",nrow(hits[hits$p.value<ewas_threshold,]),"DNAm sites associated with smoking in the unadjusted model")
  hits <- ewas.smoking$analyses$all$table
  message("There were",nrow(hits[hits$p.value<ewas_threshold,]),"DNAm sites associated with smoking in the adjusted model")
  
} else {
  # Run this command if 'Smoking_factor' column does not exist
  print("Smoking_factor not found in phenotype dataframe. If you study has smoking data,
  please ensure your phenotype dataframe has either Smoking_factor or Smoking_numeric")

}


# Age EWAS

message("Starting age EWAS")#######################################

ewas_covars_age <- c("Sex_factor",celltypes,study_specific_vars) # need to update these var names
participants <- as.character(covs$IID)
meth.temp <- norm.beta[,participants]

ewas.age <- meffil.ewas(meth.temp, variable=pheno$Smoking.factor, covariates=pheno.temp[,ewas_covars_age], sva=T, isva=F, random.seed=23) 
ewas.out <- ewas.age$analyses
save(ewas.out, file=paste0(age_ewas_stats,"_",study_name,".Robj"))
ewas.summary<-meffil.ewas.summary(ewas.age,meth.temp,parameters=ewas.parameters)                              
meffil.ewas.report(ewas.summary, output.file=paste0(age_ewas_report,"_",study_name,".html"))

hits <- ewas.age$analyses$none$table
message("There were",nrow(hits[hits$p.value<ewas_threshold,]),"DNAm sites associated with age in the unadjusted model")
hits <- ewas.age$analyses$all$table
message("There were",nrow(hits[hits$p.value<ewas_threshold,]),"DNAm sites associated with age in the adjusted model")

