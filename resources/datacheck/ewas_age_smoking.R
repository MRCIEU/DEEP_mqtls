# age and smoking EWAS
# QC phase 1
# checking for EWAS signal associated with age and smoking
# QUESTION - what are we adjusting for in these EWAS? Are we including cell counts?
# QUESTIOn - how are we coding smoking for child cohorts?
# QUESTION - are there any cohorts that may be newborns with little maternal smoking? if so we may need a different approach

# Args
arguments <- commandArgs(T)

beta_file <- arguments[1]
cov_file <- arguments[2] # [winsorized_covariates_file from check_phenotypes]
meth_array <- arguments[3]
  # QUESTION - do we need cell counts here?
#cellcount_file <- arguments[9] 
study_name <- arguments[x]
smoking_ewas_stats <- arguments[x]
smoking_ewas_report <- arguments[x]
age_ewas_stats <- arguments[x]
age_ewas_report <- arguments[x]

suppressPackageStartupMessages(library(meffil))

message("Reading in data and matching up samples across files")#######################################
covs <- load(cov_file)

# QUESTION - do we require all DNAm samples to have genetic data?
# we need to match this up to what happens in check_phenotypes
# we don't need to load the fam file if we subset the covar file to whichever participants
# we want in check_phenotypes.R (commonids_mgc line)
# but equally if not all dnam participants have genetic data, we should remove this subsetting
fam <- read.table(fam_file)[,c(1,2)]
covs <- covs[match(fam[,2], covs[,"IID"]),]

load(beta_file)
norm.beta <- norm.beta[, match(fam[,2], colnames(norm.beta))]
message(paste(nrow(fam), "samples with genetic data matched to methylation data"))

message("Setting up EWAS")#######################################

# 1. smoking ewas

# QUESTION - will all smoking be categorical? will some cohorts only have pack years?
# if so, we will need to add an extra if statement to accomodate this

# QUESTION - is there a config file option that says what array the methylation data is on?
# if so we'll pull that in here:
ewas_threshold <- ifelse(meth_array == '450k', 2.4e-7,
                ifelse(meth_array == 'EPIC', 9e-8, NA))
if(is.na(ewas_threshold)){
  message("ERROR: No EWAS threshold has been selected as the methylation array type does
          not match the expected values of 450k or EPIC. Please check your config file")
  # QUESTION - now we can either stop the script, or we can do a bonferroni correction based on the n of probes
}

ewas.parameters <- meffil.ewas.parameters(sig.threshold=ewas_threshold,  ## EWAS p-value threshold
                                          max.plots=10, ## plot at most 10 CpG sites
                                          qq.inflation.method="median",  ## measure inflation using median
                                          model="all") ## select default EWAS model; 

# QUESTION - which covars do we want for the adjusted model?
ewas_covars_smoking <- c("Age_numeric","Sex_factor") # need to updat these var names
# QUESTION : might we need to stratify smoking EWAS by sex in some cohorts?
# could that look something like if 4x as many men smoke (for example)

message("Starting smoking EWAS")#######################################

if ("Smoking_factor" %in% colnames(df)) {
  # make sure meth and covs are in same order
  participants <- as.character(covs$IID)
  meth.temp <- norm.beta[,participants]
  ewas.smoking <- meffil.ewas(meth.temp, variable=pheno$Smoking.factor, covariates=pheno.temp[,ewas_covars_smoking], sva=T, isva=F, random.seed=23) 
  # QUESTION - do we want to save out the ewas summary stats? We might only want the report?
  save(ewas.smoking, file=paste0(smoking_ewas_stats,"_",study_name,".Robj"))
  ewas.summary<-meffil.ewas.summary(ewas.smoking,meth.temp,parameters=ewas.parameters)                              
  meffil.ewas.report(ewas.summary, output.file=paste0(smoking_ewas_report,"_",study_name,".html"))
  hits <- ewas.smoking$analyses$none$table
  message("There were",nrow(hits[hits$p.value<ewas_threshold,]),"DNAm sites associated with smoking in the unadjusted model")
  hits <- ewas.smoking$analyses$all$table
  message("There were",nrow(hits[hits$p.value<ewas_threshold,]),"DNAm sites associated with smoking in the adjusted model")
  
} else {
  # Run this command if 'Smoking_factor' column does not exist
  print("Smoking_factor not found in phenotype dataframe (if this is a measure you should have, 
  please check the name of the smoking variable)")

}


# Age EWAS

message("Starting age EWAS")#######################################

ewas_covars_age <- c("Sex_factor") # need to update these var names
participants <- as.character(covs$IID)
meth.temp <- norm.beta[,participants]

ewas.age <- meffil.ewas(meth.temp, variable=pheno$Smoking.factor, covariates=pheno.temp[,ewas_covars_age], sva=T, isva=F, random.seed=23) 
# QUESTION - do we want to save out the ewas summary stats? We might only want the report?
save(ewas.age, file=paste0(age_ewas_stats,"_",study_name,".Robj"))
ewas.summary<-meffil.ewas.summary(ewas.age,meth.temp,parameters=ewas.parameters)                              
meffil.ewas.report(ewas.summary, output.file=paste0(age_ewas_report,"_",study_name,".html"))

hits <- ewas.age$analyses$none$table
message("There were",nrow(hits[hits$p.value<ewas_threshold,]),"DNAm sites associated with age in the unadjusted model")
hits <- ewas.age$analyses$all$table
message("There were",nrow(hits[hits$p.value<ewas_threshold,]),"DNAm sites associated with age in the adjusted model")

