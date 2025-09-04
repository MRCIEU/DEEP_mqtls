# age and smoking EWAS
# QC phase 1
# checking for EWAS signal associated with age and smoking

# Args
arguments <- commandArgs(T)

beta_file <- arguments[1]
updated_pheno_file <- arguments[2] # [phenotype file with predicted smoking added in]
meth_array <- arguments[3]
cellcounts_cov <- arguments[4] 
cellcount_panel <- arguments[5] 
study_name <- arguments[6]
study_specific_vars <- strsplit(arguments[7], " ")[[1]] # these will be added in the config file - batch vars and study specific factors
ewas_stats <- arguments[8]
ewas_report <- arguments[9]


suppressPackageStartupMessages(library(meffil))

message("Reading in data and matching up samples across files")#######################################
load(updated_pheno_file)
load(beta_file)
cell_counts <- read.table(cellcounts_cov, header=T)
rownames(cell_counts) <- cell_counts$IID
participants <- as.character(intersect(colnames(norm.beta),pheno$IID))
pheno <- pheno[pheno$IID%in%participants,]
norm.beta <- norm.beta[,participants]

if (cellcount_panel == "unilife") {
  celltypes <- grep("^unilife", colnames(cell_counts), value = TRUE)
} else if (cellcount_panel == "salas") {
  celltypes <- grep("^salas", colnames(cell_counts), value = TRUE)
} else {
  message("Error: cell count panel not detected")
}

# add cell counts to pheno

cellcounts_temp <- cell_counts[,celltypes]
pheno <- merge(pheno,cellcounts_temp,by="IID")

message("Setting up EWAS")#######################################

# 1. smoking ewas

# in case some cohorts only have smoking quantity,
# recode quantity into a factor if they don't have categories (0 = never, 1+= yes)
  # QUESTION - check that numerical smoking vars will indeed be called Smoking_numeric
  # rather than pack-years or something

if ("Smoking_factor" %in% colnames(pheno)) {
  message("Smoking_factor variable already exists.")
} else if ("Smoking_numeric" %in% colnames(pheno)) {
  # Create Smoking_factor based on Smoking_numeric
  pheno$Smoking_factor <- ifelse(pheno$Smoking_numeric == 0, "No", "Yes")
  message("Smoking_factor variable created based on Smoking_numeric.")
} else {
  message("Smoking_numeric column not found. Cannot create Smoking_factor. 
          This is expected if this is a child dataset; maternal smoking EWAS will run below.
          If this is an adult dataset, please check your phenotype data")
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
  # make sure meth and pheno are in same order
  participants <- as.character(pheno$IID)
  meth.temp <- norm.beta[,participants]
  ewas.smoking <- meffil.ewas(meth.temp, variable=pheno$Smoking_factor, covariates=pheno[,ewas_covars_smoking], sva=T, isva=F, random.seed=23) 
  # save out the ewas summary stats:
  ewas.out <- ewas.smoking$analyses
  save(ewas.out, file=paste0(ewas_stats,"_smoking_",study_name,"_",cellcount_panel,".Robj"))
  # generate and save html report of EWAS:
  ewas.summary<-meffil.ewas.summary(ewas.smoking,meth.temp,parameters=ewas.parameters)                              
  meffil.ewas.report(ewas.summary, output.file=paste0(ewas_report,"_smoking_",study_name,"_",cellcount_panel,".html"))
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


# 2. Age EWAS

message("Starting age EWAS")#######################################

ewas_covars_age <- c("Sex_factor","p_smoking_mcigarette",celltypes,study_specific_vars) # need to update these var names
participants <- as.character(pheno$IID)
meth.temp <- norm.beta[,participants]

ewas.age <- meffil.ewas(meth.temp, variable=pheno$Age_numeric, covariates=pheno[,ewas_covars_age], sva=T, isva=F, random.seed=23) 
ewas.out <- ewas.age$analyses
save(ewas.out, file=paste0(ewas_stats,"_age_",study_name,"_",cellcount_panel,".Robj"))
ewas.summary<-meffil.ewas.summary(ewas.age,meth.temp,parameters=ewas.parameters)                              
meffil.ewas.report(ewas.summary, output.file=paste0(ewas_report,"_age_",study_name,"_",cellcount_panel,".html"))

hits <- ewas.age$analyses$none$table
message("There were",nrow(hits[hits$p.value<ewas_threshold,]),"DNAm sites associated with age in the unadjusted model")
hits <- ewas.age$analyses$all$table
message("There were",nrow(hits[hits$p.value<ewas_threshold,]),"DNAm sites associated with age in the adjusted model")


# 3. maternal smoking ewas

# in case some cohorts only have smoking quantity,
# recode quantity into a factor if they don't have categories (0 = never, 1+= yes)
# QUESTION - check that numerical smoking vars will indeed be called Smoking_numeric
# rather than pack-years or something

if ("maternal_smoking_factor" %in% colnames(pheno)) {
  message("maternal_smoking_factor variable already exists.")
} else if ("maternal_smoking_numeric" %in% colnames(pheno)) {
  # Create maternal_smoking_factor based on Smoking_numeric
  pheno$maternal_smoking_factor <- ifelse(pheno$maternal_smoking_numeric == 0, "No", "Yes")
  message("maternal_smoking_factor variable created based on maternal_smoking_numeric")
} else {
  message("maternal_smoking_numeric column not found. Cannot create maternal_smoking_factor.
          This is expected if this is an adult dataset.
          If this is a child dataset, Please check your phenotype data")
}

ewas_covars_mat_smoking <- c("Age_numeric","Sex_factor",celltypes,study_specific_vars) # need to update these var names?

message("Starting maternal smoking EWAS")#######################################

if ("maternal_smoking_factor" %in% colnames(df)) {
  # make sure meth and pheno are in same order
  participants <- as.character(pheno$IID)
  meth.temp <- norm.beta[,participants]
  ewas.smoking <- meffil.ewas(meth.temp, variable=pheno$maternal_smoking_factor, covariates=pheno[,ewas_covars_mat_smoking], sva=T, isva=F, random.seed=23) 
  # save out the ewas summary stats:
  ewas.out <- ewas.smoking$analyses
  save(ewas.out, file=paste0(ewas_stats,"_maternal_smoking_",study_name,"_",cellcount_panel,".Robj"))
  # generate and save html report of EWAS:
  ewas.summary<-meffil.ewas.summary(ewas.smoking,meth.temp,parameters=ewas.parameters)                              
  meffil.ewas.report(ewas.summary, output.file=paste0(ewas_report,"_maternal_smoking_",study_name,"_",cellcount_panel,".html"))
  # let's have a quick glance at how many hits we get:
  hits <- ewas.smoking$analyses$none$table
  message("There were",nrow(hits[hits$p.value<ewas_threshold,]),"DNAm sites associated with smoking in the unadjusted model")
  hits <- ewas.smoking$analyses$all$table
  message("There were",nrow(hits[hits$p.value<ewas_threshold,]),"DNAm sites associated with smoking in the adjusted model")
  
} else {
  # Run this command if 'maternal_smoking_factor' column does not exist
  print("maternal_smoking_factor not found in phenotype dataframe. 
  This is expected if this is an adult dataset.
  If this is a child dataset and your study has smoking data,
  please ensure your phenotype dataframe has either maternal_smoking_factor or maternal_smoking_factor")
  
}

# 4. Sex EWAS

message("Starting sex EWAS")#######################################

ewas_covars_age <- c("Age_numeric","p_smoking_mcigarette",celltypes,study_specific_vars) # need to update these var names
participants <- as.character(pheno$IID)
meth.temp <- norm.beta[,participants]

ewas.sex <- meffil.ewas(meth.temp, variable=pheno$Sex_factor, covariates=pheno[,ewas_covars_age], sva=T, isva=F, random.seed=23) 
ewas.out <- ewas.sex$analyses
save(ewas.out, file=paste0(ewas_stats,"_sex_",study_name,"_",cellcount_panel,".Robj"))
ewas.summary<-meffil.ewas.summary(ewas.sex,meth.temp,parameters=ewas.parameters)                              
meffil.ewas.report(ewas.summary, output.file=paste0(ewas_report,"_sex_",study_name,"_",cellcount_panel,".html"))

hits <- ewas.sex$analyses$none$table
message("There were",nrow(hits[hits$p.value<ewas_threshold,]),"DNAm sites associated with age in the unadjusted model")
hits <- ewas.sex$analyses$all$table
message("There were",nrow(hits[hits$p.value<ewas_threshold,]),"DNAm sites associated with age in the adjusted model")


# 5. Scrambled Sex EWAS (negative control)

message("Starting scrambled sex (negative control) EWAS")#######################################

ewas_covars_age <- c("Age_numeric","p_smoking_mcigarette",celltypes,study_specific_vars) # need to update these var names
participants <- as.character(pheno$IID)
meth.temp <- norm.beta[,participants]
pheno$Sex_factor <- sample(pheno$Sex_factor)

ewas.sex <- meffil.ewas(meth.temp, variable=pheno$Sex_factor, covariates=pheno[,ewas_covars_age], sva=T, isva=F, random.seed=23) 
ewas.out <- ewas.sex$analyses
save(ewas.out, file=paste0(ewas_stats,"_sex_negative_control_",study_name,"_",cellcount_panel,".Robj"))
ewas.summary<-meffil.ewas.summary(ewas.sex,meth.temp,parameters=ewas.parameters)                              
meffil.ewas.report(ewas.summary, output.file=paste0(ewas_report,"_sex_negative_control_",study_name,"_",cellcount_panel,".html"))

hits <- ewas.sex$analyses$none$table
message("There were",nrow(hits[hits$p.value<ewas_threshold,]),"DNAm sites associated with age in the unadjusted model")
hits <- ewas.sex$analyses$all$table
message("There were",nrow(hits[hits$p.value<ewas_threshold,]),"DNAm sites associated with age in the adjusted model")
