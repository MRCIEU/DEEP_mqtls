# Calculate Hannum age epigenetic clock and smoking status using mcigarette
# This is a check of how well our age phenotype correlates with age as predicted by the Hannum clock
# We have chosen the Hannum clock over others as it uses relatively few cpgs; 
# takes relatively little time to calculate; 
# and is a first generation clock that tries to predict chronological age.

# QUESTION - check that numerical smoking vars will indeed be called Smoking_numeric

# Args
arguments <- commandArgs(T)

beta_file <- arguments[1]
cov_file <- arguments[2] # [winsorized_covariates_file from check_phenotypes]
fam_file <- arguments[3]
out_file <- arguments[4]
age_stats <- arguments[5]
study_name <- arguments[6]
prediction_plot <- arguments[7]

suppressPackageStartupMessages(library(meffonym))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(ggpubr))

message("Reading in data and matching up samples across files")#######################################
load(cov_file)
load(beta_file)
participants <- as.character(intersect(colnames(norm.beta),covar$IID))
covar <- covar[covar$IID%in%participants,]
norm.beta <- norm.beta[,participants]
message("Number of samples with covariate and methylation data: ", length(participants))

message(paste(nrow(fam), "samples with genetic data matched to methylation data"))

message("Predicting DNAmAge")#############################################

ret <- meffonym.score(norm.beta, "hannum")
hannum_cor <- cor(ret$score, covar$Age_numeric)
message("The correlation between predicted DNAm age and actual age is ", cor(ret$score, covar$Age_numeric))
message("The r squared between predicted DNAm age and actual age is ", cor(ret$score, covar$Age_numeric)^2)

# warning if the correlation is lower than 0.5-0.6
if(hannum_cor<0.6){
  message("The correlation between predicted DNAm age and actual age is ", cor(ret$score, covar$Age_numeric),
          "; this seems a little bit low. Is the variance in age low in your dataset? The age range is",range(covar$Age_numeric))
}
message("There were ",length(ret$sites)," DNAm sites used out of a possible 71 for prediction of Hannum age")

covar$Hannum_age <- ret$score

hannum_plot <- ggplot(data=covar, aes(x=Age_numeric,y=Hannum_age)) +
            geom_point(alpha=0.5,colour="#440154FF") +
            labs(title=paste0(cohort_name," Hannum age vs age, cor =",cor(covar$Age_numeric,ret$score),"p =",signif(cor.test(covar$Age_numeric,ret$score)$p.value),digits=2))+
            geom_smooth(method='lm')+
            theme_minimal()


message("Predicting smoking")#############################################

ret <- meffonym.score(norm.beta, "mcigarette")
covar$mcigarette <- ret$score

# in case some cohorts only have smoking quantity,
# recode quantity into a factor if they don't have categories (0 = never, 1+= yes)

if ("Smoking_factor" %in% colnames(covar)) {
  message("Smoking_factor variable already exists.")
} else if ("Smoking_numeric" %in% colnames(covar)) {
  # Create Smoking_factor based on Smoking_numeric
  covar$Smoking_factor <- ifelse(covar$Smoking_numeric == 0, "No", "Yes")
  message("Smoking_factor variable created based on Smoking_numeric.")
} else {
  message("Smoking_numeric column not found. Cannot create Smoking_factor. Please check your phenotype data")
}

# add in if statement; plot if they have smoking variable
# if they don't have a smoking variable, plot the distribution of predicted smoking

if ("Smoking_factor" %in% colnames(df)) {
  message("Smoking_factor found in phenotype dataframe: running mcigarette plot of reported smoking vs DNAm-estimated smoking")
  mcigarette_plot <- ggplot(data=covar, aes(x=Smoking_factor,y=mcigarette,color=Smoking_factor)) +
    geom_boxplot()+
    scale_colour_viridis(discrete=T,begin=0,end=0.65)+
    labs(title=paste0(cohort_name," reported vs predicted smoking, p=",signif(wilcox.test(pheno_temp[,i] ~ pheno_temp[,j], data = pheno_temp, alternative = c("two.sided"))$p.value),digits=2))+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))+
    theme(legend.position = "none")
  
} else {
  # Run this command if 'Smoking_factor' column does not exist
  print("Smoking_factor not found in phenotype dataframe (if this is a measure you should have, please check the name of the smoking variable)
        :running mcigarette plot of DNAm-estimated smoking")
  mcigarette_plot <- ggplot() +
    geom_density(data=covar, aes(x=covar$mcigarette), colour="#1F968BFF")+
    labs(title=paste0(i,":n of NAs=",sum(is.na(covar$mcigarette))))+#,color="Legend")+
    geom_vline(xintercept = mean(covar$mcigarette))+
    theme_minimal()
  
}

# print out the two plots

jpeg(filename = paste0(prediction_plot,"_",study_name,".jpg"),width = 10, height = 5, units = "in", res = 600)
plot.out <- ggarrange(hannum_plot, mcigarette_plot,
                      labels = c("A","B"),
                      ncol = 2, nrow = 1)
print(plot.out)

