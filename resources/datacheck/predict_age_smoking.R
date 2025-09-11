# Calculate Hannum age epigenetic clock and smoking status using mcigarette
# This is a check of how well our age phenotype correlates with age as predicted by the Hannum clock
# We have chosen the Hannum clock over others as it uses relatively few cpgs; 
# takes relatively little time to calculate; 
# and is a first generation clock that tries to predict chronological age.

# QUESTION - check that numerical smoking vars will indeed be called Smoking_numeric

# add in section that tabulates predicted and reported sex. Maybe also peruvian

# Args
arguments <- commandArgs(T)

beta_file <- arguments[1]
pheno_file <- arguments[2] # [winsorized_covariates_file from check_phenotypes]
fam_file <- arguments[3]
study_name <- arguments[4]
prediction_plot <- arguments[5]
updated_pheno_file <- arguments[6] 
smoking_prediction_output_file <- arguments[7] 

suppressPackageStartupMessages(library(meffonym))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(ggpubr))

message("Reading in data and matching up samples across files")#######################################
load(pheno_file)
load(beta_file)
participants <- as.character(intersect(colnames(norm.beta),pheno$IID))
pheno <- pheno[pheno$IID%in%participants,]
norm.beta <- norm.beta[,participants]
message("Number of samples with covariate and methylation data: ", length(participants))

message(paste(nrow(fam_file), "samples with genetic data matched to methylation data"))

message("Predicting DNAmAge")#############################################

ret <- meffonym.score(norm.beta, "hannum")
hannum_cor <- cor(ret$score, pheno$Age_numeric)
message("The correlation between predicted DNAm age and actual age is ", cor(ret$score, pheno$Age_numeric))
message("The r squared between predicted DNAm age and actual age is ", cor(ret$score, pheno$Age_numeric)^2)

# warning if the correlation is lower than 0.5-0.6
if(hannum_cor<0.6){
  message("The correlation between predicted DNAm age and actual age is ", cor(ret$score, pheno$Age_numeric),
          "; this seems a bit low. Is the variance in age low in your dataset? The age range is ",range(pheno$Age_numeric))
}
message("There were ",length(ret$sites)," DNAm sites used out of a possible 71 for prediction of Hannum age")

pheno$Hannum_age <- ret$score

hannum_plot <- ggplot(data=pheno, aes(x=Age_numeric,y=Hannum_age)) +
            geom_point(alpha=0.5,colour="#440154FF") +
            labs(title=paste0(study_name," Hannum age vs age, cor =",signif(cor(pheno$Age_numeric,ret$score),digits = 2),"p =",signif(cor.test(pheno$Age_numeric,ret$score)$p.value),digits=2))+
            geom_smooth(method='lm')+
            theme_minimal()

message("Recoding smoking to yes/no")#############################################

# in case some studies only have smoking quantity,
# recode quantity into a factor if they don't have categories (0 = never, 1+= yes)

if ("Smoking_factor" %in% colnames(pheno)) {
  message("Smoking_factor variable already exists.")
} else if ("Smoking_numeric" %in% colnames(pheno)) {
  # Create Smoking_factor based on Smoking_numeric
  pheno$Smoking_factor <- ifelse(pheno$Smoking_numeric == 0, "No", "Yes")
  message("Smoking_factor variable created based on Smoking_numeric.")
} else {
  message("Smoking_numeric column not found. Cannot create Smoking_factor. Please check your phenotype data")
}

message("Predicting mcigarette smoking score")#############################################

ret <- meffonym.score(norm.beta, "mcigarette")
pheno$p_smoking_mcigarette <- ret$score

# add in if statement; plot if they have smoking variable
# if they don't have a smoking variable, plot the distribution of predicted smoking

if ("Smoking_factor" %in% colnames(df)) {
  message("Smoking_factor found in phenotype dataframe: running mcigarette plot of reported smoking vs DNAm-estimated smoking")
  mcigarette_plot <- ggplot(data=pheno, aes(x=Smoking_factor,y=p_smoking_mcigarette,color=Smoking_factor)) +
    geom_boxplot()+
    scale_colour_viridis(discrete=T,begin=0,end=0.65)+
    labs(title=paste0(study_name," reported vs predicted smoking (mcigarette), p=",signif(wilcox.test(p_smoking_mcigarette ~ Smoking_factor, data = pheno, alternative = c("two.sided"))$p.value),digits=2))+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))+
    theme(legend.position = "none")
  
} else {
  # Run this command if 'Smoking_factor' column does not exist
  print("Smoking_factor not found in phenotype dataframe (if this is a measure you should have, please check the name of the smoking variable)
        :running mcigarette plot of DNAm-estimated smoking")
  mcigarette_plot <- ggplot() +
    geom_density(data=pheno, aes(x=pheno$p_smoking_mcigarette), colour="#1F968BFF")+
    labs(title=paste0(study_name," predicted smoking (mcigarette)"))+#,color="Legend")+
    geom_vline(xintercept = mean(pheno$p_smoking_mcigarette))+
    theme_minimal()
  
}


message("Predicting Elliott smoking score")#############################################

ret <- meffonym.score(norm.beta, "elliott-smoking")
pheno$p_smoking_elliott <- ret$score

# add in if statement; plot if they have smoking variable
# if they don't have a smoking variable, plot the distribution of predicted smoking

if ("Smoking_factor" %in% colnames(df)) {
  message("Smoking_factor found in phenotype dataframe: running plot of reported smoking vs Elliott DNAm-estimated smoking")
  elliott_plot <- ggplot(data=pheno, aes(x=Smoking_factor,y=p_smoking_elliott,color=Smoking_factor)) +
    geom_boxplot()+
    scale_colour_viridis(discrete=T,begin=0,end=0.65)+
    labs(title=paste0(study_name," reported vs predicted smoking (Elliott), p=",signif(wilcox.test(p_smoking_elliott ~ Smoking_factor, data = pheno, alternative = c("two.sided"))$p.value),digits=2))+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))+
    theme(legend.position = "none")
  
} else {
  # Run this command if 'Smoking_factor' column does not exist
  print("Smoking_factor not found in phenotype dataframe (if this is a measure you should have, please check the name of the smoking variable)
        :running Elliott plot of DNAm-estimated smoking")
  Elliott_plot <- ggplot() +
    geom_density(data=pheno, aes(x=pheno$p_smoking_elliott), colour="#1F968BFF")+
    labs(title=paste0(study_name," predicted smoking (Elliott)"))+#,color="Legend")+
    geom_vline(xintercept = mean(pheno$p_smoking_elliott))+
    theme_minimal()
  
}

# print out the three plots

jpeg(filename = paste0(prediction_plot,"_",study_name,".jpg"),width = 10, height = 5, units = "in", res = 600)
plot.out <- ggarrange(hannum_plot, mcigarette_plot, Elliott_plot,
                      labels = c("A","B","C"),
                      ncol = 2, nrow = 2)
print(plot.out)

# save out:
# as phenotype file
save(pheno,file=paste0(updated_pheno_file))
# as separate file for IID and predicted smoking (txt file IID and two smoking scores)

smoking_prediction_vars <- c("IID","p_smoking_mcigarette","p_smoking_elliott")
smoking_prediction <- pheno[,smoking_prediction_vars]
write.table(pheno,file=paste0(smoking_prediction_output_file))
# predicted_smoking is name of output file (txt)