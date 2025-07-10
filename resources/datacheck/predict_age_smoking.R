# Calculate Hannum age epigenetic clock and smoking status using mcigarette
# This is a check of how well our age phenotype correlates with age as predicted by the Hannum clock
# We have chosen the Hannum clock over others as it uses relatively few cpgs; 
# takes relatively little time to calculate; 
# and is a first generation clock that tries to predict chronological age.

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


message("Predicting DNAmAge")#############################################

ret <- meffonym.score(norm.beta, "hannum")
hannum_cor <- cor(ret$score, pheno_df$Age_numeric)
message("The correlation between predicted DNAm age and actual age is ", cor(ret$score, pheno_df$Age_numeric))
message("The r squared between predicted DNAm age and actual age is ", cor(ret$score, pheno_df$Age_numeric)^2)

# warning if the correlation is lower than 0.5-0.6
if(hannum_cor<0.6){
  message("The correlation between predicted DNAm age and actual age is ", cor(ret$score, pheno_df$Age_numeric),
          "; this seems a little bit low. Is the variance in age low in your dataset? The age range is",range(pheno_df$Age_numeric))
}
message("There were ",length(ret$sites)," DNAm sites used out of a possible 71 for prediction of Hannum age")

pheno_df$Hannum_age <- ret$score

hannum_plot <- ggplot(data=pheno_df, aes(x=Age_numeric,y=Hannum_age)) +
            geom_point(alpha=0.5,colour="#440154FF") +
            labs(title=paste0(cohort_name," Hannum age vs age, cor =",cor(pheno_df$Age_numeric,ret$score),"p =",signif(cor.test(pheno_df$Age_numeric,ret$score)$p.value),digits=2))+
            geom_smooth(method='lm')+
            theme_minimal()


message("Predicting smoking")#############################################

ret <- meffonym.score(norm.beta, "mcigarette")
pheno_df$mcigarette <- ret$score

# QUESTION - will all smoking be categorical? will some cohorts only have pack years?
# if so, we will need to add an extra if statement to accommodate this

# add in if statement; plot if they have smoking variable
# if they don't have a smoking variable, plot the distribution of predicted smoking

if ("Smoking_factor" %in% colnames(df)) {
  message("Smoking_factor found in phenotype dataframe: running mcigarette plot of reported smoking vs DNAm-estimated smoking")
  mcigarette_plot <- ggplot(data=pheno_df, aes(x=Smoking_factor,y=mcigarette,color=Smoking_factor)) +
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
    geom_density(data=pheno_df, aes(x=pheno_df$mcigarette), colour="#1F968BFF")+
    labs(title=paste0(i,":n of NAs=",sum(is.na(pheno_df$mcigarette))))+#,color="Legend")+
    geom_vline(xintercept = mean(pheno_df$mcigarette))+
    theme_minimal()
  
}

# print out the two plots

jpeg(filename = paste0(prediction_plot,"_",study_name,".jpg"),width = 10, height = 5, units = "in", res = 600)
plot.out <- ggarrange(hannum_plot, mcigarette_plot,
                      labels = c("A","B"),
                      ncol = 2, nrow = 1)
print(plot.out)

