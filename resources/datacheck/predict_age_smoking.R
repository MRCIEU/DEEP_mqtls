# Calculate Hannum age epigenetic clock and smoking status using mcigarette
# This is a check of how well our age phenotype correlates with age as predicted by the Hannum clock
# We have chosen the Hannum clock over others as it uses relatively few cpgs; 
# takes relatively little time to calculate; 
# and is a first generation clock that tries to predict chronological age.

# Args
arguments <- commandArgs(T)

beta_file <- arguments[1]
cov_file <- arguments[2]
fam_file <- arguments[3]
out_file <- arguments[4]
age_plot <-arguments[5]
SD <- as.numeric(arguments[6])
age_stats <- arguments[7]
#smoking_file <- arguments[8]. ## smoking should be in covariates file already??
#cellcount_file <- arguments[9]
# need to add in cohort name to the config file
cohort_name <- arguments[x]

suppressPackageStartupMessages(library(meffonym))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(ggpubr))

message("Reading in data and matching up samples across files")#######################################
covs <- read.table(cov_file, header=T)
fam <- read.table(fam_file)[,c(1,2)]
load(beta_file)
## according to original script we also load in a smoking file here;
## but why isn't smoking data in the covariates file?

# assumes that everyone in the genetic data has methylation & covariate data
covs <- covs[match(fam[,2], covs[,"IID"]),]
smoking <- read.table(smoking_file, header = T)[,c('IID', 'Smoking')]
smoking <- smoking[match(fam[,2], smoking[,"IID"]),]
load(beta_file)
norm.beta <- norm.beta[, match(fam[,2], colnames(norm.beta))]
message(paste(nrow(fam), "samples with genetic data matched to methylation data"))

# 1/7/25 - I'm not certain if most of the checking covariates is necessary. 
# it's likely this will get flagged in the initial variable check?
# and perhaps even we need to save out the edited phenotype file from the 
# new check_phenotypes script

message("Checking covariates")#######################################
# Check if the age and sex are valid for predicting biological age
age_index <- grep("^age_numeric$", names(covs), ignore.case=TRUE)
count_age <- length(unique(covs[, age_index]))
sex_index <- grep("^Sex_factor$", names(covs), ignore.case = TRUE)
count_sex <- length(unique(covs[, sex_index]))

age_valid = FALSE
sex_valid = FALSE
if(length(age_index) != 1 | length(sex_index) != 1){
  message("There should be only one column in the covariate file called 
         'Age_numeric' and 'Sex_factor'. Neither variable will be adjusted for.")
  valid_vector <- c('IID')
}else if (count_age == 1 & count_sex == 1){
  message("Age variable is constant, so will not be adjusted for.")
  message("Sex variable only has 1 level, so will not be adjusted for.")
  valid_vector <- c('IID')
}else if (count_age == 1 & count_sex != 1){
  message("Age variable is constant, so will not be adjusted for.")
  message("Sex variable has ", count_sex, " levels.")
  sex_valid = TRUE
  valid_vector <- c('IID', 'Sex_factor')
}else if (count_age != 1 & count_sex == 1){
  message("Sex variable has one level, so will not be adjusted for.")
  age_valid = TRUE
  valid_vector <- c('IID', 'Age_numeric')
}else {
  age_valid = TRUE
  sex_valid = TRUE
  valid_vector <- c('IID', 'Sex_factor', 'Age_numeric')
}

if(length(age_index) != 1){
  names(covs)[age_index] <- "Age_numeric"
}  

if(length(sex_index) != 1){
  names(covs)[sex_index] <- "Sex_factor"
}

##
# not sure this section is necessary?
pheno_df <- subset(covs, select = valid_vector) 
pheno_df <- merge(pheno_df, smoking, by.x="IID", by.y="IID", all.x=TRUE)
if (sex_valid == T){
  cortable <- subset(pheno_df, select = -c(Sex_factor))
}else{
  cortable <- pheno_df
}

if (cellcount_file == 'NULL') {
  message("No predicted cell count matrix provided.")
} else {
  cellcount <- read.table(cellcount_file, header = T, stringsAsFactors=FALSE)
  cellcount <- subset(cellcount, select=-c(Treg))
  m <- match(fam[,2], cellcount[,"IID"])
  cellcount <- cellcount[m,]
  pheno_df <- merge(pheno_df, cellcount,  by.x = 'IID', by.y = 'IID', all.x = TRUE)
}
# what does this plot? - looks like cors between dnam age measurements and age?
pdf(paste0(age_plot, '.pdf'), width=12, height=12)
name_sumstats <- c()
sumstats <- c()
##



message("Predicting DNAmAge")#############################################

ret <- meffonym.score(norm.beta, "hannum")

message("The correlation between predicted DNAm age and actual age is ", cor(ret$score, pheno_df$Age_numeric))
# should we add in a warning if the correlation is lower than 0.5-0.6?
message("There were ",length(ret$sites)," DNAm sites used out of a possible 71 for prediction of Hannum age")

pheno_df$Hannum_age <- ret$score

jpeg(filename = paste0(hannum_age_plot,".jpg"),width = 5, height = 5, units = "in", res = 600)
hannum_plot <- ggplot(data=pheno_df, aes(x=Age_numeric,y=Hannum_age)) +
            geom_point(alpha=0.5,colour="#440154FF") +
            labs(title=paste0(cohort_name," Hannum age vs age, cor =",cor(pheno_df$Age_numeric,ret$score),"p =",cor.test(pheno_df$Age_numeric,ret$score)$p.value))+
            geom_smooth(method='lm')+
            theme_minimal()
dev.off()


message("Predicting smoking")#############################################

ret <- meffonym.score(norm.beta, "mcigarette")
pheno_df$mcigarette <- ret$score

jpeg(filename = paste0(mcigarette_plot,".jpg"),width = 5, height = 5, units = "in", res = 600)
ggplot(data=pheno_df, aes(x=smoking_factor,y=mcigarette,color=mcigarette)) +
  geom_boxplot()+
  scale_colour_viridis(discrete=T,begin=0,end=0.65)+
  labs(title=paste0(cohort_name," Hannum age vs age,",x="",y=""))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))+
  theme(legend.position = "none")
dev.off()


