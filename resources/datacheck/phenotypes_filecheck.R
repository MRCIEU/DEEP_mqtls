errorlist <- list()
warninglist <- list()

library(data.table)
suppressMessages(library(matrixStats))

args <- (commandArgs(TRUE));
phenotype_file <- as.character(args[1]);
meth_ids_file <- as.character(args[2]);
sorted_methylation <- as.character(args[3]);
study_specific_vars <- strsplit(args[4], " ")[[1]] # these will be added in the config file - batch vars and study specific factors
phenotype_names <- as.character(args[5])

message("Checking phenotypes file: ", phenotype_file)
pheno <- read.table(phenotype_file,header=T)
cov1 <- dim(pheno)[1]
cov2 <- dim(pheno)[2]

meth_ids <- scan(meth_ids_file, what="character")
phenotypes <- read.csv(phenotype_names,header = F)
phenotypes <- as.character(phenotypes[,1])

#fam <- read.table(fam_file, header=FALSE, stringsAsFactors=FALSE)

#commonids_mgc <- Reduce(intersect, list(meth_ids, pheno$IID, fam[,2]))
#message("Number of samples with phenotype, methylation and genetic data: ", length(commonids_mgc))

#if(sorted_methylation == "no"){
#  if(length(commonids_mgc) < 100)
#  {
#    msg <- paste0("must have at least 100 individuals with phenotype, methylation and genetic data.")
#    errorlist <- c(errorlist, msg)
#    warning("ERROR: ", msg)
#  }
#}

w <- which(names(pheno)[1] %in% c("IID"))
if(w!=1)
{
  msg <- paste0("first column from phenotype file should be the sample identifier with the name IID")
  errorlist <- c(errorlist, msg)
  warning("ERROR: ", msg)
}

#if(cov2<3)
#{
#  msg <- paste0("are there any phenotypes missing in the phenotypes file? Sex and Age are required")
#  errorlist <- c(errorlist, msg)
#  warning("ERROR: ", msg)
#}

g1<-grep("_factor",names(pheno))
g2<-grep("_numeric",names(pheno))
g<-unique(c(g1,g2))

if(length(g)!=(cov2-1)){
  msg <- paste0("have you specified whether your phenotypes are factors or numeric in the header of the phenotypes file? Please make sure your column headers are e.g. 'IID Sex_factor Age_numeric' etc")
  errorlist <- c(errorlist, msg)
  warning("ERROR: ", msg)
}


for (i in 1:length(g1))
{
  if(length(table(na.omit(pheno[,g1[i]])))==(cov1))
  {
    msg <- paste0(g1[i], " is specified as a factor but has the same number of levels as individuals")
    errorlist <- c(errorlist, msg)
    warning("ERROR: ", msg)
  }
}

a <- apply(pheno,2,function(x) y<-length(which(is.na(x))))
if(length(which(a>0.1*length(meth_ids))))
{
  msg <- paste0("more than 10% of missingness in the following phenotypes:\n", 
                paste(names(pheno)[a>0.1*length(meth_ids)], collapse="\n")
  )
  errorlist <- c(errorlist, msg)
  warning("ERROR: ", msg)
}

if(! "Sex_factor" %in% names(pheno))
{
  msg <- paste0("There is no Sex_factor variable in the phenotype file. Please provide M/F values, even if they are all the same sex.")
  errorlist <- c(errorlist, msg)
  warning("ERROR: ", msg)
}else{
  pheno <- read.table(phenotype_file,header=T,colClasses=c('Sex_factor'='factor'))
}

if(any(is.na(pheno$Sex_factor)))
{
  msg <- paste0("There are some missing values in the Sex_factor column. Please make sure all individuals have data for this column.")
  errorlist <- c(errorlist, msg)
  warning("ERROR: ", msg)
}

index <- pheno$Sex_factor %in% c("M", "F")
if(any(!index))
{
  msg <- paste0("There are some values in the Sex_factor column that are neither M nor F. Please make sure all individuals have data for this column.")
  errorlist <- c(errorlist, msg)
  warning("ERROR: ", msg)
}

if(! "Age_numeric" %in% names(pheno))
{
  msg <- paste0("There is no Age_numeric variable in the phenotype file. Please provide age in years, even if they are all the same age.")
  errorlist <- c(errorlist, msg)
  warning("ERROR: ", msg)
}

pdf(age_distribution_plot, height=6, width=6)
hist(pheno$Age_numeric, breaks=50, xlab="Age", main=paste("age distribution (N=", length(which(!is.na(pheno$Age_numeric))),")",sep=""),cex.main=0.7)
null <- dev.off()


if(any(is.na(pheno$Age_numeric)))
{
  msg <- paste0("Some individuals don't have ages. Please make sure there are no missing values.")
  errorlist <- c(errorlist, msg)
  warning("ERROR: ", msg)
}

if(any(pheno$Age_numeric < 0))
{
  msg <- paste0("Some negative values in the age column.")
  errorlist <- c(errorlist, msg)
  warning("ERROR: ", msg)
}

if(mean(pheno$Age_numeric, na.rm=T) > 100)
{
  msg <- paste0("Average age is above 100, please make sure age is provided in years.")
  errorlist <- c(errorlist, msg)
  warning("ERROR: ", msg)
}

# check study specific vars

n_overlap <- study_specific_vars%in%colnames(pheno)
no_overlap <- setdiff(study_specific_vars, colnames(pheno))
msg <- paste0("There were ",length(n_overlap)," of ",length(study_specific_vars)," study specific variables that were specified in the config file found in the phenotype file.
              \nIf any were missing, these were: ",no_overlap,"\n please go back and add these to the config file.\n
              If you think you have included these variables in the config file, please check the spelling matches column names.")
errorlist <- c(errorlist, msg)
warning("ERROR: ", msg)

# check naming of phenotypes
 # check infection and pollution separately as they will have unknown prefixes
prefix_vars <- c("_infection_categorical","_pollution_numeric")
any_prefix_vars <- grep("(_infection_categorical|_pollution_numeric)$", 
                        colnames(pheno), 
                      value = TRUE)
phenotypes <- phenotypes[!phenotypes%in%prefix_vars]
phenotypes <- c(phenotypes,"IID")
no_overlap <- setdiff(colnames(pheno), phenotypes)
no_overlap <- no_overlap[!no_overlap%in%any_prefix_vars]
no_overlap <- no_overlap[!no_overlap%in%study_specific_vars]

msg <- paste0("There were ",no_overlap," variables found in the phenotype file that do not match the specified names.
              \nThese were: ",no_overlap,"\n Please check the spelling matches column names specified in the wiki (phenotype data section).")
errorlist <- c(errorlist, msg)
warning("ERROR: ", msg)


message("\n\nCompleted checks\n")

message("Summary of data:\n")
for(i in 1:length(cohort_summary))
{
  a <- cohort_summary[[i]]
  if(is.numeric(a)) a <- round(a, 2)
  message(names(cohort_summary)[i], ": ", paste(a, collapse=", "))
}


if(length(warninglist) > 0)
{
  message("\n\nPlease take note of the following warnings, and fix and re-run the data check if you see fit:")
  null <- sapply(warninglist, function(x)
  {
    message("- ", x)
  })
}


if(length(errorlist) > 0)
{
  message("\n\nThe following errors were encountered, and must be addressed before continuing:")
  null <- sapply(errorlist, function(x)
  {
    message("- ", x)
  })
  q(status=1)
}
message("\n\n")

