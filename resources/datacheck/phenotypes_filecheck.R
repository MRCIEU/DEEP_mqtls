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
pheno <- read.table(phenotype_file,header=T,stringsAsFactors = F, colClasses = c(Sex_factor = "character", Slide_factor = "character", Row_factor = "character", Plate_factor = "character"))
cov1 <- dim(pheno)[1]
cov2 <- dim(pheno)[2]

# Check for duplicate IIDs in the phenotype file
if (any(duplicated(pheno$IID))) {
  duplicate_ids <- unique(pheno$IID[duplicated(pheno$IID)])
  msg <- paste0("Duplicate IIDs found in phenotype file: ",
                paste(head(duplicate_ids, 5), collapse = ", "),
                if(length(duplicate_ids) > 5) "..." else "")
  errorlist <- c(errorlist, msg)
  warning("ERROR: ", msg)

}

meth_ids <- scan(meth_ids_file, what="character")
phenotypes <- read.csv(phenotype_names,header = F)
phenotypes <- as.character(phenotypes[,1])

participants <- as.character(intersect(meth_ids,pheno$IID))
pheno <- pheno[pheno$IID%in%participants,]

# set all the numeric columns to numeric variables, and all the factors to factors
# we've set stringsAsFactors = F so all the cols should currently be characters

for (col_name in colnames(pheno)[colnames(pheno) != "IID"]) {
  if (grepl("_factor$", col_name)) {
    if (length(table(na.omit(pheno[[col_name]]))) == nrow(pheno)) {
      msg <- paste0(col_name, " is specified as a factor but has the same number of levels as individuals")
      errorlist <- c(errorlist, msg)
      warning("ERROR: ", msg)
    }
    pheno[[col_name]] <- as.factor(pheno[[col_name]])
  } else if (grepl("_numeric$", col_name)) {
    unique_vals <- length(unique(na.omit(pheno[[col_name]])))
    pct <- round(100 * unique_vals / max(1, nrow(pheno)), 1)
    if (unique_vals < 0.01 * nrow(pheno)) {
      msg <- paste0(col_name, " is specified as numeric but shows low variability: ",
                    unique_vals, " unique non-missing values (", pct, "% of samples). ",
                    "This can be expected for integer or coded categorical variables. ",
                    "Please confirm the variable type; if it is categorical, consider renaming the header to '",
                    sub("_numeric$", "_factor", col_name), "'")
      errorlist <- c(errorlist, msg)
      warning("WARN: ", msg)
    }
    pheno[[col_name]] <- as.numeric(pheno[[col_name]])
  } else {
    msg <- paste0(col_name, " does not have _factor or _numeric suffix in header")
    warning("WARN: ", msg)
    warninglist <- c(warninglist, msg)
  }
}

w <- which(names(pheno)[1] %in% c("IID"))
if(w!=1)
{
  msg <- paste0("first column from phenotype file should be the sample identifier with the name IID")
  errorlist <- c(errorlist, msg)
  warning("ERROR: ", msg)
}

g1<-grep("_factor",names(pheno))
g2<-grep("_numeric",names(pheno))
g<-unique(c(g1,g2))

if(length(g)!=(cov2-1)){
  msg <- paste0("have you specified whether your phenotypes are factors or numeric in the header of the phenotypes file? Please make sure your column headers are e.g. 'IID Sex_factor Age_numeric' etc")
  errorlist <- c(errorlist, msg)
  warning("ERROR: ", msg)
}


#for (i in 1:length(g1))
#{
#  if(length(table(na.omit(pheno[,g1[i]])))==(cov1))
#  {
#    msg <- paste0(g1[i], " is specified as a factor but has the same number of levels as individuals")
#    errorlist <- c(errorlist, msg)
#    warning("ERROR: ", msg)
#  }
#}

a <- apply(pheno,2,function(x) y<-length(which(is.na(x))))
if(length(which(a>0.1*length(meth_ids))))
{
  msg <- paste0("more than 10% of missingness in the following phenotypes:\n", 
                paste(names(pheno)[a>0.1*length(meth_ids)], collapse="\n")
  )
  errorlist <- c(errorlist, msg)
  warning("ERROR: ", msg)
}

# Sex_factor: coerce existing column, report unexpected values / missing
if(! "Sex_factor" %in% names(pheno)) {
  msg <- "There is no Sex_factor variable in the phenotype file. Please provide M/F values."
  errorlist <- c(errorlist, msg); warning("ERROR: ", msg)
} else {
  pheno$Sex_factor <- as.character(pheno$Sex_factor)
  bad_sex_vals <- setdiff(unique(na.omit(pheno$Sex_factor)), c("M","F"))
  if (length(bad_sex_vals) > 0) {
    msg <- paste0("There are values in Sex_factor that are neither M nor F: ", paste(bad_sex_vals, collapse=", "))
    errorlist <- c(errorlist, msg); warning("ERROR: ", msg)
  }
  pheno$Sex_factor <- factor(pheno$Sex_factor, levels = c("M","F"))
  if (any(is.na(pheno$Sex_factor))) {
    msg <- "There are missing values in Sex_factor. Please ensure all individuals have M or F."
    errorlist <- c(errorlist, msg); warning("ERROR: ", msg)
  }
}

# Age_numeric: coerce to numeric safely and then check missing/negative/implausible
if(! "Age_numeric" %in% names(pheno)) {
  msg <- "There is no Age_numeric variable in the phenotype file."
  errorlist <- c(errorlist, msg); warning("ERROR: ", msg)
} else {
  # coerce factor/character to numeric safely
  pheno$Age_numeric <- as.numeric(as.character(pheno$Age_numeric))
  # debug helpers (optional): message("NA in Age_numeric: ", sum(is.na(pheno$Age_numeric)))
  if (any(is.na(pheno$Age_numeric))) {
    msg <- "Some individuals don't have valid numeric ages. Please check Age_numeric values."
    errorlist <- c(errorlist, msg); warning("ERROR: ", msg)
  }
  if (any(pheno$Age_numeric < 0, na.rm = TRUE)) {
    msg <- "Some negative values in the age column."
    errorlist <- c(errorlist, msg); warning("ERROR: ", msg)
  }
  if (mean(pheno$Age_numeric, na.rm = TRUE) > 100) {
    msg <- "Average age is above 100, please make sure age is provided in years."
    errorlist <- c(errorlist, msg); warning("ERROR: ", msg)
  }
}

# check study specific vars

n_overlap <- study_specific_vars%in%colnames(pheno)
no_overlap <- setdiff(study_specific_vars, colnames(pheno))

message("There were ", sum(n_overlap), " of ", length(study_specific_vars),
                      " study specific variables that were specified in the config file found in the phenotype file.")

if (length(no_overlap) > 0) {
  msg_missing <- paste0("The following study-specific variables were missing from the phenotype file: ",
                        paste(no_overlap, collapse = ", "),
                        ". Please add these to the config file or check spelling.")
  errorlist <- c(errorlist, msg_missing)
  warning("ERROR: ", msg_missing)
} else {
   message("All study-specific variables were found in the phenotype file.")
}

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

if (length(no_overlap) > 0) {
  msg <- paste0("There were ",length(no_overlap)," variables found in the phenotype file that do not match the specified names.
              \nThese were: ",no_overlap,"\n Please check the spelling matches column names specified in the wiki (phenotype data section).")
  errorlist <- c(errorlist, msg)
  warning("ERROR: ", msg)
} else {
  message("All phenotype names match those specified in the wiki (phenotype data section).")
}

message("\n\nCompleted checks\n")

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

