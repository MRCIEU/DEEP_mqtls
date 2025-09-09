# checking phenotypes
# this script will:
#       load in covariates, 
#       set them to the correct type of variable, 
#       save out variable statistics and print plots,
#       remove outliers via winsorization,
#       Then redo summary stats and plots

# Outstanding questions:
  # Do we want a raw version of the phenotypes (not winsorized)
  # There are some variables we won't want to winsorize, so how do we do this?

suppressPackageStartupMessages(library(DescTools))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(ggpubr))

arguments <- commandArgs(T)
phenotypes_file <- as.character(arguments[1]);
fam_file <- as.character(arguments[2]);
meth_ids_file <- as.character(arguments[3]);
raw_phenotype_distribution_plot <- as.character(arguments[4]);
raw_phenotype_summary_file <- as.character(arguments[5])
study_name <- arguments[6] # ${study_name} in bash
edited_phenotype_distribution_plot <- as.character(arguments[7]);
edited_phenotype_summary_file <- as.character(arguments[8])
# we'll save out the edited phenotype file as an Rdata file so we don't have to
# do anything with the variables the next time we load them in
#phenotype_outfile <- as.character(arguments[9])
winsorized_phenotype_file <- as.character(arguments[9])


################

# 1. set up data ----

################

message("Checking covariates file: ", phenotypes_file)
pheno <- read.table(phenotypes_file,header=T,stringsAsFactors = F)

# We don't require all DNAm samples to have genetic data
# we'll print out what the differential is here
meth_ids <- scan(meth_ids_file, what="character")
fam <- read.table(fam_file, header=FALSE, stringsAsFactors=FALSE)
commonids_mgc <- Reduce(intersect, list(meth_ids, pheno$IID, fam[,2]))
message("Number of samples with covariate, methylation and genetic data: ", length(commonids_mgc))
participants <- as.character(intersect(meth_ids,pheno$IID))
pheno <- pheno[pheno$IID%in%participants,]
message("Number of samples with covariate and methylation data: ", length(participants))


# set all the numeric columns to numeric variables, and all the factors to factors
# we've set stringsAsFactors = F so all the cols should currently be characters
for (col_name in colnames(pheno)[!colnames(pheno)==IID]) {
  if (grepl("_numeric", col_name)) {
    pheno[[col_name]] <- as.numeric(pheno[[col_name]])
  } else if (grepl("_factor", col_name)) {
    pheno[[col_name]] <- as.factor(pheno[[col_name]])
  }
}

# QUESTION - save out this version of pheno here so we have a raw version?

################

# 2. plot raw phenotypes and output raw summary stats ----

################

# create list to save the plots to
phenotypes <- colnames(pheno)[!colnames(pheno)=="IID"]
plot_list <- vector("list", length = length(phenotypes))
names(plot_list) <- phenotypes
summstats_list <- vector("list", length = length(phenotypes))
names(summstats_list) <- phenotypes

# run loop to generate plots (density for numeric, bar for categorical)
for(i in phenotypes){
  if(is.numeric(pheno[,i])){
    # add variable summary and N of NAs to list
    stats_out <- summary(pheno[,i])
  } else {
    stats_out <- table(pheno[,i])
  }
  summstats_list[[i]]$stats_out <- stats_out
  summstats_list[[i]]$nas_out <- sum(is.na(pheno[,i]))
  
  # Warning if there is more than 10% missingness in a variable
  if(sum(is.na(pheno[,i]))>nrow(pheno)/10){
    message("Warning: there is over 10% missingness in",i)
  }
  # remove missing cases to avoid missingness on plot
  pheno.temp <- pheno[!is.na(pheno[,i]),]

  # run plots (density for numeric variables and bar for categorical)
  if(is.numeric(pheno.temp[,i])){
    test <- ggplot() +
      geom_density(data=pheno.temp, aes_string(x=pheno.temp[,i]), colour="#1F968BFF")+
      labs(title=paste0(i,", total N = ",nrow(pheno),":n of NAs=",sum(is.na(pheno[,i]))),x=i)+#,color="Legend")+
      geom_vline(xintercept = mean(pheno.temp[,i]))+
      theme_minimal()
    
  } else if(is.factor(pheno.temp[,i])){
    
    test <- ggplot(data=pheno.temp, aes_string(x=pheno.temp[,i],fill=pheno.temp[,i])) +
      geom_bar()+
      scale_fill_viridis(discrete=T,begin=0,end=0.65)+
      labs(title=paste0(i,", total N = ",nrow(pheno),",NAs = ",sum(is.na(pheno[,i]))),x=i)+
      theme_minimal()+
      theme(legend.position="none")+
      geom_text(stat='count', aes(label=..count..), color="black", vjust=-0.1)
  }
  plot_list[[i]] <- test
  
}

# now print off all the plots from the list so that all the study variables are on one plot
# how do we automate the plot size and the ncol and nrow?
# let's start with 4 columns and we'll see how we go
# size will be something like 4xn of rows and 4 columns
# add title that has cohort name and raw distribution

n_plot_rows <- ceiling(phenotypes/4)
row_dimensions <- n_plot_rows*4

jpeg(filename = paste0(raw_phenotype_distribution_plot,"_",study_name,".jpg"),width = 12, height = row_dimensions, units = "in", res = 600)
makeplots <- ggarrange(plotlist=plot_list, ncol = 4, nrow = n_plot_rows)
annotate_figure(makeplots, top = text_grob(paste0(study_name,"; raw phenotype distributions"), 
                                      color = "black", face = "bold", size = 14))
print(makeplots)
dev.off()

# and save out the summary stats
  # QUESTION - do we like this rdata format? It's a list with two elements 
  # (variable summary and number of NAs) for each phenotype
  # we could alternatively save out as csvs
  # the only issue is that the summaries of numeric and categorical variables
  # have differing dimensions (tables vs summary) so they will be tricky to put in the same csv
save(summstats_list,file=paste0(raw_phenotype_summary_file,"_",study_name,".Rdata"))

################

# 3. remove outliers via Winsorizing numeric data ----

################

# numeric vars only

  # QUESTION - we don't want to winsorize age (as long as values are biologically plausible)
  #         so need to add in some code that checks values are reasonable and removes age and 
  #         any other such vars (perhaps smoking pack years? others?) from winsorization
  #         Is this a thing we need to code manually once we have the data harmonisation questionnaire?


numeric_phenos <- grepl("_numeric", colnames(pheno))
# add in here removal of age etc from numeric_phenos
plot_list <- vector("list", length = length(numeric_phenos))
names(plot_list) <- numeric_phenos
summstats_list <- vector("list", length = length(numeric_phenos))
names(summstats_list) <- numeric_phenos

for(i in numeric_phenos){
  if(is.numeric(pheno[,i])){
    # remove outliers
    # meeting 3/7/25 - we will winsorise as that will maintain the extreme values
    # we need to look at how that affects distribution (ie does it make it bimodal)
    
    # winsorize data
    outlier_rm <- Winsorize(pheno[,i], val = quantile(pheno[,i], probs = c(0.05, 0.95), na.rm = T))
    # print out how many observations are Winsorized
    diff_count <- sum(pheno[,i] != outlier_rm & !is.na(pheno[,i]) & !is.na(outlier_rm))
    message("There are",length(diff_count),"outliers in the variable",i,"that have been replaced by the 5%-quantile")
    # change variable to the Winsorized version
      # QUESTION - do we want to keep the original version of the variable somehow? or will this just be the
      # covariate df that we load in at the start of the script? Or do we save out a version of 
      # raw covariates above once we've converted to factor/numeric?
    if (length(diff_count)>0){pheno[,i] <- outlier_rm}
    # now re-plot and re-do summary stats
    test <- ggplot() +
      geom_density(data=pheno, aes_string(x=pheno[,i]), colour="#1F968BFF")+
      labs(title=paste0(i,", total N = ",nrow(pheno),":n of NAs=",sum(is.na(pheno[,i]))),x=i)+#,color="Legend")+
      geom_vline(xintercept = mean(pheno[,i]))+
      theme_minimal()
    plot_list[[i]] <- test
    stats_out <- summary(pheno[,i])
    summstats_list[[i]] <- stats_out
    
  } else {
    message(i,"is a categorical variable. Please check for small cell counts, as small cell counts may require some categories to be collapsed")
  }
}

################

# 4. plot cleaned phenotypes and output cleaned summary stats and winsorized data ----

################

n_plot_rows <- ceiling(numeric_phenos/4)
row_dimensions <- n_plot_rows*4

jpeg(filename = paste0(edited_phenotype_distribution_plot,"_",study_name,".jpg"),width = 12, height = row_dimensions, units = "in", res = 600)
makeplots <- ggarrange(plotlist=plot_list, ncol = 4, nrow = n_plot_rows)
annotate_figure(makeplots, top = text_grob(paste0(study_name,"; Winsorized phenotype distributions"), 
                                           color = "black", face = "bold", size = 14))
print(makeplots)
dev.off()

save(summstats_list,file=paste0(edited_phenotype_summary_file,"_",study_name,".Rdata"))

save(pheno,file=paste0(winsorized_phenotype_file))

