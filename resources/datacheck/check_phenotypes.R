# checking phenotypes
# this script will:
#       load in covariates, 
#       set them to the correct type of variable, 
#       save out variable statistics and print plots,
#       remove outliers if necessary and transform if necessary,
#       Then redo summary stats and plots

# To do:
# add in outlier removal and ? transformation
# meeting 3/7/25 - we won't do transformation
# it's hard to know what might be appropriate without seeing the data
# concern: if we transform the results will be less interpretable
# add in second pass of plots and summary stats

suppressPackageStartupMessages(library(DescTools))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(ggpubr))

args <- (commandArgs(TRUE));
covariates_file <- as.character(args[1]);
fam_file <- as.character(args[2]);
meth_ids_file <- as.character(args[3]);
sorted_methylation <- as.character(args[4]);
raw_phenotype_distribution_plot <- as.character(args[5]);
raw_phenotype_summary_file <- as.character(args[6])
# need to add in study name to the config file (should already be in the config file - ${study_name} in bash)
study_name <- arguments[7]
edited_phenotype_distribution_plot <- as.character(args[x]);
edited_phenotype_summary_file <- as.character(args[x])
# we'll save out the edited phenotype file as an Rdata file so we don't have to
# do anything with the variables the next time we load them in
phenotype_outfile <- as.character(args[x])


################

# 1. set up data ----

################

# question - what's in the covariates file?
# will we have all phenotypes, covariates, batch variables already in the file?
# that would be easiest; if not we'll load in all the files here and add them to one df
# edit 3/7/25 - yes it looks like they will already be in one file

message("Checking covariates file: ", covariates_file)
covar <- read.table(covariates_file,header=T,stringsAsFactors = F)

meth_ids <- scan(meth_ids_file, what="character")
fam <- read.table(fam_file, header=FALSE, stringsAsFactors=FALSE)

commonids_mgc <- Reduce(intersect, list(meth_ids, covar$IID, fam[,2]))
message("Number of samples with covariate, methylation and genetic data: ", length(commonids_mgc))

# set all the numeric columns to numeric variables, and all the factors to factors
# we've set stringsAsFactors = F so all the cols should currently be characters
for (col_name in colnames(covar)[!colnames(covar)==IID]) {
  if (grepl("_numeric", col_name)) {
    covar[[col_name]] <- as.numeric(covar[[col_name]])
  } else if (grepl("_factor", col_name)) {
    covar[[col_name]] <- as.factor(covar[[col_name]])
  }
}

################

# 2. plot raw phenotypes and output raw summary stats ----

################

# create list to save the plots to
phenotypes <- colnames(covar)[!colnames(covar)==IID]
plot_list <- vector("list", length = length(phenotypes))
names(plot_list) <- phenotypes
summstats_list <- vector("list", length = length(phenotypes))
names(summstats_list) <- phenotypes

# run loop to generate plots (density for numeric, bar for categorical)
for(i in phenotypes){
  if(is.numeric(covar[,i])){
    # add variable summary and N of NAs to list
    stats_out <- summary(covar[,i])
  } else {
    stats_out <- table(covar[,i])
  }
  summstats_list[[i]]$stats_out <- stats_out
  summstats_list[[i]]$nas_out <- sum(is.na(covar[,i]))
  # remove missing cases to avoid missingness on plot
  pheno.temp <- covar[!is.na(covar[,i]),]
  print(dim(pheno.temp))
  participants.temp <- as.character(pheno.temp$IID)
  
  # run plots (density for numeric variables and bar for categorical)
  if(is.numeric(pheno.temp[,i])){
    test <- ggplot() +
      geom_density(data=pheno.temp, aes_string(x=pheno.temp[,i]), colour="#1F968BFF")+
      labs(title=paste0(i,", total N = ",nrow(covar),":n of NAs=",sum(is.na(covar[,i]))),x=i)+#,color="Legend")+
      geom_vline(xintercept = mean(pheno.temp[,i]))+
      theme_minimal()
    
  } else if(is.factor(pheno.temp[,i])){
    
    test <- ggplot(data=pheno.temp, aes_string(x=pheno.temp[,i],fill=pheno.temp[,i])) +
      geom_bar()+
      scale_fill_viridis(discrete=T,begin=0,end=0.65)+
      labs(title=paste0(i,", total N = ",nrow(covar),",NAs = ",sum(is.na(covar[,i]))),x=i)+
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

# *add in the cohort name into file name*
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

# add in summ stats output
# numeric vars only

numeric_phenos <- grepl("_numeric", colnames(covar))
plot_list <- vector("list", length = length(numeric_phenos))
names(plot_list) <- numeric_phenos
summstats_list <- vector("list", length = length(numeric_phenos))
names(summstats_list) <- numeric_phenos

for(i in numeric_phenos){
  if(is.numeric(covar[,i])){
    # remove outliers
    # currently this leaves NAs?
    # do we want that or do we want to replace them with the mean?
    # meeting 3/7/25 - we will winsorise as that will maintain the extreme values
    # BUT we need to look at how that affects distribution (ie does it make it bimodal)
    # outlier <- which(data$trait<(mean(data$trait,na.rm=T)-SD*sd(data$trait,na.rm=T)) | data$trait> (mean(data$trait,na.rm=T)+SD*sd(data$trait,na.rm=T)))
    outlier_rm <- Winsorize(covar[,i], val = quantile(covar[,i], probs = c(0.05, 0.95), na.rm = T))
    diff_count <- sum(covar[,i] != outlier_rm & !is.na(covar[,i]) & !is.na(outlier_rm))
    message("There are",length(diff_count),"outliers in the variable",i,"that have been replaced by the 5%-quantile")
    if (length(diff_count)>0){covar[,i] <- outlier_rm}
    # now re-plot and re-do summary stats
    test <- ggplot() +
      geom_density(data=pheno.temp, aes_string(x=pheno.temp[,i]), colour="#1F968BFF")+
      labs(title=paste0(i,", total N = ",nrow(covar),":n of NAs=",sum(is.na(covar[,i]))),x=i)+#,color="Legend")+
      geom_vline(xintercept = mean(pheno.temp[,i]))+
      theme_minimal()
    plot_list[[i]] <- test
    stats_out <- summary(covar[,i])
    summstats_list[[i]] <- stats_out
    
  } else {
    message(i,"is a categorical variable. Please check for small cell counts, as small cell counts may require some categories to be collapsed")
  }
}

################

# 4. plot cleaned phenotypes and output cleaned summary stats ----

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

