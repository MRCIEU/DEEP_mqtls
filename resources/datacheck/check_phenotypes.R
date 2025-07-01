# checking phenotypes
# this script will:
#       load in covariates, 
#       set them to the correct type of variable, 
#       save out variable statistics and print plots,
#       remove outliers if necessary and transform if necessary,
#       Then redo summary stats and plots

# To do:
# add in outlier removal and ? transformation
# concern: if we transform the results will be less interpretable
# add in second pass of plots and summary stats

args <- (commandArgs(TRUE));
covariates_file <- as.character(args[1]);
fam_file <- as.character(args[2]);
meth_ids_file <- as.character(args[3]);
sorted_methylation <- as.character(args[4]);
raw_phenotype_distribution_plot <- as.character(args[5]);
raw_phenotype_summary_file <- as.character(args[6])
# need to add in cohort name to the config file
cohort_name <- arguments[x]

################

# 1. set up data ----

################

message("Checking covariates file: ", covariates_file)
covar <- read.table(covariates_file,header=T,stringsAsFactors = F)

meth_ids <- scan(meth_ids_file, what="character")
fam <- read.table(fam_file, header=FALSE, stringsAsFactors=FALSE)

commonids_mgc <- Reduce(intersect, list(meth_ids, covar$IID, fam[,2]))
message("Number of samples with covariate, methylation and genetic data: ", length(commonids_mgc))

# set all the numeric columns to numeric variables, and all the factors to factors
# we've set stringsAsFactors = F so all the cols should currently be characters
for (col_name in names(covar)[!names(covar)==IID]) {
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
phenotypes <- names(covar)[!names(covar)==IID]
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
  print("remove missing cases:")
  pheno.temp <- covar[!is.na(covar[,i]),]
  print(dim(pheno.temp))
  participants.temp <- as.character(pheno.temp$IID)
  
  if(is.numeric(pheno.temp[,i])){
    test <- ggplot() +
      geom_density(data=pheno.temp, aes_string(x=pheno.temp[,i]), colour="#1F968BFF")+
      labs(title=paste0(i,":n of NAs=",sum(is.na(covar[,i]))),x=i)+#,color="Legend")+
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
# let's start with 4 rows and we'll see how we go
# size will be something like 4xn of rows and 4x n of columns
# add title that has cohort name and raw distribution

n_plot_rows <- ceiling(nrow(covar)/4)
row_dimensions <- n_plot_rows*4

jpeg(filename = paste0(raw_phenotype_distribution_plot,".jpg"),width = 12, height = row_dimensions, units = "in", res = 600)
makeplots <- ggarrange(plotlist=plot_list, ncol = 4, nrow = n_plot_rows)
print(makeplots)
dev.off()

# and save out the summary stats
# QUESTION - do we like this rdata format? It's a list with two elements 
# (variable summary and number of NAs) for each phenotype
# we could alternatively save out as csvs
# the only issue is that the summaries of numeric and categorical variables
# have differing dimensions (tables vs summary) so they will be tricky to put in the same csv
save(summstats_list,file=raw_phenotype_summary_file)

################

# 3. remove outliers and ? transform data ----

################

################

# 4. plot cleaned phenotypes and output cleaned summary stats ----

################
