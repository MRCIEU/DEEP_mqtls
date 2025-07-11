# generate methylation PCs
# QC phase 1
# checking for PC distributions and associations
# QUESTION - how far do we want to test PC associations? batch? phenotypes?

# Args
arguments <- commandArgs(T)

beta_file <- arguments[1]
cov_file <- arguments[2] # [winsorized_covariates_file from check_phenotypes]
out_file <- arguments[4]
meth_array <- arguments[x]
#age_plot <-arguments[5]
# QUESTION - do we need cell counts here?
#cellcount_file <- arguments[9] 
study_name <- arguments[x]
PC1PC2_plot <- arguments[x]
PC3PC4_plot <- arguments[x]

suppressPackageStartupMessages(library(meffil))

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

message("Generating PCs")#######################################

# Create PC matrix to test PC associations (probe.range means take that number of 
# the most variable probes)
pcs <- meffil.methylation.pcs(meth,probe.range=50000)
message("there are",ncol(pcs),"PCs generated")
# QUESTION - do we want to do anything with more than the first 10 PCs
# (like an elbow plot or something)
# (I think more than 10 will be very impractical for individual PC plots)
pcs <- pcs[,1:10]
identical(rownames(pcs),rownames(covs))
pcs <- merge(x=pcs,y=covs, by.x="row.names", by.y="IID")
rownames(pcs) <- pcs$Row.names
pcs <- pcs[,-1]

test_pc_vars <- c("Age_numeric","Sex_factor") # add the vars we want to test the PCs against. Unlikely to be all
plot_pc1pc2_list <- vector("list", length = length(test_pc_vars))
names(plot_pc1pc2_list) <- test_pc_vars
plot_pc3pc4_list <- vector("list", length = length(test_pc_vars))
names(plot_pc3pc4_list) <- test_pc_vars

for(i in test_pc_vars){
  if(is.numeric(pcs[,i])){
    pc1pc2_plot <- ggplot(pcs, aes_string(x=PC1, y=PC2, color=pcs[,i])) +
      geom_point(size=1) + 
      scale_colour_viridis() +
      labs(title=paste0("PC1 vs PC2, ",i))+
      theme_bw() +
      theme(legend.position="none") 
    pc3pc4_plot <- ggplot(pcs, aes_string(x=PC3, y=PC4, color=pcs[,i])) +
      geom_point(size=1) + 
      scale_colour_viridis() +
      labs(title=paste0("PC1 vs PC2, ",i))+
      theme_bw() +
      theme(legend.position="none") 
    
  } else {
    pc1pc2_plot <- ggplot(pcs, aes_string(x=PC1, y=PC2, color=pcs[,i])) +
      geom_point(size=1) + 
      scale_colour_viridis(discrete = T) +
      labs(title=paste0("PC3 vs PC4, ",i))+
      theme_bw() #+
      #theme(legend.position="none")
    pc3pc4_plot <- ggplot(pcs, aes_string(x=PC3, y=PC4, color=pcs[,i])) +
      geom_point(size=1) + 
      scale_colour_viridis(discrete = T) +
      labs(title=paste0("PC3 vs PC4, ",i))+
      theme_bw() #+
      #theme(legend.position="none") 
    
  }
  plot_pc1pc2_list[[i]] <- pc1pc2_plot
  plot_pc3pc4_list[[i]] <- pc3pc4_plot
  
}

n_plot_rows <- ceiling(test_pc_vars/4)
row_dimensions <- n_plot_rows*4

jpeg(filename = paste0(PC1PC2_plot,"_",study_name,".jpg"),width = 12, height = row_dimensions, units = "in", res = 600)
makeplots <- ggarrange(plotlist=plot_pc1pc2_list, ncol = 4, nrow = n_plot_rows)
annotate_figure(makeplots, top = text_grob(paste0(study_name,"; PC1 vs PC2 plots"), 
                                           color = "black", face = "bold", size = 14))
print(makeplots)
dev.off()

jpeg(filename = paste0(PC3PC4_plot,"_",study_name,".jpg"),width = 12, height = row_dimensions, units = "in", res = 600)
makeplots <- ggarrange(plotlist=plot_pc3pc4_list, ncol = 4, nrow = n_plot_rows)
annotate_figure(makeplots, top = text_grob(paste0(study_name,"; PC3 vs PC4 plots"), 
                                           color = "black", face = "bold", size = 14))
print(makeplots)
dev.off()

# add lm test in here too? maybe bar plot with bars for effect size and indicator of significance?