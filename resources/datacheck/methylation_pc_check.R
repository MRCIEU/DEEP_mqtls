# generate methylation PCs
# QC phase 1
# checking for PC distributions and associations

# To do - add specificity to genetic PCs

# Args
arguments <- commandArgs(T)

beta_file <- arguments[1]
cov_file <- arguments[2] # [winsorized_covariates_file from check_phenotypes]
out_file <- arguments[3]
cellcount_file <- arguments[4] 
study_name <- arguments[5]
study_specific_vars <- arguments[6]
genetic_pc_file <- arguments[7]
scree_plot <- arguments[8]
PC1PC2_plot <- arguments[9]
PC3PC4_plot <- arguments[10]
pc_var_association_plot <- arguments[11]

suppressPackageStartupMessages(library(meffil))

message("Reading in data and matching up samples across files")#######################################
load(cov_file)
load(beta_file)
load(cellcount_file)
# TO DO: what kind of file is genetic_pc_file?
load(genetic_pc_file)

participants <- as.character(intersect(colnames(norm.beta),covar$IID))
covar <- covar[covar$IID%in%participants,]
norm.beta <- norm.beta[,participants]

# TO DO: merge genetic PCs 1:10 with covs file
  # need to know format of genetic PC file

celltypes <- colnames(cellcount_file) # replace this - what is the name of the cell counts df??

message("Generating PCs")#######################################

# Create PC matrix to test PC associations (probe.range means take that number of 
# the most variable probes)
pcs <- meffil.methylation.pcs(meth,probe.range=50000)
message("there are",ncol(pcs),"PCs generated")

# make scree plot
pca.var <- pcs$sdev^2
pca.var.explained <- pca.var / sum(pca.var)

df <- data.frame(PC = 1:length(pca.var.explained),
                 Variance = pca.var.explained)

jpeg(filename = paste0(scree_plot,"_",study_name,".jpg"),width = 4, height = 5, units = "in", res = 600)
ggplot(df, aes(x = PC, y = Variance)) +
  geom_line() +
  geom_point() +
  labs(title = paste("Scree Plot",study_name), x = "Principal Component", y = "Proportion of Variance Explained") +
  theme_minimal()
dev.off()

# reduce to top 10 PCs because that's all we will test
pcs <- pcs[,1:10]
identical(rownames(pcs),rownames(covs))
pcs <- merge(x=pcs,y=covs, by.x="row.names", by.y="IID")
rownames(pcs) <- pcs$Row.names
pcs <- pcs[,-1]

# TO DO: finish adding the vars we want to test the PCs against. Unlikely to be all. 
test_pc_vars <- c("Age_numeric","Sex_factor",study_specific_vars, celltypes) 
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

# now test how much each variable associates with each PC using a linear model
pc_analysis <- list()
for(i in 1:10){
  print(i)
  batch_vars <- paste(batch)
  model_formula <- paste0("PC",i," ~ ",paste(test_pc_vars,sep = "+"))
  temp <- lm(formula=model_formula, 
             data = pcs)
  temp <- summary(temp)
  pc_analysis[[i]] <- temp
}
summary(pc_analysis)
names(pc_analysis) <- c(1:10)
#save(pc_analysis, file=paste0(""))

pc_plotlist <- list()
for(i in 1:10){
  temp <- as.data.frame(pc_analysis[[i]]$coefficients)
  names(temp) <- c("estimate","se","t","p")
  temp <- temp[temp$p < 0.05,]
  temp <- temp[order(temp$p),]
  row.names.remove <- c("(Intercept)")
  temp <- temp[!(row.names(temp) %in% row.names.remove), ]
  temp$PC <- paste0("PC",i)
  # make barplot of p values
  temp$var <- as.character(rownames(temp))
  
  # remove slides and rows and replace with a single value
  # the plot is awful if you leave them all in
  batchvars <- as.character(temp$var)
  # TO DO: edit the below batch names:
  slidevars <- grep("^slide|^row",batchvars, value = T)
  
  temp_slide_row <- temp[temp$var %in% slidevars & temp$p < 0.05, ]
  # Calculate mean estimate for slide variables
  slide_effects <- temp_slide_row[grep("^slide", temp_slide_row$var), "estimate"]
  mean_slide_effect <- mean(slide_effects, na.rm = TRUE)
  mean_slide_p <- mean(slide_subset$p, na.rm = TRUE)
  # Calculate mean estimate for row variables
  row_effects <- temp_slide_row[grep("^row", temp_slide_row$var), "estimate"]
  mean_row_effect <- mean(row_effects, na.rm = TRUE)
  mean_row_p <- mean(row_subset$p, na.rm = TRUE)
  
  temp <- temp[!temp$var %in% slidevars,]
  temp <- temp[complete.cases(temp),]
  
  summary_slide <- data.frame(estimate = mean_slide_effect, se=NA, t=NA, p = mean_slide_p, var = "mean_slide")
  summary_row <- data.frame(estimate = mean_row_effect, se=NA, t=NA, p = mean_row_p, var = "mean_row")
  # Bind summary rows back to temp
  temp <- rbind(temp, summary_slide, summary_row)
  
  pc_plot <- ggplot(temp, aes(x=var, y=estimate, fill=var)) +
    #scale_fill_viridis(discrete = T) +  
    geom_bar(stat = "identity", alpha = 0.5) +
    scale_fill_manual(values = c("p < 0.05" = "#238A8DFF", "p ≥ 0.05" = "#FDE725FF")) +
    geom_text(aes(label = signif(-log10(p), 3)), vjust = -0.3, size = 3) +
    labs(title=paste0("PC",i),x="Variable", y="Estimate") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    ggtitle(i) 
  pc_plotlist[[i]] <- pc_plot
}

jpeg(filename = paste0(pc_var_association_plot,"_",study_name,".jpg"),width = 15, height = 20, units = "in", res = 600)
makeplots <- ggarrange(plotlist=plot_list, ncol = 3, nrow = 4)
annotate_figure(makeplots, top = text_grob(paste0(study_name,"; methylation PC associations"), 
                                           color = "black", face = "bold", size = 14))
print(makeplots)
dev.off()

