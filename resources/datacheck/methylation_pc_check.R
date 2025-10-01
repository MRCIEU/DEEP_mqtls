# generate methylation PCs
# QC phase 1
# checking for PC distributions and associations

# Args
arguments <- commandArgs(T)

beta_file <- arguments[1]
updated_pheno_file <- arguments[2] # [winsorized_covariates_file from check_phenotypes]
cellcounts_cov <- arguments[3]
study_name <- arguments[4]
genetic_pc_file <- arguments[5] # ${pcs_all} \
scree_plot <- arguments[6]
PC1PC2_plot <- arguments[7]
PC3PC4_plot <- arguments[8]
pc_var_association_plot <- arguments[9]
scripts_directory <- arguments[10]
study_specific_vars <- strsplit(arguments[11], " ")[[1]] # these will be added in the config file - batch vars and study specific factors

suppressPackageStartupMessages(library(meffil))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(ggpubr))
source(paste0(scripts_directory,"/resources/datacheck/fn_rm_constant_col.R"))

message("Reading in data and matching up samples across files")#######################################
load(updated_pheno_file)
load(beta_file)
cell_counts <- read.table(cellcounts_cov, header=T)
rownames(cell_counts) <- cell_counts$IID

# TO DO: what kind of file is genetic_pc_file?
genetic_pcs <- read.table(genetic_pc_file)[,-1]
colnames(genetic_pcs) <- c("IID", paste("genetic_pc", 1:(ncol(genetic_pcs)-1), sep=""))
rownames(genetic_pcs) <- genetic_pcs$IID

participants <- as.character(intersect(colnames(norm.beta),pheno$IID))
pheno <- pheno[pheno$IID%in%participants,]
norm.beta <- norm.beta[,participants]

# detect cell count panel prefixes
cell_count_cols <- setdiff(colnames(cell_counts), c("FID","IID"))
cellcount_panel_prefixes <- unique(sub("\\..*", "", cell_count_cols))
message("Detected cell count panel prefixes: ", paste(cellcount_panel_prefixes, collapse = ", "))

# del treg cell types if present
celltypes <- grep(paste0("^(", paste0(cellcount_panel_prefixes, collapse="|"), ")"), colnames(cell_counts), value = TRUE)
celltypes <- celltypes[!grepl("treg", celltypes, ignore.case = TRUE)]
cell_counts <- cell_counts[,c("IID",celltypes)]

# del columns with no variation
cell_counts <- remove_constant_cols(cell_counts, "cell_counts")

# del nRBC if mean age > 1
if(mean(pheno$Age_numeric) < 1){
  message("Keeping nRBC as mean age is less than 1")
  } else{
  message("Removing nRBC as mean age is greater than 1")
  celltypes <- celltypes[!grepl("nRBC", celltypes, ignore.case = TRUE)]
}

if ("Sex_factor" %in% colnames(pheno)) {
  n_sex <- length(unique(na.omit(pheno$Sex_factor)))
  if (n_sex < 2) {
    message("Sex_factor has no variation (only one value). Removing Sex_factor from pheno.")
    pheno$Sex_factor <- NULL
  }
}

if ("Age_numeric" %in% colnames(pheno)) {
  n_age <- length(unique(na.omit(pheno$Age_numeric)))
  if (n_age < 2) {
    message("Age_numeric has no variation (only one value). Removing Age_numeric from pheno.")
    pheno$Age_numeric <- NULL
  }
}

for (cellcount_panel in cellcount_panel_prefixes) {
  message("Running methylation PCs for cell count panel: ", cellcount_panel)

  cell_counts_colnames <- grep(paste0("^", cellcount_panel, "\\."), colnames(cell_counts), value = TRUE)
  cellcounts_temp <- cell_counts[,c("IID",cell_counts_colnames)]
  pheno_panel <- pheno
  pheno_panel <- merge(pheno_panel,cellcounts_temp,by="IID")

  # TO DO: merge genetic PCs 1:10 with pheno file
  # need to know format of genetic PC file
  pheno_panel <- merge(pheno_panel, genetic_pcs, by="IID")

  message("Generating PCs")#######################################

  # Create PC matrix to test PC associations (probe.range means take that number of 
  # the most variable probes)
  pcs <- meffil.methylation.pcs(norm.beta,probe.range=50000,full.obj=T)
  # message("there are ",ncol(pcs)," PCs generated")

  # make scree plot
  pca.var <- pcs$sdev^2
  pca.var.explained <- pca.var / sum(pca.var)

  df <- data.frame(PC = 1:length(pca.var.explained), Variance = pca.var.explained)

  # plot a max of 50 PCs so we can read the plot easily
  if(length(pca.var.explained) > 50){
    message("More than 50 PCs generated, scree plot will show first 50 PCs only")
    df <- df[1:50,]
    pdf(file = paste0(scree_plot,"_",study_name,"_",cellcount_panel,".pdf"), width = 4, height = 5)
    p <- ggplot(df, aes(x = PC, y = Variance)) +
      geom_line() +
      geom_point() +
      labs(title = paste("Scree Plot",study_name,"; ",length(pca.var.explained),"  total PCs"), x = "Principal Component", y = "Proportion of Variance Explained") +
      theme_minimal()
    print(p)
    dev.off()
  
  } else {
    message("50 or fewer PCs generated, scree plot will show all PCs")
    pdf(file = paste0(scree_plot,"_",study_name,"_",cellcount_panel,".pdf"), width = 4, height = 5)
    p <- ggplot(df, aes(x = PC, y = Variance)) +
      geom_line() +
      geom_point() +
      labs(title = paste("Scree Plot",study_name), x = "Principal Component", y = "Proportion of Variance Explained") +
      theme_minimal()
    print(p)
    dev.off()
  }

  # pcs <- meffil.methylation.pcs(norm.beta,probe.range=50000,full.obj=F)
  pcs <- as.data.frame(pcs$x)
  message("there are ", ncol(pcs), " PCs generated")

  # reduce to top 10 PCs because that's all we will test
  pcs <- pcs[,1:10]
  identical(rownames(pcs),rownames(pheno_panel))
  pcs <- merge(x=pcs,y=pheno_panel, by.x="row.names", by.y="IID")
  # rownames(pcs) <- pcs$Row.names
  # pcs <- pcs[,-1]

  # TO DO: finish adding the vars we want to test the PCs against. Unlikely to be all. 
  test_pc_vars <- c("Age_numeric", "Sex_factor", "Population_group_factor", study_specific_vars, celltypes, colnames(genetic_pcs)[2:11])
  test_pc_vars <- test_pc_vars[test_pc_vars %in% colnames(pcs)]
  if ("Slide_factor" %in% names(pcs)) {
  pcs$Slide_factor <- as.factor(pcs$Slide_factor)
  }
  message(paste("Test PC vars are:", paste(test_pc_vars, collapse = ", ")))
  plot_pc1pc2_list <- vector("list", length = length(test_pc_vars))
  names(plot_pc1pc2_list) <- test_pc_vars
  plot_pc3pc4_list <- vector("list", length = length(test_pc_vars))
  names(plot_pc3pc4_list) <- test_pc_vars
  
  for(i in test_pc_vars){
    if(is.numeric(pcs[,i])){
      message(paste0("Plotting PC1 vs PC2 and PC3 vs PC4 coloured by continuous variable: ",i))
      pc1pc2_plot <- ggplot(pcs, aes(x=PC1, y=PC2, color=.data[[i]])) +
        geom_point(size=1) + 
        scale_colour_viridis() +
        labs(title=paste0("PC1 vs PC2, ",i))+
        theme_bw() +
        theme(legend.position="none") 
      pc3pc4_plot <- ggplot(pcs, aes(x=PC3, y=PC4, color=.data[[i]])) +
        geom_point(size=1) + 
        scale_colour_viridis() +
        labs(title=paste0("PC3 vs PC4, ",i))+
        theme_bw() +
        theme(legend.position="none") 
    
    } else {
        message(paste0("Plotting PC1 vs PC2 and PC3 vs PC4 coloured by categorical variable: ", i))
        n_levels <- length(unique(na.omit(pcs[[i]])))

        pc1pc2_plot <- ggplot(pcs, aes(PC1, PC2, color = .data[[i]])) +
          geom_point(size = 1) +
          scale_colour_viridis(discrete = TRUE) +
          labs(title = paste0("PC1 vs PC2, ", i)) +
          theme_bw()
        if (n_levels > 10) {
          pc1pc2_plot <- pc1pc2_plot + theme(legend.position = "none")
        }

        pc3pc4_plot <- ggplot(pcs, aes(PC3, PC4, color = .data[[i]])) +
          geom_point(size = 1) +
          scale_colour_viridis(discrete = TRUE) +
          labs(title = paste0("PC3 vs PC4, ", i)) +
          theme_bw()
        if (n_levels > 10) {
          pc3pc4_plot <- pc3pc4_plot + theme(legend.position = "none")
        }}
      plot_pc1pc2_list[[i]] <- pc1pc2_plot
      plot_pc3pc4_list[[i]] <- pc3pc4_plot
  }

  n_plot_rows <- ceiling(length(test_pc_vars)/4)
  row_dimensions <- n_plot_rows*3

  pdf(file = paste0(PC1PC2_plot,"_",study_name,"_",cellcount_panel,".pdf"),width = 12, height = row_dimensions)
  makeplots <- ggarrange(plotlist=plot_pc1pc2_list, ncol = 4, nrow = n_plot_rows)
  annotate_figure(makeplots, top = text_grob(paste0(study_name,"; PC1 vs PC2 plots"), 
                                           color = "black", face = "bold", size = 14))
  print(makeplots)
  dev.off()

  pdf(file = paste0(PC3PC4_plot,"_",study_name,"_",cellcount_panel,".pdf"),width = 12, height = row_dimensions)
  makeplots <- ggarrange(plotlist=plot_pc3pc4_list, ncol = 4, nrow = n_plot_rows)
  annotate_figure(makeplots, top = text_grob(paste0(study_name,"; PC3 vs PC4 plots"), 
                                           color = "black", face = "bold", size = 14))
  print(makeplots)
  dev.off()

  # now test how much each variable associates with each PC using a linear model
  pc_analysis <- list()
  for(i in 1:10){
    model_formula <- paste0("PC",i," ~ ",paste(test_pc_vars, collapse = "+"))
    print(paste0("Linear model formula: ", model_formula))
    fit <- lm(formula = model_formula, data = pcs)
    pc_analysis[[i]] <- summary(fit)
  }
  #save(pc_analysis, file=paste0(""))
  
  pc_plotlist <- list()
  for(i in 1:10){
    temp <- as.data.frame(pc_analysis[[i]]$coefficients)
    names(temp) <- c("estimate", "se", "t", "p")
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
    slidevars <- grep("^slide", batchvars, value = TRUE, ignore.case = TRUE)
    rowvars   <- grep("^row",   batchvars, value = TRUE, ignore.case = TRUE)
    platevars <- grep("^plate", batchvars, value = TRUE, ignore.case = TRUE)

    summary_list <- list()

    if (length(slidevars) > 0) {
      temp_slide <- temp[temp$var %in% slidevars & temp$p < 0.05, ]
      mean_slide_effect <- mean(temp_slide$estimate, na.rm = TRUE)
      mean_slide_p <- mean(temp_slide$p, na.rm = TRUE)
      summary_list[["mean_slide"]] <- data.frame(estimate = mean_slide_effect, se=NA, t=NA, p = mean_slide_p, PC = paste0("PC", i), var = "mean_slide")
    }

    if (length(rowvars) > 0) {
      temp_row <- temp[temp$var %in% rowvars & temp$p < 0.05, ]
      mean_row_effect <- mean(temp_row$estimate, na.rm = TRUE)
      mean_row_p <- mean(temp_row$p, na.rm = TRUE)
      summary_list[["mean_row"]] <- data.frame(estimate = mean_row_effect, se=NA, t=NA, p = mean_row_p, PC = paste0("PC", i), var = "mean_row")
    }

    if (length(platevars) > 0) {
      temp_plate <- temp[temp$var %in% platevars & temp$p < 0.05, ]
      mean_plate_effect <- mean(temp_plate$estimate, na.rm = TRUE)
      mean_plate_p <- mean(temp_plate$p, na.rm = TRUE)
      summary_list[["mean_plate"]] <- data.frame(estimate = mean_plate_effect, se=NA, t=NA, p = mean_plate_p, PC = paste0("PC", i), var = "mean_plate")
    }

    all_batchvars <- c(slidevars, rowvars, platevars)
    temp <- temp[!temp$var %in% all_batchvars, ]
    temp <- temp[complete.cases(temp), ]

    if (length(summary_list) > 0) {
    summary_df <- do.call(rbind, summary_list)
      if (nrow(temp) == 0) {
      temp <- summary_df
      } else {
      summary_df <- summary_df[, names(temp), drop = FALSE]
      temp <- rbind(temp, summary_df)
      }
    }

    temp$ci_lower <- ifelse(is.na(temp$se), temp$estimate, temp$estimate - 1.96 * temp$se)
    temp$ci_upper <- ifelse(is.na(temp$se), temp$estimate, temp$estimate + 1.96 * temp$se)

pc_plot <- ggplot(temp, aes(x = estimate, y = var)) +
  geom_point(size = 3) +
  geom_errorbarh(
    data = subset(temp, !is.na(ci_lower) & !is.na(ci_upper) & ci_lower != ci_upper),
    aes(xmin = ci_lower, xmax = ci_upper),
    height = 0.2
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = paste0("PC", i), x = "Estimate (95% CI)", y = "Variable") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))

    pc_plotlist[[i]] <- pc_plot
    }

  pdf(file = paste0(pc_var_association_plot,"_",study_name,"_",cellcount_panel,".pdf"),width = 15, height = 20)
  makeplots <- ggarrange(plotlist=pc_plotlist, ncol = 3, nrow = 4)
  annotate_figure(makeplots, top = text_grob(paste0(study_name,"; methylation PC associations"), 
                                           color = "black", face = "bold", size = 14))
  print(makeplots)
  dev.off()
}