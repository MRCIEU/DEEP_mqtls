suppressPackageStartupMessages(library(meffil))

# on DEEP server, we will receive .RData files for control probe and then do normalization
# We will need the list of .RData files to process (list_RData)
# list_RData <- c("file1.RData", "file2.RData")
# cohort_names <- c("cohort1", "cohort2")
# qc_names_list <- lapply(seq_along(list_RData), function(i) c(list_RData[i], cohort_names[i]))
# names(qc_names_list) <- paste0("names", seq_along(list_RData))

# qc.objects <- list()
# for (rdata_file in list_RData) {
#     message("Processing ", rdata_file)
#     before <- ls()
#     load(rdata_file)
#     after <- ls()
#     # would have issue if variable names are same between any cohorts
#     new_vars <- setdiff(after, before)
#     print(new_vars)
#     for (v in new_vars) {
#         qc.objects[[v]] <- get(v)
#     }
# }
# # Before normalization, we estimate the correct number of control probe principal components to include in the normalization.
# y <- meffil.plot.pc.fit(qc.objects)
# ggsave(y$plot,filename="pc-fit-combined.pdf",height=6,width=6)

# norm.objects <- meffil.normalize.quantiles(qc.objects, number.pcs=10)

# for (i in seq_along(qc_names_list)) {
#     cohort_obj_names <- names(qc.objects[[i]]) 
#     norm_sub <- norm.objects[names(norm.objects) %in% cohort_obj_names]
#     cohort <- qc_names_list[[i]][2]
#     save(norm_sub, file = paste0("norm_objects_", cohort, ".rda"))
# }

# Then on the server of each DEEP user

main <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    cohort_name <- args[1]
    out_file <- args[2]

    load(paste0("norm_objects_", cohort_name, ".rda"))
    var = ls()
    beta_aft <- meffil.normalize.samples(var)
    save(beta_aft, file = out_file)
    message("Functional normalization between datasets [Done]")
}

main()
