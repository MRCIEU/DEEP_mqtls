# https://github.com/perishky/meffil/wiki/Functional-normalizing-separate-datasets
# https://github.com/perishky/meffil/blob/master/R/shrink.r

# Remove any information not absolutely needed for normalizing quantiles

suppressPackageStartupMessages(library(meffil))

arguments <- commandArgs(T);
idat_input_folder <- arguments[1];
study_name <- arguments[2];
home_path <- arguments[3];
stage <- arguments[4];

# qc stage 1 shrink ======= remove unnecessary information from qc objects

if (stage == "shrink") {

    vars_before <- ls()
    load(idat_input_folder)
    vars_after <- ls()
    new_vars <- setdiff(vars_after, vars_before)
    message("idat files folder: ", paste(new_vars, collapse = ", "))

    samplesheet <- meffil.create.samplesheet(idat_input_folder, recursive=TRUE)

    qc.objects <- meffil.qc(samplesheet, verbose=TRUE)

    qc_objects_shrink <- meffil.shrink.qc.object(get(new_vars[1]))

    # save original qc.objects to home directory
    save(qc.objects, file=paste0(home_path, "/processed_data/methylation_data/",study_name, "qc_objects.rda"))
    save(qc_objects_shrink, file=paste0(home_path, "/results/01/", study_name, "qc_objects_shrink.rda"))

message("Shrunk QC objects created")

}

# DEEP receives shrunk data from cohorts
# to do normalization between cohorts by dev team
# qc.objects = c(cohort1.objects, cohort2.objects...)
# norm.objects = meffil.normalize.quantiles(qc.objects, number.pcs=...)
# cohort1.norms = norm.objects[names(cohort1.objects)] ...
# send back shrunk objects after normalization in DEEP
# format: study__name + "_objects.rda"

# qc stage 2.  ========= restore information

load_new_var <- function(file) {
  vars_before <- ls()
  load(file)
  vars_after <- ls()
  new_vars <- setdiff(vars_after, vars_before)
  if (length(new_vars) == 0) stop("No new variable loaded from ", file)
  get(new_vars[1])
}

if (stage == "expand") {
  message("Loading original qc.objects")
  qc.object.ori <- load_new_var(paste0(home_path, "/processed_data/methylation_data/", study_name, "qc_objects.rda"))
  
  message("Loading normalized shrunk objects")
  message("Please put the file from DEEP server in the input folder")
  dat.norms <- load_new_var(paste0(home_path, "/input_data/", study_name, "_objects.rda"))
  
  norm.objects <- meffil.expand.norm.object(dat.norms, qc.object.ori)
  save(norm.objects, file = paste0(home_path, "/processed_data/", study_name, "_harmonized_meth.rda"))
  message("Information restored")
}