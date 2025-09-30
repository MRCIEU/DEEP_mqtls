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
  message("idat files folder: ", paste0(idat_input_folder))
  message("study name: ", study_name)
  message("Create samplesheet")
  samplesheet <- meffil.create.samplesheet(idat_input_folder, recursive=TRUE)

  qc.objects <- meffil.qc(samplesheet, verbose=TRUE)

  qc.objects.shrink = lapply(qc.objects, meffil.shrink.qc.object)

  # save original qc.objects to home directory
  save(qc.objects, file=paste0(home_path, "/processed_data/methylation_data/",study_name, ".qc_objects.rda"))
  save(qc.objects.shrink, file=paste0(home_path, "/results/01/", study_name, ".qc_objects_shrink.rda"))

  message("Create shrunk QC objects")
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
  qc.object.ori <- load_new_var(paste0(home_path, "/processed_data/methylation_data/", study_name, ".qc_objects.rda"))

  message("Loading normalized shrunk objects")
  file_path <- paste0(home_path, "/input_data/", study_name, ".norms.rda")
  if (!file.exists(file_path)) {
    message("File not found: ", file_path, ". Please place your norms.rda file in the input folder.")
    stop("Do not proceed Module 01f Step 2 until you have received the cross-cohort normalized quantiles from Bristol.")
  }
  dat.norms <- load_new_var(file_path)

  norm.objects <- mapply(meffil.expand.norm.object, dat.norms, qc.object.ori, SIMPLIFY=F)
  message("Information restored")

  meth = meffil.normalize.samples(norm.objects)

  save(meth, file = paste0(home_path, "/input_data/", study_name, "_fn_meth.RData"))

  message("Across cohorts normalized methylation data saved")
}