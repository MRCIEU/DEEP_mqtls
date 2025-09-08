# https://github.com/perishky/meffil/wiki/Functional-normalizing-separate-datasets
# https://github.com/perishky/meffil/blob/master/R/shrink.r

# Remove any information not absolutely needed for normalizing quantiles

suppressPackageStartupMessages(library(meffil))

arguments <- commandArgs(T);
betas <- arguments[1];
qc_input_file <- arguments[2];
study_name <- arguments[3];
output_path <- arguments[4];
stage <- arguments[5];

# qc stage 1 shrink ======= remove unnecessary information from qc objects

if (stage == "shrink" && betas == 0) {
required_var <- c("controls", "quantiles")

vars_before <- ls()
load(qc_input_file)
vars_after <- ls()
new_vars <- setdiff(vars_after, vars_before)
message("Loaded variable(s) from qc.object: ", paste(new_vars, collapse = ", "))

qc_objects_shrink <- meffil.shrink.qc.object(get(new_vars[1]), keep.vars = required_var)

save(qc_objects_shrink, file=paste0(output_path, "/",study_name, "qc_objects_shrink.rda"))

message("Shrunk QC objects created")

}

# DEEP receives shrunk data from cohorts
# to do normalization between cohorts by dev team
# qc.objects = c(cohort1.objects, cohort2.objects...)
# norm.objects = meffil.normalize.quantiles(qc.objects, number.pcs=...)
# cohort1.norms = norm.objects[names(cohort1.objects)] ...
# send back shrunk objects after normalization in DEEP
# format: study__name + "objects.rda"

# qc stage 2.  ========= restore information

if (stage == "expand" && betas != 0) {

    vars_before <- ls()
    load(qc_input_file)
    load(paste0(output_path, "/",study_name, "objects.rda"))
    vars_after <- ls()
    new_vars <- setdiff(vars_after, vars_before)

    qc.object.ori = get(new_vars[1])
    dat.norms = get(new_vars[2])

    norm.objects = meffil.expand.norm.object(dat.norms, qc.object.ori)
    
    save(norm.objects, file=paste0(output_path, "/processed_data/",study_name, "_harmonized_meth.rda"))

    message("Information restored")

}