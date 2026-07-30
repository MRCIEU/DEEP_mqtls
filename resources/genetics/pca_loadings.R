# Shared helpers for exporting PCA SNP loadings in a consistent format.

validate_pca_loadings <- function(loadings, n_pcs, expected_snps = NULL) {
  metadata_columns <- c(
    "CHR", "POS", "SNP", "LOADING_ALLELE", "OTHER_ALLELE"
  )
  pc_columns <- paste0("PC", seq_len(n_pcs))
  expected_columns <- c(metadata_columns, pc_columns)

  if (!identical(colnames(loadings), expected_columns)) {
    stop(
      "Unexpected SNP-loading columns. Expected: ",
      paste(expected_columns, collapse = ", ")
    )
  }
  if (nrow(loadings) == 0L) {
    stop("The SNP-loading table is empty.")
  }
  if (anyNA(loadings[, metadata_columns, drop = FALSE])) {
    stop("Missing variant metadata in the SNP-loading table.")
  }
  if (any(loadings$SNP == "") ||
      any(loadings$LOADING_ALLELE == "") ||
      any(loadings$OTHER_ALLELE == "")) {
    stop("Empty SNP or allele value in the SNP-loading table.")
  }
  if (anyDuplicated(loadings$SNP)) {
    stop("Duplicate SNP IDs in the SNP-loading table.")
  }

  loading_matrix <- as.matrix(loadings[, pc_columns, drop = FALSE])
  storage.mode(loading_matrix) <- "double"
  if (any(!is.finite(loading_matrix))) {
    stop("Non-finite values in the SNP-loading matrix.")
  }

  if (!is.null(expected_snps)) {
    expected_snps <- as.character(expected_snps)
    if (anyDuplicated(expected_snps)) {
      stop("Duplicate SNP IDs in the expected PCA SNP set.")
    }
    if (!setequal(loadings$SNP, expected_snps)) {
      missing_snps <- setdiff(expected_snps, loadings$SNP)
      unexpected_snps <- setdiff(loadings$SNP, expected_snps)
      stop(
        "SNP-loading IDs do not match the PCA SNP set (",
        length(missing_snps), " missing; ",
        length(unexpected_snps), " unexpected)."
      )
    }
  }

  invisible(TRUE)
}

write_pca_loadings <- function(loadings,
                               outfile,
                               n_pcs,
                               expected_snps = NULL) {
  validate_pca_loadings(loadings, n_pcs, expected_snps)

  output_dir <- dirname(outfile)
  if (!dir.exists(output_dir)) {
    stop("Output directory does not exist: ", output_dir)
  }

  temporary_file <- tempfile(
    pattern = paste0(".", basename(outfile), "."),
    tmpdir = output_dir
  )
  connection <- gzfile(temporary_file, open = "wt")
  write_error <- NULL
  tryCatch(
    utils::write.table(
      loadings,
      file = connection,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE
    ),
    error = function(e) {
      write_error <<- e
    },
    finally = {
      close(connection)
    }
  )
  if (!is.null(write_error)) {
    unlink(temporary_file)
    stop(write_error)
  }

  if (!file.rename(temporary_file, outfile)) {
    copied <- file.copy(temporary_file, outfile, overwrite = TRUE)
    unlink(temporary_file)
    if (!copied) {
      stop("Could not move the SNP-loading table to: ", outfile)
    }
  }

  message(
    "Saved ", nrow(loadings), " SNP loadings for ", n_pcs,
    " PCs to ", outfile
  )
  invisible(outfile)
}

save_pcair_loadings <- function(gds_file,
                                mypcair,
                                pca_snp_ids,
                                n_pcs,
                                nthreads,
                                outfile) {
  if (length(mypcair$unrels) < 2L) {
    stop("At least two unrelated samples are required for SNP loadings.")
  }

  n_pcs <- as.integer(n_pcs)
  nthreads <- as.integer(nthreads)
  if (n_pcs < 1L || n_pcs > ncol(mypcair$vectors)) {
    stop("Requested PC count is incompatible with mypcair$vectors.")
  }
  # SNPRelate permits either integer or character values in the GDS `snp.id`
  # node. DEEP commonly uses chr:pos_alleles variant IDs, so preserve the
  # original storage type when passing IDs back to snpgdsPCA(). Coercing these
  # IDs to integer would turn them into NA.
  if (length(pca_snp_ids) == 0L || anyNA(pca_snp_ids) ||
      anyDuplicated(pca_snp_ids)) {
    stop("The PC-AiR SNP set must contain unique, non-missing GDS SNP IDs.")
  }

  gds <- SNPRelate::snpgdsOpen(gds_file, allow.duplicate = TRUE)
  on.exit(SNPRelate::snpgdsClose(gds), add = TRUE)

  message(
    "Reconstructing the PCA on ", length(mypcair$unrels),
    " unrelated samples to export SNP loadings."
  )
  pca_unrelated <- SNPRelate::snpgdsPCA(
    gds,
    sample.id = mypcair$unrels,
    snp.id = pca_snp_ids,
    autosome.only = TRUE,
    eigen.cnt = n_pcs,
    num.thread = nthreads,
    verbose = FALSE
  )
  snp_loading_object <- SNPRelate::snpgdsPCASNPLoading(
    pca_unrelated,
    gdsobj = gds,
    num.thread = nthreads,
    verbose = FALSE
  )
  if (!identical(
    as.character(snp_loading_object$snp.id),
    as.character(pca_unrelated$snp.id)
  )) {
    stop("SNP order changed while calculating PC-AiR SNP loadings.")
  }

  pcair_row_index <- match(
    pca_unrelated$sample.id,
    rownames(mypcair$vectors)
  )
  if (anyNA(pcair_row_index)) {
    stop("Could not align reconstructed PCA samples with mypcair$vectors.")
  }
  pcair_unrelated <- mypcair$vectors[
    pcair_row_index,
    seq_len(n_pcs),
    drop = FALSE
  ]

  pc_correlations <- vapply(
    seq_len(n_pcs),
    function(pc) {
      stats::cor(
        pca_unrelated$eigenvect[, pc],
        pcair_unrelated[, pc],
        use = "complete.obs"
      )
    },
    numeric(1)
  )
  if (any(!is.finite(pc_correlations))) {
    stop("Non-finite PC correlation while validating reconstructed loadings.")
  }

  failed_pcs <- which(abs(pc_correlations) < 0.99)
  if (length(failed_pcs) > 0L) {
    stop(
      "Reconstructed PCs do not match mypcair$vectors for: ",
      paste0(
        "PC", failed_pcs, " (r=",
        format(round(pc_correlations[failed_pcs], 4), nsmall = 4), ")",
        collapse = ", "
      )
    )
  }
  warning_pcs <- which(abs(pc_correlations) < 0.999)
  if (length(warning_pcs) > 0L) {
    warning(
      "Reconstructed PC correlation is below 0.999 for: ",
      paste0(
        "PC", warning_pcs, " (r=",
        format(round(pc_correlations[warning_pcs], 4), nsmall = 4), ")",
        collapse = ", "
      )
    )
  }

  # SNPRelate stores this matrix as PCs x variants; output one row per SNP.
  loading_matrix <- t(snp_loading_object$snploading[
    seq_len(n_pcs),
    ,
    drop = FALSE
  ])
  # Eigenvector signs are arbitrary. Match the signs used in pca.eigenvec.
  loading_matrix <- sweep(
    loading_matrix,
    MARGIN = 2L,
    STATS = ifelse(pc_correlations < 0, -1, 1),
    FUN = "*"
  )
  colnames(loading_matrix) <- paste0("PC", seq_len(n_pcs))

  all_snp_ids <- gdsfmt::read.gdsn(
    gdsfmt::index.gdsn(gds, "snp.id")
  )
  selected_index <- match(
    as.character(snp_loading_object$snp.id),
    as.character(all_snp_ids)
  )
  if (anyNA(selected_index)) {
    stop("Could not match SNP loadings to GDS variant metadata.")
  }

  variant_ids <- as.character(gdsfmt::read.gdsn(
    gdsfmt::index.gdsn(gds, "snp.rs.id")
  )[selected_index])
  missing_variant_ids <- is.na(variant_ids) | variant_ids == ""
  variant_ids[missing_variant_ids] <- as.character(
    snp_loading_object$snp.id[missing_variant_ids]
  )

  allele_strings <- as.character(gdsfmt::read.gdsn(
    gdsfmt::index.gdsn(gds, "snp.allele")
  )[selected_index])
  split_alleles <- strsplit(allele_strings, "/", fixed = TRUE)
  # SNPRelate genotypes count the first ("A") allele in snp.allele.
  loading_allele <- vapply(
    split_alleles,
    function(x) if (length(x) >= 1L) x[[1L]] else NA_character_,
    character(1)
  )
  other_allele <- vapply(
    split_alleles,
    function(x) if (length(x) >= 2L) x[[2L]] else NA_character_,
    character(1)
  )

  output <- data.frame(
    CHR = gdsfmt::read.gdsn(
      gdsfmt::index.gdsn(gds, "snp.chromosome")
    )[selected_index],
    POS = gdsfmt::read.gdsn(
      gdsfmt::index.gdsn(gds, "snp.position")
    )[selected_index],
    SNP = variant_ids,
    LOADING_ALLELE = loading_allele,
    OTHER_ALLELE = other_allele,
    stringsAsFactors = FALSE
  )
  output <- cbind(output, loading_matrix)

  write_pca_loadings(
    output,
    outfile = outfile,
    n_pcs = n_pcs,
    expected_snps = variant_ids
  )
  message(
    "PC-AiR loading validation correlations: ",
    paste0(
      "PC", seq_len(n_pcs), "=",
      format(round(pc_correlations, 6), nsmall = 6),
      collapse = ", "
    )
  )

  invisible(pc_correlations)
}
