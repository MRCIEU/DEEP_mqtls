#!/usr/bin/env Rscript

# convert a GCTA sparse GRM into dense binary format for MLMA.
# Usage: Rscript sparse_grm_to_dense.R input_prefix output_prefix
# No SNP counts are invented: single-GRM MLMA does not require .grm.N.bin.
main <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) != 2L) stop("Usage: sparse_grm_to_dense.R input_prefix output_prefix")
    input <- args[1L]
    output <- args[2L]

    if (normalizePath(dirname(input)) == normalizePath(dirname(output)) &&
        basename(input) == basename(output)) stop("Input and output prefixes must differ")

    ids <- read.table(paste0(input, ".grm.id"), header = FALSE,
                      colClasses = "character", quote = "", comment.char = "")
    if (!nrow(ids) || !ncol(ids) %in% c(2L, 3L)) stop("Expected two or three columns in .grm.id")
    ids <- ids[, 1:2, drop = FALSE] # KING may also store an order column.
    if (anyNA(ids) || any(ids == "") || anyDuplicated(ids) || anyDuplicated(ids[[2L]])) {
        stop("Missing or duplicate sample IDs")
    }
    
    n <- nrow(ids)
    sp <- read.table(paste0(input, ".grm.sp"), header = FALSE,
                     colClasses = "numeric", quote = "", comment.char = "")
    if (ncol(sp) != 3L || !nrow(sp) || any(!is.finite(as.matrix(sp)))) {
        stop("Expected three finite numeric columns in .grm.sp")
    }
    
    i <- sp[[1L]]
    j <- sp[[2L]]
    x <- sp[[3L]]
    
    if (any(i != floor(i) | j != floor(j) | i < 0 | j < 0 | i >= n | j >= n)) {
        stop("GRM indices must be zero-based integers within the ID table")
    }
    if (any(i < j)) stop("Expected lower-triangular sparse GRM")
    if (anyDuplicated(sp[, 1:2])) stop("Duplicate GRM positions")
    if (sum(i == j) != n) stop("Missing GRM diagonal entries")

    ord <- order(i, j)
    i <- i[ord]; j <- j[ord]; x <- x[ord]
    counts <- tabulate(as.integer(i) + 1L, nbins = n)
    offsets <- c(0, cumsum(counts))
    tmp <- tempfile(pattern = paste0(basename(output), ".tmp-"), tmpdir = dirname(output))
    bin_tmp <- paste0(tmp, ".grm.bin")
    id_tmp <- paste0(tmp, ".grm.id")
    on.exit(unlink(c(bin_tmp, id_tmp)), add = TRUE)
    write_binary <- function() {
        con <- file(bin_tmp, "wb")
        on.exit(close(con))
        for (row in seq_len(n)) {
            pos <- seq.int(offsets[row] + 1, offsets[row + 1L])
            values <- numeric(row)
            values[j[pos] + 1L] <- x[pos]
            writeBin(values, con, size = 4L, endian = "little")
        }
    }
    write_binary()
    write.table(ids, id_tmp, quote = FALSE, row.names = FALSE,
                col.names = FALSE, sep = "\t")
    if (file.info(bin_tmp)$size != 4 * n * (n + 1) / 2) stop("Incorrect binary GRM size")

    # Validate every element, including implicit zeros, before publishing.
    verify_binary <- function() {
        con <- file(bin_tmp, "rb")
        on.exit(close(con))
        for (row in seq_len(n)) {
            pos <- seq.int(offsets[row] + 1, offsets[row + 1L])
            expected <- numeric(row)
            expected[j[pos] + 1L] <- x[pos]
            actual <- readBin(con, numeric(), n = row, size = 4L, endian = "little")
            tolerance <- pmax(abs(expected) * 2^-23, 2^-149)
            if (length(actual) != row || any(!is.finite(actual)) ||
                any(abs(actual - expected) > tolerance) ||
                any(actual[expected == 0] != 0)) stop("Binary GRM round-trip verification failed")
        }
    }
    verify_binary()
    check_ids <- read.table(id_tmp, header = FALSE, colClasses = "character",
                            quote = "", comment.char = "")
    if (!identical(ids, check_ids)) stop("Sample ID round-trip verification failed")
    if (!file.rename(bin_tmp, paste0(output, ".grm.bin"))) stop("Cannot publish binary GRM")
    if (!file.rename(id_tmp, paste0(output, ".grm.id"))) stop("Cannot publish GRM IDs")
    message("Saved thresholded dense GRM for ", n, " samples: ", output,
            ".grm.bin / .grm.id (original values preserved to float32 precision)")
}

main()
