#### Utility functions ####

# are we running slurm?
is_slurm <- function() nchar(Sys.which("sbatch")) > 0
is_local <- function() !is_slurm()

# are we running snakemake?
is_snakemake <- function() !interactive() && exists("snakemake")

# how many cpus do we have on the local machine?
# if we're not running on the cluster, leave one cpu free.
local_cpus <- function() {
   if (is_snakemake()) {
      snakemake@threads
   } else if (is_slurm()) {
      out <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE"))
      if(!assertthat::is.count(out)) out <- 1L
      out
   } else {
      max(parallel::detectCores() - 1, 1)
   }
}

#### convenience function for writing a file and returning its name ####

ensure_directory <- function(file) {
   d <- dirname(file)
   if (!dir.exists(d)) dir.create(d)
}

write_and_return_file <- function(x, file, ...) {
   UseMethod("write_and_return_file")
}

write_and_return_file.XStringSet <- function(x, file, ...) {
   ensure_directory(file)
   Biostrings::writeXStringSet(x, file, ...)
   file
}

write_and_return_file.data.frame <- function(x, file, type = c("rds", "tsv"), ...) {
   ensure_directory(file)
   type = match.arg(type)
   switch(
      type,
      rds = saveRDS(x, file, ...),
      tsv = readr::write_tsv(x, file, ...),
      stop("Unknown file type: ", type)
   )
   file
}

write_and_return_file.character <- function(x, file, ...) {
   ensure_directory(file)
   writeLines(x, file, ...)
   file
}

write_and_return_file.ggplot <- function(x, file, ...) {
   ensure_directory(file)
   ggplot2::ggsave(file, plot = x, ...)
   file
}

write_and_return_file.default <- function(x, file, ...) {
   ensure_directory(file)
   saveRDS(x, file, ...)
   file
}

# helper functions to work with sequences sets that may be XStringSet,
# named character, or data.frame

# guess the column name in a data frame which refers to the sequence ID
find_name_col <- function(d) {
   stopifnot(is.data.frame(d))
   if ("seq_id" %in% names(d)) return("seq_id")
   if ("name" %in% names(d)) return("name")
   if ("ASV" %in% names(d)) return("ASV")
   if ("OTU" %in% names(d)) return("OTU")
   if (ncol(d) == 2) {
      if ("seq" %in% names(d)) return(names(d)[names(d) != "seq"])
      if ("sequence" %in% names(d)) return(names(d)[names(d) != "sequence"])
      return(names(d)[1])
   }
   stop("unable to determine sequence name column:", names(d))
}

# guess the column name in a data frame which refers to the sequence
find_seq_col <- function(d) {
   stopifnot(is.data.frame(d))
   if ("seq" %in% names(d)) return("seq")
   if ("sequence" %in% names(d)) return("sequence")
   if (ncol(d) == 2) {
      if ("seq_id" %in% names(d)) return(names(d)[names(d) != "seq_id"])
      if ("name" %in% names(d)) return(names(d)[names(d) != "name"])
      if ("ASV" %in% names(d)) return(names(d)[names(d) != "ASV"])
      if ("OTU" %in% names(d)) return(names(d)[names(d) != "OTU"])
      return(names(d)[2])
   }
   stop("unable to determine sequence column:", names(d))
}

# write the sequence as a fasta file
write_sequence <- function(seq, fname, ...) {
   UseMethod("write_sequence", seq)
}

write_sequence.data.frame <- function(seq, fname, seq_col = find_seq_col(seq),
                                      name_col = find_name_col(seq), ...) {
   dplyr::select(seq, !!name_col, !!seq_col) %>%
      tibble::deframe() %>%
      Biostrings::DNAStringSet() %>%
      Biostrings::writeXStringSet(fname, ...)
}

write_sequence.character <- function(seq, fname, ...) {
   Biostrings::DNAStringSet(seq) %>%
      Biostrings::writeXStringSet(fname, ...)
}

write_sequence.XStringSet <- function(seq, fname, ...) {
   Biostrings::writeXStringSet(seq, fname, ...)
}

# select elements from the sequence set.
select_sequence <- function(seq, which, negate = FALSE, ...) {
   UseMethod("select_sequence", seq)
}

select_sequence.data.frame <- function(seq, which, name_col = find_name_col(seq), negate = FALSE, ...) {
   if (isTRUE(negate)) {
      if (is.integer(which)  ||
          (is.numeric(which) && all(which == round(which)))) return(seq[-which,])
      if (is.logical(which)) return(seq[!which,])
      if (is.character(which)) return(dplyr::filter(seq, !.data[[name_col]] %in% which))
      stop("'which' should be integer, logical, or character")
   } else if (isFALSE(negate)) {
      if (is.integer(which)  ||
          (is.numeric(which) && all(which == round(which)))) return(seq[which,])
      if (is.logical(which)) return(seq[which,])
      if (is.character(which)) return(dplyr::filter(seq, .data[[name_col]] %in% which))
      stop("'which' should be integer, logical, or character")
   }
   stop("'negate' must be TRUE or FALSE")
}

select_sequence.default <- function(seq, which, negate = FALSE, ...) {
   if (isTRUE(negate)) {
      if (is.integer(which) || (is.numeric(which) && all(which == round(which)))) return(seq[-which])
      if (is.logical(which)) return(seq[!which])
      if (is.character(which)) return(seq[setdiff(names(seq), which)])
      stop("'which' should be integer, logical, or character")
   } else if (isFALSE(negate)) {
      return(seq[which])
   }
   stop("'negate' must be true or false")
}

# get the number of sequences in the sequence set
sequence_size <- function(seq, ...) {
   UseMethod("sequence_size", seq)
}

sequence_size.XStringSet <- function(seq, ...) {
   length(seq)
}

sequence_size.default <- function(seq, ...) {
   vctrs::vec_size(seq)
}

# convert a character to an ordered factor of taxonomic ranks
TAXRANKS <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
rank2factor <- function(x) {
   factor(x, levels = rev(TAXRANKS), ordered = TRUE)
}

superranks <- function(x, ranks = TAXRANKS) {
   ranks[rank2factor(ranks) > x]
}

subranks <- function(x, ranks = TAXRANKS) {
   ranks[rank2factor(ranks) < x]
}

# Format a number with "k" or "M" suffic for thousands or millions, as appropriate
k_or_M <- function(x, ..., .sep = "", .function = format) {
   dplyr::case_when(
      abs(x) > 1e6 ~ paste(.function(x / 1e6, ...), "M", sep = .sep),
      abs(x) > 1e3 ~ paste(.function(x / 1e3, ...), "k", sep = .sep),
      is.na(x) ~ NA_character_,
      TRUE ~ paste(.function(x, ...), " ", sep = .sep)
   )
}
