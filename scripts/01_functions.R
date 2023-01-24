#' Precluster using USEARCH (As in Unite)
#'
#' This step is used to speed up blastclust by limiting the search space.
#' This can't be exactly replicated with VSEARCH yet, because it doesn't allow
#' specification of gap penalties in clustering.
#'
#' @param seqs (`character` string) a filename giving the sequences to cluster,
#' in Fasta format.
#' @param threshold (`numeric` scalar) fractional similarity threshold for
#' preclustering; a number between 0 and 1. (Default: 0.80)
#' @param exec (`character` string) path to the USEARCH executable
#'
#' @return a tibble with two columns: `seq_id` gives the names of the sequences
#' (as given in the input Fasta file), and `cluster` gives the name of the seed
#' sequence for the cluster each sequence is assigned to.
#'
#' @export
#'
usearch_cluster_fast <- function(seqs, threshold = 0.8, seqnames = names(seqs),
                                 which = TRUE, ncpu = local_cpus(),
                                 exec = Sys.which("usearch")) {
   UseMethod("usearch_cluster_fast", seqs)
}

usearch_cluster_fast.character <- function(
   seqs,
   threshold = 0.8,
   seqnames = names(seqs),
   which = TRUE,
   ncpu = local_cpus(),
   exec = Sys.which("usearch")
) {
   if (length(seqs) == 1 && file.exists(seqs)) {
      if (!missing(seqnames))
         warning("'seqnames' has no effect when 'seqs' is a file.")
      if (!missing(which))
         warning("'which' has no effect when 'seqs' is a file.")
      usearch_cluster_fast_filename(seqs, threshold, ncpu)
   }

   seqs <- Biostrings::DNAStringSet(seqs)
   usearch_cluster_fast.DNAStringSet(seqs, threshold, seqnames, which, ncpu, exec)
}

usearch_cluster_fast.DNAStringSet <- function(
   seqs,
   threshold = 0.8,
   seqnames = names(seqs),
   which = TRUE,
   ncpu = local_cpus(),
   exec = Sys.which("usearch")
) {
   # rename the sequences if necessary
   if (!isTRUE(all.equal(names(seqs), seqnames))) names(seqs) <- seqnames
   seqs <- seqs[which]
   # shortcut if only one sequence
   if (length(seqs) == 1) return(names(seqs))
   tf <- tempfile(pattern = "clust", fileext = ".fasta")
   Biostrings::writeXStringSet(seqs, tf)
   on.exit(unlink(tf))
   usearch_cluster_fast_filename(tf, threshold, ncpu, exec)
}

usearch_cluster_fast_filename <- function(
   seqs,
   threshold = 0.8,
   ncpu = local_cpus(),
   exec = Sys.which("usearch")
) {
   # temp files for IO with usearch
   tempuc <- tempfile(pattern = "clust", fileext = ".uc")
   on.exit(unlink(tempuc, force = TRUE))
   stopifnot(
      system2(
         command = exec,
         args = c(
            "--cluster_fast", seqs, # input file
            "--id", threshold, # similarity threshold
            "--gapopen", "0.0", # no gap open penalty
            "--gapext", "1.0/0.5E", # gap extension penalty 1 for internal, 0.5 for ends
            "--uc", tempuc # output file (uclust format)
         )
      ) == 0
   )
   # read and parse the uc file
   readr::read_tsv(
      tempuc,
      col_names = paste0("V", 1:10),
      col_types = "ciidfccccc",
      na = c("", "NA", "*")
   ) %>%
      dplyr::select(seq_id = V9, cluster = V10) %>%
      dplyr::mutate_all(stringr::str_replace, ";.*", "") %>%
      dplyr::mutate(cluster = dplyr::coalesce(cluster, seq_id)) %>%
      unique() %>%
      dplyr::group_by(cluster)
}

#' Title
#'
#' @param seqs
#' @param threshold
#' @param seqnames
#' @param which
#' @param ncpu
#' @param exec
#'
#' @return
#' @export
#'
#' @examples
usearch_cluster_fast_repeat <- function(
   seqs,
   threshold = c(0.975, 0.95, 0.9, 0.8, 0.6),
   seqnames = names(seqs),
   which = TRUE,
   ncpu = local_cpus(),
   exec = Sys.which("usearch")
) {
   if (is.character(seqs) && file.exists(seqs))
      stop("file inputs are not (yet) supported for usearch_cluster_fast_repeat")
   out <- list()
   for (t in threshold) {
      out <- c(
         out,
         list(usearch_cluster_fast(seqs, t, seqnames, which, ncpu, exec))
      )
      which <- unique(dplyr::last(out)$cluster)
   }
   purrr::reduce(
      out,
      ~dplyr::left_join(.x, .y, by = c("cluster" = "seq_id")) %>%
         dplyr::select(seq_id, cluster = cluster.y)
   ) %>%
      dplyr::group_by(cluster)
}

# #### its2_precluster_singletons
# its2_precluster_singletons = dplyr::filter(its2_precluster, dplyr::n() == 1)
# #### its2_precluster_clusters
# tar_group_by(
#    its2_precluster_clusters,
#    dplyr::filter(its2_precluster, dplyr::n() > 1),
#    cluster
# )

#' Cluster ASVs using blastclust
#'
#' @param seqs (`character` vector, file name, or `Biostrings::DNAStringSet`)
#' sequences to cluster
#' @param threshold (`numeric` scalar) percentage similarity threshold for
#' clustering. A number between 3 and 100
#' @param seqnames (`character`) names for the sequences (Default: `names(seqs)`)
#' @param which (`logical`, `integer`, or `character`) subset indices
#' indicating which sequences in `seqs` should be clustered. Has no effect if
#' `seqs` is a file name. (Default: `TRUE`, i.e. all the sequences)
#' @param ncpu (`integer` scalar) number of CPUs to use
#' @param outfile (`character` string giving a valid filename) filename to write
#' clustering output, if desired.
#' @param hits (`character` string giving a valid filename) filename to
#' write the hit table, if desired. For use in `blastclust_reclust()`.
#'
#' @return a `character` vector, where each line gives the sequence names
#' included in one cluster
#' @export
#'
#' @examples
blastclust <- function(seqs, threshold, seqnames = names(seqs), which = TRUE,
                       ncpu = local_cpus(), outfile = NULL, hits = NULL) {
   UseMethod("blastclust", seqs)
}

blastclust.character <- function(seqs, threshold, seqnames = names(seqs),
                                 which = TRUE, ncpu = local_cpus(),
                                 outfile = NULL, hits = NULL) {
   if (length(seqs) == 1 && file.exists(seqs)) {
      if (!missing(seqnames))
         warning("'seqnames' has no effect when 'seqs' is a file.")
      if (!missing(which))
         warning("'which' has no effect when 'seqs' is a file.")
      blastclust_filename(seqs, threshold, ncpu)
   }
   seqs <- Biostrings::DNAStringSet(seqs)
   blastclust.DNAStringSet(seqs, threshold, seqnames, which, ncpu, outfile, hits)
}

blastclust.DNAStringSet <- function(seqs, threshold, seqnames = names(seqs),
                                    which = TRUE,
                                    ncpu = local_cpus(),
                                    outfile = NULL, hits = NULL) {
   # rename the sequences if necessary
   if (!isTRUE(all.equal(names(seqs), seqnames))) names(seqs) <- seqnames
   seqs <- seqs[which]
   # shortcut if only one sequence
   if (length(seqs) == 1) return(names(seqs))
   tf <- tempfile(pattern = "clust", fileext = ".fasta")
   Biostrings::writeXStringSet(seqs, tf)
   on.exit(unlink(tf))
   blastclust_filename(tf, threshold, ncpu, outfile, hits)
}

blastclust_filename <- function(seqs, threshold, ncpu, outfile = NULL, hits = NULL) {
   if (is.null(outfile)) {
      outfile <- tempfile(pattern = "out")
      on.exit(unlink(outfile), TRUE)
   }
   stopifnot(
      system2(
         "blastclust",
         c(
            "-i", seqs, # input file
            "-S", threshold, # similarity threshold
            "-L", "0.95", # overlap threshold
            "-e", "F", # parse sequence names (FALSE)
            "-b", "F", # require coverage on both neighbors (false)
            "-W", "16", # word (kmer) size 16
            "-a", ncpu, # number of threads
            "-o", outfile, # output file
            if (!is.null(hits)) c("-s", hits) else NULL,
            "-p", "F" # sequences are proteins (FALSE)
         )
      ) == 0
   )
   readLines(outfile)
}


#' Recluster sequences already clustered with `blastclust()`
#'
#' @param hits (`character` filename) hit table generated by a call to
#' `blastclust()`
#' @param threshold (`numeric` scalar) clustering threshold
#' @param preclusters (`character` filename or vector) clustering results from
#' previous call to `blastclust()` or `blastclust_reclust()` at a lower
#' threshold (optional).
#' @param ncpu (`integer` scalar) number of threads to use
#' @param outfile (`character` filename) file to save output (e.g., to use for
#' `preclusters` in future calls to `blastclust_reclust()`)
#'
#' @return a `character` vector, where each line gives the sequence names
#' included in one cluster
#' @export
#'
#' @examples
blastclust_reclust <- function(hits, threshold, preclusters = NULL,
                               ncpu = local_cpus(), outfile = NULL) {
   if (is.list(preclusters) || length(preclusters) > 1L && (all(file.exists(preclusters)) || all(endsWith(preclusters, " ")))) {
      out <- lapply(preclusters, blastclust_reclust, hits = hits,
                    threshold = threshold, ncpu = ncpu, outfile = NULL)
      out <- unlist(out)
      if (!is.null(outfile)) {
         writeLines(out, outfile)
      }
      return(out)
   }
   if (!is.null(preclusters)) {
      stopifnot(is.character(preclusters))
      if (length(preclusters) > 1 || !file.exists(preclusters)) {
         preclusters <- trimws(preclusters)
         if (!grepl(" ", preclusters, fixed = TRUE)) {
            # shortcut if there is only one sequence
            if (!is.null(outfile)) {
               writelines(preclusters, outfile, sep = " ")
            }
            return(paste0(preclusters, " "))
         }
         tf <- tempfile("preclusters")
         writeLines(preclusters, tf, sep = " ")
         on.exit(unlink(tf))
         preclusters = tf
      }
   }
   if (is.null(outfile)) {
      outfile <- tempfile(pattern = "out")
      on.exit(unlink(outfile, force = TRUE), TRUE)
   }
   stopifnot(
      system2(
         "blastclust",
         c(
            "-r", hits, # input file
            if (!is.null(preclusters)) c("-l", preclusters) else NULL,
            "-S", threshold, # similarity threshold
            "-L", "0.95", # overlap threshold
            "-e", "F", # parse sequence names (FALSE)
            "-b", "F", # require coverage on both neighbors (false)
            "-W", "16", # word (kmer) size 16
            "-a", "8", # number of threads
            "-p", "F", # sequences are proteins (FALSE)
            "-o", outfile # output file
         )
      ) == 0
   )
   readLines(outfile)
}

blastclust2list <- function(clust) {
   if (length(clust) == 1 && file.exists(clust)) clust <- readLines(clust)
   strsplit(trimws(clust), " ", TRUE)
}

list2blastclust <- function(clust) {
   paste0(vapply(clust, paste, "", collapse = " "), " ")
}

check_blastclust_constraint <- function(preclust, clust, warnonly = FALSE) {
   if (is.character(preclust)) {
      if (length(preclust) == 1 && file.exists(preclust)) {
         preclust <- readLines(preclust)
      }
      preclust <- blastclust2list(preclust)
   }
   if (is.character(clust)) {
      if (length(clust) == 1 && file.exists(clust)) {
         clust <- readLines(clust)
      }
      clust <- blastclust2list(clust)
   }
   if (length(unlist(clust)) > length(unlist(preclust))) {
      f <- if (isTRUE(warnonly)) warning else stop
      f("output clusters contain ", length(unlist(clust)), " members,\n",
      "but there were only ", length(unlist(preclust)), " preclusters.")
   }
   if (length(unlist(preclust)) > length(unlist(clust))) {
      warning("not all preclusters are found in clust\n",
              "preclusters: ", length(unlist(preclust)), "\n",
              "clusters: ", length(unlist(clust)), "\n")
   }
   for (cl in clust) {
      inter <- vapply(preclust, optimotu::intersect_length, 1L, cl)
      w <- which(inter > 0)
      if (length(w) > 1) {
         f <- if (isTRUE(warnonly)) warning else stop
         f("inconsistent preclustering:\n cluster of size ", length(cl),
              " contains members from ", length(w), " preclusters\n",
              "with sizes ", paste(vapply(preclust[w], length, 0L), collapse = " "))
      }
   }
}

#' Do single-linkage clustering at a series of increasing similarity thresholds.
#'
#' @param seq (`character` vector, filename, or `Biostrings::DNAStringSet`)
#' sequences to cluster
#' @param seq_id (`character` vector) names for the sequences.  If they are
#' already named, this will replace the names.  Has no effect if `seq` is a
#' filename.
#' @param thresh_min (`numeric`) minimum percentage similarity threshold for
#' clustering. Number between 3 and 100.
#' @param thresh_max (`numeric`) maximum percentage similarity threshold for
#' clustering. Number between 3 and 100.
#' @param thresh_step (`numeric`) difference between successive percentage similarity thresholds for
#' clustering. Number <= thresh_max - thresh_min.
#' @param thresh_name (`character` vector) names for the thresholds.
#' @param which (`logical`, `character` or `integer` vector) subset of `seq` to
#' operate on.  Has no effect if `seq` is a filename.
#' @param ncpu (`integer` scalar) number of threads to use for calculating the
#' distance matrix and clustering
#' @param usearch (`character` scalar) path to usearch executable
#' @param ...
#'
#' @return `integer` matrix giving clustering results
#' @export
#'
#' @examples
single_link_repeat <- function(
   seq,
   seq_id = names(seq),
   thresh_min,
   thresh_max,
   thresh_step,
   thresh_name = NULL,
   which = TRUE,
   ncpu = local_cpus(),
   usearch = Sys.which("usearch"),
   ...
) {
   checkmate::assert_number(thresh_min, lower = 3, upper = 100)
   checkmate::assert_number(thresh_max, lower = thresh_min, upper = 100)
   checkmate::assert(
      checkmate::check_number(thresh_step, upper = thresh_max - thresh_min),
      thresh_min == thresh_max
   )
   nthresh <- ceiling((thresh_max - thresh_min)/thresh_step)
   checkmate::assert(
      checkmate::check_null(thresh_name),
      checkmate::check_character(thresh_name, len = nthresh)
   )
   is_file <- length(seq) == 1 && file.exists(seq)
   # if we only have one sequence, nothing to do.
   if (isFALSE(is_file)) {
      stopifnot(length(seq) == length(seq_id))
      s <- set_names(seq_is, seq_is)[which]
      if (length(s) <= 1) {
         return(
            matrix(
               rep(1L, length(s) * nthresh), nrow = nthresh,
               dimnames = list(
                  thresh = thresh_name,
                  seq = s
               )
            )
         )
      }
   }

   # normal case: multiple sequences, multiple thresholds
   out <- list()
   outfiles <- tempfile("out")
   on.exit(unlink(outfiles, force = TRUE), TRUE)
   hits <- tempfile("hits")
   on.exit(unlink(hits, force = TRUE))
   cat(sprintf("%s Generating hitlist for %d queries.\n", Sys.time(), length(seqs[which])))
   usearch_hitlist(
      seqs,
      threshold = threshold[1]/100,
      seqnames = seqnames,
      which = which,
      ncpu = ncpu,
      hits = hits,
      usearch = usearch
   )
   out <- list(
      blastclust_reclust(
         hits,
         threshold = threshold[1],
         ncpu = ncpu,
         outfile = outfiles[1]
      )
   )

   for (i in 2:length(threshold)) {
      outfiles[i] <- tempfile("out")
      out[[i]] <-
         blastclust_reclust(
            hits = hits,
            threshold = threshold[i],
            preclusters = outfiles[i-1],
            ncpu = ncpu,
         outfile = outfiles[i]
      )
   }
   if (!is.null(threshold_name)) names(out) <- threshold_name
   out
}

filter_preclust <- function(inclust = NULL, which) {
   if (is.list(which)) return(unlist(lapply(which, filter_preclust, inclust = inclust)))
   if (is.null(inclust)) {
      return(paste(c(which, ""), collapse = " "))
   }
   stopifnot(is.character(inclust))
   if (length(inclust) == 1 && file.exists(inclust)) {
      inclust <- readLines(inclust)
   }
   inclust <- trimws(inclust)
   inclust <- strsplit(inclust, split = " ")
   if (all(which %in% unlist(inclust))) {
      out <- lapply(inclust, intersect, which)
      out <- out[vapply(out, length, 0L) > 0]
      out <- vapply(out, paste, "", collapse = " ")
      out <- paste0(out, " ")
      return(out)
   } else {
      return(paste(c(which, ""), collapse = " "))
   }
}

is_list_of_character <- function(x) {
   all(vapply(x, is.character, TRUE))
}

blastclust_repeat_constrained <- function(seqs, threshold, seqnames = names(seqs),
                                          which = rep(list(TRUE), length(threshold)),
                                          threshold_name = names(threshold),
                                          ncpu = local_cpus(),
                                          hitlist_method = c("blastclust", "usearch"),
                                          usearch = Sys.which("usearch"),
                                          ...){
   stopifnot(is.null(threshold_name) || length(threshold) == length(threshold_name))
   hitlist_method <- match.arg(hitlist_method, several.ok = FALSE)
   is_file <- length(seqs) == 1 && file.exists(seqs)
   # if we only have one sequence, nothing to do.
   if (isFALSE(is_file)) {
      stopifnot(length(seqs) == length(seqnames))
      s <- set_names(seqnames, seqnames)[unique(unlist(which))]
   }

   first_threshold_is_dummy <- FALSE

   if (isFALSE(is_file) && length(s) <= 1) {
      out <- rep(list(names(s)), length(threshold))
   } else if (length(threshold) == 1) {
      # if there is only one threshold, then no need to cache the hit table
      out <- list(
         blastclust(
            seqs = seqs,
            threshold = threshold,
            seqnames = seqnames,
            which = which,
            ncpu = ncpu
         )
      )
   } else {
      # normal case: multiple sequences, multiple thresholds

      # which should either be a single specification which applies to all
      # thresholds, or a list with one element per threshold.
      # elements should be either all the logical scalar TRUE, or a character
      # giving sequence names.
      # integer indices and logical vector indices are not supported.
      if (length(which) == 1) {
         if (is.list(which)) {
            which <- rep(which, length(threshold))
         } else {
            which <- rep(list(which), length(threshold))
         }
      }
      stopifnot(length(which) == length(threshold))
      which_is_TRUE <- all(vapply(which, isTRUE, TRUE))
      which_is_character <- all(vapply(which, is.character, TRUE))
      which_is_list_of_character <- all(vapply(which, is_list_of_character, TRUE))
      stopifnot(which_is_character || which_is_list_of_character)
      if (which_is_character) which <- lapply(which, list)
      all_which <- unlist(which) %>% unique()
      if (length(which[[1]]) > 1 || !all(all_which %in% unlist(which[[1]]))) {
         which <- c(list(list(all_which)), which)
         threshold <- c(threshold[1], threshold)
         first_threshold_is_dummy <- TRUE
      }
      out <- list()
      outfiles <- tempfile("out")
      on.exit(unlink(outfiles, force = TRUE), TRUE)
      hits <- tempfile("hits")
      on.exit(unlink(hits, force = TRUE), TRUE)
      if (hitlist_method == "blastclust") {
         unpruned_hits <- tempfile("unpruned")
         on.exit(unlink(pruned_hits, force = TRUE), TRUE)
         out <- list(
            blastclust(seqs,
                       threshold = threshold[1],
                       seqnames = seqnames,
                       which = which[[1]][[1]],
                       ncpu = ncpu,
                       outfile = outfiles[1],
                       hits = unpruned_hits
            )
         )
         prune_hitlist(unpruned_hits, hits)
         unlink(unpruned_hits)
      } else if (hitlist_method == "usearch") {
         usearch_hitlist(
            seqs,
            threshold = threshold[1]/100,
            seqnames = seqnames,
            which = which[[1]][[1]],
            ncpu = ncpu,
            hits = hits,
            usearch = usearch
         )
         out <- list(
            blastclust_reclust(
               hits,
               threshold = threshold[1],
               ncpu = ncpu,
               outfile = outfiles[1]
            )
         )
      } else {
         stop("unknown hitlist method: ", hitlist_method)
      }
      for (i in 2:length(threshold)) {
         # if we have explicit selection of sequences at each level, then either
         # filter the previous round of clustering to only include our sequences,
         # or make a new, single cluster containing all the selected sequences
         # if the previous round did not include them all
         prev <- outfiles[i-1]
         if (threshold[i-1] > threshold[i]) prev <- NULL
         if (!isTRUE(which[[i]])) {
            preclust <- filter_preclust(prev, which[[i]])
            message("Checking preclusters against taxonomic clusters...\n")
            check_blastclust_constraint(which[[i]], preclust)
            message("Checking preclusters against previous clustering round...\n")
            check_blastclust_constraint(prev, preclust, warnonly = TRUE)
            prev <- preclust
         }
         outfiles[i] <- tempfile("out")
         out[[i]] <-
            blastclust_reclust(
               hits = hits,
               threshold = threshold[i],
               preclusters = prev,
               ncpu = ncpu,
               outfile = outfiles[i]
            )

         message("Checking clusters against preclusters...\n")
         check_blastclust_constraint(prev, out[[i]])
         message("Checking clusters against taxonomic clusters...\n")
         check_blastclust_constraint(which[[i]], out[[i]])
      }
   }
   if (first_threshold_is_dummy) {
      out <- out[-1]
   }
   if (!is.null(threshold_name)) names(out) <- threshold_name
   out
}

blastclust_usearch <- function(
   seqs,
   threshold,
   seqnames = names(seqs),
   which = TRUE,
   ncpu = local_cpus(),
   usearch = Sys.which("usearch")
) {
   is_file <- length(seqs) == 1 && file.exists(seqs)
   if (isFALSE(is_file)) {
      stopifnot(length(seqs) == length(seqnames))
      if (length(seqs) == 0) {
         return(character(0))
      } else if (length(seqs) == 1) {
         return(paste0(seqnames, " "))
      }
   }

   hits <- tempfile("hits")
   on.exit(unlink(hits, force = TRUE), TRUE)
   usearch_hitlist(
      seqs,
      threshold = threshold/100,
      seqnames = seqnames,
      which = which,
      ncpu = ncpu,
      hits = hits,
      usearch = usearch
   )
   blastclust_reclust(
      hits,
      threshold = threshold,
      ncpu = ncpu
   )
}

# check that two clusters are really distinct
# to keep things fast, do a small number of sequences at a time, and shortcut
# if a match is found
blast_check <- function(seq1, seq2, threshold, seqnames1 = names(seq1),
                        seqnames2 = names(seq2), ncpu = local_cpus()) {
   UseMethod("blast_check", seq1)
}

blast_check.character <- function(seq1, seq2, threshold, seqnames1 = names(seq1),
                                  seqnames2 = names(seq2), ncpu = local_cpus()) {
   stopifnot(is.character(seq2))
   isfile1 <- length(seq1) == 1 && file.exists(seq1)
   isfile2 <- length(seq2) == 1 && file.exists(seq2)
   stopifnot (isfile1 == isfile2)

   if (isfile1) {
         if (!missing(seqnames1) || !missing(seqnames2))
            warning("'seqnames' has no effect when 'seqs' is a file.")
         blast_check_filename(seq1, seq2, threshold, ncpu)
   } else {
      seq1 <- Biostrings::DNAStringSet(seq1)
      seq2 <- Biostrings::DNAStringSet(seq2)
      blast_check.DNAStringSet(seq1, seq2, threshold, seqnames1, seqnames2, ncpu)
   }
}

blast_check.DNAStringSet <- function(seq1, seq2, threshold,
                                     seqnames1 = names(seq1),
                                     seqnames2 = names(seq2),
                                     ncpu = local_cpus()) {
   # rename the sequences if necessary
   if (!isTRUE(all.equal(names(seq1), seqnames1))) names(seq1) <- seqnames1
   if (!isTRUE(all.equal(names(seq2), seqnames2))) names(seq2) <- seqnames2
   tf1 <- tempfile(pattern = "clust", fileext = ".fasta")
   tf2 <- tempfile(pattern = "clust", fileext = ".fasta")
   Biostrings::writeXStringSet(seq1, tf1)
   Biostrings::writeXStringSet(seq2, tf2)
   on.exit(unlink(c(tf1, tf2)))
   blast_check_filename(tf1, tf2, threshold, ncpu)
}

blast_check_filename <- function(seq1, seq2, threshold, ncpu) {
   if (is.null(outfile)) {
      outfile <- tempfile(pattern = "out")
      on.exit(unlink(outfile), TRUE)
   }
   stopifnot(
      system2(
         "blastclust",
         c(
            "-i", seqs, # input file
            "-S", threshold, # similarity threshold
            "-L", "0.95", # overlap threshold
            "-e", "F", # parse sequence names (FALSE)
            "-b", "F", # require coverage on both neighbors (false)
            "-W", "16", # word (kmer) size 16
            "-a", ncpu, # number of threads
            "-o", outfile, # output file
            if (!is.null(hits)) c("-s", hits) else NULL,
            "-p", "F" # sequences are proteins (FALSE)
         )
      ) == 0
   )
   readLines(outfile)
}

parse_hitlist <- function(file) {
   size <- file.size(file)
   f <- file(file, "rb")
   IDtype <- readBin(f, "integer", size = 4)
   ListSz <- readBin(f, "integer", size = 4)
   SeqID <- NULL
   nseq <- 0L
   ListFileSz <- 0L
   if (IDtype == 0L) {
      ListFileSz <- ListSz
      SeqID <- readChar(f, ListSz)
      SeqID <- trimws(SeqID)
      SeqID <- strsplit(SeqID, " ")[[1]]
      nseq <- length(SeqID)
   } else if (IDtype == 1L) {
      ListFileSz <- ListSz * 4
      SeqID <- readBin(f, "integer", n = ListSz, size = 4)
      nseq <- ListSz
   } else {
      stop("Unknown IDtype: ", IDtype)
   }
   seqlen <- readBin(f, "integer", n = nseq, size = 4)
   nhits <- (size - 8L - ListFileSz - 4L*nseq)/32L
   N1 <- integer(nhits)
   N2 <- integer(nhits)
   HSPL1 <- integer(nseq)
   HSPL2 <- integer(nseq)
   # Score <- numeric(nseq)
   PercId <- numeric(nhits)
   for (i in seq_len(nhits)) {
      N1[i] <- readBin(f, "integer", size = 4)
      N2[i] <- readBin(f, "integer", size = 4)
      HSPL1[i] <-
         readBin(f, "integer", size = 4)
      HSPL2[i] <-
         readBin(f, "integer", size = 4)
      # Score[i] <-
         readBin(f, "double", size = 8)
      PercId[i] <-
         readBin(f, "double", size = 8)
   }
   close(f)
   o <- order(N1, N2)
   out <- data.frame(
      seq1 = SeqID[N1[o] + 1],
      seq2 = SeqID[N2[o] + 1],
      cov1 = HSPL1[o]/seqlen[N1[o] + 1],
      cov2 = HSPL2[o]/seqlen[N2[o] + 1],
      dist = 1 - PercId[o]/100,
      stringsAsFactors = FALSE
   )
   structure(
      out,
      seqlen = data.frame(
         seq_id = SeqID,
         seq_len = seqlen,
         stringsAsFactors = FALSE
      )
   )
}

prune_hitlist <- function(filein, fileout, min_cover = 0.95) {
   size <- file.size(filein)
   fin <- file(filein, "rb")
   on.exit(close(fin))
   fout <- file(fileout, "wb")
   on.exit(close(fout), add = TRUE)
   IDtype <- readBin(fin, "integer", size = 4)
   writeBin(IDtype, fout, size = 4)
   ListSz <- readBin(fin, "integer", size = 4)
   writeBin(ListSz, fout, size = 4)
   ListFileSz <- 0L
   nseq <- 0L
   if (IDtype == 0L) {
      ListFileSz <- ListSz
      SeqID <- readChar(fin, ListSz)
      writeChar(SeqID, fout, eos = NULL)
      SeqID <- trimws(SeqID)
      SeqID <- strsplit(SeqID, " ")[[1]]
      nseq <- length(SeqID)
   } else if (IDtype == 1L) {
      ListFileSz <- ListSz * 4
      SeqID <- readBin(fin, "integer", n = ListSz, size = 4)
      writeBin(SeqID, fout)
      nseq <- ListSz
   } else {
      stop("Unknown IDtype: ", IDtype)
   }
   nhits <- (size - 8L - ListFileSz - 4L*nseq)/32L
   seqlen <- readBin(fin, "integer", n = nseq, size = 4)
   writeBin(seqlen, fout, size = 4)
   d <- readBin(fin, "integer", size = 4, n = 800000)
   while (length(d) > 0) {
      d <- split(d, 1:8)
      d <- as.data.frame(d)
      d <- d[d$X3/seqlen[d$X1+1]>= min_cover | d$X4/seqlen[d$X2+1] >= min_cover,]
      d <- c(apply(d, 1, c))
      writeBin(d, fout, size = 4)
      d <- readBin(fin, "integer", size = 4, n = 800000)
   }
}

split_hitlist <- function(hitlist, taxa, base = sub("\\.hits$", "_", hitlist)) {
   outnames <- paste0(base, names(taxa), ".hits")
   names(outnames) <- names(taxa)
   outfiles <- vector("list", length(taxa))
   names(outfiles) <- names(taxa)

   fin <- file(hitlist, "rb")
   on.exit(close(fin), TRUE)

   IDtype <- readBin(fin, "integer", size = 4)
   #writeBin(IDtype, fout, size = 4)
   ListSz <- readBin(fin, "integer", size = 4)
   #writeBin(ListSz, fout, size = 4)
   ListFileSz <- 0L
   nseq <- 0L
   if (IDtype == 0L) {
      ListFileSz <- ListSz
      SeqID <- readChar(fin, ListSz)
      SeqID <- trimws(SeqID)
      SeqID <- strsplit(SeqID, " ")[[1]]
      nseq <- length(SeqID)
      seqlen <- readBin(fin, "integer", n = nseq, size = 4)
      seq_key <- rep(0L, nseq)
      which_taxa <- lapply(taxa, function(x) which(SeqID %in% x))
      for (i in seq_along(which_taxa)) {
         outfiles[[i]] <- file(outnames[i], "wb")
         on.exit(close(outfiles[[i]]), TRUE)
         seq_key[which_taxa[[i]]] <- i
         writeBin(IDtype, outfiles[[i]], size = 4)
         writeBin(length(which_taxa)[[i]], outfiles[[i]], size = 4)
         writeChar(
            paste0(SeqID[which_taxa[[i]]], " ", collapse = ""),
            outfiles[[i]],
            eos = NULL
         )
         writeBin(seqLen[which_taxa[[i]]], outfiles[[i]], size = 4)
      }
   } else if (IDtype == 1L) {
      ListFileSz <- ListSz * 4
      SeqID <- readBin(fin, "integer", n = ListSz, size = 4)
      writeBin(SeqID, fout)
      nseq <- ListSz
      seqlen <- readBin(fin, "integer", n = nseq, size = 4)
      seq_key <- rep(0L, nseq)
      which_taxa <- lapply(taxa, function(x) which(SeqID %in% x))
      for (i in seq_along(which_taxa)) {
         outfiles[[i]] <- file(outnames[i], "wb")
         on.exit(close(outfiles[[i]]), TRUE)
         seq_key[which_taxa[[i]]] <- i
         writeBin(IDtype, outfiles[[i]], size = 4)
         writeBin(length(which_taxa[[i]]), outfiles[[i]], size = 4)
         writeBin(SeqID[which_taxa[[i]]], outfiles[[i]], size = 4)
         writeBin(seqLen[which_taxa[[i]]], outfiles[[i]], size = 4)
      }
   } else {
      stop("Unknown IDtype: ", IDtype)
   }
   d <- readBin(fin, "integer", size = 4, n = 800000)
   while (length(d) > 0) {
      d <- split(d, 1:8)
      d <- as.data.frame(d)
      taxon1 <- seq_key[d$X1]
      taxon2 <- seq_key[d$X2]
      keep <- taxon1 == taxon2 & taxon1 > 0
      d <- d[keep,]
      taxon1 <- taxon1[keep]
      for (t in unique(taxon1)) {
         d2 <- d[taxon1 == t,]
         d2 <- c(apply(d2, 1, c))
         writeBin(d, outfiles[[t]], size = 4)
      }
      d <- readBin(fin, "integer", size = 4, n = 800000)
   }
   outnames
}

usearch_hitlist <- function(seqs, threshold, seqnames = names(seqs), which = TRUE,
                            ncpu = local_cpus(), hits = NULL, usearch = Sys.which("usearch")) {
   UseMethod("usearch_hitlist", seqs)
}

usearch_hitlist.character <- function(seqs, threshold, seqnames = names(seqs),
                                 which = TRUE, ncpu = local_cpus(), hits = NULL, usearch = Sys.which("usearch")) {
   if (length(seqs) == 1 && file.exists(seqs)) {
      if (!missing(seqnames))
         warning("'seqnames' has no effect when 'seqs' is a file.")
      if (!missing(which))
         warning("'which' has no effect when 'seqs' is a file.")
      index <- Biostrings::fasta.seqlengths(seqs)
      do_usearch_hitlist(
         seqs,
         seqlen = index,
         names = names(index),
         threshold = threshold,
         ncpu = ncpu,
         hits = hits,
         usearch = usearch
      )
   } else {
      seqs <- Biostrings::DNAStringSet(seqs)
      usearch_hitlist.DNAStringSet(seqs, threshold, seqnames, which, ncpu,
                                   hits, usearch = usearch)
   }
}

usearch_hitlist.DNAStringSet <- function(seqs, threshold, seqnames = names(seqs),
                                    which = TRUE,
                                    ncpu = local_cpus(),
                                    hits = NULL,
                                    usearch = Sys.which("usearch")) {
   # rename the sequences if necessary
   if (!isTRUE(all.equal(names(seqs), seqnames))) names(seqs) <- seqnames
   seqs <- seqs[which]
   # shortcut if only one sequence
   if (length(seqs) == 1) return(names(seqs))
   tf <- tempfile(pattern = "clust", fileext = ".fasta")
   Biostrings::writeXStringSet(seqs, tf)
   on.exit(unlink(tf))
   do_usearch_hitlist(tf, seqlen = Biostrings::nchar(seqs), names = names(seqs),
                      threshold = threshold, ncpu = ncpu, hits = hits,
                      usearch = usearch)
}

do_usearch_hitlist <- function(seqs, seqlen, names, threshold, ncpu, hits,
                               usearch = Sys.which("usearch")) {
   if (!methods::is(hits, "connection")) {
      hits <- file(hits, open = "wb")
   }
   if (!isOpen(hits)) open(hits, "wb")

   on.exit(close(hits), TRUE)
   # list type 0 (names)
   # list size (characters in names list)
   writeBin(c(0L, sum(nchar(names)) + length(names)), hits, size = 4L)
   # list of names
   writeChar(paste0(names, " ", collapse = ""), hits, eos = NULL)
   # sequence lengths
   writeBin(seqlen, hits, size = 4L)
   seqidx <- seq_along(names) - 1L
   names(seqidx) <- names
   if (is.null(names(seqlen))) names(seqlen) <- names
   fifoname <- tempfile("fifo")
   stopifnot(system2("mkfifo", fifoname) == 0)
   on.exit(unlink(fifoname), TRUE)
   f <- fifo(fifoname)
   system2(
      usearch,
      c(
         "-calc_distmx", seqs, # input file
         "-tabbedout", fifoname, # output fifo
         "-maxdist", 1-threshold, # similarity threshold
         "-termdist", min(1, 1.5*(1-threshold)), # threshold for udist
         "-lopen", "1", # gap opening
         "-lext", "1", # gap extend
         # "-pattern", "111010010111", # pattern gives better result than kmers maybe?
         "-threads", ncpu
      ),
      wait = FALSE
   )
   open(f, mode = "r")
   on.exit(close(f), TRUE)
   d <- readLines(f, n = 100000)
   while (length(d) > 0) {
      d <- strsplit(d, "\t", fixed = TRUE)
      d <- unlist(d)
      d <- matrix(d, ncol = 3, byrow = TRUE)
      d <- as.data.frame(d)
      seq1 <- matrix(writeBin(seqidx[d$V1], raw(), size = 4), nrow = 4)
      seq2 <- matrix(writeBin(seqidx[d$V2], raw(), size = 4), nrow = 4)
      hsp1 <- seqlen[d$V1]
      hsp2 <- seqlen[d$V2]
      ident <- 1 - as.numeric(d$V3)
      score <- matrix(writeBin(ident * pmin(hsp1, hsp2), raw(), size = 8), nrow = 8)
      hsp1 <- matrix(writeBin(hsp1, raw(), size = 4), nrow = 4)
      hsp2 <- matrix(writeBin(hsp2, raw(), size = 4), nrow = 4)
      ident <- matrix(writeBin(100*ident, raw(), size = 8), nrow = 8)
      d <- c(rbind(seq1, seq2, hsp1, hsp2, score, ident))
      writeBin(d, hits)
      d <- readLines(f, n = 100000)
   }
}

#' Do single-linkage clustering at a series of increasing similarity thresholds.
#'
#' @param seq (`character` vector, filename, or `Biostrings::DNAStringSet`)
#' sequences to cluster
#' @param seq_id (`character` vector) names for the sequences.  If they are
#' already named, this will replace the names.  Has no effect if `seq` is a
#' filename.
#' @param thresh_min (`numeric`) minimum sequence dissimilarity threshold for
#' clustering. Number between 0 and 1.
#' @param thresh_max (`numeric`) maximum sequence dissimilarity threshold for
#' clustering. Number between 0 and 1.
#' @param thresh_step (`numeric`) difference between successive percentage similarity thresholds for
#' clustering. Number <= thresh_max - thresh_min.
#' @param thresh_name (`character` vector) names for the thresholds.
#' @param which (`logical`, `character` or `integer` vector, or a list of these)
#' subset(s) of `seq` to operate on.  If `seq` is a filename, then the distance
#' matrix for the full file will be calculated, but only the selected sequences
#' will be clustered.  Otherwise, distances will only be calculated for
#' sequences occurring in the subset(s)
#' @param ncpu (`integer` scalar) number of threads to use for calculating the
#' distance matrix and clustering
#' @param usearch (`character` scalar) path to usearch executable
#' @param ...
#'
#' @return `integer` matrix giving clustering results
#' @export
usearch_singlelink <- function(seq,
                               seq_id = names(seq),
                               thresh_max, thresh_min, thresh_step,
                               thresh_names = NULL,
                               which = TRUE,
                               ncpu = local_cpus(),
                               usearch = Sys.which("usearch")) {
   UseMethod("usearch_singlelink", seq)
}

usearch_singlelink.character <- function(
   seq,
   seq_id = names(seq),
   thresh_max, thresh_min, thresh_step,
   thresh_names = NULL,
   which = TRUE, ncpu = local_cpus(), usearch = Sys.which("usearch")
) {
   if (length(seq) == 1 && file.exists(seq)) {
      index <- Biostrings::fasta.seqlengths(seq)
      if (!missing(seq_id))
         warning("'seq_id' has no effect when 'seq' is a file.")
      if (!all(names(index) == as.character(seq_along(index)))) {
         # write a temp version of the file which has headers as integer indices
         tf <- tempfile(pattern = "clust", fileext = ".fasta")
         tc <- file(tf, open = "w")
         on.exit(close(tc))
         sc <- file(seq, open = "r")
         on.exit(close(sc), add = TRUE)
         nseq <- 0L
         while(length(lines <- readLines(sc, n = 10000)) > 0L) {
            headers <- grep("^>", lines)
            lines[headers] <- sprintf(">%d", seq_along(headers) + nseq - 1L)
            writeLines(lines, tc)
            nseq <- nseq + length(headers)
         }
         close(sc)
         close(tc)
         on.exit()
      } else {
         tf <- seq
      }
      do_usearch_singlelink(
         seq_file = tf,
         seq_id = names(index),
         thresh_max = thresh_max,
         thresh_min = thresh_min,
         thresh_step = thresh_step,
         thresh_names = thresh_names,
         which = which,
         ncpu = ncpu,
         usearch = usearch
      )
   } else {
      mycall <- match.call()
      mycall[[1]] <- usearch_singlelink.DNAStringSet
      newseq <- quote(Biostrings::DNAStringSet(seq))
      newseq[[2]] <- mycall$seq
      mycall$seq <- newseq
      eval(mycall, envir = parent.frame())
   }
}

usearch_singlelink.DNAStringSet <- function(seq,
                                            seq_id = names(seq),
                                            thresh_max, thresh_min, thresh_step,
                                            thresh_names = NULL,
                                            which = TRUE,
                                            ncpu = local_cpus(),
                                            usearch = Sys.which("usearch")) {
   # rename the sequences if necessary
   if (!isTRUE(all.equal(names(seq), seq_id))) names(seq) <- seq_id
   if (is.list(which)) {
      if (all(vapply(which, is.logical, TRUE))) {
         seq <- seq[Reduce(magrittr::or, which)]
      } else if (all(vapply(which, is.integer, TRUE))) {
         seq <- seq[sort(unique(unlist(which)))]
      } else if (is_list_of_character(which)) {
         seq <- seq[sort(unique(unlist(which)))]
      } else {
         stop("'which' must be a character, integer, or logical vector, or a",
              " list of one of these.")
      }
   } else {
      seq <- seq[which]
   }
   # shortcut if only one sequence
   if (length(seq) == 1) return(names(seq))
   tf <- tempfile(pattern = "clust", fileext = ".fasta")
   Biostrings::writeXStringSet(`names<-`(seq, seq_along(seq) - 1L), tf)
   on.exit(unlink(tf))
   do_usearch_singlelink(
      seq_file = tf,
      seq_id = names(seq),
      thresh_max = thresh_max,
      thresh_min = thresh_min,
      thresh_step = thresh_step,
      thresh_names = thresh_names,
      which = which,
      ncpu = ncpu,
      usearch = usearch
   )
}

do_usearch_singlelink <- function(seq_file, seq_id,
                                  thresh_max, thresh_min, thresh_step,
                                  thresh_names,
                                  which,
                                  ncpu, usearch = Sys.which("usearch")) {

   checkmate::assert_number(thresh_min, lower = 0, upper = 1)
   checkmate::assert_number(thresh_max, lower = thresh_min, upper = 1)
   checkmate::assert(
      checkmate::check_number(thresh_step, upper = thresh_max - thresh_min),
      thresh_min == thresh_max
   )
   nthresh <- length(seq(thresh_min, thresh_max, thresh_step))
   checkmate::assert(
      checkmate::check_null(thresh_names),
      checkmate::check_character(thresh_names, len = nthresh)
   )
   checkmate::assert(
      checkmate::check_integer(which, lower = 1, upper = length(seq_id), any.missing = FALSE),
      checkmate::check_subset(which, seq_id),
      checkmate::check_logical(which, any.missing = FALSE),
      checkmate::check_list(which, types = "logical", any.missing = FALSE, min.len = 1),
      checkmate::check_list(which, types = "integerish", any.missing = FALSE, min.len = 1),
      checkmate::check_list(which, types = "character", any.missing = FALSE, min.len = 1)
   )
   if (is.list(which) && !is.character(which[[1]])) {
      which <- lapply(which, `[`, x = seq_id)
   }

   fifoname <- tempfile("fifo")
   stopifnot(system2("mkfifo", fifoname) == 0)
   on.exit(unlink(fifoname), TRUE)
   f <- fifo(fifoname)
   system2(
      usearch,
      c(
         "-calc_distmx", seq_file, # input file
         "-tabbedout", fifoname, # output fifo
         "-maxdist", thresh_max, # similarity threshold
         "-termdist", min(1, 1.5*thresh_max), # threshold for udist
         "-lopen", "1", # gap opening
         "-lext", "1", # gap extend
         # "-pattern", "111010010111", # pattern gives better result than kmers maybe?
         "-threads", ncpu
      ),
      wait = FALSE
   )
   if (is.list(which)) {
      out <- optimotu::single_linkage_multi(
         file = fifoname,
         seq_id,
         dmin = thresh_min,
         dmax = thresh_max,
         dstep = thresh_step,
         preclust = which,
         threads = ncpu
      )
      out <- lapply(out, `rownames<-`, thresh_names)
   } else {
      out <- optimotu::single_linkage_pool(
         file = fifoname,
         seq_id,
         dmin = thresh_min,
         dmax = thresh_max,
         dstep = thresh_step
      )
      row.names(out) <- thresh_names
   }
   out
}

hitlist2dist <- function(hitlist, fill = 1, mincover = 0) {
   seqlens <- attr(hitlist, "seqlen")
   nseq <- nrow(seqlens)
   nhits <- nseq*(nseq-1L)/2L
   seq1 <- integer(nhits)
   seq2 <- integer(nhits)
   k <- 1
   for (i in seq_len(nseq-1L)) {
      for (j in i + seq_len(nseq-i)) {
         seq1[k] <- seqlens$seq_id[i]
         seq2[k] <- seqlens$seq_id[j]
         k <- k + 1
      }
   }
   dist <-
      dplyr::left_join(
         tibble::tibble(seq1, seq2),
         dplyr::filter(hitlist, cov1 >= mincover | cov2 >= mincover),
         by = c("seq1", "seq2")
      ) %>%
      dplyr::pull(dist) %>%
      tidyr::replace_na(fill)
   structure(
      dist,
      Size = nseq,
      Labels = seqlens$seq_id,
      Diag = FALSE,
      Upper = FALSE,
      class = "dist"
   )
}

blast <- function(query, subject, threshold, min_cover = 95, ncpu = local_cpus()) {
   subjectfile <- tempfile(fileext = ".fasta")
   dbname <- tempfile("db")
   queryfile <- tempfile(fileext = ".fasta")
   on.exit(unlink(c(subjectfile, queryfile, paste(dbname, c("ndb", "nhr", "nin", "not", "nsq", "ntf", "nto", sep = "."))), force = TRUE))
   write_sequence(subject, subjectfile)
   stopifnot(
      system2(
         "makeblastdb",
         c("-in", subjectfile, "-dbtype", "nucl", "-out", dbname)
      ) == 0
   )
   write_sequence(query, queryfile)
   blasthits <- system2(
      "blastn",
      c(
         "-query", shQuote(queryfile),
         "-db", shQuote(dbname),
         "-outfmt", shQuote("6 qseqid sseqid pident qcovs"),
         "-max_hsps", "1",
         "-max_target_seqs", "1",
         "-strand plus",
         "-perc_identity", threshold,
         "-qcov_hsp_perc", min_cover,
         "-num_threads", ncpu,
         "-mt_mode", "1"
      ),
      stdout = TRUE
   )
   stopifnot(attr(blasthits, "status") == 0)
   readr::read_tsv(
      I(blasthits),
      col_names = c("qseqid", "sseqid", "pident", "qcovs"),
      col_types = "ccdd"
   )
}

asv_table_callback <- function(prev, chunk) {
   cat("prev:\n")
   print(prev)
   cat("chunk:\n")
   print(chunk)
   dplyr::rename(chunk, ASV = 1) %>%
      #transform to long/sparse format
      tidyr::pivot_longer(
         cols = -1,
         names_to = "sample",
         values_to = "nread"
      ) %>%
      dplyr::filter(nread > 0)
}

blastclust_closed_ref <- function(query, ref, threshold, ...) {
   out <- tibble::tibble(ASV = character(0), cluster = character(0))
   while(sequence_size(query) > 0 && sequence_size(ref) > 0) {
      result <- blast(query, ref, threshold, ...) %>%
         dplyr::select(ASV = qseqid, cluster = sseqid)
      if (nrow(out) > 0) result <- dplyr::left_join(result[,"ASV"], out, by = "ASV")
      out <- dplyr::bind_rows(out, result)
      ref <- select_sequence(query, result$ASV)
      query <- select_sequence(query, result$ASV, negate = TRUE)
   }
   out
}

vsearch_usearch_global <- function(query, ref, threshold, ncpu = local_cpus()) {
   tquery <- tempfile("query", fileext = ".fasta")
   tref <- tempfile("ref", fileext = ".fasta")
   on.exit(unlink(c(tquery, tref), force = TRUE))
   write_sequence(query, tquery)
   write_sequence(ref, tref)
   uc = system(
      paste(
         "vsearch",
         "--usearch_global", tquery,
         "--db", tref,
         "--id", threshold,
         "--uc", "-",
         "--maxaccepts", "100",
         "--top_hits_only",
         "--threads", ncpu,
         "--gapopen", "1",
         "--gapext", "1",
         "--match", "1",
         "--mismatch", "-1",
         "| awk '$1==\"H\" {print $9,$10}'"
      ),
      intern = TRUE
   )
   stopifnot(attr(uc, "status") == 0)
   readr::read_delim(
      I(uc),
      col_names = c("ASV", "cluster"),
      delim = " ",
      col_types = "cc"
   )
}

vsearch_usearch_global_closed_ref <- function(query, ref, threshold, ...) {
   out <- tibble::tibble(ASV = character(0), cluster = character(0))
   while(sequence_size(query) > 0 && sequence_size(ref) > 0) {
      result <- vsearch_usearch_global(query, ref, threshold, ...)
      if (nrow(out) > 0) {
         result <- dplyr::left_join(
            result,
            out,
            by = c("cluster" = "ASV"),
            suffix = c(".orig", "")
         ) %>%
            dplyr::select(ASV, cluster)
      }
      out <- dplyr::bind_rows(out, result)
      ref <- select_sequence(query, result$ASV)
      query <- select_sequence(query, result$ASV, negate = TRUE)
   }
   out
}

vsearch_cluster_smallmem <- function(seq, threshold = 1, ncpu = local_cpus()) {
   tout <- tempfile("data", fileext = "fasta")
   tin <- tempfile("data", fileext = "uc")
   on.exit(unlink(tout, force = TRUE))
   write_sequence(seq, tout)
   uc = system(
      paste(
         "vsearch",
         "--cluster_smallmem", tout,
         "--usersort",
         "--id", threshold,
         "--uc -",
         "--threads", ncpu,
         "| awk '$1==\"H\" {print $9,$10}'"
      ),
      intern = TRUE
   )
   stopifnot(attr(uc, "status") == 0)
   readr::read_delim(
      I(uc),
      col_names = c("query", "hit"),
      delim = " ",
      col_types = "cc"
   )
}

vsearch_derep <- function(infile, outfile, ncpus = local_cpus(), sizein = FALSE) {
   result = system2(
      "vsearch",
      c(
         "--derep_fulllength", infile,
         if (isTRUE(sizein)) "--sizein",
         "--output", outfile,
         "--threads", ncpu
      )
   )
   stopifnot(result == 0)
   return(outfile)
}

# Taxonomy figure
taxon_plot <- function(
   .data, # the data
   rank, # what "rank" (column in the data) we are plotting
   ..., # conditions to filter which rows we care about
   y = reads, # should be reads or SHs
   x = Type, # what groups do we want to divide by
   weight = if ("weight" %in% names(.data)) "weight" else "1", # column to weight the read and ASV counts by
   cutoff = NULL, # groups which represent less than this fraction in all types are grouped together as "other"
   cutoff_type = c("single", "either", "both"),
   data_only = FALSE # just return the data
) {
   rank <- rlang::enquo(rank)
   rank_label <- rlang::as_label(rank)
   y <- rlang::enquo(y)
   x <- rlang::enquo(x)
   weight <- rlang::parse_expr(weight)
   ranks <- c("Phylum", "Class", "Order", "Family", "Genus")
   cutoff_type <- match.arg(cutoff_type)
   .data <- .data %>%
      dplyr::group_by(!!x) %>%
      dplyr::mutate(
         ASVs = dplyr::n_distinct(ASV),
         reads = nread/sum(nread * !!weight)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::filter(...)
   if (rank_label %in% ranks) {
      .data <- .data %>%
         dplyr::arrange_at(ranks) %>%
         dplyr::mutate_at(
            ranks,
            ~ factor(
               .,
               levels = c(NA, "other",
                          purrr::discard(unique(as.character(.)), is.na)),
               exclude = "NULL"
            )
         )
   } else if (!is.factor(dplyr::pull(.data, !!rank))) {
      .data <- .data %>%
         dplyr::mutate_at(
            rank_label,
            ~ factor(
               .,
               levels = c(NA, "other",
                          purrr::discard(unique(as.character(.)), is.na)),
               exclude = "NULL"
            )
         )
   }
   .data <- .data %>%
      dplyr::group_by(!!x, !!rank) %>%
      dplyr::summarize(reads = sum(!!weight * reads), ASVs = sum(unique(data.frame(ASV, w = !!weight))$w)/max(ASVs)) %>%
      dplyr::ungroup()

   if (data_only) return(.data)

   if (!is.null(cutoff)) {
      prelevels <- levels(dplyr::pull(.data, !!rank))
      cutoff_fun <- switch(
         cutoff_type,
         single = function(x) all(dplyr::pull(x, !!y) < cutoff),
         either = function(x) all(x$reads < cutoff, x$SHs < cutoff),
         both = function(x) all(x$reads < cutoff) | all(x$SHs < cutoff)
      )
      .data <- dplyr::group_by(.data, !!rank) %>%
         dplyr::group_map(

            ~ if (cutoff_fun(.x)) {
               dplyr::mutate(
                  .x,
                  !!rank := factor(
                     "other",
                     levels = levels(!!rank)
                  ),
                  exclude = NULL)
            } else {
               .x
            },
            .keep = TRUE
         ) %>%
         dplyr::bind_rows() %>%
         dplyr::group_by(!!x, !!rank) %>%
         dplyr::summarize(reads = sum(reads), SHs = sum(SHs)) %>%
         dplyr::ungroup() %>%
         dplyr::mutate(!!rank := factor(!!rank, levels = prelevels, exclude = NULL))
   }

   .data <- dplyr::mutate(.data, !!rank := forcats::fct_drop(!!rank))
   vals <- levels(dplyr::pull(.data, !!rank))
   # if ("other" %in% vals) vals <- c("other", vals) %>% magrittr::extract(!duplicated(.))
   # if (any(is.na(vals))) vals <- c(NA, vals) %>% magrittr::extract(!duplicated(.))

   if (rank_label == stringr::str_to_lower(rank_label)) rank_label <- stringr::str_to_title(rank_label)
   y_label <- rlang::as_label(y)
   if (y_label == stringr::str_to_lower(y_label)) y_label <- stringr::str_to_title(y_label)
   ggplot(.data, aes(x = !!x, y = !!y, fill = !!rank)) +
      geom_bar(position = position_stack(reverse = TRUE),
               stat = "identity", color = "white", size = 0.2) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            strip.background = element_blank(),
            panel.spacing = unit(3, "pt")) +
      hues::scale_fill_iwanthue(
         # type = "qual",
         # palette = 2,
         cmin = 20, cmax = 150,
         lmin = 10, lmax = 80,
         breaks = vals,
         labels = tidyr::replace_na(as.character(vals), "unidentified"),
         na.value = "gray50",
         name = rank_label,
         guide = guide_legend(ncol = 4, byrow = TRUE)
      ) +
      ylab(paste("Fraction of", y_label)) +
      coord_flip() +
      theme_bw() +
      theme(legend.position = "bottom") +
      xlab(NULL)

}

#' Update the list of unknowns in each parent taxon using the results of
#' closed-ref clustering at the parent rank.
#'
#' @param superclust (`data.frame` with columns "ASV" and "cluster") clustering
#' results from the parent rank, giving the previously 'unknown' ASV ids and
#' the id of the sequence they form a cluster with.
#' @param supertax (`data.frame` with columns "ASV" and "taxon") taxonomy of
#' the 'known' ASV sequences at the parent rank
#' @param unknowns (named `list`, where names are taxa and values are ASV ids)
#' unknown sequences, grouped by taxonomy at the parent rank. Sequences where
#' the parent rank is also unknown are in the item named `"_NA_"`.
#'
#' @return (named `list`, where names are taxa and values are ASV ids)
#' unknown sequences, grouped by taxonomy at the parent rank. Sequences where
#' the parent rank is also unknown are in the item named `"_NA_"`.
merge_unknowns <- function(superclust, supertax, unknowns) {
   out <- if (is.null(superclust)) {
      list()
      } else {
         dplyr::left_join(superclust, supertax, by = c("cluster" = "ASV")) %$%
            split(ASV, taxon)
      }

   for (n in setdiff(names(unknowns), "_NA_")) {
      out[[n]] <- c(out[[n]], unknowns[[n]])
   }
   out[["_NA_"]] <- c(out[["_NA_"]], setdiff(unknowns[["_NA_"]], unlist(out)))
   stopifnot(all.equal(length(unlist(out)), length(unlist(unknowns))))
   out
}

#' Calculate clustering thresholds for each taxon, falling back to its ancestor
#' taxa as necessary
#'
#' @param rank (`character` string) the rank within which clustering will be
#' performed
#' @param conf_level (`character` string) as in `fmeasure_optima`
#' @param taxon_table (`data.frame`) taxonomy table; column "ASV" gives the
#' sequence ID, and columns "kingdom" to "species" give the taxonomy at each
#' rank.
#' @param fmeasure_optima (`data.frame`) optimum clustering thresholds within
#' various taxa; column "rank" gives the rank which is approximated by
#' clustering; "superrank" gives the rank of the taxon within which the
#' clustering threshold was optimized; "supertaxon" gives that taxon name;
#' "conf_level" gives a string description of the confidence level threshold for
#' taxonomic assignments; "threshold" gives the optimum clustering threshold;
#' "f_measure" gives the F measure at the optimum threshold.
#' @param default (`character string`) default taxon to define threshold to use
#' when taxonomy is unknown. default: "Fungi"
#'
#'
#' @return (named `numeric`, where names are taxa and values are clustering
#' thresholds)
#' unknown sequences, grouped by taxonomy at the parent rank. Sequences where
#' the parent rank is also unknown are in the item named `"_NA_"`.
calc_taxon_thresholds <- function(rank, conf_level, taxon_table,
                                  fmeasure_optima, default = "Fungi") {
   rank_name <- rlang::sym(rank)
   dplyr::select(taxon_table, kingdom:!!rank_name) %>%
      dplyr::filter(!is.na(!!rank_name)) %>%
      unique() %>%
      purrr::reduce(
         c(superranks(rank), rank),
         function(thresholds, r) {
            dplyr::left_join(
               thresholds,
               dplyr::filter(
                  fmeasure_optima,
                  rank == subranks(!!rank)[1],
                  superrank == r,
                  conf_level == !!conf_level
               ) %>%
                  dplyr::select(
                     !!r := supertaxon,
                     !!paste0("threshold_", r) := threshold
                  ),
               by = r
            )
         },
         .init = .
      ) %>%
      dplyr::transmute(
         !!rank_name := !!rank_name,
         threshold = {.} %>%
            dplyr::select(dplyr::starts_with("threshold")) %>%
            rev() %>%
            do.call(dplyr::coalesce, .)
      ) %>%
      tibble::deframe() %>%
      c("_NA_" = dplyr::filter(
         fmeasure_optima,
         rank == subranks(!!rank)[1],
         supertaxon == default,
         conf_level == !!conf_level
      )$threshold)
}

run_protax <- function(seqs, outdir, ncpu = local_cpus()) {
   if (dir.exists(outdir)) unlink(outdir, recursive = TRUE)
   dir.create(outdir, recursive = TRUE)
   write_sequence(seqs, file.path(outdir, "all.fa"))
   status <- system2(
      "scripts/runprotax",
      c(outdir, ncpu)
   )
   stopifnot(status == 0)
   list.files(outdir, full.names = TRUE)
}

read_gsmc_names <- function(x, pos) {
   x <- x[startsWith(x, ">")]
   tidyr::extract(
      tibble::tibble(x = x),
      col = "x",
      into = c("seq_id", "sample", "OTU", "n_read"),
      regex = ">([a-f0-9]+);(?:sample=([^;]+);OTU=([a-f0-9]+);)?size=([0-9]+)",
      convert = TRUE
   )
}

summarize_by_rank <- function(rank, superrank, data) {
   dplyr::group_by(data, dplyr::across(superrank)) %>%
   dplyr::summarize(
      superrank = superrank,
      rank = rank,
      n_taxa = dplyr::n_distinct(.data[[rank]]),
      n_seq = dplyr::n_distinct(seq_id),
      seq_id = list(seq_id)
   ) %>%
   dplyr::rename(supertaxon = !!superrank)
}

assemble_GSSP_GSMc_compare_data <- function(combined_OTU, GSSP_GSMc_taxon) {
   GSMc_otus <- sum(combined_OTU$nread_GSMc > 0)
   GSSP_otus <- sum(combined_OTU$nread_GSSP > 0)
   GSMc_reads <- sum(combined_OTU$nread_GSMc)
   GSSP_reads <- sum(combined_OTU$nread_GSSP)

   dplyr::left_join(
      combined_OTU %>%
         dplyr::mutate(
            fread_GSMc = nread_GSMc/GSMc_reads,
            fread_GSSP = nread_GSSP/GSSP_reads
         ),
      dplyr::bind_rows(
         GSSP_GSMc_taxon,
         tibble::tibble(
            rank = rank2factor("kingdom"),
            OTU = unique(GSSP_GSMc_taxon$OTU),
            parent_taxonomy = NA_character_,
            taxon = "Fungi",
            prob = 1
         )
      ),
      by = c("cluster" = "OTU")
   ) %>%
      dplyr::mutate(
         taxon = ifelse(taxon=="Pezizomycotina", "Ascomycota", taxon),
         taxon = sub("dummy_(class|order|family)_.*", "incertae sedis", taxon),
         parent_taxonomy = sub("Fungi,Pezizomycotina", "Fungi,Ascomycota", parent_taxonomy),
         parent_taxonomy = gsub("dummy_(class|order|family)_[^,]+", "incertae sedis", parent_taxonomy)
      ) %>%
      dplyr::group_by(parent_taxonomy, taxon, rank) %>%
      dplyr::summarize(
         nOTU_GSSP = sum(prob*(nread_GSSP>0)),
         nOTU_GSMc = sum(prob*(nread_GSMc>0)),
         nOTU = sum(prob*(nread_GSSP > 0 | nread_GSMc > 0)),
         nread_GSSP = sum(nread_GSSP*prob),
         nread_GSMc = sum(nread_GSMc*prob),
         fread_GSSP = sum(fread_GSSP*prob),
         fread_GSMc = sum(fread_GSMc*prob),
         .groups = "drop"
      ) %>%
      dplyr::filter(!is.na(taxon)) %>%
      dplyr::mutate(
         fOTU_GSSP = nOTU_GSSP/GSSP_otus,
         fOTU_GSMc = nOTU_GSMc/GSMc_otus,
         fOTU = (fOTU_GSSP + fOTU_GSMc)/2,
         fread = (fread_GSSP + fread_GSMc)/2,
         ggof = (fOTU_GSSP - fOTU_GSMc)/(fOTU_GSSP + fOTU_GSMc),
         ggon = (nOTU_GSSP - nOTU_GSMc)/(nOTU_GSSP + nOTU_GSMc),
         ggrf = (fread_GSSP - fread_GSMc)/(fread_GSSP + fread_GSMc)
      )
}

# convert a list of data to the XML format to be sent to KronaTools
xml_format <- function(data_format) {
   lapply(data_format, vapply, sprintf, "", fmt = "<val>{%s}</val>") %>%
      vapply(paste, "", collapse = ",") %>%
      purrr::imap_chr(sprintf, fmt="<%2$s>%1$s</%2$s>") %>%
      paste(collapse = "\n")
}

krona_xml_nodes <- function(
   data,
   .rank,
   maxrank = rank2factor("species"),
   outfile,
   pre = NULL,
   post = NULL,
   taxonomy = "Fungi",
   node_data_format = NULL,
   node_xml_format = xml_format(node_data_format),
   ...
) {
   if (is.character(.rank)) .rank <- rank2factor(.rank)
   if (!methods::is(outfile, "connection")) {
      outfile <- file(outfile, open = "w")
      on.exit(close(outfile))
   }
   my_data <- data
   if (!is.null(taxonomy)) {
      my_data <- dplyr::filter(data, startsWith(parent_taxonomy, taxonomy))
   }
   xml <- dplyr::filter(my_data, rank == .rank) %>%
      dplyr::transmute(
         taxon = taxon,
         taxonomy = ifelse(is.na(parent_taxonomy), taxon, paste(parent_taxonomy, taxon, sep = ",")),
         pre = glue::glue(
            '<node name="{taxon}">',
            node_xml_format,
            .sep = "\n"
         ),
         post = "</node>"
      )
   if (!is.null(pre)) {
      writeLines(pre, outfile)
   }
   if (.rank == maxrank) {
      writeLines(paste(xml$pre, xml$post, sep = "\n"), outfile)
   } else {
      purrr::pwalk(
         xml,
         krona_xml_nodes,
         data = my_data,
         .rank = subranks(.rank)[1],
         maxrank = maxrank,
         outfile = outfile,
         GSMc_otus = GSMc_otus,
         GSSP_otus = GSSP_otus,
         GSMc_reads = GSMc_reads,
         GSSP_reads = GSSP_reads,
         node_xml_format = node_xml_format
      )
   }
   if (!is.null(post)) writeLines(post, outfile)
}

add_rank <- function(data, .rank, fulldata) {
   dplyr::left_join(
      data,
      dplyr::filter(fulldata, rank == .rank),
      by = c("OTU", "parent_taxonomy")
   ) %>%
      dplyr::bind_rows(
         dplyr::group_by(., OTU, parent_taxonomy, prob.x, rank) %>%
            dplyr::summarize(prob.y = sum(prob.y), .groups = "drop") %>%
            dplyr::filter(prob.y < prob.x) %>%
            dplyr::mutate(prob.y = prob.x - prob.y, taxon = "other")
      ) %>%
      dplyr::transmute(
         OTU = OTU,
         parent_taxonomy = dplyr::coalesce(
            stringr::str_c(parent_taxonomy, taxon, sep = ","),
            parent_taxonomy
         ),
         prob = dplyr::coalesce(prob.y, prob.x)
      ) %>%
      dplyr::group_by(OTU) %>%
      dplyr::mutate(prob = prob/sum(prob))
}

