# test the performance of different clustering thresholds for replicating
# "true" taxonomic groups.
library(tarchetypes)

protax_modeldir <- file.path("protaxFungi", "addedmodel")

   threshold_test_plan <- list(
      tar_file(
         protax_refseq_file,
         file.path(protax_modeldir, "its2.fa"),
         deployment = "main"
      ),
      tar_fst_tbl(
         protax_refseq_ids,
         tibble::tibble(
            id = names(Biostrings::fasta.seqlengths(protax_refseq_file))
         ) %>%
            tidyr::separate(
               id,
               c("SH", "accno"),
               sep = "_",
               remove = FALSE,
               extra = "merge"),
         deployment = "main"
      ),
      tar_file(
         unite_files,
         list.files("data", "UNITE.+\\.fasta\\.gz", full.names = TRUE),
         deployment = "main"
      ),
      tar_file(
         find_unite_seqs_script,
         "scripts/find_unite_seqs.awk",
         deployment = "main"
      ),
      tar_file(
         protax_refs_full,
         {
            outfile <- "data/protax_refs_full.fasta.gz"
            command <- paste(
               "zcat", paste(shQuote(unite_files), collapse = " "), "|",
               "awk", "-f", shQuote(find_unite_seqs_script),
               shQuote(protax_refseq_file), "-",
               "|", "gzip", "-c", "-",
               paste0(">", shQuote(outfile))
            )
            result <- system(command, timeout = 120)
            stopifnot(result == 0)
            outfile
         },
         deployment = "main"
      ),
      tar_file(
         protax_refs_trim1,
         cutadapt_filter_trim(
            protax_refs_full,
            primer = "GCATCGATGAAGAACGCAGC;required...ACCCGCTGAACTTAAGCATATCAATAAGCGGAGGA;required",
            trim = "data/protax_refs_trim1.fasta.gz",
            discard_untrimmed = TRUE,
            action = "retain",
            min_length = 185,
            max_err = 0.3,
            min_overlap = 8
         ),
         deployment = "main"
      ),
      tar_file(
         protax_refs_trim2,
         cutadapt_filter_trim(
            protax_refs_trim1,
            primer = "XGCATCGATGAAGAACGCAGC;required;o=8...GCATATCAATAAGCGGAGGAX;optional",
            discard_untrimmed = TRUE,
            trim = "data/protax_refs_trim2.fasta",
            min_length = 145,
            max_err = 0.3
         ),
         deployment = "main"
      ),
      tar_fst_tbl(
         protax_refs_trim2_index,
         Biostrings::fasta.index(protax_refs_trim2),
         deployment = "main"
      ),
      tar_fst_tbl(
         threshold_meta,
         tibble::tibble(
            rank_int = seq_along(TAXRANKS)[-1],
            rank = TAXRANKS[rank_int]
         ),
         deployment = "main"
      ),
      tar_file(
         protax_reftax_file,
         file.path(protax_modeldir, sprintf("ref.tax%d", threshold_meta$rank_int)),
         pattern = map(threshold_meta),
         deployment = "main"
      ),
      tar_map(
         values = tibble::tibble(
            refseq_file = rlang::syms(c(
               protax_refseq_file, # full database, ITS2 only (cut by ITSx)
               protax_refseq_file, # same
               protax_refs_trim2 # only sequences where ITS3ITS4 amplicon region could be found
            )),
            refseq_index = rlang::syms(c(
               protax_refseq_index,
               protax_refs_trim2_index,
               protax_refs_trim2_index
            )),
            refset_name = c("ITSx_all", "ITSx_match", "ITS3ITS4")
         ),
         names = refset_name,

         tar_fst_tbl(
            protax_reftax,
            readr::read_tsv(
               file = protax_reftax_file,
               col_names = c("seq_id", "taxonomy"),
               col_types = "cc"
            ) %>%
               dplyr::inner_join(refseq_index, by = c("seq_id" = "desc")) %>%
               tidyr::separate(taxonomy, into = TAXRANKS[1:threshold_meta$rank_int], sep = ","),
            pattern = map(threshold_meta, protax_reftax_file),
            iteration = "list",
            deployment = "main"
         ),
         tar_target(
            testset_select,
            purrr::map_dfr(
               superranks(threshold_meta$rank),
               summarize_by_rank,
               data = protax_reftax,
               rank = threshold_meta$rank
            ) %>%
               dplyr::filter(n_taxa >= 5 | superrank == "kingdom", n_seq >= 10),
            tidy_eval = FALSE,
            pattern = map(threshold_meta, protax_reftax),
            deployment = "main"
         ),
         tar_target(
            testset_rowwise,
            testset_select,
            deployment = "main"
         ),
         tar_target(
            threshold_testset,
            usearch_singlelink(
               seq = refseq_file,
               thresh_min = 0,
               thresh_max = 0.4,
               thresh_step = 0.001,
               thresh_names = as.character(1000 - 0:400),
               which = testset_select$seq_id,
               usearch = "bin/usearch"
            ),
            iteration = "list",
            deployment = "worker"
         ),
         tar_target(
            threshold_ntaxa,
            apply(threshold_testset, 1, dplyr::n_distinct),
            pattern = map(threshold_testset),
            iteration = "list",
            deployment = "worker"
         ),
         tar_fst_tbl(
            cluster_metrics,
            purrr::map(
               list(
                  optimotu::confusion_matrix,
                  optimotu::adjusted_mutual_information,
                  mFM = optimotu::fmeasure
               ),
               k = threshold_testset,
               c = testset_rowwise$true_taxa,
               ncpu = local_cpus()
            ) %>%
               dplyr::mutate(
                  MCC = optimotu::matthews_correlation_coefficient(.),
                  RI = optimotu::rand_index(.),
                  ARI = optimotu::adjusted_rand_index(.),
                  FMI = optimotu::fowlkes_mallow_index(.),
                  threshold = (1000 - 0:400)/10,
                  rank = testset_rowwise$rank,
                  superrank = testset_rowwise$superrank,
                  supertaxon = testset_rowwise$supertaxon
               ),
            pattern = map(testset_rowwise, threshold_testset),
            deployment = "worker"
         ),
         tar_fst_tbl(
            optima,
            cluster_metrics %>%
               tidyr::pivot_longer(
                  -c(threshold, rank, superrank, supertaxon),
                  names_to = "metric", values_to = "score"
               ) %>%
               dplyr::group_by(rank, superrank, supertaxon, metric) %>%
               dplyr::arrange(dplyr::desc(score)) %>%
               dplyr::summarize(
                  threshold = threshold[which.max(score)],
                  score = max(score),
                  .groups = "drop"
               ),
            deployment = "main"
         ),
         tar_file(
            optima_file,
            write_and_return_file(
               optima, file.path("data", sprintf("protaxFungi_%s_thresholds.tsv", refset_name)),
               "tsv"),
            deployment = "main"
         )
      )
   )
