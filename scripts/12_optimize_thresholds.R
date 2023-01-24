# test the performance of different clustering thresholds for replicating
# "true" taxonomic groups.

protax_modeldir <- file.path("protaxFungi", "addedmodel")

   threshold_test_plan <- list(
      tar_file(
         protax_refseq_file,
         file.path(protax_modeldir, "its2.fa")
      ),
      tar_fst_tbl(
         threshold_meta,
         tibble::tibble(
            rank_int = seq_along(TAXRANKS)[-1],
            rank = TAXRANKS[rank_int]
         )
      ),
      tar_file(
         protax_reftax_file,
         file.path(protax_modeldir, sprintf("ref.tax%d", threshold_meta$rank_int)),
         pattern = map(threshold_meta)
      ),

      tar_fst_tbl(
         protax_reftax,
         readr::read_tsv(
            file = protax_reftax_file,
            col_names = c("seq_id", "taxonomy"),
            col_types = "cc"
         ) %>%
            tidyr::separate(taxonomy, into = TAXRANKS[1:threshold_meta$rank_int], sep = ","),
         pattern = map(threshold_meta, protax_reftax_file),
         iteration = "list"
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
         pattern = map(threshold_meta, protax_reftax)
      ),
      tar_target(
         threshold_testset,
         usearch_singlelink(
            seq = protax_refseq_file,
            thresh_min = 0,
            thresh_max = 0.4,
            thresh_step = 0.001,
            thresh_names = as.character(1000 - 0:400),
            which = testset_select$seq_id,
            usearch = "sh_matching_pub/programs/usearch"
         ),
         iteration = "list"
      ),
      tar_target(
         threshold_ntaxa,
         apply(threshold_testset, 1, dplyr::n_distinct),
         pattern = map(threshold_testset),
         iteration = "list"
      ),
      tar_fst_tbl(
         fmeasures,
         protax_reftax[[match(testset_select$rank, TAXRANKS) - 1L]] %>%
            dplyr::select(seq_id, testset_select$rank) %>%
            dplyr::mutate(i = match(seq_id, colnames(threshold_testset))) %>%
            dplyr::filter(!is.na(i)) %>%
            {split(.[[3]], .[[2]])} %>%
            optimotu::fmeasure(
               k = apply(threshold_testset, 1, split, x = seq_len(ncol(threshold_testset))),
               c = .,
               ncpu = local_cpus()
            ) %>%
            tibble::tibble(
               f_measure = .,
               threshold = (1000 - 0:400)/10,
               rank = testset_select$rank,
               superrank = testset_select$superrank,
               supertaxon = testset_select$supertaxon
            ),
         pattern = map(testset_select, threshold_testset)
      ),
      tar_fst_tbl(
         fmeasure_optima,
         fmeasures %>%
            dplyr::group_by(rank, superrank, supertaxon) %>%
            dplyr::arrange(dplyr::desc(f_measure)) %>%
            dplyr::summarize(
               threshold = threshold[which.max(f_measure)],
               f_measure = max(f_measure),
               .groups = "drop"
            )
      ),
      tar_file(
         fmeasure_file,
         write_and_return_file(
            fmeasure_optima, file.path("data", "GSSP_thresholds.tsv"),
            "tsv")
      )
   )
