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
         protax_refseq_index,
         Biostrings::fasta.index(protax_refseq_file),
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
      tar_file_fast(
         unite_files,
         list.files("data", "UNITE.+\\.fasta\\.gz", full.names = TRUE),
         deployment = "main"
      ),
      tar_file_fast(
         find_unite_seqs_script,
         "scripts/find_unite_seqs.awk",
         deployment = "main"
      ),
      tar_file_fast(
         protax_refs_full,
         {
            outfile <- "data/protax_refs_full.fasta.gz"
            command <- paste(
               "zcat", paste(shQuote(unite_files), collapse = " "), "|",
               "awk", "-f", shQuote(find_unite_seqs_script),
               shQuote(protax_refseq_file), "-", "|",
               "gzip", "-c", "-",
               ">", shQuote(outfile)
            )
            result <- system(command, timeout = 120)
            stopifnot(result == 0)
            outfile
         },
         deployment = "main"
      ),
      tar_file_fast(
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
      tar_file_fast(
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
      tar_file_fast(
         protax_refs_itsx,
         {
            unzipped <- tempfile(fileext=".fasta")
            system(paste("gunzip -c", protax_refs_full, ">", unzipped))
            out_root <- "data/protax_refs"
            rITSx::itsx(
               in_file = unzipped,
               out_root = out_root,
               fasta = FALSE,
               preserve = TRUE,
               save_regions = "ITS2",
               complement = FALSE,
               positions = FALSE,
               table = FALSE,
               summary = FALSE,
               graphical = FALSE,
               detailed_results = FALSE,
               only_full = TRUE,
               cpu = local_cpus()
            )
            f <- paste0(out_root, ".ITS2.fasta")
            system2("gzip", f)
            paste0(f, ".gz")
         },
         deployment = "worker"
      ),
      tar_fst_tbl(
         protax_refs_itsx_index,
         Biostrings::fasta.index(protax_refs_itsx),
         deployment = "main"
      ),
      tar_fst_tbl(
         protax_refs_itsx_match_index,
         dplyr::semi_join(protax_refs_itsx_index, protax_refs_trim2_index, by = "desc"),
         deployment = "main"
      ),

      tar_file_fast(
         rfam_5.8S,
         {
            cmfile <- "data/RF00002.cm"
            rfaRm::rfamCovarianceModel("RF00002", cmfile)
            cmfile
         },
         deployment = "main"
      ),

      tar_fst_tbl(
         protax_refs_5.8S_positions,
         inferrnal::cmsearch(
            cm = rfam_5.8S,
            seq = protax_refs_trim2,
            cpu = local_cpus(),
            quiet = TRUE
         ) |>
            dplyr::rename(desc = target_name),
         deployment = "worker"
      ),

      tar_file_fast(
         protax_refs_5.8S,
         Biostrings::readDNAStringSet(protax_refs_trim2) |>
            magrittr::extract(protax_refs_5.8S_positions$desc) |>
            Biostrings::subseq(
               start = protax_refs_5.8S_positions$seq_from,
               end = protax_refs_5.8S_positions$seq_to
            ) |>
            write_and_return_file("data/protax_refs_5.8S.fasta"),
         deployment = "worker"
      ),
      tar_fst_tbl(
         protax_refs_match_index,
         dplyr::semi_join(
            protax_refs_itsx_index,
            protax_refs_trim2_index,
            by = "desc"
         ) |>
            dplyr::semi_join(
               protax_refs_5.8S_positions,
               by = "desc"
            ),
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
      tar_file_fast(
         protax_reftax_file,
         file.path(protax_modeldir, sprintf("ref.tax%d", threshold_meta$rank_int)),
         pattern = map(threshold_meta),
         deployment = "main"
      ),
      tar_map(
         values = tibble::tibble(
            refseq_file = rlang::syms(c(
               "protax_refs_itsx", # full ITS2 only (cut by ITSx)
               "protax_refs_itsx", # same
               "protax_refs_5.8S", # (partial) 5.8S
               "protax_refs_5.8S", # same
               "protax_refs_trim2",# only sequences where ITS3ITS4 amplicon region could be found
               "protax_refs_trim2" # same
            )),
            refseq_index = rlang::syms(c(
               "protax_refs_itsx_index", # all ITSx sequences
               "protax_refs_5.8S_match_index", # sequences found by ITSx, ITS3ITS4 primer trimming, and search for the 5.8S CM
               "protax_refs_5.8S_index", # all 5.8S sequences
               "protax_refs_5.8S_match_index", # sequences found by ITSx, ITS3ITS4 primer trimming, and search for the 5.8S CM
               "protax_refs_trim2_index", # ITS3ITS4 sequences
               "protax_refs_5.8S_match_index" # sequences found by ITSx, ITS3ITS4 primer trimming, and search for the 5.8S CM
            )),
            refset_name = c("ITSx_all", "ITSx_match", "5.8S_all", "5.8S_match", "ITS3ITS4_match", "ITS3ITS4")
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
            optimotu::seq_cluster_usearch(
               seq = refseq_file,
               threshold_config = optimotu::threshold_uniform(
                  from = 0,
                  to = 0.4,
                  by = 0.001,
                  thresh_names = as.character(1000 - 0:400)
               ),
               clust_config = optimotu::clust_tree(),
               parallel_config = optimotu::parallel_concurrent(local_cpus()%/%2.5),
               which = testset_select$seq_id,
               usearch = "bin/usearch",
               usearch_ncpu = local_cpus()
            ),
            iteration = "list",
            deployment = "worker"
         ),
         tar_target(
            threshold_ntaxa,
            lapply(threshold_testset, apply, 1, dplyr::n_distinct) |>
               lapply(tibble::enframe, name = "threshold", value = "ntaxa") |>
               tibble::add_column(testset_rowwise, ntaxa = _) |>
               dplyr::select(supertaxon, superrank, rank, ntaxa) |>
               tidyr::unnest(ntaxa),
            deployment = "worker"
         ),
         tar_fst_tbl(
            cluster_metrics,
            purrr::map_dfr(
               seq_along(threshold_testset),
               function(i) {
                  purrr::map_dfc(
                     .x = list(
                        optimotu::confusion_matrix,
                        optimotu::adjusted_mutual_information,
                        mFM = optimotu::fmeasure
                     ),
                     .f = purrr::exec,
                     k = threshold_testset[[i]],
                     c = testset_rowwise$true_taxa[[i]],
                     local_cpus()
                  ) |>
                     tibble::remove_rownames() %>%
                     dplyr::mutate(
                        MCC = optimotu::matthews_correlation_coefficient(.),
                        RI = optimotu::rand_index(.),
                        ARI = optimotu::adjusted_rand_index(.),
                        FMI = optimotu::fowlkes_mallow_index(.),
                        threshold = (1000 - 0:400)/10,
                        rank = testset_rowwise$rank[i],
                        superrank = testset_rowwise$superrank[i],
                        supertaxon = testset_rowwise$supertaxon[i]
                     )
               }
            ),
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
         tar_file_fast(
            optima_file,
            write_and_return_file(
               optima, file.path("data", sprintf("protaxFungi_%s_thresholds.tsv", refset_name)),
               "tsv"),
            deployment = "main"
         )
      )
   )
