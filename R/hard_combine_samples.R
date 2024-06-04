#' @title Combine samples
#' @description hard-combine results per sample.
#' @param input_dir directory with files ending in '__interactions_mvoted.rds', '__signif_interactions.rds' and '*__interactions_agg_rank.rds',
#' @param metadata path to file with corresponding metadata (RDS)
#' @param patient_var patient variable name in metadata
#' @param sample_var sample variable name in metadata (default = "Sample")
#' @param condition_var condition variable name in metadata (default = "Condition_dummy")
#' @param output_dir output directory for saving output (default = '.')
#' @export
#' @importFrom dplyr %>%
combine_samples <- function(input_dir, metadata, patient_var, condition_var = "Condition_dummy", sample_var = "Sample", output_dir = ".") {
    # Variables from metadata to add
    cols_oi <- unique(c(sample_var, condition_var, patient_var))

    # Sanity checks
    if (!file.exists(metadata)) {
        stop("Not a valid path")
    }
    all_mvoted <- list.files(input_dir,
        full.names = TRUE, pattern = glue::glue("*__interactions_mvoted.rds")
    )

    all_sign_interactions <- list.files(input_dir,
        full.names = TRUE, pattern = glue::glue("*__signif_interactions.rds")
    )

    interactions_agg_rank <- list.files(input_dir,
        full.names = TRUE, pattern = glue::glue("*__interactions_agg_rank.rds")
    )

    if (length(all_mvoted) == 0) {
        stop("No 'interactions_mvoted.rds' files found")
    }
    if (length(all_sign_interactions) == 0) {
        stop("No 'signif_interactions.rds' files found")
    }
    if (length(interactions_agg_rank) == 0) {
        stop("No 'interactions_agg_rank.rds' files found")
    }

    message("Load metadata...")
    metadata <- readRDS(metadata)

    message("Load & Hard combine 'interactions_mvoted' dataframes...")
    all_mvoted <- do.call(rbind, lapply(all_mvoted, readRDS))

    message("Load & Hard combine 'signif_interactions' dataframes...")
    all_sign_interactions <- do.call(rbind, lapply(all_sign_interactions, readRDS))


    message("Load & Hard combine 'interactions_agg_rank' dataframes...")
    interactions_agg_rank <- do.call(rbind, lapply(interactions_agg_rank, readRDS))


    # Check whether cols_oi are part of metadata
    cols <- colnames(metadata)
    cols_oi <- intersect(cols_oi, cols)
    # Only keep sample-level metadata (not single-cell)
    metadata <- metadata %>%
        dplyr::select(dplyr::all_of(cols_oi)) %>%
        tibble::remove_rownames() %>%
        dplyr::distinct()

    message("Add metadata...")
    all_mvoted <- all_mvoted %>%
        # All post-processed interaction results have the same column names, e.g. "Sample"
        dplyr::left_join(metadata, by = c("Sample" = sample_var)) %>%
        dplyr::ungroup()

    if (sample_var == patient_var) {
        all_mvoted <- all_mvoted %>% dplyr::mutate(Patient = Sample)
    }

    all_sign_interactions <- all_sign_interactions %>%
        # All post-processed interaction results have the same column names, e.g. "Sample"
        dplyr::left_join(metadata, by = c("Sample" = sample_var)) %>%
        dplyr::ungroup()

    interactions_agg_rank <- interactions_agg_rank %>%
        dplyr::left_join(metadata, by = c("Sample" = sample_var)) %>%
        dplyr::ungroup()

    # In case there is no condition given, put all samples into one group.
    if (condition_var == "Condition_dummy") {
        all_mvoted <- all_mvoted %>% dplyr::mutate(Condition_dummy == 1)
        all_sign_interactions <- all_sign_interactions %>% dplyr::mutate(Condition_dummy == 1)
        interactions_agg_rank <- interactions_agg_rank %>% dplyr::mutate(Condition_dummy == 1)
    }

    message("Save output files...")
    saveRDS(all_mvoted,
        file = glue::glue("{output_dir}/401_samples_interactions_mvoted.rds")
    )
    saveRDS(all_sign_interactions,
        file = glue::glue("{output_dir}/401_samples_sign_interactions.rds")
    )
    saveRDS(interactions_agg_rank,
        file = glue::glue("{output_dir}/401_samples_interactions_agg_rank.rds")
    )

    message("Finished!")
}
