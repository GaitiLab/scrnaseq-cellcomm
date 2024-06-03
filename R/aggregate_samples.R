#' @title Aggregate samples w/ ranked interactions
#' @description Combine p-values/scores by condition
#' @param input_file Path to file with interaction results for all samples (401_combine_samples/401_samples_interactions_agg_rank.rds)
#' @param output_dir output directory for saving output (default = '.')
#' @export
#' @importFrom dplyr %>%
aggregate_samples <- function(input_file, output_dir = ".") {
    message("Load interactions and combine p-values and scores...")
    obj <- readRDS(input_file) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
            LIANA_score = ifelse(is.na(LIANA_score), 0, LIANA_score),
            CellPhoneDB_score = ifelse(is.na(CellPhoneDB_score), 0, CellPhoneDB_score),
            CellChat_score = ifelse(is.na(CellChat_score), 0, CellChat_score)
        ) %>%
        dplyr::select(
            !!dplyr::sym(condition_var),
            complex_interaction,
            pval,
            LIANA_score,
            CellPhoneDB_score,
            CellChat_score,
            source_target
        ) %>%
        dplyr::group_by(
            !!dplyr::sym(condition_var),
            source_target,
            complex_interaction
        ) %>%
        # Combine p-values and interaction scores (CellChat, LIANA, CellphoneDB) across samples
        dplyr::mutate(
            pval = survcomp::combine.test(pval),
            LIANA_score = mean(LIANA_score),
            CellPhoneDB_score = mean(CellPhoneDB_score),
            CellChat_score = mean(CellChat_score)
        )

    message("Save results...")
    saveRDS(
        obj,
        file = glue::glue("{output_dir}/402b_aggregation_samples.rds")
    )
    message("Finished!")
}
