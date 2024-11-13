#' @title Aggregate samples w/ ranked interactions
#' @description Combine p-values/scores by condition
#' @param input_file Path to file with interaction results for all samples (401_combine_samples/401_samples_interactions_agg_rank.rds)
#' @param output_dir output directory for saving output (default = '.')
#' @param condition_var Variabile containing the 'condition' (group/category) of the samples (default = "Condition_dummy")
#' @export
#' @importFrom dplyr %>%
aggregate_samples <- function(input_file, condition_var = "Condition_dummy", output_dir = ".") {
    message("Load interactions and combine p-values and scores...")
    obj <- readRDS(input_file) %>%
        dplyr::ungroup()

    if (condition_var == "Condition_dummy") {
        obj <- obj %>% dplyr::mutate(Condition_dummy = "group1")
    }

    obj <- obj %>%
        dplyr::select(
            !!dplyr::sym(condition_var),
            complex_interaction,
            pval,
            LIANA_score,
            CellPhoneDB_score,
            CellChat_score,
            Cell2Cell_score,
            source_target
        ) %>%
        dplyr::group_by(
            !!dplyr::sym(condition_var),
            source_target,
            complex_interaction
        ) %>%
        # Combine p-values and interaction scores (CellChat, LIANA, CellphoneDB) across condition or patient (if patient = sample, then across samples)
        dplyr::reframe(
            pval = survcomp::combine.test(pval),
            LIANA_score = mean(LIANA_score, na.rm = TRUE),
            CellPhoneDB_score = mean(CellPhoneDB_score, na.rm = TRUE),
            CellChat_score = mean(CellChat_score, na.rm = TRUE),
            Cell2Cell_score = mean(Cell2Cell_score, na.rm = TRUE)
        ) %>%
        dplyr::ungroup()

    obj$pval_adj <- p.adjust(obj$pval, "BH")

    message("Save results...")
    saveRDS(
        obj,
        file = glue::glue("{output_dir}/402b_aggregation_samples.rds")
    )
    message("Finished!")
}
