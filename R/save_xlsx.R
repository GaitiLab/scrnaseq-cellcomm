#' @title Save interactions as xlsx
#' @param interactions_agg_integration path to '402c_filtering_aggregated_res.rds'
#' @param condition_var Variabile containing the 'condition' (group/category) of the samples (default = "Condition_dummy")
#' @param alpha significance threshold for filtering p-values (default = 0.05)
#' @param output_dir output directory for saving output (default = '.')
#' @param output_name name for output file (default = "interactions_summary")
#' @export
#' @importFrom dplyr %>%
save_as_xlsx <- function(
    interactions_agg_integration,
    condition_var = "Condition_dummy",
    alpha = 0.05,
    output_dir = ".",
    output_name = "interactions_summary") {
    output_filename <- glue::glue("{output_dir}/{output_name}.xlsx")

    if (file.exists(output_filename)) {
        file.remove(output_filename)
    }
    message("Load data...")
    obj <- readRDS(interactions_agg_integration)

    message("Filtering...")
    obj_filtered <- obj %>%
        filter(!is.na(lenient_condition), lenient_condition, pval < alpha) %>%
        dplyr::select(
            !!dplyr::sym(condition_var), source_target, complex_interaction,
            lenient_condition_n_patients, lenient_condition_n_samples,
            lenient_condition_samples, lenient_condition_patients,
            pval,
            # Interaction scores
            LIANA_score, CellPhoneDB_score, CellChat_score
        ) %>%
        dplyr::distinct() %>%
        as.data.frame()

    message("Write data to Excel...")
    xlsx::write.xlsx(
        obj_filtered,
        output_filename,
        col.names = TRUE, row.names = FALSE, append = TRUE
    )
    message("Finished...")
}
