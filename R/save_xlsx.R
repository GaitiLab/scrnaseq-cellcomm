#' @title Save interactions as xlsx
#' @param interactions_agg_integration path to '402c_filtering_aggregated_res.rds'
#' @param condition_var Variabile containing the 'condition' (group/category) of the samples (default = "Condition_dummy")
#' @param alpha significance threshold for filtering p-values (default = 0.05)
#' @param output_dir output directory for saving output (default = '.')
#' @param output_name name for output file (default = "interactions_summary")
#' @param is_stringent only keep stringent interactions (default = FALSE)
#' @export
#' @importFrom dplyr %>%
save_as_xlsx <- function(
    interactions_agg_integration,
    alpha = 0.05,
    output_dir = ".",
    output_name = "interactions_summary",
    is_stringent = FALSE) {
    output_filename <- file.path(output_dir, paste0(output_name, ".xlsx"))

    if (file.exists(output_filename)) {
        file.remove(output_filename)
    }
    message("Load data...")
    obj <- readRDS(interactions_agg_integration)

    message(glue::glue("Filtering (alpha = {alpha})..."))
    if (is_stringent) {
        message("Only keep interactions that passed the 'stringent' filter...")
        obj_filtered <- obj %>%
            filter(pval < alpha, stringent_condition) %>%
            dplyr::distinct() %>%
            as.data.frame()
    } else {
        obj_filtered <- obj %>%
            filter(pval < alpha) %>%
            dplyr::distinct() %>%
            as.data.frame()
    }
    message("Write data to Excel...")
    xlsx::write.xlsx(
        obj_filtered,
        output_filename,
        col.names = TRUE,
        row.names = FALSE,
        append = TRUE
    )
    message("Finished...")
}
