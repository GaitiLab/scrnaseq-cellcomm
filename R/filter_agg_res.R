#' @title Filter agg res by user constraints
#' @param input_dir directory with '402a_filtering_detect_in_multi_samples.rds' and '402b_aggregation_samples.rds', alternatively specify path separately ('interactions_agg_binarized' & 'interactions_agg_continuous') (default = "")
#' @param interactions_agg_binarized path to '402a_filtering_detect_in_multi_samples.rds' (default = "")
#' @param interactions_agg_continuous path to '402b_aggregation_samples.rds' (default = "")
#' @param condition_var Variabile containing the 'condition' (group/category) of the samples (default = "Condition_dummy")
#' @param output_dir output directory for saving output (default = '.')
#' @export
#' @importFrom dplyr %>%
filter_agg_res <- function(
    condition_var = "Condition_dummy",
    input_dir = "",
    interactions_agg_binarized = "",
    interactions_agg_continuous = "",
    output_dir = ".") {
    # Load additional libraries
    if (file.exists(input_dir)) {
        message("Load interactions after aggregation: binarized...")
        interactions_binarized <- readRDS(glue::glue("{input_dir}/402a_filtering_detect_in_multi_samples.rds"))

        message("Load interactions after aggregation: continuous...")
        interactions_continuous <- readRDS(glue::glue("{input_dir}/402b_aggregation_samples.rds"))
    } else {
        message("Load interactions after aggregation: binarized...")
        interactions_binarized <- readRDS(interactions_agg_binarized)

        message("Load interactions after aggregation: continuous...")
        interactions_continuous <- readRDS(interactions_agg_continuous)
    }

    message("Combine...")
    combi <- interactions_binarized %>%
        dplyr::left_join(interactions_continuous, by = c(condition_var, "complex_interaction", "source_target"))
    message(glue::glue("Number of rows: {nrow(combi)}"))

    message("Save results...")
    saveRDS(combi, glue::glue("{output_dir}/402c_filtering_aggregated_res.rds"))
    message("Finished!")
}
