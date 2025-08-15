#' @title Take consensus
#' @description Assess number of methods that detect an interaction in a source-target
#' @param list_of_obj list with the dataframes for CellChat, LIANA, cell2cell and cellphonedb
#' @param alpha significance threshold (default = 0.05)
#' @export
take_consensus <- function(
    interactions_df,
    alpha = 0.05
) {
    # contains all interactions, may contain NAs if interactions was missing
    scores_df <- interactions_df |> to_wide_cci_df()

    interactions_detected_in_methods_df <- interactions_df |>
        dplyr::mutate(
            is_detected_signif = pval < alpha,
        ) |>
        # Within a sample_id - source-target combination, is the interaction detected by multiple methods
        # TODO remove group_by() later in favor of using .by argument in summarise
        dplyr::group_by(sample_id, source_target, complex_interaction) |>
        dplyr::summarise(
            is_detected_signif_in_n_methods = sum(is_detected_signif),
            is_detected_signif_in_methods = ifelse(
                is_detected_signif_in_n_methods == 0,
                NA,
                paste(
                    method[is_detected_signif],
                    collapse = ", "
                )
            ),
            .by = c(sample_id, source_target, complex_interaction)
        ) |>
        # TODO remove later as will be redundant if .by is used
        # dplyr::ungroup() |>
        # Still contains the interactions that have less 3 methods (no signif filter)
        dplyr::mutate(
            lenient_voting = (is_detected_signif_in_n_methods >= 3) &
                stringr::str_detect(is_detected_signif_in_methods, "LIANA"),
            stringent_voting = is_detected_signif_in_n_methods == 4
        ) |>
        dplyr::arrange(
            dplyr::desc(is_detected_signif_in_n_methods)
        ) |>
        dplyr::left_join(scores_df) |>
        duckplyr::as_duckdb_tibble()

    message("Finished!")
    return(interactions_detected_in_methods_df)
}
