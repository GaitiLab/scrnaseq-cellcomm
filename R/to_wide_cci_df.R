#' @title Convert to side-by-side comparison format
#' @param df should contain the columns: method, interaction_score and pval
#' @return dataframe in wide-format interaction_score_<method> and pval_<method>
#' @export
to_wide_cci_df <- function(df) {
    return(
        df |>
            tidyr::pivot_wider(
                names_from = method,
                values_from = c(interaction_score, pval)
            ) |>
            duckplyr::as_duckdb_tibble()
    )
}
