#' @title Add unique identifier
#'
#' @param df dataframe
#' @param cols_to_merge columns that should be used to create the uid
#'
#' @return df
#'
#' @export
AddUniqueIdentifier <- function(
    df,
    cols_to_merge = c("sample_id", "source_target", "complex_interaction")) {
    return(
        df |>
            tidyr::unite(
                col = uid,
                dplyr::all_of(cols_to_merge),
                sep = "|"
            ) |>
            duckplyr::as_duckdb_tibble()
    )
}
