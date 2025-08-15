#' @title Format individual output from CellPhoneDB v5
#' @param df dataframe
#' @param mode character string indicating what information the dataframe contains.
#' @return dataframe with columns: interacting_pair, source_target and <mode>
#' @export
format_cpdb <- function(
    df,
    mode = c("pval", "sign_mean", "mean", "interaction_score")[1]
) {
    match.arg(mode, c("pval", "sign_mean", "mean", "interaction_score"))
    source_target_pairs <- colnames(df)[stringr::str_detect(
        colnames(df),
        "\\|"
    )]
    cols_to_keep <- c(source_target_pairs, "interacting_pair")
    if (mode == "sign_mean") {
        cols_to_keep <- c(cols_to_keep, "rank")
    }
    return(
        df |>
            duckplyr::as_duckdb_tibble() |>
            dplyr::select(dplyr::all_of(cols_to_keep)) |>
            tidyr::pivot_longer(
                cols = source_target_pairs,
                values_to = mode,
                names_to = "source_target"
            ) |>
            duckplyr::as_duckdb_tibble() |>
            dplyr::mutate(
                source_target = stringr::str_replace_all(
                    source_target,
                    "\\|",
                    "__"
                )
            )
    )
}


#' @title Format CellPhoneDB v5 results
#' @param cpdb_dfs named list of dfs (interaction_scores = 'statistical_analysis_interaction_scores', pval = 'statistical_analysis_pvalues', means = 'statistical_analysis_means, sign_mean = 'statistical_analysis_significant_means')
#' @param ref_db dataframe with reference database of interactions
#' @param sample_id character string for sample_id to add to dataframe (if non given set to NA)
#' @return dataframe with columns: source_target, interaction_score, pval, rank, sign_mean, mean, method, sample_id, complex_interaction
#' @export
format_cpdb_wrapper <- function(cpdb_dfs, ref_db, sample_id = NA) {
    df <- purrr::map2(cpdb_dfs, names(cpdb_dfs), format_cpdb) |>
        purrr::reduce(
            left_join,
            by = c("interacting_pair", "source_target")
        ) |>
        dplyr::mutate(method = "CellPhoneDBv5", sample_id = !!sample_id) |>
        dplyr::left_join(
            ref_db,
            by = c("interacting_pair" = "interaction")
        ) |>
        # in favor of 'complex_interaction'
        dplyr::select(-interacting_pair)
    message("Finished!")
    return(df)
}
