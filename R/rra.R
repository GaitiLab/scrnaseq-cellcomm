# ---- Helper functions ----
#' @title Format rankings of LIANA's internal methods
#' @param df dataframe
#' @return dataframe with columns: uid, method and rank (larger value = 'lower importance/priority')
FormatInternalMethodsLIANA <- function(df) {
    liana_internal_methods <- c(
        "connectome",
        "logfc",
        "natmi",
        "sca",
        "cytotalk"
    )
    liana_internal_methods_rank_cols <- paste0(liana_internal_methods, ".rank")

    return(
        df |>
            AddUniqueIdentifier() |>
            dplyr::rename_with(
                .cols = dplyr::all_of(liana_internal_methods_rank_cols),
                ~ paste(
                    stringr::str_remove(stringr::str_to_lower(.x), ".rank"),
                    "(LIANA)"
                )
            ) |>
            dplyr::select(
                uid,
                dplyr::ends_with("(LIANA)")
            ) |>
            tidyr::pivot_longer(
                cols = dplyr::ends_with("(LIANA)"),
                names_to = "method",
                values_to = "rank"
            ) |>
            duckplyr::as_duckdb_tibble()
    )
    # Every dataframe in list should have the following columns: uid, method and rank (larger value = 'lower importance/priority')
}
#' @title Extract dataframe with ranks
#' @description adds uid + ranks interactions and returns a dataframe with the required columns for downstream analyses
#' @param df dataframe with interactions
#' @return dataframe
ExtractDFWithRanks <- function(df, n_perm = 1e3) {
    return(
        df |>
            AddUniqueIdentifier() |>
            # Only affects CellPhoneDB
            dplyr::select(-dplyr::any_of("rank")) |>
            # As ranks are needed, add a custom rank based on -log10(pval) x interaction_score
            AddCCIRank(n_perm) |>
            dplyr::select(
                uid,
                method,
                rank
            )
    )
}

# ---- Wrapper ----

#' @title RRA interactions
#' @description Rank interactions based on CellChat, LIANA (use all methods), cell2cell and CellPhoneDB
#' @param list_of_cci_objects named list of CCI objects, should contain: CellChatv2, LIANA, Cell2Cell, CellPhoneDBv5
#' @param n_perm number of permutations to create pval_dummy variable for computing -log10(pval_dummy) * score, only applicable to CellChat, Cell2Cell and CellPhoneDB
#' @return dataframe
#' @export
AggregateCCIRanks <- function(
    list_of_cci_objects,
    n_perm = 1e3) {
    # ---- Constants ----
    methods_not_liana <- setdiff(
        names(list_of_cci_objects),
        "LIANA"
    )

    cols_common <- c(
        "source_target",
        "pval",
        "interaction_score",
        "complex_interaction",
        "method",
        "sample_id"
    )
    # Get all scores/pvalues for the methods side-by-side wide-format, for future reference + downstream analyses
    interactions_df <- list_of_cci_objects |>
        purrr::map_dfr(~ dplyr::select(.x, dplyr::all_of(cols_common))) |>
        duckplyr::as_duckdb_tibble() |>
        ToWideCCIDF()

    # ---- Defining/extracting ranks for each method ----

    # (1) Handle ranking of LIANA's internal methods
    rankings_internal_liana_methods_df <- list_of_cci_objects[["LIANA"]] |>
        FormatInternalMethodsLIANA()
    message("Extracted rankings of LIANA's internal CCI methods.")
    # (2) Handle CCI methods that don't provide a ranking

    rankings_methods_not_liana_df <- list_of_cci_objects[methods_not_liana] |>
        purrr::map_dfr(ExtractDFWithRanks, n_perm = n_perm)
    message("Added ranking for all other CCI methods.")
    # each dataframe should contain the following columns: uid, method and rank (larger value = 'lower importance/priority')

    # ---- Robust rank aggregation ----
    # Collect all individual ranks into mat
    rankings_of_all_methods_df <- list(
        rankings_internal_liana_methods_df,
        rankings_methods_not_liana_df
    ) |>
        dplyr::bind_rows() |>
        tidyr::pivot_wider(names_from = method, values_from = rank) |>
        duckplyr::as_duckdb_tibble()

    global_max_rank <- rankings_of_all_methods_df |>
        dplyr::select(dplyr::where(is.numeric)) |>
        max(
            na.rm = TRUE
        )

    rankings_of_all_methods_normalized_mat <- rankings_of_all_methods_df |>
        # Set missing interactions with max. rank
        dplyr::mutate(dplyr::across(dplyr::where(is.numeric), \(x) {
            tidyr::replace_na(x, global_max_rank)
        })) |>
        # Scale ranked matrix to a (0,1) range
        dplyr::mutate(dplyr::across(dplyr::where(is.numeric), \(x) {
            x / global_max_rank
        })) |>
        tibble::column_to_rownames("uid") |>
        as.matrix()
    message("Scaled rankings.")

    ranked_interactions_df <- rankings_of_all_methods_normalized_mat |>
        # Similar strategy as LIANA's .aggregate_rank() [.robust_rank_agg(), .rho_scores() .corr_beta_pvals()]
        RobustRankAggreg::aggregateRanks(
            rmat = _,
            method = "RRA"
        ) |>
        dplyr::rename(uid = Name) |>
        tibble::remove_rownames() |>
        tidyr::separate(
            uid,
            into = c("sample_id", "source_target", "complex_interaction"),
            sep = "\\|"
        ) |>
        duckplyr::as_duckdb_tibble() |>
        # The 'Score' can be interpreted as a pvalue as mentioned in LIANA docs.
        dplyr::rename(pval = Score) |>
        # Add information on method-level - for reference + additional downstream analyses
        dplyr::full_join(interactions_df)
    message("Completed RRA.")

    message("Finished!")
    return(ranked_interactions_df)
}
