#' @title Reshape & add interaction name
#' @description Takes the CellChat a single 2D array representing the probabilities/pvalues for an interaction in the CCI database from the cellchat object, i.e. `cc_obj@net$prob` or `cc_obj@net$pval` after running CellChat. The array (rownames=source cell type, columnnames=target cell type) is converted to the long-format.
#' @param interaction_name name of interaction generally in format: <ligand>__<receptor>
#' @param cci_ls list with 2D arrays representing the results for the interactions in the (custom) CellChat database
#' @param mode character string which CellChat output slot is formatted, i.e. 'interaction_score' or 'pval'
#' @return dataframe with 4 columns, i.e. source, target, interaction_score, interaction
ReshapeAndAddInteractionName <- function(
    interaction_name,
    cci_ls,
    mode = c("pval", "interaction_score")[1]) {
    match.arg(mode, c("pval", "interaction_score"))

    return(
        cci_ls[, , interaction_name] |>
            as.data.frame() |>
            tibble::rownames_to_column("source") |>
            tidyr::pivot_longer(
                cols = !source,
                names_to = "target",
                values_to = mode
            ) |>
            duckplyr::as_duckdb_tibble() |>
            dplyr::mutate(interaction = interaction_name)
    )
}

#' @title Format CellChat CCI results
#' @description Takes that 'raw' output from CellChat reshapes and reformats that data for downstream aggregation of the CCI results of the different tools.
#' @param cc_obj CellChat outputs
#' @param ref_db dataframe with reference database of interactions
#' @param sample_id character string for sample_id to add to dataframe (if non given set to NA)
#' @param n_cores integer indicating the number of cores to use
#' @return dataframe with columns: source_target, interaction_score, pval, complex_interaction, method, sample_id
#' @export
FormatCellChat <- function(
    cc_obj,
    ref_db,
    sample_id = NA,
    n_cores = 1) {
    cluster <- parallel::makeCluster(n_cores)

    # Needed for reshaping data, dependent on the CCI database used for running CellChat
    # net contains 'prob' and 'pval' -> both contain the same interactions
    included_interactions <- names(cc_obj@net$prob[1, 1, ])

    cc_combined_df <- purrr::map2(
        cc_obj@net[c("prob", "pval")],
        c("interaction_score", "pval"),
        \(cci_ls, mode) {
            pbapply::pblapply(
                included_interactions,
                ReshapeAndAddInteractionName,
                cci_ls = cci_ls,
                mode = mode,
                cl = cluster
            ) |>
                dplyr::bind_rows()
        }
    ) |>
        purrr::reduce(dplyr::full_join) |>
        dplyr::left_join(ref_db, by = "interaction") |>
        dplyr::select(-interaction) |>
        tidyr::unite(source_target, source, target, sep = "__") |>
        # Needed for aggregating/combining the CCI data from different tools & samples
        dplyr::mutate(method = "CellChatv2", sample_id = !!sample_id)

    message("Finished!")
    return(cc_combined_df)
}
