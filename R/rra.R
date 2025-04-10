#' @title RRA interactions
#' @description Rank interactions based on CellChat, LIANA (use all methods), cell2cell and CellPhoneDB
#' @param cellchat_obj path to file generated by 'format_cellchat()'
#' @param liana_obj path to file generated by 'format_liana()'
#' @param cell2cell_obj path to file generated by 'format_cell2cell()'
#' @param cpdb_obj path to file generated by 'format_cpdb()'
#' @param output_dir output directory for saving output (default = '.')
#' @param sample_id sample id to use for saving object
#' @param n_perm number of permutations to create pval_dummy variable for computing -log10(pval_dummy) * score, only applicable to CellChat, Cell2Cell and CellPhoneDB
#' @importFrom dplyr %>%
#' @export
rra_interactions <- function(
    cellchat_obj,
    liana_obj,
    cell2cell_obj,
    cpdb_obj,
    sample_id = NULL,
    output_dir = ".",
    n_perm = 1000) {
    #  Sanity checks
    if (
        !(file.exists(cellchat_obj) && endsWith(tolower(cellchat_obj), ".rds"))
    ) {
        stop(glue::glue(
            "'{cellchat_obj}' does not exists or is not an RDS object"
        ))
    }
    if (!(file.exists(liana_obj) && endsWith(tolower(liana_obj), ".rds"))) {
        stop(glue::glue(
            "'{liana_obj}' does not exists or is not an RDS object"
        ))
    }
    if (
        !(file.exists(cell2cell_obj) &&
            endsWith(tolower(cell2cell_obj), ".rds"))
    ) {
        stop(glue::glue(
            "'{cell2cell_obj}' does not exists or is not an RDS object"
        ))
    }
    if (!(file.exists(cpdb_obj) && endsWith(tolower(cpdb_obj), ".rds"))) {
        stop(glue::glue("'{cpdb_obj}' does not exists or is not an RDS object"))
    }
    if (!file.exists(output_dir)) {
        stop("Output directory does not exist")
    }
    if (is.null(sample_id)) {
        stop("No sample_id given")
    }

    message("Load postprocessed objects from CCIs...")
    # Pre-filtering on p-value, only impacts when downsampling is done (pval == NA, if interaction not detected in each downsampling run)
    obj_cellchat <- readRDS(cellchat_obj) %>%
        tidyr::unite(
            unique_combi,
            source_target,
            complex_interaction,
            sep = "|",
            remove = FALSE
        )

    obj_liana <- readRDS(liana_obj) %>%
        tidyr::unite(
            unique_combi,
            source_target,
            complex_interaction,
            sep = "|",
            remove = FALSE
        )

    obj_cell2cell <- readRDS(cell2cell_obj) %>%
        tidyr::unite(
            unique_combi,
            source_target,
            complex_interaction,
            sep = "|",
            remove = FALSE
        )

    obj_cpdb <- readRDS(cpdb_obj) %>%
        tidyr::unite(
            unique_combi,
            source_target,
            complex_interaction,
            sep = "|",
            remove = FALSE
        )

    message(glue::glue(
        "Number of interactions in CellChat BEFORE filtering: {nrow(obj_cellchat)}"
    ))
    message(glue::glue(
        "Number of interactions in LIANA BEFORE filtering: {nrow(obj_liana)}"
    ))
    message(glue::glue(
        "Number of interactions in Cell2Cell BEFORE filtering: {nrow(obj_cell2cell)}"
    ))
    message(glue::glue(
        "Number of interactions in CPDB BEFORE filtering: {nrow(obj_cpdb)}"
    ))

    # ---- Order from most-important (significant) to least based on interaction
    # Sorting by score/specificity (NOT p-value)
    connectome <- obj_liana %>%
        dplyr::select(
            unique_combi,
            source_target,
            complex_interaction,
            Sample,
            connectome.rank
        ) %>%
        dplyr::mutate(method = "Connectome (LIANA)")

    logfc <- obj_liana %>%
        dplyr::select(
            unique_combi,
            source_target,
            complex_interaction,
            Sample,
            logfc.rank
        ) %>%
        dplyr::mutate(method = "LogFC (LIANA)")

    natmi <- obj_liana %>%
        dplyr::select(
            unique_combi,
            source_target,
            complex_interaction,
            Sample,
            natmi.rank
        ) %>%
        dplyr::mutate(method = "NATMI (LIANA)")

    sca <- obj_liana %>%
        dplyr::select(
            unique_combi,
            source_target,
            complex_interaction,
            Sample,
            sca.rank
        ) %>%
        dplyr::mutate(method = "SCA (LIANA)")

    cytotalk <- obj_liana %>%
        dplyr::select(
            unique_combi,
            source_target,
            complex_interaction,
            Sample,
            cytotalk.rank
        ) %>%
        dplyr::mutate(method = "Cytotalk (LIANA)")

    # Sorting by p-values
    obj_cellchat <- obj_cellchat %>%
        dplyr::mutate(pval_dummy = ifelse(pval == 0, 1 / n_perm, pval)) %>%
        dplyr::mutate(log10_score = -log10(pval_dummy) * CellChat_score) %>%
        dplyr::mutate(
            CellChat.rank = dplyr::dense_rank(dplyr::desc(log10_score))
        )

    obj_cell2cell <- obj_cell2cell %>%
        dplyr::mutate(pval_dummy = ifelse(pval == 0, 1 / n_perm, pval)) %>%
        dplyr::mutate(log10_score = -log10(pval_dummy) * Cell2Cell_score) %>%
        dplyr::mutate(
            Cell2Cell.rank = dplyr::dense_rank(dplyr::desc(log10_score))
        )

    obj_cpdb <- obj_cpdb %>%
        dplyr::select(-rank) %>%
        dplyr::mutate(pval_dummy = ifelse(pval == 0, 1 / n_perm, pval)) %>%
        dplyr::mutate(log10_score = -log10(pval_dummy) * CellPhoneDB_score) %>%
        dplyr::mutate(
            CellPhoneDB.rank = dplyr::dense_rank(dplyr::desc(log10_score))
        )

    tool_res <- list(
        connectome,
        sca,
        cytotalk,
        natmi,
        logfc,
        obj_cellchat,
        obj_cell2cell,
        obj_cpdb
    )

    rankmat <- tool_res %>%
        purrr::map(function(res) {
            res %>%
                dplyr::select(unique_combi, dplyr::ends_with("rank"))
        }) %>%
        purrr::reduce(full_join, by = "unique_combi") %>%
        tibble::column_to_rownames("unique_combi") %>%
        as.matrix()
    # Set missing interactions with max. rank
    rankmat[is.na(rankmat)] <- max(rankmat, na.rm = TRUE)

    # Scale ranked matrix to a (0,1) range
    rank_mat <- rankmat %>%
        {
            . / max(.)
        }

    ranked_interactions <- RobustRankAggreg::aggregateRanks(
        rmat = rank_mat,
        method = "RRA"
    ) %>%
        dplyr::rename(unique_combi = Name) %>%
        tibble::remove_rownames() %>%
        tidyr::separate(
            unique_combi,
            into = c("source_target", "complex_interaction"),
            sep = "\\|"
        ) %>%
        dplyr::rename(pval = Score) %>%
        dplyr::mutate(
            Sample = obj_liana %>% dplyr::pull(Sample) %>% unique()
        ) %>%
        tibble::remove_rownames() %>%
        # Adding individual p-values
        dplyr::left_join(
            obj_cellchat %>%
                dplyr::rename(cellchat = pval) %>%
                dplyr::select(
                    Sample,
                    source_target,
                    complex_interaction,
                    cellchat,
                    CellChat_score
                )
        ) %>%
        dplyr::left_join(
            obj_cell2cell %>%
                dplyr::rename(cell2cell = pval) %>%
                dplyr::select(
                    Sample,
                    source_target,
                    complex_interaction,
                    cell2cell,
                    Cell2Cell_score
                )
        ) %>%
        dplyr::left_join(
            obj_liana %>%
                dplyr::rename(liana = pval) %>%
                dplyr::select(
                    Sample,
                    source_target,
                    complex_interaction,
                    liana,
                    LIANA_score
                )
        ) %>%
        dplyr::left_join(
            obj_cpdb %>%
                dplyr::rename(cpdb = pval) %>%
                dplyr::select(
                    Sample,
                    source_target,
                    complex_interaction,
                    cpdb,
                    CellPhoneDB_score
                )
        )
    message("Save RRA interactions...")
    saveRDS(
        ranked_interactions,
        glue::glue("{output_dir}/{sample_id}__interactions_agg_rank.rds")
    )
    message("Finished!")
}
