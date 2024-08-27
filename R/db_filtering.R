#' @title Filter and clean-up database
#' @param db dataframe with database of interactions, the following columns should at least be present 'source_genesymbol' and 'target_genesymbol'
#' @export
#' @importFrom dplyr %>%
db_filtering <- function(db) {
    other_cols <- setdiff(colnames(db), c("source_genesymbol", "target_genesymbol", "complex_interaction"))
    db <- db %>%
        # Split into subunits
        tidyr::separate(source_genesymbol, paste0("source_genesymbol_subunit_", seq_len(5)),
            sep = "_", remove = FALSE
        ) %>%
        tidyr::separate(target_genesymbol, paste0("target_genesymbol_subunit_", seq_len(5)),
            sep = "_", remove = FALSE
        )
    message(glue::glue("Number of interactions: {nrow(db)}..."))

    db <- db %>%
        dplyr::mutate(
            source_genesymbol_ordered := apply(db %>% dplyr::select(dplyr::starts_with("source_genesymbol_subunit")), 1, make_complex_wrapper, "_"),
            target_genesymbol_ordered := apply(db %>% dplyr::select(dplyr::starts_with("target_genesymbol_subunit")), 1, make_complex_wrapper, "_")
        ) %>%
        # Save original source and genesymbols (unordered) with prefix '_OG'
        dplyr::rename(source_genesymbol_OG = source_genesymbol, target_genesymbol_OG = target_genesymbol) %>%
        dplyr::select(-dplyr::starts_with("source_genesymbol_subunit"), -dplyr::starts_with("target_genesymbol_subunit")) %>%
        # Then rename ordered source and genesymbols
        dplyr::rename(source_genesymbol = source_genesymbol_ordered, target_genesymbol = target_genesymbol_ordered) %>%
        # Split into subunits (these are the ordered subunits)
        tidyr::separate(source_genesymbol, paste0("source_genesymbol_subunit_", seq_len(5)),
            sep = "_", remove = FALSE
        ) %>%
        tidyr::separate(target_genesymbol, paste0("target_genesymbol_subunit_", seq_len(5)),
            sep = "_", remove = FALSE
        )
    message(glue::glue("Number of interactions: {nrow(db)}..."))


    # Check if subunits have at least 2 characters (otherwise invalid gene)
    message("Filter interactions with invalid subunits...")
    target_genesymbol_mask <- apply(db %>% dplyr::select(dplyr::starts_with("target_genesymbol_subunit_")), 1, check_length_genesymbol)
    source_genesymbol_mask <- apply(db %>% dplyr::select(dplyr::starts_with("source_genesymbol_subunit_")), 1, check_length_genesymbol)
    # Remove interactions where genes only have a single character
    db <- db[source_genesymbol_mask & target_genesymbol_mask, ]
    message(glue::glue("Number of interactions: {nrow(db)}..."))

    # Adding protein ids (uniprot)
    all_genes <- db %>%
        dplyr::select(source_genesymbol, target_genesymbol) %>%
        as.list() %>%
        unlist() %>%
        stringr::str_split(., "_") %>%
        unlist() %>%
        unique() %>%
        stringr::str_to_upper()
    all_genes <- data.frame(genesymbol = all_genes[all_genes != ""] %>% unique())
    proteins <- get_proteins_wrapper(genes_df = all_genes)


    # Check mapping gene to protein (uniprot)
    source_mask <- apply(
        db %>% dplyr::select(dplyr::starts_with("source_genesymbol_subunit_")), 1, is_mapped,
        lookup_table = proteins,
        lookup_var = "genesymbol"
    )

    target_mask <- apply(db %>% dplyr::select(dplyr::starts_with("target_genesymbol_subunit_")), 1,
        is_mapped,
        lookup_table = proteins,
        lookup_var = "genesymbol"
    )
    # Remove interactions for which we don't have a matching uniprot (protein)
    db <- db[source_mask & target_mask, ]
    message(glue::glue("Number of interactions: {nrow(db)}..."))

    # Add proteins to dataframe
    db$target <- apply(
        db %>% dplyr::select(dplyr::starts_with("target_genesymbol_subunit_")),
        1,
        add_proteins,
        gene_uniprot_table = proteins
    )

    db$source <- apply(
        db %>% dplyr::select(dplyr::starts_with("source_genesymbol_subunit_")),
        1,
        add_proteins,
        gene_uniprot_table = proteins
    )

    # Ensure that the dataframe only contains interactions for which the protein IDs are available for both source AND target
    db <- db %>%
        filter(!is.na(source), !is.na(target)) %>%
        # Rename previously set 'complex_interaction' (unordered subunits)
        dplyr::rename(complex_interaction_OLD = complex_interaction) %>%
        dplyr::mutate(
            complex_interaction = paste0(source_genesymbol, "__", target_genesymbol)
        )
    message(glue::glue("Number of interactions: {nrow(db)}..."))

    db <- db %>% filter(!is.na(source_genesymbol), !is.na(target_genesymbol), !is.na(source), !is.na(target))
    message(glue::glue("Number of interactions: {nrow(db)}..."))

    collapsed_method <- db %>%
        dplyr::group_by(complex_interaction) %>%
        dplyr::reframe(method = paste0(method, collapse = ", ")) %>%
        dplyr::ungroup()

    db <- db %>%
        dplyr::select(-method) %>%
        dplyr::left_join(collapsed_method) %>%
        dplyr::arrange(dplyr::across(dplyr::all_of(other_cols), desc)) %>%
        dplyr::distinct(complex_interaction, .keep_all = TRUE) %>%
        dplyr::mutate(complex_interaction_rev = paste0(target_genesymbol, "__", source_genesymbol))
    message(glue::glue("Number of interactions: {nrow(db)}..."))

    # When reversing order i.e. {target_genesymbol}__{source_genesymbol} -> duplicates
    duplicate_interactions <- intersect(db %>%
        filter(source_genesymbol != target_genesymbol) %>%
        dplyr::pull(complex_interaction_rev), db$complex_interaction)
    message(glue::glue("Number of duplicated interactions (when undirected): {length(duplicate_interactions)}..."))

    # Not removing these interactions, just marking them with 'is_dupl_undirected'
    db <- db %>% dplyr::mutate(is_dupl_undirected = complex_interaction %in% duplicate_interactions | complex_interaction_rev %in% duplicate_interactions)
    return(list(db = db, gene_info_table = proteins))
}
