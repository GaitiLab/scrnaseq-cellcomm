#' @title Update CellPhoneDB
#' @description Update CellPhoneDB and generate required files
#' @param source_cpdb_dir directory with CellPhoneDB files
#' @param db database with at least 'source_genesymbol' and 'target_genesymbol' and their subunits starting with 'source_genesymbol_subunit' or 'target_genesymbol_subunit'
#' @param gene_info_table dataframe with uniprot, protein name, genesymbol and ensembl
#' @param output_dir path to directory for saving the generated CellPhoneDB files
#' @param return_list return list with generated dataframes (default = FALSE)
#' @return if return_list = TRUE, return list with complex_input, gene_input, interaction_input and protein_input
#' @export
#' @importFrom dplyr %>%
update_cpdb_db <- function(source_cpdb_dir, db, gene_info_table, output_dir, return_list = FALSE) {
    # ---- DB update cpdb.R ----
    # Grab column names for generating new files
    complex_input_cols <- colnames(read.csv(glue::glue("{source_cpdb_dir}/complex_input.csv")))
    gene_input_cols <- colnames(read.csv(glue::glue("{source_cpdb_dir}/gene_input.csv")))
    interaction_input_cols <- colnames(read.csv(glue::glue("{source_cpdb_dir}/interaction_input.csv")))
    protein_input_cols <- colnames(read.csv(glue::glue("{source_cpdb_dir}/protein_input.csv")))

    gene_info_table <- gene_info_table %>% dplyr::distinct(uniprot, .keep_all = TRUE)

    # Gene input
    # mandatory: gene_name, uniprot, hgnc_symbol and ensembl
    genes_input <- gene_info_table %>%
        dplyr::rename(hgnc_symbol = genesymbol) %>%
        dplyr::mutate(
            gene_name = hgnc_symbol
        ) %>%
        dplyr::select(
            gene_name, uniprot, hgnc_symbol, ensembl
        )

    # Protein input
    # mandatory: uniprot, protein_name
    protein_input <- gene_info_table %>%
        dplyr::rename(protein_name = "uniprot_entry") %>%
        dplyr::select(uniprot, protein_name) %>%
        dplyr::distinct(uniprot, .keep_all = TRUE)


    # Complex input
    # mandatory: complex_name, complex_subunits
    ligand_complexes <- db %>%
        tidyr::separate(source, paste0("uniprot_", seq_len(5)), sep = "_", remove = FALSE) %>%
        dplyr::select(source_genesymbol, dplyr::starts_with("uniprot_")) %>%
        dplyr::rename(complex_name = source_genesymbol)
    receptor_complexes <- db %>%
        tidyr::separate(target, paste0("uniprot_", seq_len(5)), sep = "_", remove = FALSE) %>%
        dplyr::select(target_genesymbol, dplyr::starts_with("uniprot_")) %>%
        dplyr::rename(complex_name = target_genesymbol)
    complex_input <- rbind(ligand_complexes, receptor_complexes) %>% dplyr::distinct()
    message(glue::glue("Number of complexes: {nrow(complex_input)}..."))
    complex_input[complex_input == ""] <- NA

    # Ensure that the dataframe only contains protein complexes, i.e. uniprot_2 should NOT be empty
    complex_input <- complex_input %>% filter(!is.na(uniprot_2))
    message(glue::glue("Number of complexes: {nrow(complex_input)}..."))

    complex_input <- cbind(complex_input, apply(complex_input %>% dplyr::select(dplyr::starts_with("uniprot_")), 1, make_complex_wrapper, collapse = "_") %>%
        data.frame() %>%
        dplyr::rename(dummy = ".")) %>%
        dplyr::distinct(dummy, .keep_all = TRUE) %>%
        dplyr::select(-dummy)
    message(glue::glue("Number of complexes: {nrow(complex_input)}..."))

    # Interaction Input
    # mandatory:  “partner_a”; “partner_b”; “annotation_strategy”; “source”
    interactions_input <- db %>%
        dplyr::select(
            # Proteins
            source, target,
            # Necessary for checks
            source_genesymbol_subunit_2, target_genesymbol_subunit_2,
            # Genesymbols
            source_genesymbol, target_genesymbol
        ) %>%
        # Rename protein columns (source, target) as partner_{a,b}
        dplyr::rename(partner_a = source, partner_b = target) %>%
        dplyr::mutate(
            annotation_strategy = "user_curated",
            version = "CellPhoneDBcore4.1",
            sources = "User curated",
            # If complex, then user 'complex_name', i.e. source_genesymbol/target_genesymbol
            partner_a = ifelse(is.na(source_genesymbol_subunit_2), partner_a, source_genesymbol),
            partner_b = ifelse(is.na(target_genesymbol_subunit_2), partner_b, target_genesymbol),
        ) %>%
        # Select
        dplyr::select(
            partner_a, partner_b, annotation_strategy, sources
        )

    # account for missing cols
    missing_cols <- setdiff(gene_input_cols, colnames(genes_input))
    genes_input[missing_cols] <- ""
    missing_cols <- setdiff(protein_input_cols, colnames(protein_input))
    protein_input[missing_cols] <- ""
    missing_cols <- setdiff(complex_input_cols, colnames(complex_input))
    complex_input[missing_cols] <- ""
    complex_input$version <- "CellPhoneDBcore4.1"
    missing_cols <- setdiff(interaction_input_cols, colnames(interactions_input))
    interactions_input[missing_cols] <- ""

    genes_input[is.na(genes_input)] <- ""
    protein_input[is.na(protein_input)] <- ""
    complex_input[is.na(complex_input)] <- ""
    interactions_input[is.na(interactions_input)] <- ""

    # Save files
    message("Save files (genes_input, protein_input, complex_input and interaction_input)...")
    write.csv(genes_input, glue::glue("{output_dir}/gene_input.csv"),
        row.names = FALSE, quote = FALSE
    )
    write.csv(protein_input, glue::glue("{output_dir}/protein_input.csv"), row.names = FALSE, quote = FALSE)
    write.csv(complex_input, glue::glue("{output_dir}/complex_input.csv"), row.names = FALSE, quote = FALSE)
    write.csv(interactions_input, glue::glue("{output_dir}/interaction_input.csv"), row.names = FALSE, quote = FALSE)

    if (return_list) {
        message("Returning outputs as list...")
        return(list(genes_input = genes_input, protein_input = protein_input, complex_input = complex_input, interactions_input = interactions_input))
    }
}
