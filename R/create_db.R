#' @title Create interactions databasee
#' @description Generates database for LIANA, CellChat, Cell2Cell and CellPhoneDB
#' @param source_cpdb_dir path to 'interaction_input.csv' from CellPhoneDB database
#' @param output_dir path to directory for saving the generated CellPhoneDB files
#' @export
create_db <- function(source_cpdb_dir, output_dir) {
    stopifnot(file.exists(source_cpdb_dir))
    stopifnot(file.exists(output_dir))
    # Interaction 'source_genesymbol' and 'target_genesymbol'
    message("(1/12) Retrieving LIANA DBs...")
    liana_db <- load_liana_dbs()

    message("(2/12) Retrieving CellChat DB...")
    cellchat_db <- load_cellchat_db()

    message("(3/12) Retrieving CPDB DB...")
    cpdb_db <- load_cpdb_db(source_cpdb_dir = source_cpdb_dir)

    message("(4/12) Combining databases...")
    combi_db <- hard_combine_dbs(cellchat_db = cellchat_db, cpdb_db = cpdb_db, liana_db = liana_db)

    message("(5/12) Filtering and unifying combined database...")
    combi_db_filtered <- db_filtering(combi_db)

    message("(6/12) Format unified database for CellPhoneDB, then save...")
    updated_cpdb_files <- update_cpdb_db(source_cpdb_dir = source_cpdb_dir, db = combi_db_filtered$db, gene_info_table = combi_db_filtered$gene_info_table, output_dir = output_dir, return_list = TRUE)

    # Update CellChat datbase + save
    message("(7/12) Format unified database for CellChat, then save...")
    update_cellchat_db(
        geneInfo = updated_cpdb_files$genes_input,
        output_dir = output_dir,
        db = combi_db_filtered$db
    )

    message("(8/12) Create simplified ref. database...")
    ref_db <- combi_db_filtered$db %>%
        dplyr::select(source_genesymbol, target_genesymbol, complex_interaction, method, is_dupl_undirected, complex_interaction_OLD) %>%
        dplyr::mutate(
            ligand_complex = stringr::str_replace_all(source_genesymbol, "_", ":"),
            receptor_complex = stringr::str_replace_all(target_genesymbol, "_", ":"),
        ) %>%
        dplyr::select(
            source_genesymbol, target_genesymbol, complex_interaction, ligand_complex, receptor_complex,
            method, is_dupl_undirected, complex_interaction_OLD
        ) %>%
        # TODO Add 'interaction' based on reordering -> testing
        dplyr::mutate(interaction = stringr::str_replace_all(complex_interaction, "__", "_")) %>%
        # Ensure that subunits are split by colon (':')
        tidyr::unite(complex_interaction, ligand_complex, receptor_complex, sep = "__", remove = FALSE)

    # Save LIANA db and Cell2Cell db
    message("(9/12) Save unified database for LIANA...")
    saveRDS(combi_db_filtered$db, glue::glue("{output_dir}/liana_db.rds"))
    message("(10/12) Save unified database for Cell2Cell (same format as LIANA)...")
    write.csv(combi_db_filtered$db, glue::glue("{output_dir}/cell2cell_db.csv"), row.names = FALSE)

    message("(11/12) Save unified & simplified database...")
    saveRDS(ref_db, glue::glue("{output_dir}/ref_db.rds"))

    message("(12/12) Zipping unified CellPhoneDB database...")
    run_script <- system.file("create_db.sh", package = "scrnaseq.cellcomm")
    python_script <- system.file("Python/update_cellphonedb.py", package = "scrnaseq.cellcomm")

    # TODO change when alternative found, make executable
    Sys.chmod(run_script, mode = "0555")
    Sys.chmod(python_script, mode = "0555")
    system(glue::glue("{run_script} {python_script} {output_dir}"))
    # Revert to previous mode
    Sys.chmod(run_script, mode = "0644")
    Sys.chmod(python_script, mode = "0644")

    cpdb_file <- list.files(output_dir, pattern = ".zip", full.names = TRUE)[1]
    file.rename(cpdb_file, glue::glue("{output_dir}/cellphonedb.zip"))

    if (file.exists(glue::glue("{output_dir}/cellphonedb.zip"))) {
        message("Clean-up redundant files used for generating 'cellphonedb.zip'...")
        files_to_remove <- c("complex_input", "gene_input", "interaction_input", "protein_input")
        for (filename in files_to_remove) {
            file.remove(glue::glue("{output_dir}/{filename}.csv"))
        }
    }
    message("Finished!")
}
