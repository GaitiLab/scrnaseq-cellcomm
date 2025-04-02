#' @title Hard combine databases
#' @description Combine databases based on 'source_genesymbol' and 'target_genesymbol'
#' @param cellchat_db cellchat database
#' @param cpdb_db cpdb_database
#' @param liana_db liana database
#' @export
#' @importFrom dplyr %>%
hard_combine_dbs <- function(cellchat_db, cpdb_db, liana_db) {
    # Load databases: CellChat and CellphoneDB
    cellchat_db <- cellchat_db %>% dplyr::mutate(method = "CellChat extracted")
    cpdb_db <- cpdb_db %>% dplyr::mutate(method = "CellphoneDB extracted")

    # Combine CellChat and CellphoneDB databases
    cpdb_cellchat_db <- rbind(cpdb_db, cellchat_db)

    # Add missing columns and set to empty string
    missing_cols <- setdiff(colnames(liana_db), colnames(cpdb_cellchat_db))
    cpdb_cellchat_db[missing_cols] <- ""

    # Combine LIANA database with CellChat and CellphoneDB databases
    db_combined <- rbind(liana_db, cpdb_cellchat_db) %>%
        tidyr::unite(
            complex_interaction,
            source_genesymbol,
            target_genesymbol,
            sep = "__",
            remove = FALSE
        )
    return(db_combined)
}
