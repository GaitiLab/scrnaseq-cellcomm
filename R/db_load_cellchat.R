#' @title Load CellChat DB
#' @description load cellcat database and format for merging
#' @return dataframe with 'source_genesymbol' and 'target_genesymbol'
#' @importFrom dplyr %>%
#' @export
load_cellchat_db <- function() {
    cellchat_db <- CellChat::CellChatDB.human
    cellchat_db_interaction <- cellchat_db$interaction

    interactions <- cellchat_db_interaction %>%
        tibble::remove_rownames() %>%
        dplyr::select(interaction_name, ligand, receptor, ligand.symbol, receptor.symbol) %>%
        dplyr::mutate(
            # Separate subunits with underscores '_'
            source_genesymbol = stringr::str_replace_all(ligand.symbol, ", ", "_"),
            target_genesymbol = stringr::str_replace_all(receptor.symbol, ", ", "_")
        ) %>%
        dplyr::mutate(
            # Ensure all genes are capitalized
            source_genesymbol = ifelse(source_genesymbol == "", toupper(ligand), source_genesymbol),
            target_genesymbol = ifelse(target_genesymbol == "", toupper(receptor), target_genesymbol)
        ) %>%
        dplyr::select(source_genesymbol, target_genesymbol)
    return(interactions)
}
