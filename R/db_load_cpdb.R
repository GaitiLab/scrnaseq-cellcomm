#' @title Load CellPhoneDB DB
#' @description load CellPhoneDB database and format for merging
#' @param source_cpdb_dir path to 'interaction_input.csv' from CellPhoneDB database
#' @return dataframe with 'source_genesymbol' and 'target_genesymbol'
#' @importFrom dplyr %>%
#' @export
load_cpdb_db <- function(source_cpdb_dir) {
    cpdb_interaction_input <- data.table::fread(
        glue::glue("{source_cpdb_dir}/interaction_input.csv"),
        sep = ","
    ) %>%
        data.frame() %>%
        # Split into subunits by "-"
        tidyr::separate(interactors, into = paste0("unit_", seq(5)), sep = "-", remove = FALSE) %>%
        dplyr::mutate(
            # Handle individual genes with a '-', i.e. merging with previous subunit if length of unit_2 is a single character, e.g. HLA-genes (HLA-A, HLA-B etc.)
            source_genesymbol_complex = dplyr::case_when(
                nchar(unit_2) == 1 ~ paste0(unit_1, "-", unit_2), TRUE ~ unit_1
            ),
            target_genesymbol_complex = dplyr::case_when(
                nchar(unit_2) == 1 ~ unit_3,
                nchar(unit_3) == 1 ~ paste0(unit_2, "-", unit_3), TRUE ~ unit_2
            )
        ) %>%
        # Format source_genesymbol and target_genesymbol
        dplyr::mutate(
            source_genesymbol = stringr::str_replace_all(source_genesymbol_complex, "\\+", "_"),
            target_genesymbol = stringr::str_replace_all(target_genesymbol_complex, "\\+", "_")
        ) %>%
        dplyr::select(source_genesymbol, target_genesymbol)

    return(cpdb_interaction_input)
}
