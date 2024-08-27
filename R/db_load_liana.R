#' @title Load LIANA databases
#' @description Loads and combines the LIANA Consensus and Ramilowski 2015 database (hard combine, missing columns set to NA)
#' @return dataframe with the two databases
#' @export
#' @importFrom dplyr %>%
load_liana_dbs <- function() {
    # Loading Consensus database
    liana_consensus <- liana::select_resource("Consensus")[[1]] %>% dplyr::mutate(method = "LIANA consensus")

    # load Ramilowski 2015 database
    liana_ramilowski <- liana::select_resource("Ramilowski2015")[[1]] %>% dplyr::mutate(method = "LIANA Ramilowski2015")

    # Account for missing columns
    missing_liana <- setdiff(colnames(liana_consensus), colnames(liana_ramilowski))
    missing_ramilowski <- setdiff(colnames(liana_ramilowski), colnames(liana_consensus))

    liana_consensus[missing_ramilowski] <- ""
    liana_ramilowski[missing_liana] <- ""

    # Combining the two LIANA databases
    liana_db <- rbind(liana_consensus, liana_ramilowski)

    return(liana_db)
}
