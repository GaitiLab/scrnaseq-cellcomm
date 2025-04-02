#' @title Update CellChat DB using (updated) CellPhoneDB
#' @param geneInfo dataframe object containing contents from gene_input.csv from CellPhoneDB
#' @param output_dir path to directory for saving the generated CellPhoneDB database (.rds)
#' @param db database
#' @export
#' @importFrom dplyr %>%
update_cellchat_db <- function(geneInfo, db, output_dir) {
    # Ref: https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Update-CellChatDB.html
    cellchat_db <- CellChat::CellChatDB.human
    cellchat_db_interaction <- cellchat_db$interaction
    cellchat_db_complex <- cellchat_db$complex
    cellchat_db_cofactor <- cellchat_db$cofactor
    cellchat_db_geneinfo <- cellchat_db$geneInfo

    geneInfo$Symbol <- geneInfo$hgnc_symbol
    geneInfo$Symbol <- geneInfo$hgnc_symbol
    geneInfo <- dplyr::select(geneInfo, -c("ensembl"))
    geneInfo <- unique(geneInfo)

    geneInfo <- geneInfo %>%
        dplyr::rename(EntryID.uniprot = uniprot) %>%
        dplyr::select(EntryID.uniprot, Symbol) %>%
        dplyr::mutate()

    # Get missing columns
    missing_cols <- setdiff(colnames(cellchat_db_geneinfo), colnames(geneInfo))
    geneInfo[, missing_cols] <- ""

    cellchat_db_geneinfo <- rbind(cellchat_db_geneinfo, geneInfo) %>%
        dplyr::distinct(Symbol, .keep_all = TRUE)

    # Create new CellChatDB
    cellchatDB_Omni <- liana:::cellchat_formatDB(
        ccDB = list(
            complex = cellchat_db_complex,
            geneInfo = cellchat_db_geneinfo,
            interaction = cellchat_db_interaction,
            cofactor = cellchat_db_cofactor
        ),
        op_resource = db,
        exclude_anns = c()
    )

    # Save CellChatDB
    saveRDS(cellchatDB_Omni, file = glue::glue("{output_dir}/cellchat_db.rds"))
}
