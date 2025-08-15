#' @title Run CellChat
#' @description Uses a Seurat object to run CellChat, follows the standard CellChat workflow.
#' @param gene_expr_path seurat object with gene expression (rds file)
#' @param annot Column in metadata containing the cell type labels
#' @param interactions_db_path path to (custom) CellChat database which should be an '.rds' file.
#' @param min_cells Minimum number of cells required in each cell group for cell-cell communication (default = 5)
#' @param n_perm Number of permutations for permutation testing (default = 1000)
#' @return cellchat object
#' @export
run_cellchat <- function(
    gene_expr_path,
    interactions_db_path,
    annot,
    min_cells = 5,
    n_perm = 1000
) {
    options(future.globals.maxSize = 8000 * 1024**2)

    if (min_cells < 5) {
        # This is a CellChat specific constraint
        stop("Min cells has to be >= 5...")
    }
    message("Load data...")
    seurat_obj <- readRDS(gene_expr_path)

    message("Extract gene expression and convert to matrix...")
    gene_expr_mat <- as.matrix(seurat_obj@assays$RNA@data)
    metadata_df <- seurat_obj@meta.data
    metadata_df[, annot] <- factor(metadata_df[, annot])

    cc_object <- CellChat::createCellChat(
        object = gene_expr_mat,
        meta = metadata_df,
        group.by = annot
    )
    cc_object <- CellChat::addMeta(cc_object, meta = metadata_df)
    # set 'labels' as default cell identity
    cc_object <- CellChat::setIdent(cc_object, ident.use = annot)
    message("Created CellChat object...")

    cc_object@DB <- readRDS(interactions_db_path)
    message("Loaded custom database with interactions...")

    # This step is necessary even if using the whole database
    cc_object <- CellChat::subsetData(cc_object)

    cc_object <- CellChat::identifyOverExpressedGenes(cc_object)
    cc_object <- CellChat::identifyOverExpressedInteractions(cc_object)
    message("Preprocessed the expression data...")

    cc_object <- CellChat::computeCommunProb(
        cc_object,
        nboot = n_perm,
        population.size = TRUE
    )

    cc_object <- CellChat::filterCommunication(cc_object, min.cells = min_cells)

    # b. aggregated cell-cell communication network
    cc_object <- CellChat::aggregateNet(cc_object)
    message("Inferred cell-cell interactions...")

    message("Finished...")
    return(cc_object)
}
