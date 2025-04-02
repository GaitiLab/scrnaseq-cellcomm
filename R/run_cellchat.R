#' @title Run CellChat
#' @param gene_expr seurat object with gene expression (rds file)
#' @param annot Column in metadata containing the cell type labels
#' @param interactions_db (custom) CellChat database
#' @param output_dir output directory for saving output (default = '.')
#' @param min_cells Minimum number of cells required in each cell group for cell-cell communication (default = 5)
#' @param n_perm Number of permutations for permutation testing (default = 1000)
#' @export
run_cellchat <- function(
    gene_expr,
    annot,
    interactions_db,
    output_dir = ".",
    min_cells = 5,
    n_perm = 1000) {
    options(future.globals.maxSize = 8000 * 1024**2)
    #  Sanity checks
    if (!(file.exists(gene_expr) && endsWith(tolower(gene_expr), ".rds"))) {
        stop(
            "Seurat object ('gene_expr') does not exists or is not an RDS object"
        )
    }
    if (
        !(file.exists(interactions_db) &&
            endsWith(tolower(interactions_db), ".rds"))
    ) {
        stop(
            "Interactions database ('interactions_db') does not exists or is not an RDS object"
        )
    }
    if (!file.exists(output_dir)) {
        stop("Output directory does not exist")
    }
    if (min_cells < 5) {
        stop("Min cells has to be >= 5...")
    }
    message("Load data...")
    seurat_obj <- readRDS(gene_expr)

    message("Extract gene expression and convert to matrix...")
    mat <- as.matrix(seurat_obj@assays$RNA@data)
    meta <- seurat_obj@meta.data
    meta[, annot] <- factor(meta[, annot])

    message("Create CellChat object...")
    cellchat <- CellChat::createCellChat(
        object = mat,
        meta = meta,
        group.by = annot
    )
    cellchat <- CellChat::addMeta(cellchat, meta = meta)
    cellchat <- CellChat::setIdent(cellchat, ident.use = annot) # set 'labels' as default cell identity

    message("Load custom database with interactions...")
    cellchat@DB <- readRDS(interactions_db)

    message("Preprocessing the expression data...")
    cellchat <- CellChat::subsetData(cellchat) # This step is necessary even if using the whole database

    cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
    cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)

    message("Infer cell-cell interactions...")
    cellchat <- CellChat::computeCommunProb(
        cellchat,
        nboot = n_perm,
        population.size = TRUE
    )
    cellchat <- CellChat::filterCommunication(cellchat, min.cells = min_cells)

    # b. aggregated cell-cell communication network
    cellchat <- CellChat::aggregateNet(cellchat)

    out_filename <- GaitiLabUtils::get_name(gene_expr)
    message("Save CellChat object...")
    saveRDS(
        cellchat,
        file = glue::glue("{output_dir}/cellchat__{out_filename}__raw_obj.rds")
    )

    message("Post-processing...")
    interactions <- names(cellchat@net$prob[1, 1, ])
    res <- pbapply::pblapply(interactions, function(interaction) {
        # Handle probabilities
        cci <- reshape2::melt(cellchat@net$prob[, , interaction], )
        colnames(cci) <- c("source", "target", "proba")
        cci["interaction"] <- interaction

        # Handle pvalues
        pval_long <- reshape2::melt(cellchat@net$pval[, , interaction])
        colnames(pval_long) <- c("source", "target", "pval")
        cci["pval"] <- pval_long$pval
        return(cci)
    })
    message("Concatenate results...")
    res_concat <- do.call("rbind", res)

    message("Save formatted CellChat results...")
    saveRDS(
        res_concat,
        glue::glue("{output_dir}/cellchat__{out_filename}.rds"),
    )
    message("Finished...")
}
