#' @title Prepare data for CCI inference
#' @description Checks if there are enough cell types with at least `min_cells`, then continues with normalization. Groups (cell types) with insufficient no. cells are removed.
#' @param seurat_obj Seurat object
#' @param annot variable in metadata containing the cell annotation
#' @param min_cells Minimum number of cells required in each cell group for cell-cell communication (default = 5)
#' @return Seurat object
#' @export
prepare_data <- function(
    seurat_obj,
    annot,
    min_cells = 5
) {
    if (min_cells < 5) {
        stop("Min. cells has to be at least 5.")
    }
    if (is.null(annot) && !(annot %in% colnames(seurat_obj@meta.data))) {
        stop("Given annotation not in Seurat object")
    }

    seurat_obj <- seurat_obj |>
        remove_groups_with_low_ncells(
            annot = annot,
            min_cells = min_cells
        )
    message(paste("Filtered out groups with less than", min_cells, "cells."))

    # Determine number of cell types present with at least min_cells
    n_groups_with_enough_cells <- length(unique(seurat_obj@meta.data[[annot]]))
    message(paste(
        "No. groups (cell types) with enough cells:",
        n_groups_with_enough_cells
    ))
    if (n_groups_with_enough_cells < 2) {
        stop("Not enough cell types present (at least 2 required)...")
    } else {
        seurat_obj <- Seurat::NormalizeData(seurat_obj)
        message("Normalized data.")

        # Ensure factor only contain cell types that are included in the object (aka passed the filtering)
        metadata <- seurat_obj@meta.data
        metadata[, annot] <- factor(
            metadata[, annot],
            levels = unique(metadata[, annot])
        )
        seurat_obj <- Seurat::AddMetaData(seurat_obj, metadata = metadata)
    }
    message("Finished.")
    return(seurat_obj)
}
