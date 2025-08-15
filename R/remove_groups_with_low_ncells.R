#' @title Filtering the Seurat object to reduce size of object
#' @description Only keep cell types with at least N cells (user-defined)
#' @param seurat_obj Seurat object
#' @param annot Annotation to use for filtering
#' @param min_cells Minimum number of cells per annotation (default = 5)
#' @return seurat_obj filtered seurat object
#' @examples dontrun{seurat_obj <- filtering(seurat_obj, "custom_annot", 300, genes_oi)}
#' @export
remove_groups_with_low_ncells <- function(seurat_obj, annot, min_cells = 5) {
    if (min_cells < 5) {
        stop("Min_cells should be >= 5")
    }
    # Get annotations of cells (rownames = cell_id)
    cells_annotated <- seurat_obj@meta.data[annot]
    counts_per_label <- table(cells_annotated)
    labels_to_keep <- names(counts_per_label)[counts_per_label >= min_cells]
    cells_to_keep <- rownames(cells_annotated)[
        cells_annotated[[annot]] %in% labels_to_keep
    ]
    seurat_obj <- subset(seurat_obj, cells = cells_to_keep)
    return(seurat_obj)
}
