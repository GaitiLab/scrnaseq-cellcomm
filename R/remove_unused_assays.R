#' @title Reduce seurat object size
#' @description remove redundant layers or assays, i.e. only keeping the 'RNA' assay
#' @param seurat_obj Seurat object
#' @export
remove_unused_assays <- function(seurat_obj) {
    Seurat::DefaultAssay(seurat_obj) <- "RNA"

    # Removing all assays except RNA
    for (assay_name in names(seurat_obj)) {
        if (assay_name == "RNA") {
            next
        } else {
            message(glue::glue("Removing assay: {assay_name}..."))
            try(seurat_obj[[assay_name]] <- NULL, FALSE)
        }
    }
    return(seurat_obj)
}
