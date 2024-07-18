#' @title Get metadata 
#' @description Extracts the metadata and saves it as and '.rds' and '.csv' file in user-specified output directory.
#' @param input_file path to seurat object (rds file)
#' @param output_dir output directory for saving output (default = '.')
#' @return NULL
#' @examples 
#' \dontrun{
#' filepath <- "/johndoe/data/seurat_obj.rds"
#' output_dir <- "/johndoe/output"
#' get_metadata(filepath, output_dir)
#' }
#' @export
get_metadata <- function (input_file, output_dir = ".") {
    if (!(file.exists(input_file) && endsWith(input_file, ".rds"))) { 
        stop("File does not exists or is not an RDS object")
    }
    if (!file.exists(output_dir)) {
        stop("Output directory does not exist")
    }
    message("Load input file...")
    seurat_obj <- readRDS(input_file)
    out_filename <- GaitiLabUtils::get_name(input_file)
    message('Save metadata as rds file...')
    saveRDS(seurat_obj@meta.data, glue::glue("{output_dir}/{out_filename}__metadata.rds"))
    message("Save metadata as csv file...")
    write.csv(seurat_obj@meta.data, glue::glue("{output_dir}/{out_filename}__metadata.csv"))
    message("Finished...")
}

#' @title Reduce seurat object size 
#' @description remove redundant layers or assays, i.e. only keeping the 'RNA' assay
#' @param input_file path to seurat object (rds file)
#' @param output_dir output directory for saving output (default = '.')
#' @export 
reduce_seurat_obj <- function (input_file, output_dir = ".") {
    if (!(file.exists(input_file) && endsWith(input_file, ".rds"))) {
          stop("File does not exists or is not an RDS object")
      }
    if (!file.exists(output_dir)) {
          stop("Output directory does not exist")
      }
    message("Load input file...")
    seurat_obj <- readRDS(input_file)
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
    # Set a name for the output file
    out_filename <- GaitiLabUtils::get_name(
        stringr::str_split(GaitiLabUtils::get_name(input_file), "__", simplify = TRUE)[1])

    message("Saving reduced Seurat object...")
    saveRDS(seurat_obj, glue::glue("{output_dir}/{out_filename}_reduced_size.rds"))
    message("Finished...")
}

#' @title Split seurat object by variable
#' @description Splits seurat object by a user-specified variable, e.g. Sample
#' @param input_file path to seurat object (rds file)
#' @param output_dir output directory for saving output (default = '.')
#' @param sample_var variable to split the object by, must be part of metadata (default = "Sample")
#' @export
split_seurat_into_samples <- function (input_file, output_dir = ".", sample_var = "Sample") { 
    message("Loading Seurat object...")
    seurat_obj <- readRDS(input_file)
    print(seurat_obj)
    obj_list <- Seurat::SplitObject(seurat_obj, split.by = sample_var)
    obj_ids <- names(obj_list)

    message(glue::glue("Split object by {sample_var}..."))
    for (obj_id in obj_ids) {
        message(glue::glue("Sample: {obj_id}..."))
        obj <- obj_list[[obj_id]]

        message(glue::glue("Saving {sample_var}: {obj_id}"))
            saveRDS(
                obj,
                glue::glue("{output_dir}/{obj_id}.rds")
            )
    }

}

#' @title Filtering the Seurat object to reduce size of object
#' @description  1. Only keep genes that are part of the interactions db, 2. Only keep cell types with at least N cells (user-defined)
#' @param seurat_obj Seurat object
#' @param annot Annotation to use for filtering
#' @param min_cells Minimum number of cells per annotation (default = 5)
#' @param genes_oi list of genes of interest from database
#' @return seurat_obj Filtered seurat object
#' @examples dontrun{seurat_obj <- filtering(seurat_obj, "custom_annot", 300, genes_oi)}
#' @export
filtering <- function(seurat_obj, annot, min_cells = 5) {
    if(min_cells < 5) {
        stop("Min_cells should be >= 5")
    }
    # Get annotations of cells (rownames = cell_id)
    cells_annotated <- seurat_obj@meta.data[annot]
    counts_per_label <- table(cells_annotated)
    labels_to_keep <- names(counts_per_label)[counts_per_label >= min_cells]
    cells_to_keep <- rownames(cells_annotated)[cells_annotated[[annot]] %in% labels_to_keep]
    seurat_obj <- subset(seurat_obj, cells = cells_to_keep)
    return(seurat_obj)
}