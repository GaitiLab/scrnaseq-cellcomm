#' @title Prepare data for CCI inference
#' @description filters data and then normalizes the gene expression
#' @param input_file path to seurat object (rds file)
#' @param output_dir output directory for saving output (default = '.')
#' @param annot variable in metadata containing the cell annotation
#' @param min_cells Minimum number of cells required in each cell group for cell-cell communication (default = 5)
#' @param sample_id sample id to use for saving formatted CCI results (default = NULL; will determine sample id based on 'input_interactions' for this input_interactions has to be of the format 'cpdb__{sample_id}.csv')
#' @export
prepare_data <- function(
    input_file, annot,
    sample_id = NULL, output_dir = ".",
    is_confident = FALSE, min_cells = 5) {
    #  Sanity checks
    if (!(file.exists(input_file) && endsWith(input_file, ".rds"))) {
        stop("File does not exists or is not an RDS object")
    }
    if (min_cells < 5) {
        stop("Min cells has to be >= 5...")
    }

    message("Loading Seurat object...")
    seurat_obj <- readRDS(input_file)

    if (is.null(annot) && !(annot %in% colnames(seurat_obj@meta.data))) {
        stop("Given annotation not in Seurat object")
    }

    message("Create output directories...")
    output_seurat <- glue::glue("{output_dir}/seurat")
    output_mtx <- glue::glue("{output_dir}/mtx")
    GaitiLabUtils::create_dir(output_seurat)
    GaitiLabUtils::create_dir(output_mtx)


    # Only relevant for internal project
    if (is_confident) {
        seurat_obj <- subset(seurat_obj, subset = Confident_Annotation)
    }
    message(glue::glue("Only keep cell type groups with at least {min_cells} cells"))
    seurat_obj <- filtering(
        seurat_obj,
        annot = annot,
        min_cells = min_cells
    )

    message(glue::glue(
        "Check number of cell types after filtering (>= 2)..."
    ))
    # Determine number of cell types present with at least min_cells
    n_cell_types <- length(unique(seurat_obj@meta.data[[annot]]))

    if (is.null(sample_id)) {
        output_name <- stringr::str_split(
            GaitiLabUtils::get_name(input_file), "__",
            simplify = TRUE
        )
    } else {
        output_name <- sample_id
    }
    if (n_cell_types < 2) {
        message("Not enough cell types present (at least 2 necessary)...")
    } else {
        message("Normalizing data...")
        seurat_obj <- Seurat::NormalizeData(seurat_obj)

        # Ensure factor only contain cell types that are included in the object (aka passed the filtering)
        metadata <- seurat_obj@meta.data
        metadata[, annot] <- factor(metadata[, annot], levels = unique(metadata[, annot]))
        seurat_obj <- Seurat::AddMetaData(seurat_obj, metadata = metadata)


        message("Saving Seurat object...")
        saveRDS(
            seurat_obj,
            glue::glue("{output_seurat}/{output_name}.rds")
        )
        message("Convert to mtx format (for Cell2Cell) and save...")
        mtx_output_dir <- glue::glue("{output_mtx}/{output_name}")
        if (file.exists(mtx_output_dir)) {
            file.remove(mtx_output_dir)
        }

        mat <- seurat_obj[["RNA"]]@data
        DropletUtils::write10xCounts(
            mtx_output_dir,
            mat
        )
    }
    message("Finished...")
}
