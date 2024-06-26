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

#' @title Ranking interactions for consensus using all tools
#' @param cell_type_pair cell type of interest for filtering
#' @param connectome dataframe with connectome results (interactions ordered by descending score)
#' @param sca dataframe with sca results (interactions ordered by descending score)
#' @param cytotalk dataframe with cytotalk results (interactions ordered by descending score)
#' @param natmi dataframe with natmi results (interactions ordered by descending score)
#' @param logfc dataframe with logfc results (interactions ordered by descending score)
#' @param obj_cell2cell dataframe with cell2cell results (interactions ordered by (ascending) p-values)
#' @param obj_cellchat dataframe with cellchat results (interactions ordered by (ascending) p-values)
#' @param obj_cpdb dataframe with cellphonedb results (interactions ordered by (ascending) p-values)
#' @importFrom dplyr %>%
#' @export
rank_interactions <- function(cell_type_pair, connectome, sca, cytotalk, natmi, logfc, obj_cell2cell, obj_cellchat, obj_cpdb) {
    ordered_interactions <- list(
        connectome = connectome %>%
            filter(source_target == cell_type_pair) %>%
            dplyr::pull(complex_interaction),
        sca = sca %>%
            filter(source_target == cell_type_pair) %>%
            dplyr::pull(complex_interaction),
        cytotalk = cytotalk %>%
            filter(source_target == cell_type_pair) %>%
            dplyr::pull(complex_interaction),
        natmi = natmi %>%
            filter(source_target == cell_type_pair) %>%
            dplyr::pull(complex_interaction),
        logfc = logfc %>%
            filter(source_target == cell_type_pair) %>%
            dplyr::pull(complex_interaction),
        cell2cell = obj_cell2cell %>%
            filter(source_target == cell_type_pair) %>%
            dplyr::pull(complex_interaction),
        cellchat = obj_cellchat %>%
            filter(source_target == cell_type_pair) %>%
            dplyr::pull(complex_interaction),
        cpdb = obj_cpdb %>%
            filter(source_target == cell_type_pair) %>%
            dplyr::pull(complex_interaction)
    )

    r <- RobustRankAggreg::rankMatrix(ordered_interactions, full = TRUE)
    #                connectome         sca   cytotalk       natmi       logfc cell2cell  cellchat      cpdb
    # LGALS9__PTPRK 0.001605136 0.096308186 0.13643660 0.011235955 0.027287319 0.3103652 0.9557242 0.5319594
    # VSIR__IGSF11  0.003210273 0.104333868 0.08025682 0.009630819 0.014446228        NA 0.9993675 0.9832736
    # TGFB1__ITGB8  0.004815409 0.117174960 0.14125201 0.008025682 0.008025682 0.4568978 0.9658444 0.3416965
    # PDGFB__PDGFRA 0.006420546 0.101123596 0.07062600 0.004815409 0.003210273 0.8428212 0.9917774 0.4447431
    # C3__LRP1      0.008025682 0.006420546 0.13322632 0.165329053 0.006420546 0.7042733 0.9807084 0.1412784
    # SPP1__ITGA9   0.009630819 0.075441413 0.09630819 0.062600321 0.001605136 0.8926512 0.9946237 0.3578256

    ranked_interactions <- RobustRankAggreg::aggregateRanks(rmat = r, method = "RRA")
    #                                  Name        Score
    # NRG2__ERBB2:ERBB3   NRG2__ERBB2:ERBB3 2.000351e-07
    # NRG2__ERBB2:ERBB4   NRG2__ERBB2:ERBB4 8.001405e-07
    # ANXA1__FPR2:FPR3     ANXA1__FPR2:FPR3 5.000878e-06
    # THY1__ITGAM:ITGB2   THY1__ITGAM:ITGB2 9.801721e-06
    # THY1__ITGAX:ITGB2   THY1__ITGAX:ITGB2 1.280225e-05
    # ICAM1__ITGAM:ITGB2 ICAM1__ITGAM:ITGB2 1.620285e-05
    return(ranked_interactions %>% dplyr::mutate(source_target = cell_type_pair))
}
