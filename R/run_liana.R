#' @title Run CellChat
#' @param gene_expr seurat object with gene expression (rds file)
#' @param annot Column in metadata containing the cell type labels
#' @param interactions_db (custom) CellChat database (rds file)
#' @param output_dir output directory for saving output (default = '.')
#' @param min_cells Minimum number of cells required in each cell group for cell-cell communication (default = 5)
#' @param n_perm Number of permutations for permutation testing (default = 1000)
#' @param min_pct Minimum fraction of cells expressing a gene (default = 0.1 = 10%), max = 1
#' @export
run_liana <- function(
    gene_expr,
    interactions_db,
    annot,
    output_dir = ".",
    min_cells = 5,
    min_pct = 0.1,
    n_perm = 1000) {
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
    if (min_pct > 1) {
        stop("min_pct > 1...")
    }

    # ---- Constants ----
    methods <- c("natmi", "connectome", "logfc", "sca", "cytotalk")
    supp_columns <- c("ligand.expr", "receptor.expr")
    permutation_params <- list(
        nperms = n_perm
    )
    assay <- "RNA"

    # ---- Loading data ----
    message("Loading Seurat object...")
    seurat_obj <- readRDS(gene_expr)

    message("Loading database with interactions...")
    custom_resource <- readRDS(interactions_db)

    # ---- Run LIANA ----
    liana_obj <- liana::liana_wrap(
        seurat_obj,
        method = methods,
        resource = "custom",
        external_resource = custom_resource,
        idents_col = annot,
        supp_columns = supp_columns,
        return_all = TRUE,
        permutation.params = permutation_params,
        assay = assay,
        min_cells = min_cells,
        expr_prop = min_pct
    )
    message("Save LIANA results...")
    out_filename <- GaitiLabUtils::get_name(gene_expr)
    saveRDS(
        liana_obj,
        file = glue::glue("{output_dir}/liana__{out_filename}.rds")
    )

    message("Finished...")
}
