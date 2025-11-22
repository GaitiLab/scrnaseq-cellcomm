#' @title Run LIANA
#' @description Uses a Seurat object to run LIANA, follows the standard LIANA workflow.
#' @param gene_expr_path seurat object with gene expression (rds file). This object should contain the RNA assay.
#' @param annot Column in metadata containing the cell type labels
#' @param interactions_db_path path to (custom) LIANA database which should be an '.rds' file
#' @param min_cells Minimum number of cells required in each cell group for cell-cell communication (default = 5)
#' @param n_perm Number of permutations for permutation testing (default = 1000)
#' @param min_pct Minimum fraction of cells expressing a gene (default = 0.1 = 10%), should be a value between 0 and 1
#' @return liana_obj returns output of `liana::liana_wrap`
#' @export
RunLIANA <- function(
    gene_expr_path,
    interactions_db_path,
    annot,
    min_cells = 5,
    min_pct = 0.1,
    n_perm = 1000) {
    # ---- Define constants
    # Run all method
    methods <- c("natmi", "connectome", "logfc", "sca", "cytotalk")
    supp_columns <- c("ligand.expr", "receptor.expr")
    # Define the no. of permutations for permutation testing when applicable
    permutation_params <- list(
        nperms = n_perm
    )
    # ---- Perform sanity checks ----
    # is_valid_filepath(gene_expr_path)
    # is_valid_filepath(interactions_db_path)
    # Minimum of 5 cells enforced/required by LIANA
    if (min_cells < 5) {
        stop("Min cells has to be >= 5...")
    }
    # In documentation of LIANA `min_pct` actually represents a fraction/proportion, not a percentage. Therefore value should not be greater than 1.
    if (min_pct > 1) {
        stop("min_pct > 1...")
    }

    assay <- "RNA"
    if (!assay %in% Seurat::Assays(assay)) {
        stop("`RNA` assay is not present...")
    }

    # ---- Loading data
    message("Loading Seurat object...")
    seurat_obj <- readRDS(gene_expr_path)

    message("Loading database with interactions...")
    custom_resource <- readRDS(interactions_db_path)

    # ---- Run LIANA
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
    message("Finished...")
    return(liana_obj)
}
