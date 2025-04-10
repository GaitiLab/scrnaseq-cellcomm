#' @title Format CellChat CCI results
#' @param input_interactions path to file with the inferred interactions from CellChat (output file generated by run_cellchat()) (rds file)
#' @param output_dir output directory for saving output (default = '.')
#' @param sample_id sample id to use for saving formatted CCI results (default = NULL; will determine sample id based on 'input_interactions' for this input_interactions has to be of the format 'cellchat__{sample_id}.rds')
#' @param ref_db Path to interactions database (default = "data/interactions_db/ref_db.rds")
#' @export
#' @importFrom dplyr %>%
format_cellchat <- function(
    input_interactions,
    ref_db = "data/interactions_db/ref_db.rds",
    output_dir = ".",
    sample_id = NULL) {
    #  Sanity checks
    if (
        !(file.exists(input_interactions) &&
            endsWith(tolower(input_interactions), ".rds"))
    ) {
        stop(
            "Interactions file ('input_interactions') does not exists or is not an RDS object"
        )
    }
    # TODO in future release remove, when switching to new database
    if (!file.exists(ref_db) && endsWith(tolower(input_interactions), ".rds")) {
        stop(glue::glue("{ref_db} path does not exist"))
    }

    if (!file.exists(output_dir)) {
        stop("Output directory does not exist")
    }

    if (is.null(sample_id)) {
        # Setting sample + run id depending on downsampling or not
        split_sample_id <- stringr::str_split(sample_id, "__", simplify = TRUE)
        sample_id <- split_sample_id[1]
    }
    message(glue::glue("Formatting sample={sample_id}..."))

    message("Load database...")
    # TODO in future release remove, when switching to new database
    ref_db <- readRDS(ref_db) %>%
        dplyr::select(
            # simple_interaction,
            complex_interaction,
            interaction
        )

    message("Load interactions")
    interactions <- readRDS(input_interactions) %>%
        dplyr::rename(CellChat_score = proba) %>%
        dplyr::left_join(ref_db, by = "interaction") %>%
        dplyr::select(-interaction) %>%
        tidyr::unite(source_target, source, target, sep = "__") %>%
        dplyr::mutate(method = "CellChatv2", Sample = sample_id)

    message("Save output...")
    saveRDS(
        interactions,
        glue::glue("{output_dir}/cellchat__{sample_id}__postproc.rds")
    )
    message("Finished!")
}
