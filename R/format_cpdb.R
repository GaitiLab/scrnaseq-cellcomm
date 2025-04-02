#' @title Format cpdb CCI results
#' @param interaction_scores path to 'statistical_analysis_interaction_scores' txt file
#' @param pval path to 'statistical_analysis_pvalues' txt file
#' @param sign_means path to 'statistical_analysis_significant_means' txt file
#' @param means path to 'statistical_analysis_means' txt file
#' @param output_dir output directory for saving output (default = '.')
#' @param sample_id sample id to use for saving formatted CCI results (default = NULL; will determine sample id based on 'input_interactions' for this input_interactions has to be of the format 'cpdb__{sample_id}.csv')
#' @param ref_db Path to interactions database (default = "data/interactions_db/ref_db.rds")
#' @export
#' @importFrom dplyr %>%
format_cpdb <- function(
    interaction_scores,
    pval,
    sign_means,
    means,
    ref_db = "data/interactions_db/ref_db.rds",
    output_dir = ".",
    sample_id = NULL) {
    #  Sanity checks
    if (
        !(file.exists(interaction_scores) &&
            endsWith(tolower(interaction_scores), ".txt") &&
            stringr::str_detect(
                interaction_scores,
                "statistical_analysis_interaction_scores"
            ))
    ) {
        stop(
            "'statistical_analysis_interaction_scores' file ('interaction_scores') does not exists or is not a txt file"
        )
    }
    if (
        !(file.exists(pval) &&
            endsWith(tolower(pval), ".txt") &&
            stringr::str_detect(pval, "statistical_analysis_pvalues"))
    ) {
        stop(
            "'statistical_analysis_pvalues' file ('pval') does not exists or is not a txt file"
        )
    }
    if (
        !(file.exists(sign_means) &&
            endsWith(tolower(sign_means), ".txt") &&
            stringr::str_detect(
                sign_means,
                "statistical_analysis_significant_means"
            ))
    ) {
        stop(
            "'statistical_analysis_significant_means' file ('sign_means') does not exists or is not a txt file"
        )
    }
    if (
        !(file.exists(means) &&
            endsWith(tolower(means), ".txt") &&
            stringr::str_detect(means, "statistical_analysis_means"))
    ) {
        stop(
            "'statistical_analysis_means' file ('means') does not exists or is not a txt file"
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

    message("Load CellPhoneDB output...")
    pval <- read.table(pval, sep = "\t", header = TRUE, check.names = FALSE)
    sign_means <- read.table(
        sign_means,
        sep = "\t",
        header = TRUE,
        check.names = FALSE
    )
    means <- read.table(means, sep = "\t", header = TRUE, check.names = FALSE)
    interaction_scores <- read.table(
        interaction_scores,
        sep = "\t",
        header = TRUE,
        check.names = FALSE
    )

    message("Process 'pval'...")
    colnames_oi <- colnames(pval)[stringr::str_detect(colnames(pval), "\\|")]
    pval <- pval %>%
        dplyr::select(interacting_pair, dplyr::all_of(colnames_oi)) %>%
        reshape2::melt(
            "interacting_pair",
            value.name = "pval",
            variable.name = "source_target"
        ) %>%
        dplyr::mutate(
            source_target = stringr::str_replace_all(source_target, "\\|", "__")
        )

    message("Process 'sign_means'...")
    colnames_oi <- colnames(sign_means)[stringr::str_detect(
        colnames(sign_means),
        "\\|"
    )]
    sign_means <- sign_means %>%
        dplyr::select(interacting_pair, rank, dplyr::all_of(colnames_oi)) %>%
        reshape2::melt(
            c("interacting_pair", "rank"),
            value.name = "sign_mean",
            variable.name = "source_target"
        ) %>%
        dplyr::mutate(
            source_target = stringr::str_replace_all(source_target, "\\|", "__")
        )

    message("Process 'means'...")
    colnames_oi <- colnames(means)[stringr::str_detect(colnames(means), "\\|")]
    means <- means %>%
        dplyr::select(interacting_pair, dplyr::all_of(colnames_oi)) %>%
        reshape2::melt(
            c("interacting_pair"),
            value.name = "mean",
            variable.name = "source_target"
        ) %>%
        dplyr::mutate(
            source_target = stringr::str_replace_all(source_target, "\\|", "__")
        )

    message("Process 'interaction_scores'...")
    colnames_oi <- colnames(interaction_scores)[stringr::str_detect(
        colnames(interaction_scores),
        "\\|"
    )]
    interaction_scores <- interaction_scores %>%
        dplyr::select(interacting_pair, dplyr::all_of(colnames_oi)) %>%
        reshape2::melt(
            c("interacting_pair"),
            value.name = "interaction_score",
            variable.name = "source_target"
        ) %>%
        dplyr::mutate(
            source_target = stringr::str_replace_all(source_target, "\\|", "__")
        )

    message("Merge CellPhoneDB output...")
    interactions <- pval %>%
        dplyr::left_join(
            sign_means,
            by = c("interacting_pair", "source_target")
        ) %>%
        dplyr::left_join(means, by = c("interacting_pair", "source_target")) %>%
        dplyr::left_join(
            interaction_scores,
            by = c("interacting_pair", "source_target")
        ) %>%
        dplyr::left_join(ref_db, by = c("interacting_pair" = "interaction")) %>%
        dplyr::mutate(method = "CellPhoneDBv5", Sample = sample_id) %>%
        dplyr::rename(CellPhoneDB_score = interaction_score)

    message("Save output...")
    saveRDS(
        interactions,
        glue::glue("{output_dir}/cpdb__{sample_id}__postproc.rds")
    )

    message("Finished!")
}
