#' @title Combine samples
#' @description hard-combine results per sample.
#' @param input_dir df
#' @param metadata path to file with corresponding metadata (RDS)
#' @param patient_var patient variable name in metadata
#' @param sample_var sample variable name in metadata (default = "Sample")
#' @param condition_var condition variable name in metadata (default = "Condition_dummy")
#' @export
combine_samples <- function(
    df,
    metadata,
    patient_var = "sample_id",
    condition_var = "condition_dummy",
    sample_var = "sample_id"
) {
    # Variables from meta_df to add
    cols_oi <- unique(c(
        sample_var,
        condition_var,
        patient_var
    ))

    # Only keep sample-level meta_df (not single-cell)
    meta_df <- meta_df |>
        dplyr::select(dplyr::any_of(cols_oi)) |>
        tibble::remove_rownames() |>
        dplyr::distinct()

    if (condition_var == "condition_dummy") {
        meta_df <- meta_df |> dplyr::mutate(condition = NA)
        condition_var <- "condition"
    }

    if (patient_var == sample_var) {
        meta_df <- meta_df |>
            dplyr::mutate(patient_id = !!dplyr::sym(sample_var)) |>
            dplyr::rename(
                sample_id = sample_var,
                condition = condition_var
            )
    }

    df <- df |> dplyr::left_join(meta_df)
    message("Finished!")
    return(df)
}
