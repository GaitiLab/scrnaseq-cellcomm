#' @title Combine samples
#' @description hard-combine results per sample.
#' @param df dataframe
#' @param metadata dataframe with metadata
#' @param patient_var patient variable name in metadata
#' @param sample_var sample variable name in metadata (default = "Sample")
#' @param condition_var condition variable name in metadata (default = "Condition_dummy")
#' @export
AddMetaData <- function(
    df,
    meta_df,
    patient_var = "sample_id",
    condition_var = "condition_dummy",
    sample_var = "sample_id") {
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
        lookup <- setNames(c(sample_var, condition_var), c("sample_id", "condition"))
        meta_df <- meta_df |>
            dplyr::mutate(patient_id = !!dplyr::sym(sample_var)) |>
            dplyr::rename(dplyr::all_of(lookup))
    } else {
        lookup <- setNames(c(patient_var, sample_var, condition_var), c("patient_id", "sample_id", "condition"))
        meta_df <- meta_df |>
            dplyr::rename(dplyr::all_of(lookup))
    }
    return(df)
}
