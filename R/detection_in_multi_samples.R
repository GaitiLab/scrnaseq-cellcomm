#' @title Add detection flags
#' @description adds columns indicating how many patients/samples detected a source-target x interaction combination, within a voting strategy and condition/group. Also adds columns with the exact patients/samples where the combination was found.
#' @param df dataframe
#' @return df dataframe
add_detection_flags <- function(df) {
    return(
        df |>
            duckplyr::as_duckdb_tibble() |>
<<<<<<< HEAD
            dplyr::group_by(
                type_of_vote,
                condition,
                source_target,
                complex_interaction
            ) |>
=======
            # TODO remove later, in favor of using .by =
            # dplyr::group_by(
            #     type_of_vote,
            #     condition,
            #     source_target,
            #     complex_interaction
            # ) |>
>>>>>>> 56557e8e7e4633569c0cf9e182e7c14d86197058
            dplyr::summarise(
                # Is an interaction x source-target pair identified in multiple sample_ids from the same patient_id?
                detected_in_n_sample_ids = dplyr::n_distinct(sample_id),
                detected_in_sample_ids = paste(
                    sort(unique(sample_id)),
                    collapse = ", "
                ),

                # Is an interaction x source-target pair identified in multiple patient_ids within the condition for lenient/stringent vote
                detected_in_n_patient_ids = dplyr::n_distinct(patient_id),
                detected_in_patient_ids = paste(
                    sort(unique(patient_id)),
                    collapse = ", "
                ),
<<<<<<< HEAD
            ) |>
=======
                .by = c(
                    type_of_vote,
                    condition,
                    source_target,
                    complex_interaction
                )
            ) |>
            # TODO remove later, will be redundant if .by argument is used
>>>>>>> 56557e8e7e4633569c0cf9e182e7c14d86197058
            dplyr::ungroup()
    )
}
#' @title Filter by detection in multiple samples
#' @description Filtering detected interactions by looking at recurrence in multiple samples/patients. Keeping source-target x interaction pairs that were found in multiple patients (default >2 ) within a condition/group. Performing filtering for the two voting strategies 'lenient' and 'stringent' (based on detection in multiple tools)
#' @param df dataframe from 'samples_interactions_mvoted.rds'
#' @param min_patients Minimum number of patients for an interaction to be kept (default = 2)
#' @return dataframe
#' @export
filter_by_detection_in_multi_samples <- function(
    df,
    min_patients = 2
) {
    cols_oi <- c(
        "patient_id",
        "sample_id",
        "condition",
        "source_target",
        "complex_interaction"
    )

    df_long <- df |>
        dplyr::select(dplyr::all_of(cols_oi), dplyr::ends_with("_voting")) |>
        dtplyr::lazy_dt() |>
        # Conversion for easy data wrangling
        tidyr::pivot_longer(
            cols = dplyr::ends_with("_voting"),
            values_to = "is_voted",
            names_to = "type_of_vote",
        ) |>
        # Should only pass interactions that have been identified in multiple tools as defined in take_consensus_across_tools.R, either by 'lenient' or 'stringent' filter
        dplyr::filter(is_voted) |>
        # Column remains unused after
        dplyr::select(-is_voted) |>
        dplyr::collect() |>
        duckplyr::as_duckdb_tibble() |>
        # for simplification, just keep 'lenient' or 'stringent' as label
        dplyr::mutate(
            type_of_vote = stringr::str_remove(type_of_vote, "_voting")
        )

    df_long_w_detection_flags <- df_long |>
        add_detection_flags() |>
        dplyr::mutate(
            detected_in_enough_patients = detected_in_n_patient_ids >=
                min_patients
        )
    message(
        "Added columns indicating in how many and in which sample_ids and patient_ids a source-target x interaction combination was detected for each condition/group."
    )

    df_wide_filtered_for_multi_patient_detection <- df_long_w_detection_flags |>
        # For boosting confidence, only keep source-target pair x interaction combination that were found in multiple patients WITHIN conditions/groups.
        # In case patients can have multiple sample_ids, then a source-target pair x interaction combination may be found in multiple samples of the same patient, but not captures in a different patient, then this finding may not be robust.
        dplyr::filter(detected_in_enough_patients) |>
        tidyr::pivot_wider(
            names_from = type_of_vote,
            names_glue = "{type_of_vote}_voting_{.value}",
            # These are the only columns affected by voting
            values_from = dplyr::starts_with("detected_in"),
            # Need to use original value column in list
            values_fill = list(detected_in_enough_patients = FALSE)
            # NOTE in the resulting df, only two cominations for `detected_in_enough_patients` should be
            # (1) lenient_voting = TRUE, stringent_voting = FALSE and
            # (2) lenient_voting = TRUE, stringent_voting = FALSE (i.e. for some sample_id the source-target x interaction pair may not be detected in all methods/cci tools).
            # If the option (1), then the columns starting with 'stringent_voting' will be NA.
        )
    message(paste(
        "Filtered for detection of a source-target x interaction combination in at least",
        min_patients,
        "within each group/condition."
    ))

    if (identical(df$sample_id, df$patient_id)) {
        # If you don't have a nested structure, i.e. patients having multiple samples, then the sample_id and patient_id columns are the same, therefore remove the sample_id related columns.
        df_wide_filtered_for_multi_patient_detection <- df_wide_filtered_for_multi_patient_detection |>
            dplyr::select(-ends_with("sample_ids"))
        message(
            "The columns 'sample_id' and 'patient_id' are the same, removed columns related to sample_ids."
        )
    }
    message("Finished!")
    return(df_wide_filtered_for_multi_patient_detection)
}
