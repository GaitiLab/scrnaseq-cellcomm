#' @title Filter by detection in multiple samples
#' @param input_file Path to '401_samples_interactions_mvoted.rds' file
#' @param annot variable that contains the cell type labels (default = '')
#' @param min_patients Minimum number of patients for an interaction to be kept (default = 2)
#' @param condition_var Variabile containing the 'condition' (group/category) of the samples (default = "Condition_dummy")
#' @param output_dir output directory for saving output (default = '.')
#' @export
#' @importFrom dplyr %>%
filter_by_detection_in_multi_samples <- function(
    input_file, annot = "",
    min_patients = 2,
    condition_var = "Condition_dummy",
    output_dir = ".") {
    cols_oi <- c("Patient", "Sample", condition_var, "source_target", "complex_interaction")

    if (!file.exists(input_file)) {
        stop(glue::glue("{input_file} is not a valid path"))
    }
    message("Loading input file...")
    input_file <- readRDS(input_file)

    if (sum(stringr::str_detect(colnames(input_file), annot)) == 0) {
        stop(glue::glue("'{annot}' is not a column in the input_file"))
    }

    if (sum(stringr::str_detect(colnames(input_file), "Patient")) == 0) {
        message("Patient column missing, assuming Patient = Sample...")
        input_file <- input_file %>%
            dplyr::mutate(Patient = Sample)
    }

    # ---- Lenient voting ---- #
    # Check detection of the same interaction in the same patient across samples for the same region
    lenient_voting <- input_file %>%
        filter(lenient_voting) %>%
        dplyr::select(dplyr::all_of(cols_oi)) %>%
        dplyr::group_by_at(
            dplyr::vars((c(condition_var, "Patient", "source_target", "complex_interaction")))
        ) %>%
        dplyr::reframe(
            lenient_detected_same_patient = paste0(Sample, collapse = ", "), lenient_N_samples_same_patient = dplyr::n()
        )

    lenient_voting_by_condition <- lenient_voting %>%
        dplyr::group_by_at(
            dplyr::vars((c(condition_var, "source_target", "complex_interaction")))
        ) %>%
        dplyr::reframe(
            lenient_condition_patients = paste0(Patient, collapse = ", "), lenient_condition_n_patients = dplyr::n(),
            lenient_condition_n_samples = sum(lenient_N_samples_same_patient),
            lenient_condition_samples = paste0(lenient_detected_same_patient, collapse = ", ")
        )
    n_before <- nrow(lenient_voting_by_condition)

    message(glue::glue("Only keep interactions that are found in at least {min_patients} patients..."))
    lenient_voting_by_condition <- lenient_voting_by_condition %>% filter(lenient_condition_n_patients >= min_patients)
    n_after <- nrow(lenient_voting_by_condition)

    message(glue::glue("Before filtering: {n_before}"))
    message(glue::glue("After filtering: {n_after}"))

    # ---- Stringent voting ---- #
    stringent_voting <- input_file %>%
        filter(stringent_voting) %>%
        dplyr::select(dplyr::all_of(cols_oi)) %>%
        dplyr::group_by_at(dplyr::vars(
            (c(condition_var, "Patient", "source_target", "complex_interaction"))
        )) %>%
        dplyr::reframe(
            stringent_detected_same_patient = paste0(Sample, collapse = ", "), stringent_N_samples_same_patient = dplyr::n()
        )

    stringent_voting_by_condition <- stringent_voting %>%
        dplyr::group_by_at(dplyr::vars((c(condition_var, "source_target", "complex_interaction")))) %>%
        dplyr::reframe(
            stringent_condition_patients = paste0(Patient, collapse = ", "), stringent_condition_n_patients = dplyr::n(),
            stringent_condition_n_samples = sum(stringent_N_samples_same_patient), stringent_condition_samples = paste0(stringent_detected_same_patient, collapse = ", ")
        )

    n_before <- nrow(stringent_voting_by_condition)

    message(glue::glue("Only keep interactions that are found in at least {min_patients} patients..."))
    stringent_voting_by_condition <- stringent_voting_by_condition %>%
        filter(stringent_condition_n_patients >= min_patients)

    n_after <- nrow(stringent_voting_by_condition)

    message(glue::glue("Before filtering: {n_before}"))
    message(glue::glue("After filtering: {n_after}"))

    # ---- Combine stringent and lenient voting ---- #
    message("Combine stringent and lenient voting...")
    # Add column for visualization scripts later on
    lenient_voting_by_condition <- lenient_voting_by_condition %>%
        dplyr::mutate(lenient_condition = TRUE)

    stringent_voting_by_condition <- stringent_voting_by_condition %>%
        dplyr::mutate(stringent_condition = TRUE)

    combined_voting <- merge(lenient_voting_by_condition,
        stringent_voting_by_condition,
        by = c(condition_var, "source_target", "complex_interaction"),
        all = TRUE
    )

    message("Save results...")
    saveRDS(combined_voting, glue::glue("{output_dir}/402a_filtering_detect_in_multi_samples.rds"))
    message("Finished!")
}
