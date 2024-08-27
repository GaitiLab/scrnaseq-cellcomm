#' @title Filter by detection in multiple samples
#' @param input_file Path to '401_samples_interactions_mvoted.rds' file
#' @param min_patients Minimum number of patients for an interaction to be kept (default = 2)
#' @param condition_var Variabile containing the 'condition' (group/category) of the samples (default = "Condition_dummy")
#' @param output_dir output directory for saving output (default = '.')
#' @export
#' @importFrom dplyr %>%
filter_by_detection_in_multi_samples <- function(
    input_file,
    min_patients = 2,
    condition_var = "Condition_dummy",
    output_dir = ".") {
    cols_oi <- unique(c("Patient", "Sample", condition_var, "source_target", "complex_interaction"))

    if (!file.exists(input_file)) {
        stop(glue::glue("{input_file} is not a valid path"))
    }
    message("Loading input file...")
    input_file <- readRDS(input_file)

    if (sum(stringr::str_detect(colnames(input_file), "Patient")) == 0) {
        message("Patient column missing, assuming Patient = Sample...")
        input_file <- input_file %>%
            dplyr::mutate(Patient = Sample)
    }

    if (condition_var == "Condition_dummy") {
        input_file <- input_file %>% dplyr::mutate(Condition_dummy = "group1")
    }

    # ---- Lenient voting ---- #
    # Check detection of the same interaction in the same patient across samples for the same region
    message("Voting mode: 'lenient' (interaction needs to be detected in LIANA + 2 other tools)...")
    lenient_voting <- input_file %>%
        # Only look at interactions that pass the lenient-filter
        filter(lenient_voting) %>%
        dplyr::select(dplyr::all_of(cols_oi)) %>%
        # Count in how many samples detected an interaction per patient
        # lenient_N_samples_same_patient will be 1 if there is only a single sample per patient
        dplyr::group_by_at(
            dplyr::vars((c(condition_var, "Patient", "source_target", "complex_interaction")))
        ) %>%
        dplyr::reframe(
            lenient_detected_same_patient = paste0(Sample, collapse = ", "),
            lenient_N_samples_same_patient = dplyr::n()
        ) %>%
        dplyr::ungroup()

    lenient_voting_by_condition <- lenient_voting %>%
        dplyr::group_by_at(
            dplyr::vars((c(condition_var, "source_target", "complex_interaction")))
        ) %>%
        # Count in how many patients an interaction was detected per condition or group
        # NOTE: if sample == patient (i.e. 1 sample per patient), this will count how many samples detected an interaction per condition or group
        dplyr::reframe(
            lenient_condition_n_patients = dplyr::n(),
            lenient_condition_patients = paste0(Patient, collapse = ", "),
            lenient_condition_n_samples = sum(lenient_N_samples_same_patient),
            lenient_condition_samples = paste0(lenient_detected_same_patient, collapse = ", ")
        ) %>%
        dplyr::ungroup()
    # NOTE if sample == patient: lenient_N_samples_same_patient == lenient_condition_n_patients


    n_before <- nrow(lenient_voting_by_condition)

    message(glue::glue("Only keep interactions that are found in at least {min_patients} patients (if patient = sample, then samples)..."))
    lenient_voting_by_condition <- lenient_voting_by_condition %>% filter(lenient_condition_n_patients >= min_patients)
    n_after <- nrow(lenient_voting_by_condition)

    message(glue::glue("Before filtering: {n_before}"))
    message(glue::glue("After filtering: {n_after}"))

    # ---- Stringent voting ---- #
    message("Voting mode: 'stringent' (interaction needs to be detected in all tools)...")

    stringent_voting <- input_file %>%
        # Only look at interactions that pass the stringent-filter
        filter(stringent_voting) %>%
        dplyr::select(dplyr::all_of(cols_oi)) %>%
        # Count in how many samples detected an interaction per patient
        # lenient_N_samples_same_patient will be 1 if there is only a single sample per patient
        dplyr::group_by_at(dplyr::vars(
            (c(condition_var, "Patient", "source_target", "complex_interaction"))
        )) %>%
        dplyr::reframe(
            stringent_N_samples_same_patient = dplyr::n(),
            stringent_detected_same_patient = paste0(Sample, collapse = ", ")
        ) %>%
        dplyr::ungroup()

    stringent_voting_by_condition <- stringent_voting %>%
        # Count in how many patients an interaction was detected per condition or group
        # NOTE: if sample == patient (i.e. 1 sample per patient), this will count how many samples detected an interaction per condition or group
        dplyr::group_by_at(dplyr::vars((c(condition_var, "source_target", "complex_interaction")))) %>%
        dplyr::reframe(
            stringent_condition_n_samples = sum(stringent_N_samples_same_patient),
            stringent_condition_samples = paste0(stringent_detected_same_patient, collapse = ", "),
            stringent_condition_n_patients = dplyr::n(),
            stringent_condition_patients = paste0(Patient, collapse = ", ")
        ) %>%
        dplyr::ungroup()
    # NOTE if sample == patient: stringent_N_samples_same_patient == stringent_condition_n_patients

    n_before <- nrow(stringent_voting_by_condition)

    message(glue::glue("Only keep interactions that are found in at least {min_patients} patients..."))
    stringent_voting_by_condition <- stringent_voting_by_condition %>%
        filter(stringent_condition_n_patients >= min_patients)

    n_after <- nrow(stringent_voting_by_condition)

    message(glue::glue("Before filtering: {n_before}"))
    message(glue::glue("After filtering: {n_after}"))

    # ---- Combine stringent and lenient voting ---- #
    message("Combine stringent and lenient voting...")
    stringent_voting_by_condition <- stringent_voting_by_condition %>%
        dplyr::mutate(stringent_condition = TRUE)

    # lenient-filtered + stringent (stringent is always a subset of lenient)
    combined_voting <- lenient_voting_by_condition %>%
        dplyr::left_join(
            stringent_voting_by_condition,
        ) %>%
        # If not found by stringent, then automatically 'lenient', set in that case to 'FALSE'
        dplyr::mutate(stringent_condition = ifelse(is.na(stringent_condition), FALSE, stringent_condition)) %>%
        dplyr::ungroup()

    # Handling sample == patient (i.e. only 1 sample per patient) -> remove redundant columns
    if (sum(combined_voting$lenient_condition_n_patients != combined_voting$lenient_condition_n_samples) == 0 & sum(combined_voting$lenient_condition_samples != combined_voting$lenient_condition_samples) == 0) {
        combined_voting <- combined_voting %>% dplyr::select(
            -lenient_condition_n_samples, -lenient_condition_samples,
            -stringent_condition_n_samples, -stringent_condition_samples
        )
    }

    message("Save results...")
    saveRDS(combined_voting, glue::glue("{output_dir}/402a_filtering_detect_in_multi_samples.rds"))
    message("Finished!")
}
