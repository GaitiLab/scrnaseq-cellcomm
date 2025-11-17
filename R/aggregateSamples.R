#' @title Aggregate samples w/ ranked interactions
#'
#' @description Combine p-values/scores by condition
#'
#' @param df dataframe with interaction results for all samples (samples_interactions_agg_rank.rds)
#' @param method correction method
#'
#' @return dataframe
#'
#' @export
AggregateSamples <- function(df, method = "BH") {
    match.arg(
        method,
        c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
    )
    df <- df |>
        # Do not need methods' individual pvalues; only using the pvalue from RRA
        dplyr::select(-dplyr::starts_with("pval_")) |>
        duckplyr::as_duckdb_tibble() |>
        dplyr::group_by(condition, source_target, complex_interaction) |>
        dplyr::summarise(
            dplyr::across(dplyr::starts_with("interaction_score"), \(x) {
                mean(x, na.rm = TRUE)
            }),
            pval = survcomp::combine.test(pval),
        ) |>
        dplyr::ungroup()
    df$pval_adj <- p.adjust(df$pval, method)
    message("Finished!")
    return(df)
}
