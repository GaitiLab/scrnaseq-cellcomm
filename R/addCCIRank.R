#' @title Add CCI rank
#' @description Add cci rank based on p-value x interaction score. Needed for CCI methods that do not provide a rank, i.e. CellChat, CellPhoneDB and cell2cell.
#' @param df dataframe with at least the columns pval and interaction_score
#' @param n_perm no. permutations to define smallest possible pvalue
#' @return df with columns pval_corr, log10_score and rank
#' @export
AddCCIRank <- function(df, n_perm = 1e3) {
    return(
        df |>
            dplyr::mutate(
                # with n_perm permutation, smallest p =  1/n_perm
                pval_corr = ifelse(R.utils::isZero(pval), 1 / n_perm, pval),
                # Higher interaction_score is better, smaller pval is 'better' -> take log
                log10_score = -log10(pval_corr) * interaction_score
            ) |>
            dplyr::mutate(
                # Lower/smaller value = better (more important)
                rank = dplyr::dense_rank(dplyr::desc(log10_score))
            )
    )
}
