#' @title Check genesymbol length
#' @param vec vector with genesymbols
#' @export
check_length_genesymbol <- function(vec) {
    mask <- (nchar(vec[!is.na(vec)]) > 1)
    sum(mask) == length(mask)
}

#' @title order vector
#' @param vec vector
#' @description order elements in vector alphabetically and return vector with same length
#' @export
order_vector <- function(vec) {
    sorted_vec <- sort(vec)
    placeholder <- rep(NA, length(vec))
    placeholder[seq(length(sorted_vec))] <- sorted_vec
    return(placeholder)
}

#' @title make complex
#' @param vec vector
#' @param collapse character string to separate items default = ':'
#' @description Merge items in vector by collapse character specified by user, default = ':'
#' @export
make_complex <- function(vec, collapse = collapse) {
    paste0(vec[!is.na(vec)], collapse = collapse)
}

#' @title make complex wrapper
#' @description (1) orders vector (alphabetically, NAs last), (2) combines into single string
#' @param vec vector
#' @param collapse character string to separate items default = ':'
#' @export
make_complex_wrapper <- function(x, collapse = ":") {
    make_complex(order_vector(x), collapse = collapse)
}

#' @title Look up uniprot based on genesymbol
#' @param gene genesymbol
#' @param gene_uniprot_table dataframe containing at least 'genesymbol' and 'uniprot' columns
#' @export
#' @importFrom dplyr %>%
lookup_uniprot <- function(gene, gene_uniprot_table) {
    gene_uniprot_table %>%
        filter(genesymbol == gene) %>%
        dplyr::select(uniprot)
}

#' @title Check mapping
#' @description Check whether all entries in vector are found in a dataframe
#' @param vec vector
#' @param lookup_table dataframe for looking up valuees
#' @param lookup_var variable in lookup_table for checking overlap with vec
#' @return returns whether all entries of the vector where found in the lookup_table (boolean)
is_mapped <- function(vec, lookup_table, lookup_var = "genesymbol") {
    vec_not_na <- vec[!is.na(vec)]
    # Number of entries found should match with the length of vector without missing values
    return(
        length(
            # Check number of entries are found in table
            intersect(
                vec_not_na,
                lookup_table %>% dplyr::pull(!!dplyr::sym(lookup_var))
            )
        ) == length(vec_not_na)
    )
}

#' @title Add proteins
#' @param vec vector
#' @param gene_uniprot_table dataframe containing at least 'genesymbol' and 'uniprot' columns
#' @return string with all protein complexes that are involved separate by an underscore, if no protein found then returns vector of NAs (same length as input vector)
add_proteins <- function(vec, gene_uniprot_table) {
    tmp <- paste0(lapply(vec[!is.na(vec)], lookup_uniprot, gene_uniprot_table = gene_uniprot_table) %>% unlist(), collapse = "_")
    if (tmp == "") {
        return(NA)
    } else {
        tmp
    }
}

#' @title get proteins from  OmniPathR based on list of genes
#' @param genes_df dataframe with a single column called 'genesymbol'
#' @export
#' @importFrom dplyr %>%
get_proteins_wrapper <- function(genes_df) {
    # CellPhoneDB needs: gene_name, uniprot, hgnc_symbol, ensemble, uniprot, protein_name
    proteins <- OmnipathR::translate_ids(
        # Query
        genes_df,
        # Match ID
        genesymbol = genesymbol,
        # Variables to retrieve from OmnipathR
        uniprot, uniprot_entry
    )
    proteins <- OmnipathR::translate_ids(
        # Query
        proteins,
        # Match ID
        uniprot = uniprot,
        ensembl
    )
    # Remove duplicates based on genesymbol
    proteins <- proteins %>%
        na.omit() %>%
        dplyr::distinct(genesymbol, .keep_all = TRUE)
    return(proteins)
}
