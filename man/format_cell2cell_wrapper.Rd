% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format_cell2cell.R
\name{format_cell2cell_wrapper}
\alias{format_cell2cell_wrapper}
\title{Format cell2cell CCI results}
\usage{
format_cell2cell_wrapper(
  input_interactions_pval,
  input_interactions_scores,
  ref_db = "data/interactions_db/ref_db.rds",
  output_dir = ".",
  sample_id = NULL
)
}
\arguments{
\item{input_interactions_pval}{path to file with the inferred interactions pvalues from cell2cell (output file generated by Python/202_cci_cell2cell.py, a csv file)}

\item{input_interactions_scores}{path to file with the inferred interactions pvalues from cell2cell (output file generated by Python/202_cci_cell2cell.py, a csv file)}

\item{ref_db}{Path to interactions database (default = "data/interactions_db/ref_db.rds")}

\item{output_dir}{output directory for saving output (default = '.')}

\item{sample_id}{sample id to use for saving formatted CCI results (default = NULL; will determine sample id based on 'input_interactions' for this input_interactions has to be of the format 'cell2cell__{sample_id}.csv')}
}
\description{
Format cell2cell CCI results
}
