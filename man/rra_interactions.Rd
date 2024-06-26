% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rra.R
\name{rra_interactions}
\alias{rra_interactions}
\title{RRA interactions}
\usage{
rra_interactions(
  cellchat_obj,
  liana_obj,
  cell2cell_obj,
  cpdb_obj,
  sample_id = NULL,
  output_dir = "."
)
}
\arguments{
\item{cellchat_obj}{path to file generated by 'format_cellchat()'}

\item{liana_obj}{path to file generated by 'format_liana()'}

\item{cell2cell_obj}{path to file generated by 'format_cell2cell()'}

\item{cpdb_obj}{path to file generated by 'format_cpdb()'}

\item{sample_id}{sample id to use for saving object}

\item{output_dir}{output directory for saving output (default = '.')}
}
\description{
Rank interactions based on CellChat, LIANA (use all methods), cell2cell and CellPhoneDB
}
