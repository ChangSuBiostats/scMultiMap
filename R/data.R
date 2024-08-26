#' A simulated independent peak-gene pair
#'
#' This data set is created only for illustrative purposes and is used to test
#' scMultiMap_IRLS() and scMultiMap_WLS().
#' The source code for generating this data set is in data-raw/ on Github.
#'
#' @format ## `example_pair`
#' Two lists (one named 'gene' and the other named 'peak'), each with 2 elements:
#' \describe{
#'   \item{counts}{A length 2000 vector of gene/peak counts}
#'   \item{seq_depths}{A length 2000 vector of sequencing depths of scRNA-seq/scATAC-seq}
#' }
"example_pair"


#' A simulated perfectly correlated peak-gene pair with multiple biological samples present
#'
#' This data set is created only for illustrative purposes and is used to test
#' scMultiMap_IRLS() and scMultiMap_WLS().
#' The source code for generating this data set is in data-raw/ on Github.
#'
#' @format ## `example_pair_bsample`
#' Three lists:
#' The first two list are named 'gene' and 'peak' respectively, each with 2 elements:
#' \describe{
#'   \item{counts}{A length 2000 vector of gene/peak counts}
#'   \item{seq_depths}{A length 2000 vector of sequencing depths of scRNA-seq/scATAC-seq}
#' },
#' the third list is named 'bsample' and hosts a length n vector of sample labels
"example_pair_bsample"


#' A small pedagogical Seurat object
#'
#' This Seurat object is created only for illustrative purposes and is used to test
#' scMultiMap().
#' The source code for generating this data set is in data-raw/ on Github.
#'
#' @format ## `small_obj`
#' A Seurat object with 100 cells and two assay
#' \describe{
#'   \item{RNA}{A gene count matrix with 100 genes and 100 cells}
#'   \item{peak}{A peak count matrix with 250 genes and 100 cells}
#' }
"small_obj"

#' A data frame of peak-gene pairs, to be used with `small_obj`
#'
#' This data frame holds the candidate peak-gene pairs to be inferred with
#' scMultiMap() and `small_obj`.
#' The source code for generating this data set is in data-raw/ on Github.
#'
#' @format ## `small_pairs_df`
#' A data frame with 250 rows (peak-gene pairs) and 2 columns
#' \describe{
#'   \item{gene}{gene name}
#'   \item{peak}{peak name}
#' }
"small_pairs_df"


#' A toy example of estimated results from scMultiMap_IRLS() on `small_obj`.
#'
#'
#' The source code for generating this data set is in data-raw/ on Github.
#'
#' @format ## `small_irls_list`
#' A list of two lists generated from the 1st step of scMultiMap
#' \describe{
#'   \item{gene}{\emph{scMultiMap_IRLS} results on genes}
#'   \item{peak}{\emph{scMultiMap_IRLS} results on peaks}
#' }
"small_irls_list"
