#' WLS utility function.
#'
#' This function performs WLS for one gene and all peaks in its candidate pairs.
#'
#' @param x_gene A length n vector of gene counts. n denotes the number of cells.
#' @param X_peak A n by p peak count matrix. p denotes the number of peaks.
#' @param m_gene A length n vector of gene counts' mean across cells.
#' @param M_peak A n by p matrix of peak counts' mean across cells for all peaks.
#' @param w_gene A length n vector used to calculate regression weights.
#' @param W_peak A n by p matrix used to calculate regression weights.
#' @param var_gene A length n vector of gene counts' variance across cells.
#' @param Var_peak A n by p matrix of peak counts' variance across cells for all peaks.
#' @param gene_seq_depth A length n vector of sequencing depths from the RNA modality.
#' @param peak_seq_depth A length n vector of sequencing depths from the ATAC/peak modality.
#'
#' @return A data.frame with p rows and three columns
#' \describe{
#'   \item{pval}{p-value}
#'   \item{test_stat}{test statistics}
#'   \item{covar}{estimated covariance}
#' }
#' Each row correspond to scMultiMap's WLS result on one gene-peak pair.
#' @export
#'
#' @source
#' Cell-type-specific mapping of enhancer and target genes from single-cell multimodal data.
#' Chang Su, Dongsoo Lee, Peng Jin, Jingfei Zhang;
scMultiMap_WLS_by_gene <- function(x_gene, X_peak,
                                   m_gene, M_peak,
                                   w_gene, W_peak,
                                   var_gene, Var_peak,
                                   gene_seq_depth, peak_seq_depth){
  # obtain centered counts
  x_gene <- x_gene - m_gene
  X_peak <- X_peak - M_peak
  # obtain predictor, i.e. products of seq depths
  s_gp <- gene_seq_depth * peak_seq_depth
  # calculate regression weights
  g <- w_gene * W_peak
  # calculate covariance estimates
  nume <- colSums(x_gene * X_peak / g * s_gp)
  est <- nume / colSums(s_gp^2 / g)
  # calculate test statistics
  deno <- colSums(var_gene * Var_peak / g^2 * s_gp^2)
  ts <- nume / sqrt(deno)
  # calculate p values
  pval <- stats::pnorm(abs(ts), lower.tail = F) * 2
  # format output
  df <- data.frame(pval = pval, test_stat = ts, covar = est)
  return(df)
}
