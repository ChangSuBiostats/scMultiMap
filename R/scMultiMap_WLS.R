#' Weighted least squares (WLS) in scMultiMap
#'
#' This function implements the 2nd step of scMultiMap:
#' estimate and infer peak-gene association with WLS.
#' It takes the output from the 1st step, \emph{scMultiMap_IRLS} as input.
#'
#' @param count_list A list of two p by n count matrices.
#'        The first element is for gene and the second for peak.
#' @param seq_depth_list A list of two length p vectors of sequencing depths.
#'        The first element is for scRNA-seq and the second for scATAC-seq.
#' @param pairs_df A data frame of candidate peak-gene pairs to evaluate.
#'        The first column is gene name, and the second is peak name.
#' @param irls_list A list of two lists generated from the 1st step of scMultiMap,
#'        i.e. \emph{scMultiMap_IRLS}. The first list contains the results on genes,
#'        the second on peaks.
#' @param verbose Whether to print detailed messages. Default to FALSE.
#'
#' @return A data.frame with five columns
#' \describe{
#'   \item{gene}{gene in the peak-gene pair}
#'   \item{peak}{peak in the peak-gene pair}
#'   \item{pval}{p-value}
#'   \item{test_stat}{test statistics}
#'   \item{covar}{estimated covariance}
#' }
#' Each row corresponds to one peak-gene pair's results from scMultiMap.
#' \emph{pval} and \emph{test_stat} denote the statistical significance of peak-gene association,
#' and \emph{covar} denote the magnitude of association.
#' The number of rows is equal to the number of pairs in \emph{pairs_df} present in the Seurat object.
#' @export
#'
#' @examples
#' count_list <- list(gene = t(as.matrix(small_obj[['RNA']]$counts)),
#'  peak = t(as.matrix(small_obj[['peak']]$counts)))
#' seq_depth_list <- list(gene = Matrix::colSums(small_obj[['RNA']]$counts),
#'  peak = Matrix::colSums(small_obj[['peak']]$counts))
#' wls_res <- scMultiMap_WLS(count_list, seq_depth_list, small_pairs_df, small_irls_list)
#'
#' @source
#' Cell-type-specific mapping of enhancer and target genes from single-cell multimodal data.
#' Chang Su, Dongsoo Lee, Peng Jin, Jingfei Zhang;
#' https://www.biorxiv.org/content/10.1101/2024.09.24.614814v1
scMultiMap_WLS <- function(count_list,
                           seq_depth_list,
                           pairs_df,
                           irls_list,
                           verbose = F){
  unique_genes <- unique(pairs_df$gene)
  message(sprintf('There are %i unique genes in the peak-gene pairs.', length(unique_genes)))
  # for each gene, simultaneously infer its association with all neighboring peaks in pairs_df
  wls_res <- list()
    
  for(g in unique_genes){
    g_peaks <- pairs_df$peak[pairs_df$gene == g]
    if(verbose) print(sprintf('%i peaks around gene %s', length(g_peaks), g))
    ##
    # prepare data input to WLS for each gene
    ##
    # extract counts for this gene and for neighboring peaks
    x_gene <- count_list[[1]][,g]
    X_peak <- count_list[[2]][,g_peaks,drop=F]
    # extract fitted mean matrix
    if('M' %in% names(irls_list[[1]])){
      m_gene <- irls_list[[1]]$M[,g]
      M_peak <- irls_list[[2]]$M[,g_peaks,drop=F]
    }else{
      m_gene <- seq_depth_list[[1]] * irls_list[[1]]$mu[g]
      M_peak <- outer(seq_depth_list[[2]], irls_list[[2]]$mu[g_peaks])
    }
    # extract weights
    w_gene <- m_gene + m_gene^2 / irls_list[[1]]$theta
    W_peak <- M_peak + M_peak^2 / irls_list[[2]]$theta
    w_gene[is.na(w_gene)|w_gene <= 0] = 1
    W_peak[is.na(W_peak)|W_peak <= 0] = 1
    # extract variance
    var_gene <- m_gene + seq_depth_list[[1]]^2 * irls_list[[1]]$sigma_sq[g]
    Var_peak <- M_peak + outer(seq_depth_list[[2]]^2, irls_list[[2]]$sigma_sq[g_peaks])
    ##
    # run WLS for this gene
    ##
    wls_df <- scMultiMap_WLS_by_gene(x_gene, X_peak,
                                     m_gene, M_peak,
                                     w_gene, W_peak,
                                     var_gene, Var_peak,
                                     seq_depth_list[[1]], seq_depth_list[[2]])

    wls_df$gene <- g
    wls_df$peak <- g_peaks
    wls_df <- wls_df[, c('gene', 'peak', 'pval', 'test_stat', 'covar')]
    wls_res[[g]] <- wls_df
  }
  return(do.call(rbind, wls_res))
}
