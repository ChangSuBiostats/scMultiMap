#' scMultiMap
#'
#' Run scMultiMap to infer cell-type-specific peak-gene associations.
#'
#' @param obj A Seurat object of single-cell multimodal data.
#' @param pairs_df A data frame of peak-gene pairs to evaluate. The first column is gene name,
#' and the second column is peak name.
#' @param bsample Default to NULL when all cells are from the same biological sample.
#' When multiple biological samples are present, \emph{bsample} should be either the column name for
#' the sample label in Seurat object meta data, or a length n vector of labels for biological samples.
#' @param gene_assay Name of the gene assay in Seurat object. Default to 'RNA'.
#' @param peak_assay Name of the peak assay in Seurat object. Default to 'peak'.
#' @param seq_depth_list A list of two length n vectors of sequencing depths.
#' If NULL, they will be calculated as the total number of gene/peak counts in each cell
#' for two modalities, respectively.
#' @param p.adjust.method Method for multiple testing adjustment, passed to stats::p.adjust.
#' Default to 'BH' for FDR control.
#' @param irls Whether to run IRLS to estimate mean and variance parameters in step 1.
#' Default to TRUE.
#' @param verbose Whether to print detailed messages. Default to FALSE.
#'
#' @return
#' A data frame of scMultiMap results with six columns.
#' \describe{
#'   \item{gene}{gene in the peak-gene pair}
#'   \item{peak}{peak in the peak-gene pair}
#'   \item{pval}{p-value}
#'   \item{test_stat}{test statistics}
#'   \item{covar}{estimated covariance}
#'   \item{cor}{estimated correlation}
#' }
#' Each row corresponds to one peak-gene pair's results from scMultiMap.
#' \emph{pval} and \emph{test_stat} denote the statistical significance of peak-gene association,
#' and \emph{covar} denote the magnitude of association.
#' The number of rows is equal to the number of pairs in \emph{pairs_df} present in the Seurat object.
#' @export
#'
#' @examples
#' # For illustrative purposes, we provide
#' # a toy Seurat object `small_obj` and
#' # a data frame of candidate peak-gene pairs `small_pairs_df`
#' # in this R package.
#' res <- scMultiMap(small_obj, small_pairs_df)
#' head(res,2)
#' #    gene                    peak      pval test_stat        covar
#' # 1 MALAT1 chr11-65013281-65014841 0.1521555 1.4319591 3.100282e-05
#' # 2 MALAT1 chr11-65026934-65027962 0.7938414 0.2613256 5.204042e-06
#'
#' @source
#' Cell-type-specific mapping of enhancer and target genes from single-cell multimodal data.
#' Chang Su, Dongsoo Lee, Peng Jin, Jingfei Zhang;
#' https://www.biorxiv.org/content/10.1101/2024.09.24.614814v1
scMultiMap <- function(obj, pairs_df,
                       bsample = NULL,
                       gene_assay = 'RNA', peak_assay = 'peak',
                       seq_depth_list = NULL,
                       p.adjust.method = 'BH',
                       irls = T, verbose = F){
  start_time <- Sys.time()
  ##
  # prepare input data for scMultiMap: count and sequencing depth
  ##
  # extract the UMI count matrix from the single cell object
  count_list <- list()
  assay_names <- c(gene_assay, peak_assay)
  if(is.null(seq_depth_list)){
    seq_depth_list <- list()
  }else{
    names(seq_depth_list) <- assay_names
  }
  for(i in 1:2){
    mod <- assay_names[i]
    # extract 'full' counts for all features
    fcounts <- obj[[mod]]$counts
    # calculate sequencing depths if not given
    if(is.null(seq_depth_list[[mod]])){
      seq_depth_list[[mod]] <- Matrix::colSums(fcounts)
    }else{
      if(length(seq_depth_list[[mod]]) != ncol(fcounts)){
        stop(sprintf("The length of the %s sequencing depth does not match the number of cells.", mod))
      }
    }
    # check if all genes/peaks are in the Seurat object
    unique_features <- unique(pairs_df[[i]])
    int_features <- rownames(fcounts)[rownames(fcounts) %in% unique_features]
    if(length(int_features) < length(unique_features)){
      warning(sprintf('%i %s counts were not found in the Seurat object',
                      length(unique_features) - length(int_features), colnames(pairs_df)[i]))
      # remove pairs that do not overlap with the Seurat object
      pairs_df <- pairs_df[pairs_df[[i]] %in% int_features,]
    }
    # Convert to a dense matrix since most operations below do not benefit from sparsity 
    # (e.g., `X_centered = X - M` in scMultiMap_IRLS.R). 
    # However, this approach may require substantial memory, especially for large single-cell datasets 
    # with millions of cells. 
    # TODO: we leave optimizing memory efficiency an area for future development.
    count_list[[mod]] <- as.matrix(Matrix::t(fcounts[int_features,,drop=F]))
  }

  ##
  # scMultiMap step 1: IRLS to estimate mean and variance
  ##
  irls_res <- list()
  # preprocess the vector for biological samples, if specified
  if(!is.null(bsample) & length(bsample) != ncol(obj)){
    if(!bsample %in% colnames(obj@meta.data)){
      stop(sprintf('bsample label %s was not found in the metadata of Seurat object', bsample))
    }else{
      bsample <- obj[[bsample]][[1]]
    }
  }
  # run IRLS for each modality, respectively
  message('Start step 1: IRLS')
  for(i in 1:2){
    mod <- assay_names[i]
    message(sprintf('Start IRLS for %s', mod))
    irls_res[[mod]] <- scMultiMap_IRLS(count_list[[mod]],
                                       seq_depth_list[[mod]],
                                       bsample=bsample, irls=irls,
                                       verbose=verbose)
  }

  ##
  # scMultiMap step 2: WLS to infer peak-gene association
  ##
  message('Start step 2: WLS')
  wls_res <- scMultiMap_WLS(count_list,
                            seq_depth_list,
                            pairs_df,
                            irls_res,
                            verbose = verbose)
  # match with the ordering in pairs_df
  pairs <- paste(pairs_df$gene, pairs_df$peak, sep = '.')
  wls_pairs <- paste(wls_res$gene, wls_res$peak, sep = '.')
  wls_res <- wls_res[match(pairs, wls_pairs),]
  # compute adjusted p values
  wls_res$padj <- stats::p.adjust(wls_res$pval, method = p.adjust.method)
  # append correlation estimates
  gene_var <- irls_res[[1]]$sigma_sq[match(pairs_df$gene, names(irls_res[[1]]$sigma_sq))]
  peak_var <- irls_res[[2]]$sigma_sq[match(pairs_df$peak, names(irls_res[[2]]$sigma_sq))]
  # due to mom estimates, some variance estimates may be zero
  # these genes/peaks tend to have lower abundance
  # we set those involved correlations to 0
  gene_var[gene_var <= 0] <- NA
  peak_var[peak_var <= 0] <- NA
  wls_res$cor <- wls_res$covar / sqrt(gene_var * peak_var)
  wls_res$cor[is.na(wls_res$cor)] <- 0
  # postprocess out-of-bound correlation estimates
  wls_res$cor[wls_res$cor > 1] <- 1
  wls_res$cor[wls_res$cor < -1] <- -1
  rownames(wls_res) <- NULL

  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "secs")
  message(sprintf("scMultiMap elapsed time: %02d:%02d (mm:ss)\n",
              (as.integer(elapsed) %% 3600) %/% 60, # Minutes
              as.integer(elapsed) %% 60))        # Seconds
  return(wls_res)
}
