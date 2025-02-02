#' Iteratively reweighted least squares (IRLS) in scMultiMap
#'
#' This function implements the 1st step of scMultiMap:
#' estimate mean and variance for peaks or genes with IRLS.
#' The output will be used as input to the 2nd step, scMultiMap_WLS(), to
#' estimate and infer the association between genes and peaks.
#'
#' @param X A n by p count matrix, where n denotes the number of cells and p denotes the number of genes/peaks.
#'          It can be a length p vector when n=1.
#' @param seq_depth A length n vector of sequencing depths
#' @param bsample A length n vector of the names of biological sample the cell is from. When all cells are from the same
#'                sample, set to NULL (Default).
#' @param irls Whether to use IRLS. Default to TRUE. If FALSE, use ordinary least squares.
#' @param verbose Whether to print detailed messages. Default to FALSE.
#'
#' @return A list of two length p vectors:
#' \describe{
#'   \item{mu}{estimated mean}
#'   \item{sigma_sq}{estimated variance}
#' }
#' @export
#'
#' @examples
#' # The following codes illustrate the use of `scMultiMap_IRLS` with toy datasets
#' # `small_obj` provided in this R package.
#' irls_list <- list()
#' irls_list[['gene']] <- scMultiMap_IRLS(t(as.matrix(small_obj[['RNA']]$counts)),
#'  Matrix::colSums(small_obj[['RNA']]$counts))
#' irls_list[['peak']] <- scMultiMap_IRLS(t(as.matrix(small_obj[['peak']]$counts)),
#'  Matrix::colSums(small_obj[['peak']]$counts))
#' print(str(irls_list[['gene']])) # estimated mean and variance for genes
#' print(str(irls_list[['peak']])) # estimated mean and variance for peaks
#'
#' @source
#' Cell-type-specific mapping of enhancer and target genes from single-cell multimodal data.
#' Chang Su, Dongsoo Lee, Peng Jin, Jingfei Zhang;
#' https://www.biorxiv.org/content/10.1101/2024.09.24.614814v1
scMultiMap_IRLS <- function(X, seq_depth, bsample=NULL, irls=T, verbose=F){
  if(is.null(dim(X))) X <- matrix(X, ncol = 1)
  if (is.null(seq_depth)) {
    seq_depth = apply(X, 1, sum, na.rm = T)
  }
  if(nrow(X) != length(seq_depth)){
    print(nrow(X))
    print(length(seq_depth))
    stop('The length of the sequencing depth must match the number of cells.')
  }
  n_cell = nrow(X)
  if(is.null(bsample)){
    bsample_mtx <- matrix(1, nrow = n_cell, ncol = 1)
  }else{
    if(!is.factor(bsample)) bsample <- factor(bsample, levels = unique(bsample), labels = unique(bsample))
    bsample_mtx <- stats::model.matrix(~bsample-1)
    colnames(bsample_mtx) <- levels(bsample)
  }
  n_bsample <- ncol(bsample_mtx)
  bind_mtx <- matrix(as.logical(bsample_mtx), nrow = n_cell, ncol = n_bsample)
  colnames(bind_mtx) <- colnames(bsample_mtx)
  n_gene = ncol(X)
  seq_depth_sq = seq_depth^2
  seq_2 = sum(seq_depth_sq)
  seq_4 = sum(seq_depth^4)
  mu_mat <- matrix(NA, nrow = n_gene, ncol = n_bsample)
  M <- matrix(NA, nrow = n_cell, ncol = n_gene)
  tmp <- X * seq_depth
  mu_mat <- t(crossprod(bind_mtx, tmp) / colSums(seq_depth_sq * bind_mtx))
  M <- tcrossprod(seq_depth * bind_mtx, mu_mat)
  X_centered = X - M
  sigma2 = colSums(((X_centered^2 - M) * seq_depth_sq))/seq_4

  if(irls){
    theta = mu_mat^2/sigma2
    j = 0
    delta = Inf

    while( delta > 0.05 & j <= 10 ){
      theta_previous = theta
      # obtain common median for each biological sample
      theta_median <- numeric(n_cell)
      for(i_b in 1:n_bsample){
        theta_median[bind_mtx[,i_b]] <- stats::quantile(theta[, i_b][theta[, i_b] > 0], na.rm = T, probs = 0.5)
      }
      theta[theta < 0] = Inf
      names(theta_median) <- NULL
      w <- M + M^2 / theta_median
      w[is.na(w)|w <= 0] = 1
      mu_mat <- matrix(NA, nrow = n_gene, ncol = n_bsample)
      num_tmp <- (X/w) * seq_depth
      deno_tmp <- seq_depth_sq/w
      mu_mat <- t(crossprod(bind_mtx, num_tmp)[, , drop = FALSE]) / crossprod(deno_tmp, bind_mtx)[, , drop = FALSE]
      M <- tcrossprod(seq_depth * bind_mtx, mu_mat)
      X_centered = X - M
      h = (M^2/theta_median + M)^2
      h[h <= 0] = 1
      sigma2 = colSums(((X_centered^2 - M)/h * seq_depth_sq))/ colSums(seq_depth_sq^2/h)
      theta = mu_mat^2/sigma2
      j = j+1
      delta = max(abs(log((theta/theta_previous)[theta > 0 & theta_previous > 0])), na.rm = T)
    }
    if(j == 10 & delta > 0.05 & verbose){
      print('IRLS failed to converge after 10 iterations. Please check your data.')
    }else if(verbose){
      print(sprintf('IRLS converged after %i iterations.', j))
    }
    theta_median <- numeric(n_cell)
    for(i_b in 1:n_bsample){
      theta_median[bind_mtx[,i_b]] <- stats::quantile(theta[, i_b][theta[, i_b] > 0], na.rm = T, probs = 0.5)
    }
    theta[theta < 0] = Inf
    w = M + M^2 / theta_median
    w[is.na(w)|w <= 0] = 1
  }else{
    w <- NULL
  }
  rownames(mu_mat) <- names(sigma2) <- colnames(X)
  names(theta_median) <- rownames(X)
  if(n_bsample == 1){
    mu_mat <- mu_mat[,1,drop=TRUE]
  }else if(n_gene == 1){
    mu_mat <- mu_mat[1,,drop=TRUE]
  }
  if(n_bsample == 1){
    return(list(mu = mu_mat, sigma_sq = sigma2,
                theta = unique(theta_median)))
  }else{
    rownames(M) <- rownames(X)
    colnames(M) <- colnames(X)
    return(list(mu = mu_mat, sigma_sq = sigma2,
                theta = theta_median,
                M = M))
  }
}
