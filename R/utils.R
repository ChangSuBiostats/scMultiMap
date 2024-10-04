#' get_top_peak_gene_pairs
#'
#' For peaks and genes with high abundance, construct candidate pairs where
#' peaks are in close proximity to the gene.
#'
#' @param obj A Seurat object of single-cell multimodal data.
#' @param gene_top Number of top highly expressed genes. Default to 2000.
#' @param peak_top Number of top highly accessible peaks. Default to 20000.
#' @param distance Distance threshold for peaks to be considered as possibly associated with genes. Default to 5e+05.
#' @param gene_assay Name of the gene assay in Seurat object. Default to 'RNA'.
#' @param peak_assay Name of the peak assay in Seurat object. Default to 'peak'.
#'
#' @return
#' A data frame of candidate peak-gene pairs.
#' \describe{
#'   \item{gene}{gene in the peak-gene pair}
#'   \item{peak}{peak in the peak-gene pair}
#' }
#' @export
#'
get_top_peak_gene_pairs <- function(obj, gene_top=2000, peak_top=20000,
                                    distance = 5e+05,
                                    gene_assay = 'RNA', peak_assay = 'peak'){
  # focus on top genes
  top_genes <- order(Matrix::rowSums(obj[[gene_assay]]$counts), decreasing = T)[1:gene_top]
  # and top peaks
  top_peaks <- order(Matrix::rowSums(obj[[peak_assay]]$counts), decreasing = T)[1:peak_top]
  # obtain gene locations
  # https://github.com/stuart-lab/signac/blob/HEAD/R/links.R#L281C1-L287C6
  annot <- Signac::Annotation(object = obj[[peak_assay]])
  gene.coords <- Signac:::CollapseToLongestTranscript(
    ranges = annot
  )
  peaks <- Signac::granges(x = obj[[peak_assay]])
  gene_names <- rownames(obj[[gene_assay]]$counts)
  gene_names_top <- gene_names[1:gene_top]
  # genes with high expression levels and known locations from genome annotation
  int_gene_names <- intersect(gene_names_top, gene.coords$gene_name)
  # construct a peak-gene pair sparse matrix with where peaks and genes are within a certain distance
  # https://github.com/stuart-lab/signac/blob/HEAD/R/links.R#L351C3-L355C4
  suppressWarnings({
  peak_distance_matrix <- Signac:::DistanceToTSS(
    peaks = peaks[top_peaks],
    genes = gene.coords[match(int_gene_names, gene.coords$gene_name)],
    distance = distance
  )
  })
  summ <- Matrix::summary(peak_distance_matrix)
  df <- data.frame(gene = colnames(peak_distance_matrix)[summ$j],
                   peak = rownames(peak_distance_matrix)[summ$i])
  return(df)
}
