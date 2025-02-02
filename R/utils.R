utils::globalVariables(c(".","seqnames", "start", "end", "strand", "gene_biotype", "gene_name"))


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
  gene.coords <- CollapseToLongestTranscript(
    ranges = annot
  )
  peaks <- Signac::granges(x = obj[[peak_assay]])
  gene_names <- rownames(obj[[gene_assay]]$counts)
  gene_names_top <- gene_names[top_genes]
  # genes with high expression levels and known locations from genome annotation
  int_gene_names <- gene_names_top[gene_names_top %in% gene.coords$gene_name]
  # construct a peak-gene pair sparse matrix with where peaks and genes are within a certain distance
  # https://github.com/stuart-lab/signac/blob/HEAD/R/links.R#L351C3-L355C4
  suppressWarnings({
  peak_distance_matrix <- DistanceToTSS(
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


#' validate_with_assay
#'
#' Validate peak-gene associations computed with single-cell multimodal data
#' by evaluating the consistency with orthogonal assay.
#'
#' @param res_df A data frame of peak-gene associations results, containing columns of `gene`, `peak`, and pvar (see below).
#' @param sig_assay_gr A GRanges object for associated peak-gene pairs from other orthogonal assays, containing `gene` as a column in metadata
#' @param pvar Which variable in `res_df` to use for claiming significant associations. Default to 'padj'.
#' @param p_cutoff Cutoff for significance based on `pvar`. Default to 0.2.
#'
#' @return
#' A named numeric vector that summarizes the overlap between `res_df` results and `sig_assay_gr` results
#' \describe{
#'   \item{pval}{p value from one-sided Fisher exact test}
#'   \item{n_overlap}{numebr of overlapped peak-gene pairs}
#'   \item{oddsratio}{odds ratio in the two-by-two table}
#'   \item{ratio}{overlapped pairs among all significant pairs from `res_df`}
#'   \item{enr}{enrichment of overlapped pairs among significant pairs compared to all background pairs}
#' }
#' @export
#'
validate_with_assay <- function(res_df, sig_assay_gr,
                                pvar = 'padj', p_cutoff = 0.2){
  grange_peaks <- Signac::StringToGRanges(res_df$peak, sep = c('-', '-'))
  sig_res_inds <- res_df[[pvar]] < p_cutoff

  # evaluate the overlap between peaks
  bovp_gr = GenomicRanges::findOverlaps(grange_peaks, sig_assay_gr)
  bovp = data.frame(
      d1_idx = S4Vectors::queryHits(bovp_gr),
      d2_idx = S4Vectors::subjectHits(bovp_gr)
  )
  sovp_gr <- GenomicRanges::findOverlaps(grange_peaks[sig_res_inds], sig_assay_gr)
  sovp = data.frame(
      d1_idx = S4Vectors::queryHits(sovp_gr),
      d2_idx = S4Vectors::subjectHits(sovp_gr)
  )
  # evaluate the overlap between peak-gene pairs
  assay_sig_inds <- res_df$gene[sig_res_inds][sovp$d1_idx] == sig_assay_gr$gene[sovp$d2_idx]
  n_assay_sig <- sum(assay_sig_inds, na.rm=T)
  n_assay <- sum(res_df$gene[bovp$d1_idx] == sig_assay_gr$gene[bovp$d2_idx], na.rm=T)
  res_gene_inds <- (res_df$gene %in% sig_assay_gr$gene)

  # compute statistics on the significance of overlap
  n_not_assay <- nrow(res_df[res_gene_inds,]) - n_assay
  n_sig <- sum(sig_res_inds & res_gene_inds)

  pval <- stats::phyper(n_assay_sig, n_assay, n_not_assay, n_sig, lower.tail=F)
  n_overlap <- n_assay_sig
  ratio <- n_assay_sig/n_sig
  enr <- (n_assay_sig/n_sig) / (n_assay/(n_assay+n_not_assay))
  A <- n_assay_sig
  B <- n_sig - A
  C <- n_assay - A
  D <- (n_assay+n_not_assay) - (A + B + C)
  oddsratio <-  (A*D)/(B*C)
  print(sprintf('p value from one-sided Fisher exact test: %.2e, #overlap: %i, odds ratio: %.2f',
                pval, n_overlap, oddsratio))
  return(c(pval = pval, n_overlap = n_overlap, oddsratio = oddsratio, ratio = ratio, enr = enr))
}

#' make_gene_to_peak_link
#'
#' make a data frame that links peak to a given gene as input to Signac::LinkPlot()
#'
#' @param df A data frame of peak-gene associations results, containing columns of `gene`, `peak`, and pvar (see below).
#' @param gene Character, gene of interest.
#' @param obj Seurat object that contains the data on the gene and the peaks.
#' @param peak_assay Name for peak assay. Default to 'peak'.
#'
#' @return
#' A grange object with genomic ranges corresponding to the gene and the peak, and with the following metadata:
#' \describe{
#'   \item{score}{correlation from the peak-gene association results}
#'   \item{gene}{gene name}
#'   \item{peak}{peak name}
#'   \item{test_stat}{test statistics from the peak-gene association results`}
#'   \item{pval}{p-value from the peak-gene association results}
#' }
#' @export
#'
make_gene_to_peak_link <- function(df, gene, obj, peak_assay='peak'){
  # subset to the gene of interest
  df <- df[df$gene == gene,]
  # extract the start of the peak
  tmp <- strsplit(df$peak, split='-')
  chr <- sapply(tmp, function(pe) pe[1])
  peak_start <- sapply(tmp, function(pe) pe[2])
  # obtain gene TSS
  # https://github.com/stuart-lab/signac/blob/HEAD/R/links.R#L281C1-L287C6
  annot <- Signac::Annotation(object = obj[[peak_assay]])
  gene.coords <- CollapseToLongestTranscript(
    ranges = annot
  )
  tmp <- as.data.frame(gene.coords[gene.coords$gene_name == gene])[1,]
  gene.TSS <- ifelse(tmp$strand == '+', tmp$start, tmp$end)
  # construct a dataframe of peak-gene links
  gene_to_peak_df <- data.frame(chr = chr,
                                start = peak_start,
                                end = gene.TSS)
  rev_inds <- gene_to_peak_df$start > gene_to_peak_df$end
  tmp <- gene_to_peak_df$start[rev_inds]
  gene_to_peak_df$start[rev_inds] <- gene_to_peak_df$end[rev_inds]
  gene_to_peak_df$end[rev_inds] <- tmp
  # convert the dataframe into a grange object
  gene_to_peak_gr <- GenomicRanges::makeGRangesFromDataFrame(gene_to_peak_df,
                                              seqnames.field='chr',
                                              start.field="start",
                                              end.field='end')
  gene_to_peak_gr$score <- df$cor
  gene_to_peak_gr$gene <- df$gene
  gene_to_peak_gr$peak <- df$peak
  gene_to_peak_gr$test_stat <- df$test_stat
  gene_to_peak_gr$pval <- df$pval
  return(gene_to_peak_gr)
}


#' CollapseToLongestTranscript
#'
#' Obtain gene coordinates, a duplicate of Signac:::CollapseToLongestTranscript from Signac 1.14.0
#'
#' @param ranges Signac annotation
#'
#' @return A grange object of gene coordinates
#'
#' @import data.table
CollapseToLongestTranscript <- function(ranges) {
    range.df <- as.data.table(x = ranges)
    range.df$strand <- as.character(x = range.df$strand)
    range.df$strand <- ifelse(test = range.df$strand == "*", 
        yes = "+", no = range.df$strand)
    collapsed <- range.df[, .(unique(seqnames), min(start), max(end), 
        strand[[1]], gene_biotype[[1]], gene_name[[1]]), "gene_id"]
    colnames(x = collapsed) <- c("gene_id", "seqnames", "start", 
        "end", "strand", "gene_biotype", "gene_name")
    collapsed$gene_name <- make.unique(names = collapsed$gene_name)
    gene.ranges <- GenomicRanges::makeGRangesFromDataFrame(df = collapsed, keep.extra.columns = TRUE)
    return(gene.ranges)
}

#' Find peaks near genes
#'
#' Find peaks that are within a given distance threshold to each gene
#' a duplicate of Signac:::DistanceToTSS from Signac 1.14.0
#'
#' @param peaks A GRanges object containing peak coordinates
#' @param genes A GRanges object containing gene coordinates
#' @param distance Distance threshold. Peaks within this distance from the gene
#' will be recorded.
#' @param sep Separator for peak names when creating results matrix
#'
#' @return Returns a sparse matrix
DistanceToTSS <- function (peaks, genes, distance = 2e+05, sep = c("-", "-")) {
    tss <- GenomicRanges::resize(x = genes, width = 1, fix = "start")
    genes.extended <- suppressWarnings(expr = Signac::Extend(x = tss, 
        upstream = distance, downstream = distance))
    overlaps <- GenomicRanges::findOverlaps(query = peaks, subject = genes.extended, 
        type = "any", select = "all")
    hit_matrix <- Matrix::sparseMatrix(i = S4Vectors::queryHits(x = overlaps), j = S4Vectors::subjectHits(x = overlaps), 
        x = 1, dims = c(length(x = peaks), length(x = genes.extended)))
    rownames(x = hit_matrix) <- Signac::GRangesToString(grange = peaks, 
        sep = sep)
    colnames(x = hit_matrix) <- genes.extended$gene_name
    return(hit_matrix)
}
