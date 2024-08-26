## code to prepare `small_obj` dataset goes here

data_dir <- '/Users/CSU30/Documents/projects/GRN/R package/validation_w_cluster/no_bsample'
obj <- readRDS(sprintf('%s/CD14 Mono_seurat_obj.rds', data_dir))
pairs_df <- read.table(sprintf('%s/CD14 Mono_20000_peak_2000_gene_pairs.txt', data_dir),
                       header = T)
# create an artificial small Seurat object to illustrate scMultiMap usage
# for the gene counts, select 100 genes and 100 cells
uni_genes <- unique(pairs_df$gene)
small_obj <- Seurat::CreateSeuratObject(
  counts = obj[['RNA']]$counts[uni_genes[1:100],1:100],
  assay = 'RNA')
# for the peak counts, select the same 100 cells and peaks that are in candidate pairs with 10 selected genes
small_obj[['peak']] <- SeuratObject::CreateAssayObject(
  counts = obj[['peaks_ct']]$counts[pairs_df$peak[pairs_df$gene %in% uni_genes[1:10]],1:100])

usethis::use_data(small_obj, overwrite = TRUE)
