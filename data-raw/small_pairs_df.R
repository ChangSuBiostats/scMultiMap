## code to prepare `example_pairs` dataset goes here

pairs_df <- read.table(sprintf('%s/CD14 Mono_20000_peak_2000_gene_pairs.txt', data_dir),
                       header = T)
# genes and peaks randomly selected in small_obj
uni_genes <- unique(pairs_df$gene)
genes <- uni_genes[1:100]
peaks <- pairs_df$peak[pairs_df$gene %in% uni_genes[1:10]]
small_pairs_df <- pairs_df[pairs_df$gene %in% uni_genes[1:10], c('gene', 'peak')]

usethis::use_data(small_pairs_df, overwrite = TRUE)
