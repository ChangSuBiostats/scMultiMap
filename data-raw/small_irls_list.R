## code to prepare `small_irls_list` dataset goes here

# `small_obj` is also a dataset in this package.
small_irls_list <- list()
small_irls_list[['gene']] <- scMultiMap_IRLS(Matrix::t(small_obj[['RNA']]$counts), Matrix::colSums(small_obj[['RNA']]$counts))
small_irls_list[['peak']] <- scMultiMap_IRLS(Matrix::t(small_obj[['peak']]$counts), Matrix::colSums(small_obj[['peak']]$counts))
usethis::use_data(small_irls_list, overwrite = TRUE)
