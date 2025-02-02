test_that("IRLS works as expected when bsample=1", {
  irls_example <- scMultiMap_IRLS(matrix(example_pair$gene$counts, ncol=1),
                                  example_pair$gene$seq_depths)
  expect_equal(irls_example$mu, 0.0003696779, tolerance=1e-7)
  expect_equal(irls_example$sigma_sq, 3.941435e-08, tolerance=1e-7)

  irls_example <- scMultiMap_IRLS(matrix(example_pair$peak$counts, ncol=1),
                                  example_pair$peak$seq_depths)
  expect_equal(irls_example$mu, 9.382313e-05, tolerance=1e-7)
  expect_equal(irls_example$sigma_sq, 1.909467e-09, tolerance=1e-7)
})


test_that("IRLS works as expected when bsample>2", {
  irls_example <- scMultiMap_IRLS(example_pair_bsample$gene$counts,
                                  example_pair_bsample$gene$seq_depths,
                                  example_pair_bsample$bsample)
  expected_mu <- c(0.0003515204, 0.0007114155)
  names(expected_mu) <- c('bsample1', 'bsample2')
  expect_equal(irls_example$mu, expected_mu, tolerance=1e-7)
  expect_equal(irls_example$sigma_sq, 4.225381e-08, tolerance=1e-7)

  irls_example <- scMultiMap_IRLS(example_pair_bsample$peak$counts,
                                  example_pair_bsample$peak$seq_depths,
                                  example_pair_bsample$bsample)
  expected_mu <- c(8.851425e-05, 0.0001859865)
  names(expected_mu) <- c('bsample1', 'bsample2')
  expect_equal(irls_example$mu, expected_mu, tolerance=1e-7)
  expect_equal(irls_example$sigma_sq, 2.146999e-09, tolerance=1e-7)
})
