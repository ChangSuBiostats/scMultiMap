## code to prepare `example_pair` dataset goes here

# Generate peak and gene seq depths from two log normal distributions
n <- 2000
set.seed(2024)
gene_seq_depths <- exp(rnorm(n, 7.8, sd=0.6))
gene_seq_depths[gene_seq_depths < 400] <- 400
gene_seq_depths <- round(gene_seq_depths)
summary(gene_seq_depths)

peak_seq_depths <- exp(rnorm(n, 9.5, sd=0.3))
peak_seq_depths[peak_seq_depths < 2000] <- 2000
peak_seq_depths <- round(peak_seq_depths)
summary(peak_seq_depths)

seq_depths_mat <- matrix(c(gene_seq_depths, peak_seq_depths), nrow = n, ncol = 2)

# Generate a pair of independent peak and gene
z_mat <- matrix(c(rgamma(n, shape=4.83, scale=7.56e-05),
                  rgamma(n, shape=4.83, scale=7.56e-05/4)), # smaller mean and var for hte peak
                nrow = n, ncol = 2)
x_mat <- matrix(rpois(n = 2*n, lambda=z_mat * seq_depths_mat), nrow = n, ncol = 2)
summary(x_mat[,1])
summary(x_mat[,2])

example_pair <- list(gene = list(counts = x_mat[,1], seq_depths = seq_depths_mat[,1]),
                     peak = list(counts = x_mat[,2], seq_depths = seq_depths_mat[,2]))

usethis::use_data(example_pair, overwrite = TRUE)
