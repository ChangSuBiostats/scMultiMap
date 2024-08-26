## code to prepare `example_pair_bsample` dataset goes here

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

# Generate a pair of perfectly correlated peak and gene
mu <- 4.83 * 7.56e-05
sigma_sq <- 4.83 * 7.56e-05^2
# Simulate dataset with two subjcts, where each sample has 1000 cells,
# and cells from different subjects have subject-specific mean for gene expression / peak accessibility
mu_vec <- c(mu, mu * 2)
sigma_sq_vec <- rep(sigma_sq, 2)

z_gene <- c(rgamma(n/2, shape = mu_vec[1]^2 / sigma_sq_vec[1], scale = sigma_sq_vec[1] / mu_vec[1]),
            rgamma(n/2, shape = mu_vec[2]^2 / sigma_sq_vec[2], scale = sigma_sq_vec[2] / mu_vec[2]))
z_peak <- z_gene / 4 # correlation=1
z_mat <- matrix(c(z_gene, z_peak),
                nrow = n, ncol = 2)
x_mat <- matrix(rpois(n = 2*n, lambda=z_mat * seq_depths_mat), nrow = n, ncol = 2)
summary(x_mat[,1])
summary(x_mat[,2])

example_pair_bsample <- list(gene = list(counts = x_mat[,1], seq_depths = seq_depths_mat[,1]),
                             peak = list(counts = x_mat[,2], seq_depths = seq_depths_mat[,2]),
                             bsample = rep(paste0('bsample', 1:2), each = 1000))

usethis::use_data(example_pair_bsample, overwrite = TRUE)
