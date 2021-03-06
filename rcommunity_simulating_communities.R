#' Generate random ecological community abundance distributions
#'
#' Sample a family of Multinomial distributions with probabilities generated by
#' a multivariate continuous-state Markov chain. Use to simulate an ecological
#' assemblage with temporally and spatially varying abundance distributions.
#'
#' @param n Number of samples for each of `sites` by `iterations` Multinomial
#'   distributions
#' @param size Number of individuals (total abundance) in each Multinomial
#'   sample
#' @param alpha "Local" species richness (richness in each sample)
#' @param gamma "Global" species richness (maximum richness of pooled samples)
#' @param theta Overall similarity in the Multinomial sampling probabilities
#' @param sites Number of times each species is repeated in the Markov chain
#'   (potentially representing spatially distinct samples)
#' @param beta Dissimilarity between `sites` in the Multinomial sampling
#'   probabilities
#' @param iterations Number of steps in the Markov chain (potentially
#'   representing temporally distinct samples)
#' @param sigma Similarity between `iterations` in the Multinomial sampling
#'   probabilities
#' @param shift Integer added to the abundance of each species in each
#'   multinomial sample; the default shift adds one individual to enforce
#'   richness equal to `alpha`
#'
#' Valid parameter values are
#' \itemize{
#' \item{`n`, `size`, `alpha`, `gamma`, `sites`, `iterations`, and `shift`:
#'   positive integer}
#' \item{`theta`: real number in (0, Inf)}
#' \item{`sigma`, `beta`: real number in [0, 1]}
#' }
#'
#' @examples
#' # Generate 5 samples from a single Multinomial distribution
#' samp <- rcommunity(n = 5, size = 100, sites = 1, iterations = 1,
#'                    alpha = 10, gamma = 20, beta = 0, sigma = 1)
#' 
#' # Generate a multivariate time series across sites, each of three
#' # site-by-iteration samples returning abundances for two species drawn from
#' # a common three-species pool
#' demo <- rcommunity(n = 3, size = 100, sites = 3, iterations = 20,
#'                    alpha = 2, gamma = 3,
#'                    beta = 0, sigma = 1, theta = 1)
#' \dontrun{
#' library(ggplot2)
#' ggplot(demo, aes(x = iteration, y = abundance, color = species)) +
#'   geom_point() +
#'   facet_wrap(~ site, nrow = 3, labeller = label_both)
#' }
#'
rcommunity <- function(n = 1, size, alpha, gamma,
                       theta = 1, beta = 1, sigma = 1,
                       sites = 1, iterations = 1, shift = TRUE) {
  
  # VARX(1) process will have jj variables (species * sites)
  # and kk steps (iterations)
  jj <- gamma * sites
  kk <- iterations
  
  # Covariance of exogenous forcing (noise)
  
  ## (1) species within site are independent with variance chosen
  ## in conjunction with the autoregresive term (see below)
  species_Sigma <- diag(1 - (1 - sigma)^2, nrow = gamma)
  species_Sigma_L <- chol.safe(species_Sigma)
  ## (2) species between sites get compound symmetry with positive correlation
  site_Sigma <- matrix(1 - beta, nrow = sites, ncol = sites)
  diag(site_Sigma) <- 1
  site_Sigma_L <- chol.safe(site_Sigma)
  ## the Kronecker product of the factorizations combines site and species
  ## covariance
  Sigma_L <- site_Sigma_L %x% species_Sigma_L

  # autoregressive process [VARX(1)]: x_{t+1} = a %*% x_t + b %*% w_t
    
  ## the autoregressive term, affecting turnover, is a scalar
  ## chosen with site_Sigma in mind to eliminate affect of
  ## sigma parameter on the stationary solution
  a <- 1 - sigma
  b <- t(Sigma_L)
  
  ## stationary covariance (assuming -1 < a < 1)
  ## < xx* > = (1 - a^2)^{-1} bb* = (1 - (1-sigma)^2)^{-1} Sigma

  ## allocate matrix for time series
  x <- matrix(0, nrow=jj, ncol=kk)
  
  ## initializing from the stationary distribution of x is done with
  x[, 1] <- (1 - a^2)^(-1/2) * b %*% rnorm(jj)

  ## simulation
  a <- diag(a, ncol = jj, nrow = jj) 
  for (i in seq_len(kk - 1) + 1) {
    x[, i] <- a %*% x[, i-1] + b %*% rnorm(jj)
  }
  
  # select at most alpha species in x
  dim(x) <- c(gamma, sites, iterations)
  I <- apply(x, c(2, 3), function(z) which(rank(-z) <= alpha))
  Ijk <- arrayInd(1:length(I), dim(I))
  Ijk[, 1] <- c(I)
  y <- scale(x[Ijk])
  
  # apply evenness power transform and normalize to probabilities
  dim(y) <- c(alpha, sites, iterations)
  y <- exp(y - max(y))^(1/theta)
  p <- y / rep(apply(y, c(2, 3), sum), each = alpha)
  
  # sample n multinomials for each site and iteration
  if (shift) {
    s <- as.integer(shift)
    abund <- apply(p, c(2, 3), rmultinom, n = n, size = (size - s * alpha)) + s
  } else {
    abund <- apply(p, c(2, 3), rmultinom, n = n, size = size)
  }
  dim(abund) <- c(alpha, n, sites, iterations)

  # assign to species number (out of 1:gamma)
  ijkl <- which(abund > 0, arr.ind = TRUE)
  Ijkl <- ijkl
  Ijkl[ , 1] <- I[matrix(ijkl[, c(1, 3, 4)], ncol = 3)]
  
  # return result as data frame
  result <- data.frame(Ijkl)
  result <- as.data.frame(lapply(result, factor))
  colnames(result) <- c('species', 'sample', 'site', 'iteration')
  result <- result[c('site', 'iteration', 'sample', 'species')]
  result$abundance <- abund[ijkl]
  return(result)
}

# Choleski decomposition within a tryCatch that safely
# handles positive semi-definite covariance matrices
chol.safe <- function(Sigma) {
  tryCatch(chol(Sigma), error = function(e) {
    q <- chol(Sigma, TRUE)
    r <- attr(q, 'rank')
    q[-(1:r), -(1:r)] <- 0
    q[, order(attr(q, 'pivot'))]
  })
}
             
