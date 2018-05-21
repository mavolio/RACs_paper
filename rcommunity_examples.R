library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(codyn)

## Compare the relationship between evenness and richness
## across a range of parameter settings.
# set parameter ranges
param <- expand.grid(
  rep = 1:5,
  alpha = 2^(1:6),                        
  theta = c(0.1, 0.5, 1, 2, 10),          
  sigma = c(0.1, 0.5, 0.9),
  beta = c(0, 0.1, 0.4, 0.6)
)
# run simulations
sims <- mapply(rcommunity, n = 1, size = 1000,
               gamma = max(param$alpha),
               alpha = param$alpha,
               theta = param$theta,
               sigma = param$sigma,
               beta = param$beta,
               SIMPLIFY = FALSE
)
# calculate avg richness and evenness across parameter rep
mats <- lapply(sims, function(sim) {
  sim %>%
    spread('species', 'abundance') %>%
    select(-sample, -site, -iteration)
})
param$richness <- sapply(mats, function(mat) specnumber(mat))
param$evenness <- sapply(mats, function(mat) (diversity(mat, 'invsimpson') - 1) / (specnumber(mat) - 1))
param <- param %>%
  group_by(alpha, theta, sigma, beta) %>%
  summarise(avg_richness = mean(richness), avg_evenness = mean(evenness))
# visualize result
param$theta <- factor(param$theta)
ggplot(param, aes(x = avg_richness, y = avg_evenness, color = theta)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(trans = 'log2') +
  facet_grid(beta ~ sigma, labeller = label_both)
# Comments:
# The plots confirm the model design: `sigma` and `beta` have no effect on richness or evenness at a single site and iteration.
# The parameter `alpha` sets the richness level exactly (when `shift` takes the default value of TRUE), and
# parameter `theta` corresponds to evenness (here inverse Simpson shifted and scaled to span [0, 1]) of the abundance distribution.


## Compare pairwise similarities for comparing samples across
## sites versus across iterations.
# set parameter ranges
param <- expand.grid(
  rep = 1:25,
  theta = c(0.1, 0.5, 1, 2, 10),          
  sigma = c(0.1, 0.5, 0.9),
  beta = c(0, 0.1, 0.4, 0.6)
)
# run simulations
sims <- mapply(rcommunity, n=1, size=500, sites=3, iterations=3,
               alpha=16, gamma=32,
               theta=param$theta,
               sigma=param$sigma,
               beta=param$beta,
               SIMPLIFY=FALSE
)
# calculate avg community dissimilarity across reps, separately for across-site pairs and across-iteration pairs.
sims <- lapply(sims, function(sim) {
  sim <- spread(sim, 'species', 'abundance', fill = 0, sep='_')
  dist_across_iteration <- sim %>%
    group_by(site) %>%
    do({
      x <- select(., starts_with('species'))
      bc <- as.matrix(vegdist(x))[-1,]
      data.frame(bray=diag(bc))
    }) %>%
    ungroup() %>%
    summarise(across='iteration', avg_bray = mean(bray))
  dist_across_site <- sim %>%
    group_by(iteration) %>%
    do({
      x <- select(., starts_with('species'))
      bc <- as.matrix(vegdist(x))[-1,]
      data.frame(bray=diag(bc))
    }) %>%
    ungroup() %>%
    summarise(across='site', avg_bray = mean(bray))
  rbind(dist_across_site, dist_across_iteration)
})
param <- param[rep(seq_len(nrow(param)), each = 2), ]
param <- cbind(param, do.call('rbind', sims))
# visualize result
param$theta <- factor(param$theta)
ggplot(param, aes(x=across, y=avg_bray, fill=theta)) +
  facet_grid(beta ~ sigma, labeller = label_both) +
  geom_boxplot() +
  scale_y_continuous(limits=c(0, 1)) +
  xlab('Comparisons across ...') +
  ylab('Bray-Curtis Dissimilarity')
# Observations: The main effect of increasing `sigma` is to increase
# dissimilarity across iterations, although it has small and inconsistent efects
# on dissimilarity across sites. The main effect of increasing `beta` is to 
# increase dissimilarity across sites, and it does not effect similarity across 
# iterations. The effect of increasing theta is context dependent, although 
# overall it causes all pairwise dissimilarities to approach a constant value 
# (here about 0.5). For high values of theta, variation among species in 
# sampling probability that exists between sites or iterations is washed out, as
# anticipated. These simulations do not show that `alpha` (here 16) and `gamma` 
# (here 32) control the asymptotic (for high theta) value of disimilarity.
