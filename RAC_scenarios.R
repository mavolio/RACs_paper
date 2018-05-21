library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(codyn)
# targets
# richness: 5, 20, 50
# evenness (inv Simpsons as E1/D from Smith & Wilson 1996): 0.2, 0.5, 0.8
# beta diversity (gamma / average alpha): 1.5, 3 to 4 ## FIXME: should this be normalized by number of sites?
# codyn turnover: 0.2, 0.7

# first for single site and single iteration (to set alpha, gamma, theta)
param <- expand.grid(
  rep = 1:50,
  alpha = c(5, 20, 50),
  gamma = 120,
  theta = c(0.7, 1.4, 2.5)
)
# run simulations
sims <- mapply(rcommunity, n = 1, size = 1000,
               gamma = param$gamma,
               alpha = param$alpha,
               theta = param$theta,
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
  group_by(alpha, theta) %>%
  summarise(avg_richness = mean(richness), avg_evenness = mean(evenness))
# visualize result
param$theta <- factor(param$theta)
ggplot(param, aes(x = avg_richness, y = avg_evenness, color = theta)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(limits = c(0, 1), breaks = c(0.2, 0.5, 0.8)) +
  scale_x_continuous(limits = c(0, 55), breaks = c(5, 20, 50))

# now add on sites and iterations, setting beta, sigma
grid_param <- expand.grid(
  rep = 1:5,
  alpha = c(5, 20, 50),
  theta = c(0.7, 1.4, 2.5)
)
param <- rbind(
  mutate(grid_param, scenario = 'a', # high beta diversity, high turnover
         gamma = 10 * alpha,
         beta = 1,
         sigma = 0.3),
  mutate(grid_param, scenario = 'b', # low beta diversity, low turnover
         gamma = 3 * alpha,
         beta = 0.1,
         sigma = 0.03),
  mutate(grid_param, scenario = 'c', # high beta diversity, low turnover
         gamma = 10 * alpha,
         beta = 1,
         sigma = 0.02),
  mutate(grid_param, scenario = 'd', # low beta diversity, high turnover
         gamma = 3 * alpha,
         beta = 0.1,
         sigma = 0.7)
)
# run simulations
sims <- mapply(rcommunity, n = 1, size = 1000, sites = 3, iterations = 3,
               gamma = param$gamma,
               alpha = param$alpha,
               theta = param$theta,
               beta = param$beta,
               sigma = param$sigma,
               SIMPLIFY = FALSE
)
# calculate turnover
param$avg_turnover <- sapply(sims, function(sim) {
  sim_turnover <- sim %>%
    select(-sample) %>%
    group_by(site) %>%
    mutate(time = as.numeric(iteration)) %>%
    do({
      turnover(., time.var = 'time', species.var = 'species',
               abundance.var = 'abundance')
    }) %>%
    summarize(avg_turnover = mean(total))
  return(mean(sim_turnover$avg_turnover))
})
# calculate beta diversity
param$avg_pooled_S <- sapply(sims, function(sim) {
  sim_pooled_S <- sim %>%
    select(-sample) %>%
    group_by(iteration) %>%
    summarize(pooled_S = n_distinct(species))
  return(mean(sim_pooled_S$pooled_S))
})
param <- param %>%
  mutate(avg_beta_diversity = avg_pooled_S / alpha)
# visualize result
ggplot(param, aes(x=avg_turnover, y=avg_beta_diversity, color = scenario)) +
  facet_grid(alpha ~ theta, labeller = label_both) +
  geom_point() +
  scale_x_continuous(limits = c(0, 1), breaks = c(0.2, 0.7), minor_breaks = NULL) +
  scale_y_continuous(limits = c(0, 4), breaks = c(1.5, 3.5), minor_breaks = NULL) +
  xlab('codyn turnover') +
  ylab('pooled S / alpha')

####TO USE
#combining params and sim dataframes
#make sites and iterations 10
#
grid_param <- expand.grid(
  rep = 1:10,
  alpha = c(5, 20, 50),
  theta = c(0.7, 1.4, 2.5)
)
param <- rbind(
  mutate(grid_param, scenario = 'a', # high beta diversity, high turnover
         gamma = 10 * alpha,
         beta = 1,
         sigma = 0.3),
  mutate(grid_param, scenario = 'b', # low beta diversity, low turnover
         gamma = 3 * alpha,
         beta = 0.1,
         sigma = 0.03),
  mutate(grid_param, scenario = 'c', # high beta diversity, low turnover
         gamma = 10 * alpha,
         beta = 1,
         sigma = 0.02),
  mutate(grid_param, scenario = 'd', # low beta diversity, high turnover
         gamma = 3 * alpha,
         beta = 0.1,
         sigma = 0.7)
)
param$id<-paste(param$alpha, param$theta, param$scenario, param$rep, sep="_")

sims <- mapply(rcommunity, n = 1, size = 1000, sites = 10, iterations = 10,
               gamma = param$gamma,
               alpha = param$alpha,
               theta = param$theta,
               beta = param$beta,
               sigma = param$sigma,
               shift=T,# prevents richness from varying
               SIMPLIFY = FALSE
)

df<-data.frame()
for (i in 1:nrow(param)){
  sim<-sims[[i]]
  sim$id<-param[i,"id"]##take for row i the id column
  df<-rbind(df, sim)
}

write.csv(df, "C:\\Users\\megha\\Dropbox\\SESYNC\\SESYNC_RACs\\R Files\\SimCom_Sept28.csv")
write.csv(df, "C:\\Users\\megha\\Dropbox\\SESYNC\\SESYNC_RACs\\R Files\\SimCom_May15_ShiftT.csv")
