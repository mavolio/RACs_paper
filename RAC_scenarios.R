library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(codyn)
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

