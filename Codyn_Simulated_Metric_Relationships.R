library(tidyverse)
library(vegan)
library(gridExtra)
library(grid)
library(gtable)
library(reldist)
#library(gtools)
library(devtools)
install_github("mavolio/codyn", ref = "RACs_cleaner")
library(codyn)

# Read in Data ------------------------------------------------------------
# for the sim dataset, communites are differentiated by alpha (richness), theta (evenness) and scenario (rate of turnover and spatial heterogeniety: four scenarios: a: high turnover, high spatial heterogeniety; b: low turnover, low spatial heterogeniety; c: low turnover, high spatial heterogeniety; d: high turnover, low spatial heterogeniety"). For each richness-evennes combination (9 combinations) there are each community type. Each of these 10 community types have 10 replicates, called "sites" at a given point in time. Each community type, time, and site is then replicated 10 time.

#home
sim<-read.csv("~/Dropbox/SESYNC/SESYNC_RACs/R Files/SimCom_Sept28.csv")%>%
  mutate(time=as.numeric(iteration),
         id2=paste(id, site, sep="::"))%>%
  select(-X, -sample, -iteration)%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))
  
codyndat<-read.csv("~/Dropbox/CoDyn/R Files/11_06_2015_v7/relative cover_nceas and converge_12012015_cleaned.csv")%>%
  gather(species, abundance, sp1:sp99)%>%
  filter(site_code!="MISS")

codyndat_info<-read.csv("~/Dropbox/CoDyn/R Files/11_06_2015_v7/siteinfo_key.csv")%>%
  filter(site_project_comm!="")%>%
  select(-site_project_comm)%>%
  mutate(site_project_comm = paste(site_code, project_name, community_type, sep="_"))

#work
sim<-read.csv('C:\\Users\\megha\\Dropbox\\SESYNC\\SESYNC_RACs\\R Files\\SimCom_May15_ShiftT.csv')%>%
  mutate(time=as.numeric(iteration),
         id2=paste(id, site, sep="::"))%>%
  select(-X, -sample, -iteration)%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))

codyndat<-read.csv('C:\\Users\\megha\\Dropbox\\CoDyn\\R Files\\11_06_2015_v7\\relative cover_nceas and converge_12012015_cleaned.csv')%>%
  gather(species, abundance, sp1:sp99)%>%
  filter(site_code!="MISS")

codyndat_info<-read.csv("C:\\Users\\megha\\Dropbox\\CoDyn\\R Files\\11_06_2015_v7\\siteinfo_key.csv")%>%
  filter(site_project_comm!="")

###CLEANING CODYN DATASET
#restrict to species that are present in an experiment
splist<-codyndat%>%
  group_by(site_code, project_name, community_type, species)%>%
  summarize(present=sum(abundance))%>%
  filter(present!=0)%>%
  select(-present)

#merge back and will drop species that do not exist in a dataset
codyndat_clean<-merge(codyndat, splist, by=c("site_code","project_name","community_type","species"))%>%
  select(-site_project_comm)%>%
  mutate(site_project_comm = paste(site_code, project_name, community_type, sep = "_"))%>%
  select(-X, -sitesubplot, -site_code, -project_name, -community_type)%>%
  mutate(id=paste(site_project_comm, plot_id, sep="::"))



# Richness Evenness Metrics -----------------------------------------------

#codyn dataset
spc<-unique(codyndat_clean$site_project_comm)
codyn_div_eq<-data.frame()

for (i in 1:length(spc)){
  subset<-codyndat_clean%>%
    filter(site_project_comm==spc[i])
  
  out<-community_structure(subset, time.var = 'experiment_year', abundance.var = 'abundance', replicate.var = 'plot_id', metric = "EQ")
  out$site_project_comm<-spc[i]
  
  codyn_div_eq<-rbind(codyn_div_eq, out)
}

codyn_div_esimp<-data.frame()

for (i in 1:length(spc)){
  subset<-codyndat_clean%>%
    filter(site_project_comm==spc[i])
  
  out<-community_structure(subset, time.var = 'experiment_year', abundance.var = 'abundance', replicate.var = 'plot_id', metric = "SimpsonEvenness")
  out$site_project_comm<-spc[i]
  
  codyn_div_esimp<-rbind(codyn_div_esimp, out)
}

codyn_div_evar<-data.frame()

for (i in 1:length(spc)){
  subset<-codyndat_clean%>%
    filter(site_project_comm==spc[i])
  
  out<-community_structure(subset, time.var = 'experiment_year', abundance.var = 'abundance', replicate.var = 'plot_id', metric = "Evar")
  out$site_project_comm<-spc[i]
  
  codyn_div_evar<-rbind(codyn_div_evar, out)
}


codyn_div1<-merge(codyn_div_eq, codyn_div_esimp, by=c("site_project_comm",'plot_id','richness','experiment_year'))
codyn_div<-merge(codyn_div1, codyn_div_evar, by=c("site_project_comm",'plot_id','richness','experiment_year'))

codyndat_diversity_mean <- codyn_div%>%
  group_by(site_project_comm, experiment_year)%>%
  summarize(Sp=mean(richness),
            EQ=mean(EQ, na.rm=T),
            ESimp=mean(SimpsonEvenness),
            Evar=mean(Evar, na.rm=T))


#calculating gini coefficeint using the gini function in the reldist package
#' @x the vector of abundances of each species
Gini<-function(x){
  x1<-x[x!=0]
  1-reldist::gini(x1)
}

codyndat_div_gini <- group_by(codyndat_clean, site_project_comm, experiment_year, plot_id) %>% 
  summarize(Gini=Gini(abundance))%>%
  tbl_df()%>%
  group_by(site_project_comm, experiment_year)%>%
  summarize(EGini=mean(Gini))

codyn_div_all<-merge(codyndat_diversity_mean, codyndat_div_gini, by=c("experiment_year",'site_project_comm'))

#sim dataset
com_rep<-unique(sim$id)

sim_div_eq<-data.frame()
for (i in 1:length(com_rep)){
  
  subset<-sim%>%
    filter(id==com_rep[i])
  
  out <- community_structure(df = subset, time.var = "time", abundance.var = "abundance", replicate.var = "site", metric = "EQ")
  out$id<-com_rep[i]
  
  sim_div_eq<-rbind(sim_div_eq, out)  
}

sim_div_esimp<-data.frame()
for (i in 1:length(com_rep)){
  
  subset<-sim%>%
    filter(id==com_rep[i])
  
  out <- community_structure(df = subset, time.var = "time", abundance.var = "abundance", replicate.var = "site", metric = "SimpsonEvenness")
  out$id<-com_rep[i]
  
  sim_div_esimp<-rbind(sim_div_esimp, out)  
}

sim_div_evar<-data.frame()
for (i in 1:length(com_rep)){
  
  subset<-sim%>%
    filter(id==com_rep[i])
  
  out <- community_structure(df = subset, time.var = "time", abundance.var = "abundance", replicate.var = "site", metric = "Evar")
  out$id<-com_rep[i]
  
  sim_div_evar<-rbind(sim_div_evar, out)  
}

sim_div1<-merge(sim_div_eq, sim_div_esimp, by=c("id",'site','richness','time'))
sim_div<-merge(sim_div1, sim_div_evar, by=c("id",'site','richness','time'))

sim_diversity_mean<-sim_div%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))%>%
  group_by(id3, time, site)%>%#average over replicates
  summarize(Sp=mean(richness),
            EQ=mean(EQ, na.rm=T),
            ESimp=mean(SimpsonEvenness),
            Evar=mean(Evar, na.rm=T))%>%
  ungroup()%>%
  group_by(id3, time)%>%#average over sites
  summarize(Sp=mean(Sp),
            EQ=mean(EQ, na.rm=T),
            ESimp=mean(ESimp),
            Evar=mean(Evar, na.rm=T))

sim_div_gini <- group_by(sim, id, time, site) %>% 
  summarize(Gini=Gini(abundance))%>%
  tbl_df()%>%
  group_by(id, time)%>%
  summarize(EGini=mean(Gini))%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))%>%
  group_by(id3, time)%>%
  summarize(EGini=mean(EGini))
  

sim_div_all<-merge(sim_diversity_mean, sim_div_gini, by=c("time",'id3'))


# RAC changes -------------------------------------------------------------


##Codyn Dataset

codyndat_rac_change<-data.frame()
spc<-unique(codyndat_clean$site_project_comm)

for (i in 1:length(spc)){
  
  subset<-codyndat_clean%>%
    filter(site_project_comm==spc[i])
  
  out <- RAC_change(df = subset, time.var = "experiment_year", species.var = "species", abundance.var = "abundance", replicate.var = "plot_id")
  
  out$site_project_comm<-spc[i]
  
  codyndat_rac_change<-rbind(codyndat_rac_change, out)  
}

codyndat_rac_change_average<-codyndat_rac_change%>%
  group_by(site_project_comm, experiment_year, experiment_year2)%>%
  summarise(S=mean(richness_change),
            E=mean(evenness_change,na.rm=T),
            R=mean(rank_change),
            G=mean(gains),
            L=mean(losses))




##SIM dataset
sim_rac_change<-data.frame()

com_rep<-unique(sim$id)

for (i in 1:length(com_rep)){
  
  subset<-sim%>%
    filter(id==com_rep[i])
  
  out <- RAC_change(df = subset, time.var = "time", species.var = "species", abundance.var = "abundance", replicate.var = "site")
  
  out$id<-com_rep[i]
  
  sim_rac_change<-rbind(sim_rac_change, out)  
}

sim_rac_change_mean<-sim_rac_change%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))%>%
  group_by(id3, time, time2, site)%>%
  summarize(S=mean(richness_change),
            E=mean(evenness_change,na.rm=T),
            R=mean(rank_change),
            G=mean(gains),
            L=mean(losses))%>%
  ungroup()%>%
  group_by(id3, time, time2)%>%
  summarize(S=mean(S),
            E=mean(E,na.rm=T),
            R=mean(R),
            G=mean(G),
            L=mean(L))

# Mean Change and Dispersion ----------------------------------------------
#codyn dataset

codyn_multchange<-data.frame()
spc<-unique(codyndat_clean$site_project_comm)

for (i in 1:length(spc)){
  
  subset<-codyndat_clean%>%
    filter(site_project_comm==spc[i])
  
  out <- multivariate_change(df = subset, time.var = "experiment_year", species.var = "species", abundance.var = "abundance", replicate.var = "plot_id")
  
  out$site_project_comm<-spc[i]

  codyn_multchange<-rbind(codyn_multchange, out)  
}

#Sim dataset
sim_mult_change<-data.frame()

com_rep<-unique(sim$id)

for (i in 1:length(com_rep)){
  
  subset<-sim%>%
    filter(id==com_rep[i])
  
  out <- multivariate_change(df = subset, time.var = "time", species.var = "species", abundance.var = "abundance", replicate.var = "site")
  
  out$id<-com_rep[i]
  
  sim_mult_change<-rbind(sim_mult_change, out)  
}

sim_multchange_mean<-sim_mult_change%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))%>%
  group_by(id3, time, time2)%>%
  summarize(composition_change=mean(composition_change),
            dispersion_change=mean(dispersion_change))

# Curve change ------------------------------------------------------------

####codyn first
codyn_curvechange<-data.frame()
spc<-unique(codyndat_clean$site_project_comm)

for (i in 1:length(spc)){
  
  subset<-codyndat_clean%>%
    filter(site_project_comm==spc[i])
  
  out <- curve_change(df = subset, time.var = "experiment_year", species.var = "species", abundance.var = "abundance", replicate.var = "plot_id")
  
  out$site_project_comm<-spc[i]
  
  codyn_curvechange<-rbind(codyn_curvechange, out)  
}

codyn_curvechange_mean<-codyn_curvechange%>% 
  group_by(site_project_comm, experiment_year, experiment_year2)%>%
  summarise(curve_change=mean(curve_change))

#sim dataset
sim_curve_change<-data.frame()

com_rep<-unique(sim$id)

for (i in 1:length(com_rep)){
  
  subset<-sim%>%
    filter(id==com_rep[i])
  
  out <- curve_change(df = subset, time.var = "time", species.var = "species", abundance.var = "abundance", replicate.var = "site")
  
  out$id<-com_rep[i]
  
  sim_curve_change<-rbind(sim_curve_change, out)  
}

sim_cc_ave<-sim_curve_change%>% 
  group_by(id, time, time2)%>%
  summarise(curve_change=mean(curve_change))%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))%>%
  group_by(id3, time, time2)%>%
  summarize(curve_change=mean(curve_change))


sim_shannon <- data.frame()
for (i in 1:length(com_rep)){
  
  subset<-sim%>%
    filter(id==com_rep[i])
  
  out <- community_diversity(df = subset, time.var = "time", abundance.var = "abundance", replicate.var = "site")
  
  out$id<-com_rep[i]
  
  sim_shannon<-rbind(sim_shannon, out)  
}

sim_shannon_ave<-sim_shannon%>% 
  group_by(id, time)%>%
  summarise(H=mean(Shannon))%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))%>%
  group_by(id3, time)%>%
  summarize(H=mean(H))

corr.test<-sim_cc_ave%>%
  left_join(sim_shannon_ave)

plot(corr.test$H, corr.test$curve_change)
cor.test(corr.test$H, corr.test$curve_change)
# Merging all metrics to single datasets ----------------------------------

####MERGING TO A SINGE DATASET

codyndat_allmetrics<-codyn_multchange%>%
  left_join(codyndat_rac_change_average)%>%
  left_join(codyn_curvechange_mean)%>%
  left_join(codyn_div_all)%>%
  left_join(codyndat_info)

#sim
sim_allmetrics<-sim_multchange_mean%>%
  left_join(sim_rac_change_mean)%>%
  left_join(sim_cc_ave)%>%
  left_join(sim_div_all)%>%
  separate(id3, into=c("alpha","even","comtype"), sep="_")%>%
  mutate(comtype2 = as.factor(comtype))

write.csv(codyndat_allmetrics,'~/Dropbox/SESYNC/SESYNC_RACs/R Files/codyn_allmetrics_April2018.csv', row.names = F)
write.csv(sim_all_metrics, 'C:\\Users\\megha\\Dropbox\\SESYNC\\SESYNC_RACs\\R Files\\sim_allmetrics_May2018.csv', row.names = F)

# pair plot graphs --------------------------------------------------------

theme_set(theme_bw(12))

codyndat_allmetrics<-read.csv('~/Dropbox/SESYNC/SESYNC_RACs/R Files/codyn_allmetrics_April2018.csv')%>%
  mutate(absS = abs(S),
         absE = abs(E))

# codyndat_allmetrics_old<-read.csv('~/Dropbox/SESYNC/SESYNC_RACs/R Files/codyn_allmetrics_Jan2018.csv')%>%
#   separate(experiment_year_pair, into=c("experiment_year", "experiment_year2"), sep = "-")
# 
# t1<-codyndat_allmetrics%>%
#   select(experiment_year, curve_change, site_project_comm)
# 
# t2<-codyndat_allmetrics_old%>%
#   select(experiment_year, curve_change, site_project_comm)%>%
#   mutate(experiment_year=as.numeric(experiment_year),
#          curve_change_old=curve_change,
#          exyear=experiment_year,
#          spc=site_project_comm)%>%
#   select(-curve_change, -experiment_year, -site_project_comm)
# 
# merge<-cbind(t1, t2)%>%
#   mutate(diff=curve_change_old-curve_change)%>%
#   filter(diff!=0)
# write.csv(merge, "~/Dropbox/SESYNC/SESYNC_RACs/R Files/curve_change_inconsistentcies.csv", row.names = F)
# 
sim_allmetrics<-read.csv('C:\\Users\\megha\\Dropbox\\SESYNC\\SESYNC_RACs\\R Files\\sim_allmetrics_May2018.csv')

#graphing this
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  test <- cor.test(x,y) 
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 1),
                   symbols = c("*", " "))
  
  
  text(0.5, 0.5, txt, cex = 2)
  text(0.8, 0.5, Signif, cex=5, col="red")
}

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, ...)
}

##codyn graphs
#comparing with raw S and E changes
colnames(codyndat_allmetrics)
pairs(codyndat_allmetrics[,c(6:10,3,4,11)], col=codyndat_allmetrics$taxa, labels=c("Richness\nChange", "Evenness\nChange","Rank\nChanges","Species\nGains","Species\nLosses","Compositional\nChange","Dispersion\nChange", "Curve\nChange"), font.labels=2, cex.labels=2, upper.panel = panel.cor,oma=c(4,4,4,10))
par(xpd=T)

#comparing absolute S and E changes
pairs(codyndat_allmetrics[,c(38,39, 8:10, 3, 4, 11)], col=codyndat_allmetrics$taxa, labels=c(" Abs. Richness\nChange", "Abs. Evenness\nChange","Rank\nChanges","Species\nGains","Species\nLosses","Compositional\nChange","Dispersion\nChange","Curve\nChange"), font.labels=2, cex.labels=2, upper.panel = panel.cor,oma=c(4,4,4,10))
par(xpd=T)

#comparing absolute S and E changes
pairs(codyndat_allmetrics[,c(38,39, 8:10, 3, 4, 11)], col=codyndat_allmetrics$taxa,upper.panel = panel.cor,diag.panel = panel.hist, cex.axis = 2)
par(xpd=T)


#how do these correlate with experiment parameters.
codyndat_allmetrics3<-codyndat_allmetrics%>%
  group_by(site_project_comm)%>%
  summarize_at(vars(absS, absE, R, G, L, curve_change, composition_change, dispersion_change, spatial_extent, plot_size, num_plots), funs(mean), rm.na=T)%>%
  mutate(spatialExtent=log(spatial_extent),plotSize=log(plot_size))
                 
colnames(codyndat_allmetrics3)

pairs(codyndat_allmetrics3[,c(2:9, 12:14)], upper.panel = panel.cor)
cor(codyndat_allmetrics3[,c(2:9, 12:14)], method = "pearson")


#How are evenness metrics affected by richness?
cor.test(sim_div_all$Sp, sim_div_all$EQ)
cor.test(sim_div_all$Sp, sim_div_all$EGini)
cor.test(sim_div_all$Sp, sim_div_all$ESimp)
cor.test(sim_div_all$Sp, sim_div_all$Evar)
cor.test(codyn_div_all$Sp, codyn_div_all$EQ)
cor.test(codyn_div_all$Sp, codyn_div_all$EGini)
cor.test(codyn_div_all$Sp, codyn_div_all$ESimp)
cor.test(codyn_div_all$Sp, codyn_div_all$Evar)

pairs(codyn_div_all[,4:7])
pairs(sim_div_all[,4:7])


# Effect of rich and even on CHANGE sim data ------------------------------

#How are CHANGE metrics affected by the richness and evenness of the community?
#Richness

#doing for all communities

#rich with delta rank
with(sim_allmetrics,  cor.test(Sp, R))
with(subset(sim_allmetrics, comtype =="a"), cor.test(Sp, R))
with(subset(sim_allmetrics, comtype =="b"), cor.test(Sp, R))
with(subset(sim_allmetrics, comtype =="c"), cor.test(Sp, R))
with(subset(sim_allmetrics, comtype =="d"), cor.test(Sp, R))

##removing extra division by Sp to demonstate to why it is necessary to divide again by the size of the species pool.
sim_allmetrics2<-sim_allmetrics%>%
  mutate(MRS = R*Sp)
with(sim_allmetrics2, cor.test(Sp, MRS))
with(subset(sim_allmetrics2, comtype =="a"), cor.test(Sp, MRS))
with(subset(sim_allmetrics2, comtype =="b"), cor.test(Sp, MRS))
with(subset(sim_allmetrics2, comtype =="c"), cor.test(Sp, MRS))
with(subset(sim_allmetrics2, comtype =="d"), cor.test(Sp, MRS))

ggplot(data=sim_allmetrics2, aes(x=Sp, y=MRS, color = Evar))+
  geom_point()+
  xlab("Simulated Community Richness")+
  ylab("Rank Change")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~comtype, labeller=labeller(comtype = labels), ncol =4)

#rich with gains or lossess same thing
with(sim_allmetrics, cor.test(Sp, G))
with(subset(sim_allmetrics, comtype =="a"), cor.test(Sp, G))
with(subset(sim_allmetrics, comtype =="b"), cor.test(Sp, G))
with(subset(sim_allmetrics, comtype =="c"), cor.test(Sp, G))
with(subset(sim_allmetrics, comtype =="d"), cor.test(Sp, G))

#rich with composition
with(sim_allmetrics, cor.test(Sp, composition_change))
with(subset(sim_allmetrics, comtype =="a"), cor.test(Sp, composition_change))
with(subset(sim_allmetrics, comtype =="b"), cor.test(Sp, composition_change))
with(subset(sim_allmetrics, comtype =="c"), cor.test(Sp, composition_change))
with(subset(sim_allmetrics, comtype =="d"), cor.test(Sp, composition_change))

#rich with dispersion
with(sim_allmetrics, cor.test(Sp, dispersion_change))
with(subset(sim_allmetrics, comtype =="a"), cor.test(Sp, dispersion_change))
with(subset(sim_allmetrics, comtype =="b"), cor.test(Sp, dispersion_change))
with(subset(sim_allmetrics, comtype =="c"), cor.test(Sp, dispersion_change))
with(subset(sim_allmetrics, comtype =="d"), cor.test(Sp, dispersion_change))

#rich with curve
with(sim_allmetrics, cor.test(Sp, curve_change))
with(subset(sim_allmetrics, comtype =="a"), cor.test(Sp, curve_change))
with(subset(sim_allmetrics, comtype =="b"), cor.test(Sp, curve_change))
with(subset(sim_allmetrics, comtype =="c"), cor.test(Sp, curve_change))
with(subset(sim_allmetrics, comtype =="d"), cor.test(Sp, curve_change))

#even with delta rank
with(sim_allmetrics, cor.test(Evar, R))
with(subset(sim_allmetrics, comtype =="a"), cor.test(Evar, R))
with(subset(sim_allmetrics, comtype =="b"), cor.test(Evar, R))
with(subset(sim_allmetrics, comtype =="c"), cor.test(Evar, R))
with(subset(sim_allmetrics, comtype =="d"), cor.test(Evar, R))

#even with gains/losses
with(sim_allmetrics, cor.test(Evar, G))
with(subset(sim_allmetrics, comtype =="a"), cor.test(Evar, G))
with(subset(sim_allmetrics, comtype =="b"), cor.test(Evar, G))
with(subset(sim_allmetrics, comtype =="c"), cor.test(Evar, G))
with(subset(sim_allmetrics, comtype =="d"), cor.test(Evar, G))


#even with composition
with(sim_allmetrics, cor.test(Evar, composition_change))
with(subset(sim_allmetrics, comtype =="a"), cor.test(Evar, composition_change))
with(subset(sim_allmetrics, comtype =="b"), cor.test(Evar, composition_change))
with(subset(sim_allmetrics, comtype =="c"), cor.test(Evar, composition_change))
with(subset(sim_allmetrics, comtype =="d"), cor.test(Evar, composition_change))

#even with dispersion
with(sim_allmetrics, cor.test(Evar, dispersion_change))
with(subset(sim_allmetrics, comtype =="a"), cor.test(Evar, dispersion_change))
with(subset(sim_allmetrics, comtype =="b"), cor.test(Evar, dispersion_change))
with(subset(sim_allmetrics, comtype =="c"), cor.test(Evar, dispersion_change))
with(subset(sim_allmetrics, comtype =="d"), cor.test(Evar, dispersion_change))

#even with curve
with(sim_allmetrics, cor.test(Evar, curve_change))
with(subset(sim_allmetrics, comtype =="a"), cor.test(Evar, curve_change))
with(subset(sim_allmetrics, comtype =="b"), cor.test(Evar, curve_change))
with(subset(sim_allmetrics, comtype =="c"), cor.test(Evar, curve_change))
with(subset(sim_allmetrics, comtype =="d"), cor.test(Evar, curve_change))


###Figures of significant relationships
labels <-c(a = "High Spatial, High Temporal",
           b = "Low Spatial, Low Temporal",
           c = "High Spatail, Low Temporal",
           d = "Low Spatial, High Temporal")

ggplot(data=sim_allmetrics, aes(x=Sp, y=Evar, shape = comtype))+
  geom_point(size = 2)+
  xlab("Simulated Community Richness")+
  ylab("Simulated Community Evenness")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_shape_manual(name = "Community Type", labels = c("High Spatial, High Temporal","Low Spatial, Low Temporal","High Spatail, Low Temporal","Low Spatial, High Temporal"), values = c(15,16,17,18))

sim_tograph<-sim_allmetrics%>%
  gather(metric, value, composition_change:curve_change)%>%
  filter(metric != "S")%>%
  filter(metric != "E")%>%
  filter(metric != "L")


labels_change <-c(composition_change = "Compositon Change",
           curve_change = "Curve Change",
           dispersion_change = "Dispersion Change",
           G = "Gains/Losses",
           R = "Rank Change")

#change with richness
ggplot(data=sim_tograph, aes(x=as.factor(Sp), y=value, color = comtype))+
  geom_boxplot()+
  scale_color_manual(name = "Community Type", labels = c("High Spatial, High Temporal","Low Spatial, Low Temporal","High Spatail, Low Temporal","Low Spatial, High Temporal"), values= c("gold","orange","darkorange4","red"))+
  xlab("Simulated Community Richness")+
  ylab("Metric Value")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~metric, ncol = 3, scales = "free", labeller = labeller(metric = labels_change))+
  theme(strip.background = element_rect(fill = 0))


#change with evenness
ggplot(data=sim_tograph, aes(x=as.factor(even), y=value, color = comtype))+
  geom_boxplot()+
  scale_color_manual(name = "Community Type", labels = c("High Spatial, High Temporal","Low Spatial, Low Temporal","High Spatail, Low Temporal","Low Spatial, High Temporal"), values= c("gold","orange","darkorange4","red"))+
  xlab("Simulated Community Evenness")+
  ylab("Metric Value")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~metric, ncol = 3, scales = "free", labeller = labeller(metric = labels_change))+
  theme(strip.background = element_rect(fill = 0))+
  scale_x_discrete(labels = c("Low","Mid","High"))

sim_formeans<-sim_allmetrics%>%
  mutate(abs_dc = abs(dispersion_change))%>%
  gather(metric, value, c(composition_change:curve_change, abs_dc))

sim_means<-sim_formeans%>%
  group_by(comtype, metric)%>%
  summarize(mean = mean (value),
            se = sd(value)/sqrt(81))

##make all the plots - for richness change
rrc<-ggplot(data=sim_allmetrics, aes(x=Sp, y=R, color = Evar))+
  geom_point()+
  xlab("Simulated Community Richness")+
  ylab("Rank Change")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~comtype, labeller=labeller(comtype = labels), ncol=4)
rgl<-ggplot(data=sim_allmetrics, aes(x=Sp, y=G, color = Evar))+
  geom_point()+
  xlab("Simulated Community Richness")+
  ylab("Gains/Losses")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~comtype, labeller=labeller(comtype = labels), ncol=4)
rcc<-ggplot(data=sim_allmetrics, aes(x=Sp, y=composition_change, color = Evar))+
  geom_point()+
  xlab("Simulated Community Richness")+
  ylab("Composition Change")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~comtype, labeller=labeller(comtype = labels), ncol=4)
rcuc<-ggplot(data=sim_allmetrics, aes(x=Sp, y=curve_change, color = Evar))+
  geom_point()+
  xlab("Simulated Community Richness")+
  ylab("Curve Change")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~comtype, labeller=labeller(comtype = labels), ncol=4)
grid.arrange(rrc, rgl, rcc, rcuc, ncol=1)

##make all the plots - for evenness change
erc<-ggplot(data=sim_allmetrics, aes(x=Evar, y=R, color = as.factor(Sp)))+
  geom_point()+
  xlab("Simulated Community Evenness")+
  ylab("Rank Change")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~comtype, labeller=labeller(comtype = labels), ncol = 4)
ecc<-ggplot(data=sim_allmetrics, aes(x=Evar, y=curve_change, color = as.factor(Sp)))+
  geom_point()+
  xlab("Simulated Community Evenness")+
  ylab("Composition Change")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~comtype, labeller=labeller(comtype = labels), ncol = 4)
ecuc<-ggplot(data=sim_allmetrics, aes(x=Evar, y=curve_change, color = as.factor(Sp)))+
  geom_point()+
  xlab("Simulated Community Evenness")+
  ylab("Curve Change")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~comtype, labeller=labeller(comtype = labels), ncol = 4)
grid.arrange(erc, ecc, ecuc, ncol=1)




# Effet of rich and even on DIFF sim data ---------------------------------

sim_diff<-sim%>%
  mutate(treatment = as.factor(ifelse(site <5, "T1", "T2")))

#RAC diff
sim_rac_diff<-data.frame()
com_rep<-unique(sim_diff$id)

for (i in 1:length(com_rep)){
  
  subset<-sim_diff%>%
    filter(id==com_rep[i])
  
  out <- RAC_difference(df = subset, time.var = "time", species.var = "species", abundance.var = "abundance", replicate.var = "site", pool = TRUE, treatment.var = "treatment")
  
  out$id<-com_rep[i]
  
  sim_rac_diff<-rbind(sim_rac_diff, out)  
}

sim_rac_diff_mean<-sim_rac_diff%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))%>%
  group_by(id3, time)%>%
  summarize(S=mean(richness_diff),
            E=mean(evenness_diff,na.rm=T),
            R=mean(rank_diff),
            D=mean(species_diff))

##multivariate_diff
sim_mult_diff<-data.frame()
com_rep<-unique(sim_diff$id)

for (i in 1:length(com_rep)){
  
  subset<-sim_diff%>%
    filter(id==com_rep[i])
  
  out <- multivariate_difference(df = subset, time.var = "time", species.var = "species", abundance.var = "abundance", replicate.var = "site", treatment.var = "treatment")
  
  out$id<-com_rep[i]
    sim_mult_diff<-rbind(sim_mult_diff, out)  
}

sim_mult_diff_mean<-sim_mult_diff%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))%>%
  group_by(id3, time)%>%
  summarize(composition_diff=mean(composition_diff),
            abs_dispersion_diff=mean(abs_dispersion_diff))

###curve diff
sim_curve_diff<-data.frame()
com_rep<-unique(sim_diff$id)

for (i in 1:length(com_rep)){
  subset<-sim_diff%>%
    filter(id==com_rep[i])
  
  out <- curve_difference(df = subset, time.var = "time", species.var = "species", abundance.var = "abundance", replicate.var = "site", pool = TRUE, treatment.var = "treatment")
 
  out$id<-com_rep[i]
  sim_curve_diff<-rbind(sim_curve_diff, out)  
}

sim_cc_diff<-sim_curve_diff%>% 
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))%>%
  group_by(id3, time)%>%
  summarize(curve_diff=mean(curve_diff))

sim_diff_div<-sim_div_evar%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))%>%
  group_by(id3, time, site)%>%#average over replicates
  summarize(Sp=mean(richness),
            Evar=mean(Evar, na.rm=T))%>%
  ungroup()%>%
  group_by(id3, time)%>%#average over sites
  summarize(Sp=mean(Sp),
            Evar=mean(Evar, na.rm=T))

sim_diff_allmetrics<-sim_rac_diff_mean%>%
  left_join(sim_cc_diff)%>%
   left_join(sim_mult_diff_mean)%>%
  left_join(sim_diff_div)%>%
  separate(id3, into=c("alpha","even","comtype"), sep="_")%>%
  mutate(comtype2 = as.factor(comtype))

###doing the correlations with rich and even for all community types


#rich with delta rank
with(sim_diff_allmetrics, cor.test(Sp, R))
with(subset(sim_diff_allmetrics, comtype =="a"), cor.test(Sp, R))
with(subset(sim_diff_allmetrics, comtype =="b"), cor.test(Sp, R))
with(subset(sim_diff_allmetrics, comtype =="c"), cor.test(Sp, R))
with(subset(sim_diff_allmetrics, comtype =="d"), cor.test(Sp, R))

#rich with species_diff
with(sim_diff_allmetrics,cor.test(Sp, D))
with(subset(sim_diff_allmetrics, comtype =="a"), cor.test(Sp, D))
with(subset(sim_diff_allmetrics, comtype =="b"), cor.test(Sp, D))
with(subset(sim_diff_allmetrics, comtype =="c"), cor.test(Sp, D))
with(subset(sim_diff_allmetrics, comtype =="d"), cor.test(Sp, D))

#rich with composition
with(sim_diff_allmetrics, cor.test(Sp, composition_diff))
with(subset(sim_diff_allmetrics, comtype =="a"), cor.test(Sp, composition_diff))
with(subset(sim_diff_allmetrics, comtype =="b"), cor.test(Sp, composition_diff))
with(subset(sim_diff_allmetrics, comtype =="c"), cor.test(Sp, composition_diff))
with(subset(sim_diff_allmetrics, comtype =="d"), cor.test(Sp, composition_diff))

#rich with dispersion
with(sim_diff_allmetrics, cor.test(Sp, abs_dispersion_diff))
with(subset(sim_diff_allmetrics, comtype =="a"), cor.test(Sp, abs_dispersion_diff))
with(subset(sim_diff_allmetrics, comtype =="b"), cor.test(Sp, abs_dispersion_diff))
with(subset(sim_diff_allmetrics, comtype =="c"), cor.test(Sp, abs_dispersion_diff))
with(subset(sim_diff_allmetrics, comtype =="d"), cor.test(Sp, abs_dispersion_diff))

#rich with curve
with(sim_diff_allmetrics, cor.test(Sp, curve_diff))
with(subset(sim_diff_allmetrics, comtype =="a"), cor.test(Sp, curve_diff))
with(subset(sim_diff_allmetrics, comtype =="b"), cor.test(Sp, curve_diff))
with(subset(sim_diff_allmetrics, comtype =="c"), cor.test(Sp, curve_diff))
with(subset(sim_diff_allmetrics, comtype =="d"), cor.test(Sp, curve_diff))

#even with delta rank
with(sim_diff_allmetrics, cor.test(Evar, R))
with(subset(sim_diff_allmetrics, comtype =="a"), cor.test(Evar, R))
with(subset(sim_diff_allmetrics, comtype =="b"), cor.test(Evar, R))
with(subset(sim_diff_allmetrics, comtype =="c"), cor.test(Evar, R))
with(subset(sim_diff_allmetrics, comtype =="d"), cor.test(Evar, R))

#even with species difference
with(sim_diff_allmetrics, cor.test(Evar, D))
with(subset(sim_diff_allmetrics, comtype =="a"), cor.test(Evar, D))
with(subset(sim_diff_allmetrics, comtype =="b"), cor.test(Evar, D))
with(subset(sim_diff_allmetrics, comtype =="c"), cor.test(Evar, D))
with(subset(sim_diff_allmetrics, comtype =="d"), cor.test(Evar, D))

#even with composition
with(sim_diff_allmetrics, cor.test(Evar, composition_diff))
with(subset(sim_diff_allmetrics, comtype =="a"), cor.test(Evar, composition_diff))
with(subset(sim_diff_allmetrics, comtype =="b"), cor.test(Evar, composition_diff))
with(subset(sim_diff_allmetrics, comtype =="c"), cor.test(Evar, composition_diff))
with(subset(sim_diff_allmetrics, comtype =="d"), cor.test(Evar, composition_diff))

#even with dispersion
with(sim_diff_allmetrics, cor.test(Evar, abs_dispersion_diff))
with(subset(sim_diff_allmetrics, comtype =="a"), cor.test(Evar, abs_dispersion_diff))
with(subset(sim_diff_allmetrics, comtype =="b"), cor.test(Evar, abs_dispersion_diff))
with(subset(sim_diff_allmetrics, comtype =="c"), cor.test(Evar, abs_dispersion_diff))
with(subset(sim_diff_allmetrics, comtype =="d"), cor.test(Evar, abs_dispersion_diff))

#even with curve
with(sim_diff_allmetrics, cor.test(Evar, curve_diff))
with(subset(sim_diff_allmetrics, comtype =="a"), cor.test(Evar, curve_diff))
with(subset(sim_diff_allmetrics, comtype =="b"), cor.test(Evar, curve_diff))
with(subset(sim_diff_allmetrics, comtype =="c"), cor.test(Evar, curve_diff))
with(subset(sim_diff_allmetrics, comtype =="d"), cor.test(Evar, curve_diff))

###looking at some figures

sim_diff_tograph<-sim_diff_allmetrics%>%
  gather(metric, value, S:abs_dispersion_diff)%>%
  filter(metric != "S")%>%
  filter(metric != "E")

labels_diff <-c(composition_diff = "Compositon Difference",
                  curve_diff = "Curve Difference",
                  abs_dispersion_diff = "Dispersion Difference",
                  D = "Species Differences",
                  R = "Rank Difference")
#diff with richness
ggplot(data=sim_diff_tograph, aes(x=as.factor(Sp), y=value, color = comtype))+
  geom_boxplot()+
  scale_color_manual(name = "Community Type", labels = c("High Spatial, High Temporal","Low Spatial, Low Temporal","High Spatail, Low Temporal","Low Spatial, High Temporal"), values= c("gold","orange","darkorange4","red"))+
  xlab("Simulated Community Richness")+
  ylab("Metric Value")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~metric, ncol = 3, scales = "free", labeller = labeller(metric = labels_diff))+
  theme(strip.background = element_rect(fill = 0))


#diff with evenness
ggplot(data=sim_diff_tograph, aes(x=as.factor(even), y=value, color = comtype))+
  geom_boxplot()+
  scale_color_manual(name = "Community Type", labels = c("High Spatial, High Temporal","Low Spatial, Low Temporal","High Spatail, Low Temporal","Low Spatial, High Temporal"), values= c("gold","orange","darkorange4","red"))+
  xlab("Simulated Community Evenness")+
  ylab("Metric Value")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~metric, ncol = 3, scales = "free", labeller = labeller(metric = labels_diff))+
  theme(strip.background = element_rect(fill = 0))+
  scale_x_discrete(labels = c("Low","Mid","High"))

diff_means<-sim_diff_tograph%>%
  group_by(comtype, metric)%>%
  summarize(mean = mean (value),
            se = sd(value)/sqrt(90))


ggplot(data=sim_diff_allmetrics, aes(x=Evar, y=composition_diff))+
  geom_point()+
  facet_wrap(~comtype, scales ="free")

##make all the plots - for richness change
rrd<-ggplot(data=sim_diff_allmetrics, aes(x=Sp, y=R, color = Evar))+
  geom_point()+
  xlab("Simulated Community Richness")+
  ylab("Rank Difference")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~comtype, labeller=labeller(comtype = labels), ncol=4)
rdd<-ggplot(data=sim_diff_allmetrics, aes(x=Sp, y=abs_dispersion_diff, color = Evar))+
  geom_point()+
  xlab("Simulated Community Richness")+
  ylab("Dispersion Difference")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~comtype, labeller=labeller(comtype = labels), ncol=4)
rcd<-ggplot(data=sim_diff_allmetrics, aes(x=Sp, y=curve_diff, color = Evar))+
  geom_point()+
  xlab("Simulated Community Richness")+
  ylab("Curve Difference")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~comtype, labeller=labeller(comtype = labels), ncol=4)
grid.arrange(rrd,rdd, rcd, ncol=1)

##make all the plots - for evenness change
edd<-ggplot(data=sim_diff_allmetrics, aes(x=Evar, y=abs_dispersion_diff, color = as.factor(Sp)))+
  geom_point()+
  xlab("Simulated Community Evenness")+
  ylab("Dispersion Difference")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~comtype, labeller=labeller(comtype = labels), ncol = 4)
ecd<-ggplot(data=sim_diff_allmetrics, aes(x=Evar, y=curve_diff, color = as.factor(Sp)))+
  geom_point()+
  xlab("Simulated Community Evenness")+
  ylab("Composition Difference")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~comtype, labeller=labeller(comtype = labels), ncol = 4)
grid.arrange(ecd, edd, ncol=1)

# example of a curve comparision for the paper ----------------------------


#######example in paper for curve comparision.
time<-c(1,1,1,1,1,1,2,2,2,2,2,2)
sp<-c(1,2,3,4,5,6,1,3,4,6,7,8)
relrank<-c(0.333,0.5, 0.667, 0.167, 1, 0.833, 0.167, 0.583, 0.333, 1, 0.833, 0.583)
cumabund<-c(90,110,125,50,132,131,70,130,110,163,161,150)

df<-data.frame(time, sp, relrank, cumabund)%>%
  arrange(relrank)

result <- df %>%
  do({
    y <- unique(.$time)###assumption this is a length 2 list
    df1 <- filter(., time==y[[1]])
    df2 <- filter(., time==y[[2]])
    sf1 <- stepfun(df1$relrank, c(0, df1$cumabund))
    sf2 <- stepfun(df2$relrank, c(0, df2$cumabund))
    r <- sort(unique(c(0, df1$relrank, df2$relrank)))
    h <- abs(sf1(r) - sf2(r))
    w <- c(diff(r), 0)
    data.frame(Dstar=sum(w*h))#do has to output a dataframe
  })

ggplot(df, aes(x=relrank, y=cumabund, group=time))+
  geom_step(color=time)+
  xlab("Relative Rank")+
  ylab("Cumulative Abundance")+
  theme_bw()

result <- df %>%
  do({
    y <- unique(.$time)###assumption this is a length 2 list
    df1 <- filter(., time==y[[1]])
    df2 <- filter(., time==y[[2]])
    sf1 <- stepfun(df1$relrank, c(0, df1$cumabund))
    sf2 <- stepfun(df2$relrank, c(0, df2$cumabund))
    r <- sort(unique(c(0, df1$relrank, df2$relrank)))
    h <- abs(sf1(r) - sf2(r))
    w <- c(diff(r), 0)
    data.frame(Dstar=sum(w*h))#do has to output a dataframe
  })

ggplot(df, aes(x=relrank, y=cumabund, group=time))+
  geom_step(color=time)+
  xlab("Relative Rank")+
  ylab("Cumulative Abundance")+
  theme_bw()


df <- df[order(df$time, df$cumabund),]

timestep2 <- unique(df$time)#assumes this is a length of 2

df1 <- df[df$time == timestep2[1],]
df2 <- df[df$time == timestep2[2],]
sf1 <- stepfun(df1$relrank, c(0, df1$cumabund))
sf2 <- stepfun(df2$relrank, c(0, df2$cumabund))
r <- sort(unique(c(0, df1$relrank, df2$relrank)))
h <- abs(sf1(r) - sf2(r))
w <- c(diff(r), 0)
CC=sum(w*h)

# How well did the simualtions do -----------------------------------------

###looking at turnover
sim_turnover<-data.frame()
com_rep<-unique(sim$id)

for (i in 1:length(com_rep)){
  
  subset<-sim%>%
    filter(id==com_rep[i])
  
  out <- turnover(df = subset, time.var = "time", species.var = "species", abundance.var = "abundance", replicate.var = "site")
  
  out$id<-com_rep[i]
  
  sim_turnover<-rbind(sim_turnover, out)  
}

sim_turnvoer_ave<-sim_turnover%>% 
  group_by(id, time)%>%
  summarise(total=mean(total))%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))%>%
  group_by(id3, time)%>%
  summarize(turnover=mean(total))%>%
  ungroup()%>%
  separate(id3, into=c("alpha","even","comtype"), sep="_")%>%
  group_by(comtype)%>%
  summarize(turnover=mean(turnover))
  
#Whittacker beta diversity (gamma / average alpha)
gammadiv<-sim%>%
  filter(abundance!=0)%>%
  group_by(id3, time, rep, species)%>%
  summarize(splist=mean(abundance))%>%
  ungroup()%>%
  group_by(id3, time, rep)%>%
  summarize(gamma=length(species))

##richness
S<-function(x){
  x1<-x[x!=0]
  length(x1)
}

rich <- group_by(sim, id3, time, rep, site) %>% 
  summarize(S=S(abundance))%>%
  ungroup()%>%
  group_by(id3, time, rep)%>%
  summarize(alpha=mean(S))

beta_div<-gammadiv%>%
  left_join(rich)%>%
  mutate(wbeta=gamma/alpha)%>%
  group_by(id3, time)%>%
  summarize(wbeta=mean(wbeta))%>%
  ungroup()%>%
  separate(id3, into=c("alpha","even","comtype"), sep="_")%>%
  group_by(comtype)%>%
  summarize(betadiv=mean(wbeta))

