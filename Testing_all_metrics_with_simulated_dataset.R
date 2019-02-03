library(tidyverse)
library(vegan)
library(gridExtra)
library(grid)
library(gtable)
library(gtools)
library(devtools)
install_github("NCEAS/codyn", ref = "anderson")
library(codyn)


theme_set(theme_bw(12))

# Read in Data ------------------------------------------------------------

#simualted dataset
sim<-df%>%#run rcommunity_simulating_communities.R then run RAC_scenarios, which creates df
  mutate(time=as.numeric(iteration), 
         id2=paste(id, site, sep="::"))%>%
  select(-sample, -iteration)%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))

# in the sim dataset, communites are differentiated by alpha (richness), theta (evenness) and scenario (rate of temporal and spatial variability: four scenarios: a: high temporal and high spatial variability; b: low temporal, low spatial variability; c: low temporal, high spatial variability; d: high temporal, low spatial variability"). For each richness-evennes combination (9 combinations) there are each community type. Each of these 10 community types have 10 replicates, called "sites" at a given point in time. Each community type, time, and site is then replicated 10 times.


# Making tables for Appendix 5 -----------------------------------------

###looking at temporal variability
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

### looking at spatial variability using Whittacker's beta diversity (gamma / average alpha)
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



# Richness Evenness Metrics -----------------------------------------------
com_rep<-unique(sim$id)

sim_div_evar<-data.frame()
for (i in 1:length(com_rep)){
  
  subset<-sim%>%
    filter(id==com_rep[i])
  
  out <- community_structure(df = subset, time.var = "time", abundance.var = "abundance", replicate.var = "site", metric = "Evar")
  out$id<-com_rep[i]
  
  sim_div_evar<-rbind(sim_div_evar, out)  
}

sim_diversity_mean<-sim_div_evar%>%
  separate(id, into=c("alpha","theta","scenario","rep"), sep="_", remove=F)%>%
  mutate(id3=paste(alpha, theta, scenario, sep="_"))%>%
  group_by(id3, time, site)%>%#average over replicates
  summarize(Sp=mean(richness),
            Evar=mean(Evar, na.rm=T))%>%
  ungroup()%>%
  group_by(id3, time)%>%#average over sites
  summarize(Sp=mean(Sp),
            Evar=mean(Evar, na.rm=T))


###making figure for appendix 5
sim_div_tograph<-sim_diversity_mean%>%
  separate(id3, into=c("alpha","theta","scenario"), sep="_", remove=F)

ggplot(data = sim_div_tograph, aes(x = Sp, y = Evar, color = scenario))+
  geom_point(size = 2)+
  scale_color_brewer(type = "div", name = "Community Type", labels = c("High Spatial, High Temporal", "Low Spatial, Low Temporal", "High Spatial, Low Temporal", "Low Spatial, High Temporal"))+
  xlab("Simulated Community Richness")+
  ylab("Simulated Community Evenness (Evar)")+
  theme(panel.grid = element_blank())



# Change Metrics -------------------------------------------------------------

### RAC Change
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

##Multivariate Change
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
  summarize(composition_change=mean(composition_change, na.rm = T),
            dispersion_change=mean(dispersion_change, na.rm =T))

## Curve change

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


# Difference Metrics ------------------------------------------------------

sim_diff<-sim%>%
  mutate(treatment = as.factor(ifelse(as.integer(site) < 5, "T1", "T2")))

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
  summarize(composition_diff=mean(composition_diff, na.rm = T),
            abs_dispersion_diff=mean(abs_dispersion_diff, na.rm=T))

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

# Merging all metrics to single datasets ----------------------------------

sim_allmetrics<-sim_multchange_mean%>%
  left_join(sim_rac_change_mean)%>%
  left_join(sim_cc_ave)%>%
  left_join(sim_diversity_mean)%>%
  separate(id3, into=c("alpha","even","comtype"), sep="_")%>%
  mutate(comtype2 = as.factor(comtype))


sim_diff_allmetrics<-sim_rac_diff_mean%>%
  left_join(sim_cc_diff)%>%
  left_join(sim_mult_diff_mean)%>%
  left_join(sim_diff_div)%>%
  separate(id3, into=c("alpha","even","comtype"), sep="_")%>%
  mutate(comtype2 = as.factor(comtype))

# Effect of rich and even on CHANGE and DIFFERNCE sim data Table 2 ------------------------------

#How are CHANGE metrics affected by the richness and evenness of the community?
#Richness
#rich with delta rank
with(sim_allmetrics,  cor.test(Sp, R))

##removing extra division by Sp to demonstate to why it is necessary to divide again by the size of the species pool.
sim_allmetrics2<-sim_allmetrics%>%
  mutate(MRS = R*Sp)
with(sim_allmetrics2, cor.test(Sp, MRS))

ggplot(data=sim_allmetrics2, aes(x=Sp, y=MRS, color = Evar))+
  geom_point()+
  xlab("Simulated Community Richness")+
  ylab("Rank Change")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~comtype, labeller=labeller(comtype = labels), ncol =4)

#rich with gains or lossess same thing
with(sim_allmetrics, cor.test(Sp, G))

#rich with composition
with(sim_allmetrics, cor.test(Sp, composition_change))

#rich with dispersion
with(sim_allmetrics, cor.test(Sp, dispersion_change))

#rich with curve
with(sim_allmetrics, cor.test(Sp, curve_change))

#even with delta rank
with(sim_allmetrics, cor.test(Evar, R))

#even with gains/losses
with(sim_allmetrics, cor.test(Evar, G))

#even with composition
with(sim_allmetrics, cor.test(Evar, composition_change))

#even with dispersion
with(sim_allmetrics, cor.test(Evar, dispersion_change))

#even with curve
with(sim_allmetrics, cor.test(Evar, curve_change))

###doing the correlations for the difference metrics

#rich with delta rank
with(sim_diff_allmetrics, cor.test(Sp, R))

#rich with species_diff
with(sim_diff_allmetrics,cor.test(Sp, D))

#rich with composition
with(sim_diff_allmetrics, cor.test(Sp, composition_diff))

#rich with dispersion
with(sim_diff_allmetrics, cor.test(Sp, abs_dispersion_diff))

#rich with curve
with(sim_diff_allmetrics, cor.test(Sp, curve_diff))

#even with delta rank
with(sim_diff_allmetrics, cor.test(Evar, R))

#even with species difference
with(sim_diff_allmetrics, cor.test(Evar, D))

#even with composition
with(sim_diff_allmetrics, cor.test(Evar, composition_diff))

#even with dispersion
with(sim_diff_allmetrics, cor.test(Evar, abs_dispersion_diff))

#even with curve
with(sim_diff_allmetrics, cor.test(Evar, curve_diff))


# Looking at figures of rich/even on metrics Appendix 6 figures and Table -----------


#Change metrics
sim_tograph<-sim_allmetrics%>%
  gather(metric, value, composition_change:curve_change)%>%
  filter(metric != "S")%>%
  filter(metric != "E")%>%
  filter(metric != "L")


labels_change <-c(composition_change = "Composition Change",
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

sim_means<-sim_allmetrics%>%
  mutate(abs_dc = abs(dispersion_change))%>%
  gather(metric, value, c(composition_change:curve_change, abs_dc))%>%
  group_by(comtype, metric)%>%
  summarize(mean = mean (value),
            se = sd(value)/sqrt(81))



###Difference Metrics

sim_diff_tograph<-sim_diff_allmetrics%>%
  gather(metric, value, S:abs_dispersion_diff)%>%
  filter(metric != "S")%>%
  filter(metric != "E")

labels_diff <-c(composition_diff = "Composition Difference",
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

