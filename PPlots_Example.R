##worked example with pplots
setwd("~/Dropbox/converge_diverge/datasets/Longform")
setwd("C:\\Users\\megha\\Dropbox\\SESYNC\\SESYNC_RACs\\For R package")

library(tidyverse)
library(vegan)
library(gridExtra)
library(grid)
library(gtable)
library(devtools)
install_github("mavolio/codyn", ref = "RACs_cleaner")
library(codyn)


theme_set(theme_bw(12))

# read in data ------------------------------------------------------------


dat<-read.csv("~/Dropbox/converge_diverge/datasets/Longform/SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)

dat<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\Longform\\SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)
traits<-read.csv("C:\\Users\\megha\\Dropbox\\pplots\\site_review_2017\\traits_final.csv")
traits<-read.csv("~/Dropbox/pplots/site_review_2017/traits_final.csv")

##pplots
pplots<-dat%>%
  filter(project_name=="pplots",calendar_year==2002|calendar_year==2011, treatment=="N1P0"|treatment=="N2P3"|treatment=="N2P0")%>%
  tbl_df()%>%
  group_by(calendar_year, plot_id)%>%
    mutate(rank=rank(-relcov, ties.method="average"))%>%
  ungroup()%>%
  mutate(treatment=factor(treatment, levels=c("N1P0","N2P0", "N2P3")))###need to do this b/c otherwise R will remember every treatment

trts<-pplots%>%
  ungroup()%>%
  select(plot_id, treatment)%>%
  unique()

# NMDS --------------------------------------------------------------------


##step 1. do NMDS of pretreatment and last year of data
pplots_wide<-dat%>%
  filter(project_name=="pplots",calendar_year==2002|calendar_year==2011,treatment=="N1P0"|treatment=="N2P3"|treatment=="N2P0")%>%
  select(treatment, calendar_year, plot_id, genus_species, relcov)%>%
  spread(genus_species, relcov, fill=0)

plots<-pplots_wide[,1:3]
mds<-metaMDS(pplots_wide[,4:53], autotransform=FALSE, shrink=FALSE)
mds

adonis(pplots_wide[,4:53]~treatment, pplots_wide)

#differences in dispersion?
dist<-vegdist(pplots_wide[,4:53])
betadisp<-betadisper(dist,pplots_wide$treatment,type="centroid")
betadisp
permutest(betadisp)

scores <- data.frame(scores(mds, display="sites"))  # Extracts NMDS scores for year "i" #
scores2<- cbind(plots, scores) # binds the NMDS scores of year i to all years previously run


toplot<-scores2
##plot NMDS
ggplot(subset(toplot, treatment=="N1P0"), aes(x=NMDS1, y=NMDS2, color=as.factor(calendar_year)))+
  geom_point(size=5)+
  scale_color_manual(name="", values=c("black","red"), labels=c("Pre-Treatment","Experiment Year 12"))+
  scale_x_continuous(limits=c(-1,1.2))+
  scale_y_continuous(limits=c(-1.1,1.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none")

ggplot(subset(toplot, treatment=="N2P3"), aes(x=NMDS1, y=NMDS2, color=as.factor(calendar_year)))+
  geom_point(size=5)+
  scale_color_manual(name="", values=c("black","red"), labels=c("Pre-Treatment","Experiment Year 12"))+
  scale_x_continuous(limits=c(-1,1.2))+
  scale_y_continuous(limits=c(-1.1,1.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")

ggplot(subset(toplot, treatment=="N2P0"), aes(x=NMDS1, y=NMDS2, color=as.factor(calendar_year)))+
  geom_point(size=5)+
  scale_color_manual(name="", values=c("black","red"), labels=c("Pre-Treatment","Experiment Year 12"))+
  scale_x_continuous(limits=c(-1,1.2))+
  scale_y_continuous(limits=c(-1.1,1.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")


# Make RACs ---------------------------------------------------------------


##plot RACS

traits2<-traits%>%
  mutate(genus=tolower(Genus),
         genus_species=paste(genus, Species, sep=" "))

ractoplot<-merge(pplots, traits2, by="genus_species", all=T)%>%
  select(-Spnum, -Genus, -Species, -Family, -cot, -canopy, -clonality, -native, -bloom, -Seed_weight__g_, -genus)%>%
  mutate(cat=ifelse(life=="P"&form=="G"&C3_C4=="C4","G1", ifelse(life=="P"&form=="G"&C3_C4=="C3","G2", ifelse(life=="A"&form=="G", "G3", ifelse(life=="P"&form=="F"&n_fixer=="N", "F1",ifelse(life=="P"&form=="S"&n_fixer=="N", "F1", ifelse(life=="P"&form=="F"&n_fixer=="Y", "F2", ifelse(life=="P"&form=="S"&n_fixer=="Y","F2", "F3"))))))))%>%
  na.omit%>%
  mutate(colorfill=ifelse(genus_species=="andropogon gerardii","ag",ifelse(genus_species=="andropogon scoparius","as", ifelse(genus_species=="sorghastrum nutans", "sn", ifelse(genus_species=="solidago canadensis","sc",ifelse(genus_species=="solidago missouriensis", "sm", ifelse(genus_species=="oxalis stricta", "os", ifelse(genus_species=="dichanthelium oligosanthes","do","other"))))))))

pplots_ave<-pplots%>%
  group_by(calendar_year, treatment, genus_species)%>%
  summarise(relcov=mean(relcov))

#label top three species in control plots blue, top 3 in N or N+P purple and top 3 in NP red.

#top 3 control: ag, as, sn, top 3 n/np: os, sc, sm, do

#ag=blue, as=lightblue, do=darkred, os = pink, other = green3, sc = red, sn = cornflowerblue, sm = orange

#contorls
ggplot(data=subset(ractoplot, treatment=="N1P0"&calendar_year==2002) , aes(x=rank, y=relcov))+
   geom_line(aes(group=plot_id), color="darkgray", size=1)+
   geom_point(aes(color=colorfill), size=5)+
  scale_color_manual(values = c("blue","lightblue","darkred", "pink","green3","orange","cornflowerblue"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
ggplot(data=subset(ractoplot, treatment=="N1P0"&calendar_year==2011) , aes(x=rank, y=relcov))+
  geom_line(aes(group=plot_id), color="darkgray", size=1)+
  geom_point(aes(color=colorfill), size=5)+
  scale_color_manual(values = c("blue","lightblue","darkred","pink","green3","orange","cornflowerblue"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
#N+P
ggplot(data=subset(ractoplot, treatment=="N2P3"&calendar_year==2002) , aes(x=rank, y=relcov))+
  geom_line(aes(group=plot_id), color="darkgray", size=1)+
  geom_point(aes(color=colorfill), size=5)+
  scale_color_manual(values = c("blue","lightblue","darkred", "green3","red","orange","cornflowerblue"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
ggplot(data=subset(ractoplot, treatment=="N2P3"&calendar_year==2011) , aes(x=rank, y=relcov))+
  geom_line(aes(group=plot_id), color="darkgray", size=1)+
  geom_point(aes(color=colorfill), size=5)+
  scale_color_manual(values = c("blue","lightblue","darkred", "pink","green3","red","orange","cornflowerblue"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
#N only
ggplot(data=subset(ractoplot, treatment=="N2P0"&calendar_year==2002) , aes(x=rank, y=relcov))+
  geom_line(aes(group=plot_id), color="darkgray", size=1)+
  geom_point(aes(color=colorfill), size=5)+
  scale_color_manual(values = c("blue","lightblue","darkred", "green3","orange","cornflowerblue"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
ggplot(data=subset(ractoplot, treatment=="N2P0"&calendar_year==2011) , aes(x=rank, y=relcov))+
  geom_line(aes(group=plot_id), color="darkgray", size=1)+
  geom_point(aes(color=colorfill), size=5)+
  scale_color_manual(values = c("blue","lightblue",'darkred',"pink","green3","orange","cornflowerblue"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())

# make CC examples --------------------------------------------------------

###figuring out cumulative curve
average<-ractoplot%>%
  group_by(calendar_year, treatment, genus_species)%>%
  summarize(relcov=mean(relcov))%>%
  mutate(sumabund=sum(relcov),
         relcov2=relcov/sumabund)%>%
  mutate(rank=rank(-relcov2, ties.method="average"),
         maxrank=max(rank),
         relrank=rank/maxrank)%>%
  arrange(relrank)%>%
  mutate(cumabund=cumsum(relcov2))%>%
  mutate(year=as.character(calendar_year))
  
ccplot<-ractoplot%>%
  mutate(sumabund=sum(relcov),
         relcov2=relcov/sumabund)%>%
  mutate(rank=rank(-relcov2, ties.method="average"),
         maxrank=max(rank),
         relrank=rank/maxrank)%>%
  arrange(relrank)%>%
  mutate(cumabund=cumsum(relcov2))%>%
  mutate(year=as.character(calendar_year))


ggplot(data=subset(ccplot, treatment=="N1P0"), aes(x=relrank, y=cumabund, color=year))+
  geom_step(size=2)+
  scale_color_manual(values=c("black","red"))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
    ylab("Cumulative Relative Abundance")+
    xlab("Relative Rank")+
  scale_x_continuous(breaks=c(0.25, 0.75))+
  facet_wrap(~plot_id, ncol=3)+
  theme(strip.background = element_blank(),strip.text.x = element_blank())


ggplot(data=subset(ccplot, treatment=="N2P3"), aes(x=relrank, y=cumabund, color=year))+
  geom_step(size=2)+
  scale_color_manual(values=c("black","red"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  ylab("Cumulative Relative Abundance")+
  xlab("Relative Rank")+
  scale_x_continuous(breaks=c(0.25, 0.75))+
  facet_wrap(~plot_id, ncol=3)+
  theme(strip.background = element_blank(),strip.text.x = element_blank())

ggplot(data=subset(ccplot, treatment=="N2P0"), aes(x=relrank, y=cumabund, color=year))+
  geom_step(size=2)+
  scale_color_manual(values=c("black","red"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  ylab("Cumulative Relative Abundance")+
  xlab("Relative Rank")+
  scale_x_continuous(breaks=c(0.25, 0.75))+
  facet_wrap(~plot_id, ncol=3)+
  theme(strip.background = element_blank(),strip.text.x = element_blank())


# doing RAC change and curve change -------------------------------------------------------------

rac <- RAC_change(pplots, time.var = "calendar_year", species.var = "genus_species", abundance.var = "relcov", replicate.var = "plot_id")

##doing curve change
cc <- curve_change(pplots, time.var = "calendar_year", species.var = "genus_species", abundance.var = "relcov", replicate.var = "plot_id")

allmetrics<-rac%>%
  left_join(cc)%>%
  left_join(trts)%>%
  gather(metric, value, richness_change:curve_change)%>%
  group_by(treatment, calendar_year, calendar_year2, metric)%>%
    summarize(vmean=mean(value),
              vn=length(plot_id),
              vsd=sd(value))%>%
  mutate(vse=vsd/sqrt(vn))
            
theme_set(theme_bw(12))
S<-ggplot(data=subset(allmetrics, metric=="richness_change"), aes(x=treatment, y=vmean))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=vmean-vse, ymax=vmean+vse), width=.2)+
  scale_x_discrete(name="Treatment", labels=c("Control", "N", "N+P"))+
  ylab("Richness Change")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
E<-ggplot(data=subset(allmetrics, metric=="evenness_change"), aes(x=treatment, y=vmean))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=vmean-vse, ymax=vmean+vse), width=.2)+
  scale_x_discrete(name="Treatment", labels=c("Control", "N","N+P"))+
  ylab("Evenness Change")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
R<-ggplot(data=subset(allmetrics, metric=="rank_change"), aes(x=treatment, y=vmean))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=vmean-vse, ymax=vmean+vse), width=.2)+
  scale_x_discrete(name="Treatment", labels=c("Control", "N", "N+P"))+
  ylab("Rank Change")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  annotate('text', label="C", x=1, y = 0.175, size=5)+
  annotate('text', label="B", x=2, y = 0.24, size=5)+
  annotate('text', label="A", x=3, y = 0.28, size=5)
G<-ggplot(data=subset(allmetrics, metric=="gains"), aes(x=treatment, y=vmean))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=vmean-vse, ymax=vmean+vse), width=.2)+
  scale_x_discrete(name="Treatment", labels=c("Control", "N","N+P"))+
  ylab("Speices Gains")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
L<-ggplot(data=subset(allmetrics, metric=="losses"), aes(x=treatment, y=vmean))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=vmean-vse, ymax=vmean+vse), width=.2)+
  scale_x_discrete(name="Treatment", labels=c("Control", "N","N+P"))+
  ylab("Species Losses")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
C<-ggplot(data=subset(allmetrics, metric=="curve_change"), aes(x=treatment, y=vmean))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=vmean-vse, ymax=vmean+vse), width=.2)+
  scale_x_discrete(name="Treatment", labels=c("Control", "N","N+P"))+
  ylab("Curve Change")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

grid.arrange(S, E,R, G, L, C, ncol=3)

allmetrics_full<-rac%>%
  left_join(cc)%>%
  left_join(trts)%>%
  gather(metric, value, richness_change:curve_change)

summary(aov(value~treatment, data=subset(allmetrics_full, metric=="richness_change")))# p = 0.906
summary(aov(value~treatment, data=subset(allmetrics_full, metric=="evenness_change")))  # p = 0.053
summary(aov(value~treatment, data=subset(allmetrics_full, metric=="rank_change")))  # p < 0.001
TukeyHSD(aov(value~treatment, data=subset(allmetrics_full, metric=="rank_change")))

summary(aov(value~treatment, data=subset(allmetrics_full, metric=="gains"))) # p = 0.353
summary(aov(value~treatment, data=subset(allmetrics_full, metric=="losses"))) # p = 0.174
summary(aov(value~treatment, data=subset(allmetrics_full, metric=="curve_change"))) # p = 0.681


# appendix fig of change difference through time --------------------------
pplots_allyears<-dat%>%
  filter(project_name=="pplots", treatment=="N1P0"|treatment=="N2P3"|treatment=="N2P0")

rac_allyears<-RAC_change(pplots_allyears, time.var = "calendar_year", species.var = "genus_species", abundance.var = "relcov", replicate.var = "plot_id")
cc_allyears<- curve_change(pplots_allyears, time.var = "calendar_year", species.var = "genus_species", abundance.var = "relcov", replicate.var = "plot_id")
mult_change_allyears<-multivariate_change(pplots_allyears, time.var = "calendar_year", species.var = "genus_species", abundance.var = "relcov", replicate.var = "plot_id", treatment.var = "treatment")

rac_cc_mean<-rac_allyears%>%
  left_join(cc_allyears)%>%
  left_join(trts)%>%
  gather(metric, value, richness_change:curve_change)%>%
  group_by(treatment, calendar_year, calendar_year2, metric)%>%
  summarize(vmean=mean(value),
            vn=length(plot_id),
            vsd=sd(value))%>%
  mutate(vse=vsd/sqrt(vn))

theme_set(theme_bw(12))
bc<-
  ggplot(data=mult_change_allyears, aes(x=calendar_year2, y=composition_change, group=treatment, color=treatment))+
  geom_point(size=3)+
  scale_color_manual(name="Treatment", values=c("black","red","purple"))+
  geom_line(size=1)+
  ylab("Compositional Change")+
  xlab("Year-Pair")+
  scale_x_continuous(breaks=c(2004, 2008, 2012))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.title.x = element_blank())
disp<-
  ggplot(data=mult_change_allyears, aes(x=calendar_year2, y=dispersion_change, group=treatment, color=treatment))+
  geom_point(size=3)+
  scale_color_manual(name="Treatment", values=c("black","red","purple"))+
  geom_line(size=1)+
  ylab("Dispersion Change")+
  xlab("Year-Pair")+
  scale_x_continuous(breaks=c(2004, 2008, 2012))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.title.x = element_blank())


S<-
ggplot(data=subset(rac_cc_mean,metric=="richness_change"), aes(x=calendar_year2, y=vmean, color=treatment))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=vmean-vse, ymax=vmean+vse), width=.2)+
  scale_color_manual(values=c("black","red","purple"))+
  ylab("Richness Change")+
  geom_line(size=1, aes(group=treatment))+
  scale_x_continuous(breaks=c(2004, 2008, 2012))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.title.x=element_blank())
E<-
  ggplot(data=subset(rac_cc_mean,metric=="evenness_change"), aes(x=calendar_year2, y=vmean, color=treatment))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=vmean-vse, ymax=vmean+vse), width=.2)+
  scale_color_manual(values=c("black","red","purple"))+
  ylab("Evenness Change")+
  geom_line(size=1, aes(group=treatment))+
  scale_x_continuous(breaks=c(2004, 2008, 2012))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.title.x=element_blank())
R<-
  ggplot(data=subset(rac_cc_mean,metric=="rank_change"), aes(x=calendar_year2, y=vmean, color=treatment))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=vmean-vse, ymax=vmean+vse), width=.2)+
  scale_color_manual(values=c("black","red","purple"))+
  ylab("Rank Change")+
  geom_line(size=1, aes(group=treatment))+
  scale_x_continuous(breaks=c(2004, 2008, 2012))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.title.x=element_blank())
G<-
  ggplot(data=subset(rac_cc_mean,metric=="gains"), aes(x=calendar_year2, y=vmean, color=treatment))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=vmean-vse, ymax=vmean+vse), width=.2)+
  scale_color_manual(values=c("black","red","purple"))+
  ylab("Species Gains")+
  geom_line(size=1, aes(group=treatment))+
  scale_x_continuous(breaks=c(2004, 2008, 2012))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.title.x=element_blank())
L<-
  ggplot(data=subset(rac_cc_mean,metric=="losses"), aes(x=calendar_year2, y=vmean, color=treatment))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=vmean-vse, ymax=vmean+vse), width=.2)+
  scale_color_manual(name = "Treatment", label=c("Control", "N", "N+P"),values=c("black","red","purple"))+
  ylab("Species Losses")+
  geom_line(size=1, aes(group=treatment))+
  scale_x_continuous(breaks=c(2004, 2008, 2012))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x=element_blank())

legend=gtable_filter(ggplot_gtable(ggplot_build(L)), "guide-box") 
grid.draw(legend)


grid.arrange(arrangeGrob(bc+theme(legend.position="none"),
                         disp+theme(legend.position="none"),
                         S+theme(legend.position="none"),
                         E+theme(legend.position="none"),
                         R+theme(legend.position="none"),
                         G+theme(legend.position="none"),
                         L+theme(legend.position="none"),
                         ncol=2),legend, 
             widths=unit.c(unit(1, "npc") - legend$width, legend$width),nrow=1)



# Differences -------------------------------------------------------------

rac_diff<-RAC_difference(df = pplots, time.var="calendar_year", species.var = "genus_species", abundance.var = "relcov", treatment.var = "treatment", replicate.var = "plot_id", pool = TRUE)

cc_diff<-curve_difference(df = pplots, time.var="calendar_year", species.var = "genus_species", abundance.var = "relcov", treatment.var = "treatment", replicate.var = "plot_id", pool = TRUE)

##graph all differences for all years
rac_diff_allyears<-RAC_difference(pplots_allyears, time.var="calendar_year", species.var = "genus_species", abundance.var = "relcov", treatment.var = "treatment", replicate.var = "plot_id", pool = TRUE)

cc_diff_allyears<-curve_difference(pplots_allyears, time.var="calendar_year", species.var = "genus_species", abundance.var = "relcov", treatment.var = "treatment", replicate.var = "plot_id", pool = TRUE)

mult_diff_allyears<-multivariate_difference(pplots_allyears, time.var="calendar_year", species.var = "genus_species", abundance.var = "relcov", treatment.var = "treatment", replicate.var = "plot_id")%>%
  mutate(group1=paste(treatment, treatment2, sep="_"))

bc_d<-
ggplot(data=mult_diff_allyears, aes(x=as.numeric(calendar_year), y=composition_diff, color=group1, group=group1))+
  geom_point(size=3)+
  geom_line(size=1)+
  scale_color_manual(name = "Treatment\nComparision", label=c("Control - N", "Control - N+P", "N - N+P"),values=c("orange","darkgoldenrod4","darkred"))+
  ylab("Compositional Difference")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank())
disp_d<-
  ggplot(data=mult_diff_allyears, aes(x=as.numeric(calendar_year), y=abs_dispersion_diff, color=group1, group=group1))+
  geom_point(size=3)+
  geom_line(size=1)+
  scale_color_manual(name = "Treatment\nComparision", label=c("Control - N", "Control - N+P", "N - N+P"),values=c("orange","darkgoldenrod4","darkred"))+
  ylab("Absolute Disperion Difference")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank())
s_d<-
  ggplot(data=rac_diff_allyears, aes(x=as.numeric(calendar_year), y=richness_diff, color=group1, group=group1))+
  geom_point(size=3)+
  geom_line(size=1)+
  scale_color_manual(name = "Treatment\nComparision", label=c("Control - N", "Control - N+P", "N - N+P"),values=c("orange","darkgoldenrod4","darkred"))+
  ylab("Compositional Difference")+
  ylab("Richness Difference")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank())
e_d<-
  ggplot(data=rac_diff_allyears, aes(x=as.numeric(calendar_year), y=evenness_diff, color=group1, group=group1))+
  geom_point(size=3)+
  geom_line(size=1)+
  scale_color_manual(name = "Treatment\nComparision", label=c("Control - N", "Control - N+P", "N - N+P"),values=c("orange","darkgoldenrod4","darkred"))+
  ylab("Compositional Difference")+
  ylab("Evenness Difference")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank())
r_d<-
  ggplot(data=rac_diff_allyears, aes(x=as.numeric(calendar_year), y=rank_diff, color=group1, group=group1))+
  geom_point(size=3)+
  geom_line(size=1)+
  scale_color_manual(name = "Treatment\nComparision", label=c("Control - N", "Control - N+P", "N - N+P"),values=c("orange","darkgoldenrod4","darkred"))+
  ylab("Compositional Difference")+
  ylab("Rank Difference")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank())
sp_d<-
  ggplot(data=rac_diff_allyears, aes(x=as.numeric(calendar_year), y=species_diff, color=group1, group=group1))+
  geom_point(size=3)+
  geom_line(size=1)+
  scale_color_manual(name = "Treatment\nComparision", label=c("Control - N", "Control - N+P", "N - N+P"),values=c("orange","darkgoldenrod4","darkred"))+
  ylab("Compositional Difference")+
  ylab("Species Difference")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank())
  

legend=gtable_filter(ggplot_gtable(ggplot_build(sp_d)), "guide-box") 
grid.draw(legend)


grid.arrange(arrangeGrob(bc_d+theme(legend.position="none"),
                         disp_d+theme(legend.position="none"),
                         s_d+theme(legend.position="none"),
                         e_d+theme(legend.position="none"),
                         r_d+theme(legend.position="none"),
                         sp_d+theme(legend.position="none"),
                         ncol=2),legend, 
             widths=unit.c(unit(1, "npc") - legend$width, legend$width),nrow=1)


###examples for the appendix
Evar <- function(x, S = length(x)) {
  x1 <- x[x!=0]
  lnx <- log(x1)
  theta <- (S - 1) / S * var(lnx)
  return(1 - 2 / pi * atan(theta))
} 

t1=c(40,20,15,50,1,6,0,0)
Evar(t1)
t2=c(70,0,20,40,0,2,11,20)
Evar(t2)