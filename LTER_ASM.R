library(codyn)

#make dataframe
species<-c("A","A","B","B","C","C")
abundance<-c(2,18,0,30,67,75)
replicate<-c(1,2,1,2,1,2)
time<-c(1,2,1,2,1,2)

cols<-cbind(species, abundance, replicate, time)
dat<-as.data.frame(cols)
dat$abundance<-as.numeric(as.character(dat$abundance))
str(dat)

#richness and evenness
community_structure(df = dat, abundance.var = "abundance", replicate.var = "replicate")

#diversity
community_diversity(df = dat, abundance.var = "abundance", replicate.var = "replicate")

#rank changes
RAC_change(df = dat, time.var = "replicate", species.var = "species", abundance.var = "abundance" )

#rank diff
RAC_difference(df = dat, species.var = "species",abundance.var = "abundance",replicate.var = "time")

#abundance change
abundance_change(df = dat, time.var = "replicate", species.var = "species", abundance.var = "abundance")

abundance_difference(df = dat, species.var = "species",abundance.var = "abundance",replicate.var = "time")

#curve change
curve_change(df = dat, time.var = "replicate", species.var = "species", abundance.var = "abundance")

#multivariate difference
species2<-c("A","A","A","A", "B","B","B","B","C","C","C","C")
abundance2<-c(2,18,19,15,0,30, 20, 25, 67,75, 70, 60)
replicate2<-c(1,2,3,4,1,2,3,4,1,2,3,4)
treat<-c("C","T","C","T","C","T","C","T")
cols2<-cbind(species2, abundance2, replicate2, treat)
dat2<-as.data.frame(cols2)
dat2$abundance2<-as.numeric(as.character(dat2$abundance2))
dat2$replicate2<-as.numeric(as.character(dat2$replicate2))
str(dat2)

multivariate_difference(df = dat2, species.var = "species2", abundance.var = "abundance2", replicate.var = "replicate2", treatment.var = "treat")

