library(tidyverse)
library(vegan)
library(gridExtra)
library(grid)
library(gtable)
#library(gtools)
library(devtools)
install_github("NCEAS/codyn", ref = "anderson")
library(codyn)

# Read in Data ------------------------------------------------------------

#community dynamics dataset

##read in full dataset
infile1  <- "https://pasta.lternet.edu/package/data/eml/edi/15/5/f69c8fe563067164191d61b6e33eff03" 
infile1 <- sub("^https","http",infile1) 
dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "rowID",     
                 "sitesubplot",     
                 "experiment_year",     
                 "site_code",     
                 "project_name",     
                 "community_type",     
                 "site_project_comm",     
                 "plot_id",     
                 "species",     
                 "relcover"    ), check.names=TRUE)


# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt1$rowID)!="factor") dt1$rowID<- as.factor(dt1$rowID)
if (class(dt1$sitesubplot)!="factor") dt1$sitesubplot<- as.factor(dt1$sitesubplot)
if (class(dt1$experiment_year)=="factor") dt1$experiment_year <-as.numeric(levels(dt1$experiment_year))[as.integer(dt1$experiment_year) ]
if (class(dt1$site_code)!="factor") dt1$site_code<- as.factor(dt1$site_code)
if (class(dt1$project_name)!="factor") dt1$project_name<- as.factor(dt1$project_name)
if (class(dt1$community_type)!="factor") dt1$community_type<- as.factor(dt1$community_type)
if (class(dt1$site_project_comm)!="factor") dt1$site_project_comm<- as.factor(dt1$site_project_comm)
if (class(dt1$plot_id)!="factor") dt1$plot_id<- as.factor(dt1$plot_id)
if (class(dt1$species)!="factor") dt1$species<- as.factor(dt1$species)
if (class(dt1$relcover)=="factor") dt1$relcover <-as.numeric(levels(dt1$relcover))[as.integer(dt1$relcover) ]


codyndat<-dt1%>%
  filter(site_code!="MISS")


##read in site and experiment information
infile2  <- "https://pasta.lternet.edu/package/data/eml/edi/15/5/8284876afe3a1cb0a919d37e1164357f" 
infile2 <- sub("^https","http",infile2) 
dt2 <-read.csv(infile2,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "site_code",     
                 "site_project_comm",     
                 "project_name",     
                 "community_type",     
                 "location",     
                 "Country",     
                 "Continent",     
                 "Lat",     
                 "Long",     
                 "MAP_mm",     
                 "plot_size",     
                 "spatial_extent",     
                 "succession",     
                 "lifespan",     
                 "trophic_level",     
                 "taxa",     
                 "ANPP",     
                 "broad_ecosystem_type",     
                 "num_plots",     
                 "temp_C",     
                 "dataset_length",     
                 "time_step"    ), check.names=TRUE)


# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt2$site_code)!="factor") dt2$site_code<- as.factor(dt2$site_code)
if (class(dt2$site_project_comm)!="factor") dt2$site_project_comm<- as.factor(dt2$site_project_comm)
if (class(dt2$project_name)!="factor") dt2$project_name<- as.factor(dt2$project_name)
if (class(dt2$community_type)!="factor") dt2$community_type<- as.factor(dt2$community_type)
if (class(dt2$location)!="factor") dt2$location<- as.factor(dt2$location)
if (class(dt2$Country)!="factor") dt2$Country<- as.factor(dt2$Country)
if (class(dt2$Continent)!="factor") dt2$Continent<- as.factor(dt2$Continent)
if (class(dt2$Lat)=="factor") dt2$Lat <-as.numeric(levels(dt2$Lat))[as.integer(dt2$Lat) ]
if (class(dt2$Long)=="factor") dt2$Long <-as.numeric(levels(dt2$Long))[as.integer(dt2$Long) ]
if (class(dt2$MAP_mm)=="factor") dt2$MAP_mm <-as.numeric(levels(dt2$MAP_mm))[as.integer(dt2$MAP_mm) ]
if (class(dt2$plot_size)=="factor") dt2$plot_size <-as.numeric(levels(dt2$plot_size))[as.integer(dt2$plot_size) ]
if (class(dt2$spatial_extent)=="factor") dt2$spatial_extent <-as.numeric(levels(dt2$spatial_extent))[as.integer(dt2$spatial_extent) ]
if (class(dt2$succession)!="factor") dt2$succession<- as.factor(dt2$succession)
if (class(dt2$lifespan)!="factor") dt2$lifespan<- as.factor(dt2$lifespan)
if (class(dt2$trophic_level)!="factor") dt2$trophic_level<- as.factor(dt2$trophic_level)
if (class(dt2$taxa)!="factor") dt2$taxa<- as.factor(dt2$taxa)
if (class(dt2$ANPP)=="factor") dt2$ANPP <-as.numeric(levels(dt2$ANPP))[as.integer(dt2$ANPP) ]
if (class(dt2$broad_ecosystem_type)!="factor") dt2$broad_ecosystem_type<- as.factor(dt2$broad_ecosystem_type)
if (class(dt2$num_plots)=="factor") dt2$num_plots <-as.numeric(levels(dt2$num_plots))[as.integer(dt2$num_plots) ]
if (class(dt2$temp_C)=="factor") dt2$temp_C <-as.numeric(levels(dt2$temp_C))[as.integer(dt2$temp_C) ]
if (class(dt2$dataset_length)=="factor") dt2$dataset_length <-as.numeric(levels(dt2$dataset_length))[as.integer(dt2$dataset_length) ]
if (class(dt2$time_step)=="factor") dt2$time_step <-as.numeric(levels(dt2$time_step))[as.integer(dt2$time_step) ]

codyndat_info<-dt2

###CLEANING CODYN DATASET
#restrict to species that are present in an experiment
splist<-codyndat%>%
  group_by(site_code, project_name, community_type, species)%>%
  summarize(present=sum(relcover))%>%
  filter(present!=0)%>%
  select(-present)

#merge back and will drop species that do not exist in a dataset
codyndat_clean<-merge(codyndat, splist, by=c("site_code","project_name","community_type","species"))%>%
  select(-site_project_comm)%>%
  mutate(site_project_comm = paste(site_code, project_name, community_type, sep = "_"))%>%
  select(-sitesubplot, -site_code, -project_name, -community_type)%>%
  mutate(id=paste(site_project_comm, plot_id, sep="::"))



# Richness Evenness Metrics -----------------------------------------------
spc<-unique(codyndat_clean$site_project_comm)
codyn_div_evar<-data.frame()

for (i in 1:length(spc)){
  subset<-codyndat_clean%>%
    filter(site_project_comm==spc[i])
  
  out<-community_structure(subset, time.var = 'experiment_year', abundance.var = 'relcover', replicate.var = 'plot_id', metric = "Evar")
  out$site_project_comm<-spc[i]
  
  codyn_div_evar<-rbind(codyn_div_evar, out)
}


codyndat_diversity_mean <- codyn_div_evar%>%
  group_by(site_project_comm, experiment_year)%>%
  summarize(Sp=mean(richness),
            Evar=mean(Evar, na.rm=T))

# RAC changes -------------------------------------------------------------

codyndat_rac_change<-data.frame()
spc<-unique(codyndat_clean$site_project_comm)

for (i in 1:length(spc)){
  
  subset<-codyndat_clean%>%
    filter(site_project_comm==spc[i])
  
  out <- RAC_change(df = subset, time.var = "experiment_year", species.var = "species", abundance.var = "relcover", replicate.var = "plot_id")
  
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

# Mean Change and Dispersion ----------------------------------------------
codyn_multchange<-data.frame()
spc<-unique(codyndat_clean$site_project_comm)

for (i in 1:length(spc)){
  
  subset<-codyndat_clean%>%
    filter(site_project_comm==spc[i])
  
  out <- multivariate_change(df = subset, time.var = "experiment_year", species.var = "species", abundance.var = "relcover", replicate.var = "plot_id")
  
  out$site_project_comm<-spc[i]

  codyn_multchange<-rbind(codyn_multchange, out)  
}

# Curve change ------------------------------------------------------------

codyn_curvechange<-data.frame()
spc<-unique(codyndat_clean$site_project_comm)

for (i in 1:length(spc)){
  
  subset<-codyndat_clean%>%
    filter(site_project_comm==spc[i])
  
  out <- curve_change(df = subset, time.var = "experiment_year", species.var = "species", abundance.var = "relcover", replicate.var = "plot_id")
  
  out$site_project_comm<-spc[i]
  
  codyn_curvechange<-rbind(codyn_curvechange, out)  
}

codyn_curvechange_mean<-codyn_curvechange%>% 
  group_by(site_project_comm, experiment_year, experiment_year2)%>%
  summarise(curve_change=mean(curve_change))

# Merging all metrics to single datasets ----------------------------------

####MERGING TO A SINGE DATASET

codyndat_allmetrics<-codyn_multchange%>%
  left_join(codyndat_rac_change_average)%>%
  left_join(codyn_curvechange_mean)%>%
  left_join(codyndat_diversity_mean)%>%
  left_join(codyndat_info)%>%
  select(-ANPP)%>%
  na.omit #remove 19 points for which multivariate compostion or dipsersion could not be calcualted


# Making figure 4 --------------------------------------------------------

theme_set(theme_bw(12))

codyndat_allmetrics<-codyndat_allmetrics%>%
  mutate(absS = abs(S),
         absE = abs(E))

#Code to make pairs graph
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

colnames(codyndat_allmetrics)# double check columns are correct
#comparing absolute S and E changes
pairs(codyndat_allmetrics[,c(34,35, 8:10, 3, 4, 11)], col=codyndat_allmetrics$taxa,upper.panel = panel.cor,diag.panel = panel.hist, cex.axis = 2)
par(xpd=T)

# Make table 3 --------

#how do these correlate with experiment parameters.
codyndat_allmetrics3<-codyndat_allmetrics%>%
  group_by(site_project_comm)%>%
  summarize_at(vars(absS, absE, R, G, L, curve_change, composition_change, dispersion_change, spatial_extent, plot_size, num_plots), funs(mean), rm.na=T)%>%
  mutate(spatialExtent=log(spatial_extent),plotSize=log(plot_size))
                 
colnames(codyndat_allmetrics3)
pairs(codyndat_allmetrics3[,c(2:9, 12:14)], upper.panel = panel.cor)
cor(codyndat_allmetrics3[,c(2:9, 12:14)], method = "pearson")
