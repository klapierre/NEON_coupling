################################################################################
##  KNZ_prelim.R: Initial data exploration and figures for the Konza Prairie NEON site only.
##  Author: Kimberly La Pierre
##  Date created: September 6, 2018
################################################################################

# install.packages("devtools")
# install.packages("raster")
# source("http://bioconductor.org/biocLite.R")
# library(devtools)
# install_github("NEONScience/NEON-utilities/neonUtilities")
# install_github("NEONScience/NEON-geolocation/geoNEON")
# biocLite("rhdf5")

library(neonUtilities) #includes stackByTable, zipsByProduct
library(geoNEON) #uses neon spatial data
library(rhdf5)
library(raster)
library(codyn)
library(Hmisc)
library(network)
library(igraph)
library(tidyverse)
options(stringsAsFactors = F) #turns off default that makes any strings into factors

#function to load many csv files into one dataframe by binding rows
load_data <- function(path) { 
  files <- dir(path, pattern = '\\.csv', full.names = TRUE)
  tables <- lapply(files, read.csv)
  do.call(rbind, tables)
}

# setwd('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\grants\\NSF_FY2019\\NSF_ecosystems_FY2018\\NEON_KNZ_data')
# setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\grants\\NSF_FY2019\\NSF_ecosystems_FY2018\\NEON_KNZ_data')


##############################
### beetle data

# stackByTable('NEON_count-beetles.zip') #unzips and merges files; downloaded data had zip files within zip files, and then tons of csv files within that; this function unzips the files and puts them all together into four data files (one for each data step); this should work on most NEON data, but not all (will give an intelligible error message when it can't work; doesn't work on non-tabular data) #note, didn't work for the beetle data, with error that some files are not found (although I see them right there)--unzipping by hand
beetleRaw <- load_data("NEON_count-beetles")%>% #load in all sorted beetle data and bind together using load data function and get sum by species and plot
  filter(sampleType=='carabid')%>%
  group_by(plotID, collectDate, taxonID)%>%
  summarise(beetle_count=sum(individualCount))%>%
  ungroup()

#calculate average beetle count across all sampling dates for each plot
beetleCount <- beetleRaw%>%
  group_by(plotID, collectDate)%>%
  summarise(count=sum(beetle_count))%>%
  ungroup()%>%
  group_by(plotID)%>%
  summarise(beetle_count=mean(count))%>%
  ungroup()

#calculate average beetle richness across all sampling dates for each plot
beetleRichness <- beetleRaw%>%
  group_by(plotID, collectDate)%>%
  summarise(num_taxa=length(taxonID))%>%
  ungroup()%>%
  group_by(plotID)%>%
  summarise(beetle_richness=mean(num_taxa))%>%
  ungroup()

#merge beetle data
beetleSummary <- beetleCount%>%left_join(beetleRichness)
rm(beetleRaw, beetleRichness, beetleCount)


##############################
### bird data

birdRaw <- read.csv('NEON_count-landbird\\NEON.D06.KONZ.DP1.10003.001.brd_countdata.2017-06.expanded.20180418T200232Z.csv')%>% #load in all bird count data and get sum by species and plot
  group_by(plotID, startDate, taxonID)%>%
  summarise(bird_count=sum(clusterSize))%>%
  ungroup()%>%
  filter(taxonID!='')

#calculate average bird count across all sampling dates for each plot
birdCount <- birdRaw%>%
  group_by(plotID, startDate)%>%
  summarise(count=sum(bird_count))%>%
  ungroup()%>%
  group_by(plotID)%>%
  summarise(bird_count=mean(count))%>%
  ungroup()

#calculate average bird richness across all sampling dates for each plot
birdRichness <- birdRaw%>%
  group_by(plotID, startDate)%>%
  summarise(num_taxa=length(taxonID))%>%
  ungroup()%>%
  group_by(plotID)%>%
  summarise(bird_richness=mean(num_taxa))%>%
  ungroup()

#merge bird data
birdSummary <- birdCount%>%left_join(birdRichness)
rm(birdRaw, birdRichness, birdCount)


##############################
### mammal data

mammalRaw <- load_data("NEON_count-small-mammals")%>%  #load in all mammal count data and get sum by species and plot
  filter(trapStatus=='4 - more than 1 capture in one trap'|trapStatus=='5 - capture')%>%
  group_by(plotID, collectDate, taxonID)%>%
  summarise(mammal_count=length(taxonID))%>%
  ungroup()

#calculate average mammal count across all sampling dates for each plot
mammalCount <- mammalRaw%>%
  group_by(plotID, collectDate)%>%
  summarise(count=sum(mammal_count))%>%
  ungroup()%>%
  group_by(plotID)%>%
  summarise(mammal_count=mean(count))%>%
  ungroup()

#calculate average mammal richness across all sampling dates for each plot
mammalRichness <- mammalRaw%>%
  group_by(plotID, collectDate)%>%
  summarise(num_taxa=length(taxonID))%>%
  ungroup()%>%
  group_by(plotID)%>%
  summarise(mammal_richness=mean(num_taxa))%>%
  ungroup()

#merge mammal data
mammalSummary <- mammalCount%>%left_join(mammalRichness)%>%
  mutate(plotID=ifelse(plotID=='KONZ_027', 'KONZ_025', plotID)) #changes plotID for the one plot that is most closely adjacent to where the corresponding bird and base plots are
rm(mammalRaw, mammalRichness, mammalCount)


##############################
### foliar chemistry data

foliarCNRaw <- read.csv('NEON_chem-phys-foliar\\NEON.D06.KONZ.DP1.10026.001.cfc_carbonNitrogen.2017-05.basic.20180619T234037Z.csv')%>%  #load in all foliar C and N data and get sum by plot
  group_by(plotID)%>%
  summarise(foliar_C=mean(carbonPercent), foliar_N=mean(nitrogenPercent))%>%
  ungroup()

foliarElementsRaw <- read.csv('NEON_chem-phys-foliar\\NEON.D06.KONZ.DP1.10026.001.cfc_elements.2017-05.expanded.20180619T234037Z.csv')%>%  #load in all foliar elemental data and get sum by plot
  group_by(plotID)%>%
  summarise(foliar_P=mean(foliarPhosphorusConc), foliar_K=mean(foliarPotassiumConc), foliar_Ca=mean(foliarCalciumConc), foliar_Mg=mean(foliarMagnesiumConc), foliar_S=mean(foliarSulfurConc), foliar_Mn=mean(foliarManganeseConc), foliar_Fe=mean(foliarIronConc), foliar_Cu=mean(foliarCopperConc), foliar_B=mean(foliarBoronConc), foliar_Zn=mean(foliarZincConc))%>%
  ungroup()

ligninRaw <- read.csv('NEON_chem-phys-foliar\\NEON.D06.KONZ.DP1.10026.001.cfc_lignin.2017-05.basic.20180619T234037Z.csv')%>%  #load in all foliar lignin data and get sum by plot
  group_by(plotID)%>%
  summarise(foliar_lignin=mean(ligninPercent), foliar_cellulose=mean(cellulosePercent))%>%
  ungroup()

lmaRaw <- read.csv('NEON_chem-phys-foliar\\NEON.D06.KONZ.DP1.10026.001.cfc_LMA.2017-05.basic.20180619T234037Z.csv')%>%  #load in all foliar Leaf Mass per Area (LMA) data and get sum by plot
  group_by(plotID)%>%
  summarise(LMA=mean(leafMassPerArea))%>%
  ungroup()

#merge foliar chemistry data
foliarChemistrySummary <- foliarCNRaw%>%left_join(foliarElementsRaw)%>%left_join(ligninRaw)%>%left_join(lmaRaw)%>%
  mutate(plotID=ifelse(plotID=='KONZ_010', 'KONZ_006', ifelse(plotID=='KONZ_007', 'KONZ_005', plotID))) #changes plotID for the two plots that are most closely adjacent to where the corresponding bird and base plots are
rm(foliarCNRaw, foliarElementsRaw, ligninRaw, lmaRaw, all=T)


##############################
### plant data
plantRaw <- load_data("NEON_presence-cover-plant")%>%  #load in all plant and get sum by species and plot
  filter(divDataType=='plantSpecies')%>% #filter out all non-plant stuff (like water, rock, etc)
  separate(col=endDate, into=c('year', 'month', 'day'))%>% #split the sampling date into components to get years
  group_by(plotID, subplotID, year, taxonID)%>%
  summarise(cover=max(percentCover))%>% #get max percent cover for each taxon in each subplot for each year
  ungroup()

#calculate plant relative cover
plantTotalCover <- plantRaw%>% 
  group_by(plotID, subplotID, year)%>%
  summarise(total_cover=sum(cover))%>% #get total percent cover for each subplot for each year
  ungroup()

plantRelativeCover <- plantRaw%>%
  left_join(plantTotalCover)%>%
  mutate(plant_cover=cover/total_cover)%>% #calculate relative cover for each species
  group_by(plotID, year, taxonID)%>%
  summarise(year_cover=mean(plant_cover))%>%
  ungroup()%>%
  group_by(plotID, taxonID)%>%
  summarise(plant_cover=mean(year_cover))%>%
  ungroup()%>%
  filter(!is.na(plant_cover))

#calculate plant richness and evenness in each plot
plantRichnessEvenness <- plantRelativeCover%>% #this might need rethinking, because a couple plots have multiple measures across years, so could be artifically inflating richness
  community_structure(abundance.var='plant_cover', replicate.var='plotID')%>%
  mutate(plant_richness=richness, plant_evenness=Evar)%>%
  select(-richness, -Evar)

#plant PCA to get first two axes for all plots
plantWide <- plantRelativeCover%>%
  spread(key=taxonID, value=plant_cover, fill=0)
plantPCA <- prcomp(plantWide[,2:289])
plantAxes <- predict(plantPCA, newdata=plantWide)%>%
  cbind(plantWide[,1])%>%
  select(plotID, PC1, PC2)


#merge plant data
plantSummary <- plantRichnessEvenness%>%left_join(plantAxes)
rm(plantRaw, plantTotalCover, plantRelativeCover, plantRichnessEvenness, plantWide, plantAxes, plantPCA)


##############################
### foliar isotope data

foliarIsotope <- load_data("NEON_isotope-foliar")%>%  #load in all foliar C and N data and get sum by plot
  group_by(plotID)%>%
  summarise(foliar_d15N=mean(d15N), foliar_d13C=mean(d13C))%>%
  ungroup()


##############################
### soil abiotic data

#soil chemistry
soilChemistry <- read.csv('NEON_chem-soil-distrib-initial\\NEON.D06.KONZ.DP1.10008.001.spc_biogeochem.2015-10.expanded.20171027T193900Z.csv')%>%  #load in all foliar C and N data and get sum by plot
  group_by(plotID)%>%
  summarise(soil_ph=mean(phH2o), soil_C=mean(carbonTot), soil_N=mean(nitrogenTot), soil_organic_content=mean(estimatedOC), soil_Su=mean(sulfurTot), soil_Al=mean(alOxalate), soil_Fe=mean(alOxalate), soil_Mn=mean(mnOxalate), soil_P=mean(pOxalate), soil_Si=mean(siOxalate))%>% #means across soil horizons
  ungroup()%>%
  mutate(plotID=ifelse(plotID=='KONZ_020', 'KONZ_003', ifelse(plotID=='KONZ_010', 'KONZ_005', ifelse(plotID=='KONZ_024', 'KONZ_008', plotID)))) #match missing plots to nearest neighbor

#soil particle size
soilPhysical <- read.csv('NEON_phys-soil-distrib-initial\\NEON.D06.KONZ.DP1.10047.001.spc_particlesize.2015-10.basic.20171027T165241Z.csv')%>%  #load in all foliar C and N data and get sum by plot
  group_by(plotID)%>%
  summarise(soil_sand=mean(sandTotal), soil_silt=mean(siltTotal), soil_clay=mean(clayTotal), soil_frag_2to5=mean(coarseFrag2To5), soil_frag_5to20=mean(coarseFrag5To20))%>% #means across soil horizons
  ungroup()%>%
  mutate(plotID=ifelse(plotID=='KONZ_020', 'KONZ_003', ifelse(plotID=='KONZ_010', 'KONZ_005', ifelse(plotID=='KONZ_024', 'KONZ_008', plotID)))) #match missing plots to nearest neighbor


##############################
### merging all data
allData <- beetleSummary%>%full_join(birdSummary)%>%full_join(mammalSummary)%>%full_join(foliarChemistrySummary)%>%full_join(plantSummary)%>%full_join(foliarIsotope)%>%full_join(soilChemistry)%>%full_join(soilPhysical)%>%
  na.omit()



##############################
### get correlation coefficients

#get p values
pMatrix <- rcorr(as.matrix(allData[,-1]), type='pearson') #calculate all possible correlations to get p values
p <- pMatrix$P
p[upper.tri(p, diag=T)] <- 'NA'

allP <- p%>%
  as.data.frame%>%
  rownames_to_column(var = 'var1')%>%
  gather(var2, value, -var1)%>%
  rename(p=value)
allMatrix <- as.matrix(allData[,-1])%>% 
  cor(method='pearson') #calculate all possible correlation
allMatrix[upper.tri(allMatrix, diag=T)] <- 'NA'

#get pearson correlation coefficients
allCorrelation <- allMatrix%>%
  as.data.frame%>%
  rownames_to_column(var = 'var1')%>%
  gather(var2, value, -var1)%>%
  rename(pearson=value)%>%
  left_join(allP)%>% #join with p-values
  filter(pearson!='NA')%>%
  mutate(pearson_sig=ifelse(p>0.05, 0, pearson))

#calculate average correlation coefficient
avgCorrelationSig <- mean(abs(as.numeric(allCorrelation$pearson_sig)))
avgCorrelation <- mean(abs(as.numeric(allCorrelation$pearson)))

##############################
### make figures

#get list of nodes
nodes <- allCorrelation%>%
  mutate(label=as.character(var1))%>%
  select(label)%>%
  unique()%>%
  rowid_to_column('id')%>%
  mutate(type=ifelse(label=='beetle_richness', 'beetle', ifelse(label=='bird_count'|label=='bird_richness', 'bird', ifelse(label=='mammal_count'|label=='mammal_richness', 'mammal', ifelse(label=='soil_ph'|label=='soil_C'|label=='soil_N'|label=='soil_organic_content'|label=='soil_Su'|label=='soil_Al'|label=='soil_Fe'|label=='soil_Mn'|label=='soil_P'|label=='soil_Si'|label=='soil_sand'|label=='soil_silt'|label=='soil_clay'|label=='soil_frag_2to5'|label=='soil_frag_5to20', 'soil', 'plant')))))%>%
  mutate(color=ifelse(type=='soil', '#c00000', ifelse(type=='plant', '#00b050', ifelse(type=='beetle', '#ed7d31', ifelse(type=='mammal', '#4472c4', '#a325a3')))))
nodes[nrow(nodes)+1,] = list(42,'beetle_count', 'beetle', '#ed7d31') #add in beetle_count, which is the one variable not in var1 column (but is in var 2 column)

#get list of edge strengths
edges <- allCorrelation%>%
  left_join(nodes, by=c('var1'='label'))%>%
  select(-var1)%>%
  mutate(var1=as.integer(id))%>%
  select(-id)%>%
  left_join(nodes, by=c('var2'='label'))%>%
  select(-var2)%>%
  mutate(var2=as.integer(ifelse(!is.na(id), id, 19)))%>%
  select(var1, var2, pearson_sig)%>%
  filter(pearson_sig!=0)



###plot network
network <- network(edges, vertex.attr=nodes, matrix.type='edgelist', ignore.eval=T)

plot.network(network, mode='circle', usearrows=F, vertex.cex=2, vertex.col=nodes$color, edge.lwd=1, edge.col='#626363')


###prelim regression plots
theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=30, vjust=-0.8), axis.text.x=element_text(size=25),
             axis.title.y=element_text(size=30, angle=90, vjust=0.7, margin=margin(l=0,r=10,t=0,b=0)), axis.text.y=element_text(size=25),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

ggplot(data=allData, aes(x=LMA, y=beetle_richness)) + #sig
  geom_point(color='#c17313', size=5) + 
  geom_smooth(method='lm', se=F, color='#c17313', size=2) +
  xlab('Leaf Mass Area') + ylab('Ground Beetle\nRichness') +
  annotate("text", x=52, y=3.42, label='r = 0.919\np = 0.003', size=6)
#export at 600x400

ggplot(data=allData, aes(x=LMA, y=beetle_count)) + #non-sig
  geom_point(color='#60605e', size=5) + 
  xlab('Leaf Mass Area') + ylab('Ground Beetle\nAbundance')
#export at 600x400

ggplot(data=allData, aes(x=soil_ph, y=soil_P)) + #sig
  geom_point(color='#c00000', size=5) + 
  geom_smooth(method='lm', se=F, color='#c00000', size=2) +
  xlab('Soil pH') + ylab('\nSoil P') +
  annotate("text", x=6.2, y=275, label='r = 0.763\np = 0.046', size=6)
#export at 600x400



#hypothesized pattern figure
coupling = rnorm(42, mean=0.35, sd=0.15)

#net ecosystem exchange
sigma2 = n*800
eps = rnorm(coupling,mean=0,sd=sqrt(sigma2))
net_ecosystem_exchange = -60+150*coupling + eps
simData <- data.frame(coupling, net_ecosystem_exchange)

ggplot(data=simData, aes(x=coupling, y=net_ecosystem_exchange)) +
  geom_point(color='black', size=5) + 
  xlab('Ecosystem Coupling') + ylab('Net Ecosystem\nExchange') +
  geom_smooth(method='lm', se=F, color='black')


#net ecosystem exchange
sigma2 = n*0.04
eps = rnorm(coupling,mean=0,sd=sqrt(sigma2))
ecosystem_stability = 0+coupling + eps
simData <- data.frame(coupling, ecosystem_stability)

ggplot(data=simData, aes(x=coupling, y=abs(ecosystem_stability))) +
  geom_point(color='black', size=5) + 
  xlab('Ecosystem Coupling') + ylab('Ecosystem\nStability') +
  geom_smooth(method='lm', se=F, color='black')
