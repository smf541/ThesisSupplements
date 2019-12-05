#script to plot Bayes Factors of partitioned phylogenetic analyses. Based on CombiBF_error.R
#plots the bfs for one dataset 

#packages
require(ggplot2)
require(tidyr)
require(stringr)
require(readr)
require(tibble)
require(utils)
require(ape)
require(Quartet)
require(gtools)
require(tidyverse)
require(naniar)
require(grDevices)

floor_any = function(x, accuracy, f=floor){f(x/ accuracy) * accuracy}

#rootdir <- mutations
#files are located in subdirs CEA etc.
#variables for filename and path:     dataset
data <- 'THER'    #HYO - HYOclock - HYOexp, SCO - SCOgamma - SCOexp, THER -THERgamma - THERclock, CEA - CEAexp - CEAclock, OZL - OZLgamma - OZLclock,   
#set working directory 
rootDir <-"/Volumes/MasterData/dxsb43/Partitioning_Strategies/mutations"

#choose group of analyses based on dataset
group <- if (data == 'CEA') {c('CEA', 'CEAexp', 'CEAclock')
  } else if (data == 'OZL') {c('OZL', 'OZLgamma', 'OZLclock')
  } else if (data == 'SCO') {c('SCO', 'SCOexp', 'SCOgamma')
  } else if (data == 'HYO') {c('HYO', 'HYOexp', 'HYOclock')
  } else if (data == 'THER') {c('THER', 'THERgamma', 'THERclock')}

#initialise data frame to rbind() onto
allMLs <- data.frame(partstrat = character(),
                     ML = double(),
                     BF = double(),
                     se = double(),
                     sd = double(),
                     brlenspr = character())

#grab code from CombiBF_error for making data frames of MLs and calculating Bayes Factors

#extract MLs from .ss files         e.g. HYO '_' anatomy '_ss.nex.ss'
# add up columns of .ss file
# if a column has any positive values, discard that column
# calculate mean of the remaining run MLs -> ML for whole analysis
for (dataSet in group) {
  setwd(paste0(rootDir, '/', dataSet))#  dir is e.g. .../mutations/SCO
#read in optimal tree .nex file - the published tree - use this to calculate QD of startTree from optTree
optTree <- read.nexus(paste0('../Rest of files/optimal_trees/', dataSet,'_optimal_tree.nex'))
optTree <- multi2di(optTree)

#read in unpartitioned tree .ss file - use this for calculating Bayes factor
unpartMat <- read_tsv(list.files(pattern = paste0(dataSet, '_unpart_ss.nex.ss'))[1],comment = '[')
unpartMat <- unpartMat[,c(3:10)]
unpartSD <- sd(colSums(unpartMat)[which(colSums(unpartMat) < 0)],na.rm=TRUE)
unpartSE <- unpartSD/sqrt(8)  #calculate standard error
#sum run columns over all steps, DISCARD POSITIVE RUNS, then take the mean
unpartML <- mean(colSums(unpartMat)[which(colSums(unpartMat) < 0)])


#read in neotrans .ss file - use this for calculating Bayes factor
neoMat <- read_tsv(list.files(pattern = paste0(dataSet, '_neotrans_ss.nex.ss'))[1],comment = '[')
neoMat <- neoMat[,c(3:10)]
neoSD <- sd(colSums(neoMat)[which(colSums(neoMat) < 0)],na.rm=TRUE)
neoSE <- neoSD/sqrt(8)  #calculate standard error
#sum run columns over all steps, DISCARD POSITIVE RUNS, then take the mean
neoML <- mean(colSums(neoMat)[which(colSums(neoMat) < 0)])

#read in neotrans2 .ss file - use this for calculating Bayes factor
neo2Mat <- read_tsv(list.files(pattern = paste0(dataSet, '_neotrans2_ss.nex.ss'))[1],comment = '[')
neo2Mat <- neo2Mat[,c(3:10)]
neo2SD <- sd(colSums(neo2Mat)[which(colSums(neo2Mat) < 0)],na.rm=TRUE)
neo2SE <- neo2SD/sqrt(8)  #calculate standard error
#sum run columns over all steps, DISCARD POSITIVE RUNS, then take the mean
neo2ML <- mean(colSums(neo2Mat)[which(colSums(neo2Mat) < 0)])

#read in anatomy .ss file - use this for calculating Bayes factor
anaMat <- read_tsv(list.files(pattern = paste0(dataSet, '_anatomy_ss.nex.ss'))[1],comment = '[')
anaMat <- anaMat[,c(3:10)]
anaSD <- sd(colSums(anaMat)[which(colSums(anaMat) < 0)],na.rm=TRUE)
anaSE <- anaSD/sqrt(8)  #calculate standard error
#sum run columns over all steps, DISCARD POSITIVE RUNS, then take the mean
anaML <- mean(colSums(anaMat)[which(colSums(anaMat) < 0)])

#read in IC .ss file - use this for calculating Bayes factor
ICMat <- read_tsv(list.files(pattern = paste0(dataSet, '_IC_ss.nex.ss'))[1],comment = '[')
ICMat <- ICMat[,c(3:10)]
ICSD <- sd(colSums(ICMat)[which(colSums(ICMat) < 0)],na.rm=TRUE)
ICSE <- ICSD/sqrt(8)  #calculate standard error
#sum run columns over all steps, DISCARD POSITIVE RUNS, then take the mean
ICML <- mean(colSums(ICMat)[which(colSums(ICMat) < 0)])

#read in F .ss file - use this for calculating Bayes factor
FMat <- read_tsv(list.files(pattern = paste0(dataSet, '_F_ss.nex.ss'))[1],comment = '[')
FMat <- FMat[,c(3:10)]
FSD <- sd(colSums(FMat)[which(colSums(FMat) < 0)],na.rm=TRUE)
FSE <- FSD/sqrt(8)  #calculate standard error
#sum run columns over all steps, DISCARD POSITIVE RUNS, then take the mean
FML <- mean(colSums(FMat)[which(colSums(FMat) < 0)])

#read in random partition .ss file - use this for calculating Bayes factor
ranMat <- read_tsv(list.files(pattern = paste0(dataSet, '_random_ss.nex.ss'))[1],comment = '[')
ranMat <- ranMat[,c(3:10)]
ranSD <- sd(colSums(ranMat)[which(colSums(ranMat) < 0)],na.rm=TRUE)
ranSE <- ranSD/sqrt(8)  #calculate standard error
#sum run columns over all steps, DISCARD POSITIVE RUNS, then take the mean
ranML <- mean(colSums(ranMat)[which(colSums(ranMat) < 0)])

#read in optimal_tree .ss file - use this for calculating Bayes factor
optMat <- read_tsv(list.files(pattern = paste0(dataSet, '_optTree_ss.nex.ss'))[1],comment = '[')
optMat <- optMat[,c(3:10)]
optSD <- sd(colSums(optMat)[which(colSums(optMat) < 0)],na.rm=TRUE)
optSE <- optSD/sqrt(8)  #calculate standard error
#sum run columns over all steps, discard positive runs, then take the mean
optML <- mean(colSums(optMat)[which(colSums(optMat) < 0)])

#make matrix of neotransML etc, calculate Bayes Factor for each relative to unpartitioned ML
otherMLs <- data.frame(neoML, neo2ML, anaML,
                       #ana2ML,  #OZL and HYO ONLY
                       ICML, FML, ranML,optML, unpartML)
otherMLs <- gather(otherMLs, key='partstrat', value='ML', factor_key = TRUE)
otherMLs$BF <- otherMLs$ML - unpartML #these are the BFs of the hand-partitioned analyses
otherMLs$sd <- c(neoSD, neo2SD, anaSD, ICSD, FSD, ranSD, optSD, unpartSD)
otherMLs$se <- c(neoSE, neo2SE, anaSE, ICSE, FSE, ranSE, optSE, unpartSE)
#otherMLs$brlenspr <- rep(dataSet)

if (dataSet == data) {
  otherMLs$brlenspr <- rep('default')
} else
  if (dataSet == paste0(data, 'exp')) {
    otherMLs$brlenspr <- rep('exponential')
  } else 
    if (dataSet == paste0(data, 'gamma')) {
      otherMLs$brlenspr <- rep('gamma')
    } else 
      if (dataSet == paste0(data, 'clock')) {
        otherMLs$brlenspr <- rep('clock')
      }

allMLs <- rbind(allMLs, otherMLs)
}


allMLs$brlenspr <- as.factor(allMLs$brlenspr)
allMLs$ML <- round(allMLs$ML,2)
allMLs$BF <- round(allMLs$BF,2)
allMLs$se <- round(allMLs$se,2)
allMLs$sd <- round(allMLs$se,2)

####################################################################################################################
####################################################################################################################
write_delim(allMLs, 
            path = paste0("~/Dropbox/MScR Thesis/Results/MLs_BFs/", data, '_MLsBFs.csv'),
            delim = ';'
            )

####################################################################################################################
####################################################################################################################

#plot Bayes Factors - difference of ML of each Bayesian result tree against ML of result tree generated from unpartitioned analysis
#this will tell us if model fit is impacted if the tree that is used for weighting is not ideal!
#we already know that the ability of MC3 to find a good tree is not impacted negatively 

#plot means and error bars representing standard error
ggplot(data = allMLs) +
  facet_grid(. ~ brlenspr)+ ## only when presenting data for all brlens priors
  geom_hline(aes(yintercept=BF, colour = partstrat)) +
  geom_crossbar(aes(ymin = BF - 2*sd, 
                    ymax = BF + 2*sd,
                    x=1,
                    y=BF,
                    colour = partstrat,
                    fill = partstrat,
                    alpha = 0.3), size = 0.1, show.legend=FALSE)+

  scale_y_continuous(name = "Bayes Factor",
                     breaks = seq(floor_any(min(allMLs$BF),10),ceiling(max(allMLs$BF)+5),10)) + ####limits = c(min(0,min(bfMat$BayesFactor)),max(bfMat$BayesFactor))
  scale_x_continuous(name= data) +
  theme_light() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()) +
  scale_colour_discrete(name='Partitioning Strategy',
                        breaks=c('neoML',   'neo2ML',   'anaML', 'ICML','FML','ranML', 'optML','unpartML'),
                        labels=c('neotrans','neotrans 2','anatomy','IC','Homoplasy (random)',  'random','Homoplasy (preferred tree)','unpartitioned'))


  
####
  
ggsave(filename = paste0('plot_', data,'_sd_facet.pdf'),
       device = cairo_pdf, 
       path = paste0(rootDir, '/Rest of files/BayesFactorPlots'),
       width = 7,
       height = 5,
       units = 'in'
)

