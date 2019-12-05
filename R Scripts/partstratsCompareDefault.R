
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
allSets <- c('HYO','CEA','OZL','SCO','THER')
#specify root directory
rootDir <-"/Volumes/MasterData/dxsb43/Partitioning_Strategies/mutations" #path


#initialise data frame to rbind() onto - contains MLs for all datasets
allMLs <- data.frame(data=character(),
                     partstrat = character(),
                     ML = double(),
                     BF = double(),
                     se = double(),
                     brlenspr = character())


for (data in allSets) {
          #extract MLs from .ss files         e.g. HYO '_' anatomy '_ss.nex.ss'
        # add up columns of .ss file
        # if a column has any positive values, discard that column
        # calculate mean of the remaining run MLs -> ML for whole analysis
  #choose group of analyses based on dataset
  group <- if (data == 'CEA') {c('CEA', 'CEAexp', 'CEAclock')
  } else if (data == 'OZL') {c('OZL', 'OZLgamma', 'OZLclock')
  } else if (data == 'SCO') {c('SCO', 'SCOexp', 'SCOgamma')
  } else if (data == 'HYO') {c('HYO', 'HYOexp', 'HYOclock')
  } else if (data == 'THER') {c('THER', 'THERgamma', 'THERclock')}
  
  #initialise data frame to rbind() onto
  MLs <- data.frame(data = character(),
                    partstrat = character(),
                    ML = double(),
                    BF = double(),
                    se = double(),
                    brlenspr = character())
  
        for (dataSet in group) {
          setwd(paste0(rootDir, '/', dataSet))#  dir is e.g. .../mutations/SCO
          #read in optimal tree .nex file - the published tree - use this to calculate QD of startTree from optTree
          optTree <- read.nexus(paste0('../Rest of files/optimal_trees/', dataSet,'_optimal_tree.nex'))
          optTree <- multi2di(optTree)
          
          #read in unpartitioned tree .ss file - use this for calculating Bayes factor
          unpartMat <- read_tsv(list.files(pattern = paste0(dataSet, '_unpart_ss.nex.ss'))[1],comment = '[')
          unpartMat <- unpartMat[,c(3:10)]
          unpartSE <- sd(colSums(unpartMat)[which(colSums(unpartMat) < 0)],na.rm=TRUE)/sqrt(8)  #calculate standard error
          #sum run columns over all steps, DISCARD POSITIVE RUNS, then take the mean
          unpart <- mean(colSums(unpartMat)[which(colSums(unpartMat) < 0)])
          
          
          #read in neotrans .ss file - use this for calculating Bayes factor
          neoMat <- read_tsv(list.files(pattern = paste0(dataSet, '_neotrans_ss.nex.ss'))[1],comment = '[')
          neoMat <- neoMat[,c(3:10)]
          neoSE <- sd(colSums(neoMat)[which(colSums(neoMat) < 0)],na.rm=TRUE)/sqrt(8)  #calculate standard error
          #sum run columns over all steps, DISCARD POSITIVE RUNS, then take the mean
          neo <- mean(colSums(neoMat)[which(colSums(neoMat) < 0)])
          
          #read in neotrans2 .ss file - use this for calculating Bayes factor
          neo2Mat <- read_tsv(list.files(pattern = paste0(dataSet, '_neotrans2_ss.nex.ss'))[1],comment = '[')
          neo2Mat <- neo2Mat[,c(3:10)]
          neo2SE <- sd(colSums(neo2Mat)[which(colSums(neo2Mat) < 0)],na.rm=TRUE)/sqrt(8)  #calculate standard error
          #sum run columns over all steps, DISCARD POSITIVE RUNS, then take the mean
          neo2 <- mean(colSums(neo2Mat)[which(colSums(neo2Mat) < 0)])
          
          #read in anatomy .ss file - use this for calculating Bayes factor
          anaMat <- read_tsv(list.files(pattern = paste0(dataSet, '_anatomy_ss.nex.ss'))[1],comment = '[')
          anaMat <- anaMat[,c(3:10)]
          anaSE <- sd(colSums(anaMat)[which(colSums(anaMat) < 0)],na.rm=TRUE)/sqrt(8)  #calculate standard error
          #sum run columns over all steps, DISCARD POSITIVE RUNS, then take the mean
          ana <- mean(colSums(anaMat)[which(colSums(anaMat) < 0)])
          
          # # OZL & HYO ONLY - read in anatomy2 .ss file - use this for calculating Bayes factor
          # ana2Mat <- read_tsv(list.files(pattern = paste0(dataSet, '_anatomy2_ss.nex.ss'))[1],comment = '[')
          # ana2Mat <- ana2Mat[,c(3:10)]
          # ana2SE <- sd(colSums(ana2Mat)[which(colSums(neoMat) < 0)],na.rm=TRUE)/sqrt(8)  #calculate standard error
          # #sum run columns over all steps, DISCARD POSITIVE RUNS, then take the mean
          # ana2 <- mean(colSums(ana2Mat)[which(colSums(ana2Mat) < 0)])
          
          ## these files do exist      
          #read in IC .ss file - use this for calculating Bayes factor
          ICMat <- read_tsv(list.files(pattern = paste0(dataSet, '_IC_ss.nex.ss'))[1],comment = '[')
          ICMat <- ICMat[,c(3:10)]
          ICSE <- sd(colSums(ICMat)[which(colSums(ICMat) < 0)],na.rm=TRUE)/sqrt(8)  #calculate standard error
          #sum run columns over all steps, DISCARD POSITIVE RUNS, then take the mean
          IC <- mean(colSums(ICMat)[which(colSums(ICMat) < 0)])
          
          #read in F .ss file - use this for calculating Bayes factor
          FMat <- read_tsv(list.files(pattern = paste0(dataSet, '_F_ss.nex.ss'))[1],comment = '[')
          FMat <- FMat[,c(3:10)]
          FSE <- sd(colSums(FMat)[which(colSums(FMat) < 0)],na.rm=TRUE)/sqrt(8)  #calculate standard error
          #sum run columns over all steps, DISCARD POSITIVE RUNS, then take the mean
          F <- mean(colSums(FMat)[which(colSums(FMat) < 0)])
          
          #read in random partition .ss file - use this for calculating Bayes factor
          ranMat <- read_tsv(list.files(pattern = paste0(dataSet, '_random_ss.nex.ss'))[1],comment = '[')
          ranMat <- ranMat[,c(3:10)]
          ranSE <- sd(colSums(ranMat)[which(colSums(ranMat) < 0)],na.rm=TRUE)/sqrt(8)  #calculate standard error
          #sum run columns over all steps, DISCARD POSITIVE RUNS, then take the mean
          ran <- mean(colSums(ranMat)[which(colSums(ranMat) < 0)])
          
          #read in optimal_tree .ss file - use this for calculating Bayes factor
          optMat <- read_tsv(list.files(pattern = paste0(dataSet, '_optTree_ss.nex.ss'))[1],comment = '[')
          optMat <- optMat[,c(3:10)]
          optSE <- sd(colSums(optMat)[which(colSums(optMat) < 0)],na.rm=TRUE)/sqrt(8)  #calculate standard error
          #sum run columns over all steps, discard positive runs, then take the mean
          opt <- mean(colSums(optMat)[which(colSums(optMat) < 0)])
          
          #make matrix of neotransML etc, calculate Bayes Factor for each relative to unpartitioned ML
          otherMLs <- data.frame(neo, neo2, ana,
                                 #ana2,  #OZL and HYO ONLY
                                 IC, F, ran,opt, unpart)
          otherMLs <- gather(otherMLs, key='partstrat', value='ML', factor_key = TRUE)
          otherMLs$BF <- otherMLs$ML - unpart #these are the BFs of the hand-partitioned analyses
          
          otherMLs$se <- c(neoSE, neo2SE, anaSE, 
                           #ana2SE, 
                           ICSE, FSE, ranSE, optSE, unpartSE)
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
          otherMLs$data <- data
          MLs <- rbind(MLs, otherMLs)
        }
allMLs <- rbind(allMLs, MLs)
  }


allMLs$data <- as.factor(allMLs$data)
allMLs$brlenspr <- as.factor(allMLs$brlenspr)
allMLs$ML <- round(allMLs$ML,2)
allMLs$BF <- round(allMLs$BF,2)
allMLs$se <- round(allMLs$se,2)



####################################################################################################################
####################################################################################################################


#plot Bayes Factors - difference of ML of each Bayesian result tree against ML of result tree generated from unpartitioned analysis


##bumps chart:
require(ggrepel)
ggplot(data = allMLs[which(allMLs$brlenspr=='default'),], aes(x=data, y=BF, group=partstrat)) +
  geom_line(aes(colour = partstrat),alpha=0.6, size=1.5) +
  geom_point(aes(colour = partstrat, size=se, alpha=-se))+
  scale_size(guide =FALSE) +
  scale_alpha(guide=FALSE)+
  geom_text_repel(data = allMLs[which(allMLs$brlenspr=='default'),] %>% filter(data == "CEA"),
                   aes(label=partstrat, colour=partstrat),
                  direction="y",
                  nudge_x=-1.2)+
  geom_text_repel(data = allMLs[which(allMLs$brlenspr=='default'),] %>% filter(data == "THER"),
                  aes(label=partstrat, colour=partstrat),
                  direction="y",
                  nudge_x=1)+
  scale_y_continuous(name = "Bayes Factor",
                     breaks = seq(floor_any(min(allMLs$BF),10),max(allMLs$BF),10)) + ####limits = c(min(0,min(bfMat$BayesFactor)),max(bfMat$BayesFactor))
  scale_x_discrete(name = "Dataset",
                   limits = c("CEA", "HYO","OZL","SCO","THER")) +
  theme_light() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()) +
  scale_colour_discrete(name='Partitioning Strategy',
                        breaks=c('neoML',   'neo2ML',   'anaML', 'ICML','FML','ranML', 'optML','unpartML'),
                        labels=c('neotrans','neotrans 2','anatomy','IC','Homoplasy (random)',  'random','Homoplasy (preferred tree)','unpartitioned'))


ggsave(filename = 'summary_allDatasets_error.pdf',
       device = cairo_pdf, 
       path = paste0(rootDir, '/Rest of files/BayesFactorPlots'),
       width = 8,
       height = 6,
       units = 'in'
)
