#script to extract MLs and calculate standard dev of MLs from analyses in groups HB and MG. 
#Based on SS_BF.R
#written during corrections phase

#packages
require(tidyr)
require(tidyverse)

floor_any = function(x, accuracy, f=floor){f(x/ accuracy) * accuracy}

#variables for filename and path: dataSet, PS = partstrat, group = mod group, high burnin or more gens
dataSet <- 'oldCEAexp'
PS <- 'neotrans2'
group <- 'HB'

#set working directory 
rootDir <-"/Volumes/MasterData/dxsb43/CorrectionsData"
setwd(paste0(rootDir, '/', group))#  dir is e.g. .../CorrectionsData/HB


#grab code from CombiBF_error for making data frames of MLs and calculating Bayes Factors

#extract MLs from .ss files         e.g. HYO '_' HB '_' anatomy '_ss.nex.ss'
# add up columns of .ss file
# if a column has any positive values, discard that column
# calculate mean of the remaining run MLs -> ML for whole analysis

#read in .ss file - use this for calculating Bayes factor
preMLMat <- read_tsv(list.files(pattern = paste0(dataSet,'_',group,'_',PS,'_ss.nex.ss'))[1],comment = '[')

#check average standard deviation of split frequencies (ideal <0.01, good <0.05)
split <- max(preMLMat[, 11])
print(split)

MLMat <- preMLMat[,c(3:10)]
#check standard deviation of step MLs (if max is low, all is probably fine. If dodgy, view matrix)
max(apply(MLMat, 1, sd))

#check standard deviation of run MLs
sd(colSums(MLMat))

#overall logML: sum run columns over all steps, DISCARD POSITIVE RUNS, then take the mean
ML <- mean(colSums(MLMat)[which(colSums(MLMat) < 0)])
print(ML)

####################################################################################################################
####################################################################################################################
