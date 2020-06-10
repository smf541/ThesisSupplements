#check for relationship between preferred tree, partitioning tree and reconstructed tree
#make matrix of quartet similarities, run linear regressions
#     check:      QS(reconstructed vs preferred tree) vs QS(partitioning vs preferred tree)
#                 QS(reconstructed vs partitioning tree) vs BF  
#written during corrections phase
################################################################################

#packages
require(ggplot2)
require(tidyr)
require(stringr)
require(readr)
require(ape)
require(Quartet)
require(tidyverse)
require(nlcor)
require(devtools)

floor_any = function(x, accuracy, f=floor){f(x/ accuracy) * accuracy}


#set dataset and perturbation move
dataSet <- 'HYO'         #HYO, SCO, THER

#set working directory 
rootDir <-"/Volumes/MasterData/dxsb43/Partitioning_Strategies/mutations"
setwd(paste0(rootDir, '/', dataSet))####################################################   dir is e.g. .../mutations/SCO
if (!dir.exists('SingleStartTrees')) dir.create('SingleStartTrees')


#read in trees from files of 100 perturbed trees
for (ourFile in list.files('StartingTrees', pattern='*.nex')) {
  trees <- read.nexus(paste0('StartingTrees/', ourFile))
  
  for (i in seq_along(trees)) {
    oneTree <- trees[[i]]
    write.nexus(oneTree, file=paste0('SingleStartTrees/', ourFile, '.',i, '_tree.nex'))
  }
}


#read in optimal tree .con.tre file - use this to calculate QD of startTree from optTree
optTree <- read.nexus(paste0('../Rest of files/optimal_trees/',dataSet,'_optimal_tree.nex'))
optTree <- multi2di(optTree)


#initialise tibble for Bayes Factors
bfMat <- as_tibble(matrix(data=NA,nrow=length(list.files('StartingTrees',pattern='*.nex'))*100,ncol=15))
colnames(bfMat) <- c('file','QS_start_to_published','QS_result_to_published','QS_start_to_result','resTreeLength', 'BF1','BF2','BF3','BF4','BF5','BF6','BF7','BF8','BayesFactor','se')
bfMat$file <- str_sort(list.files('MrBayes', pattern = '.lstat'), numeric=TRUE)

#read in lstat file and populate amMat with ML values, indexing by the filename

for (file in str_sort(list.files('MrBayes', pattern = '.lstat'), numeric=TRUE)) {
  
  #####for (i in c(1:100)) {    #e.g. SCO_NNI_chain.nex.1.nex.lstat
  #file <- 'SCO_SPR_chain.nex.4.nex.lstat'       ## for testing
  #E:\dxsb43\Partitioning_Strategies\mutations\SCO\MrBayes
  
  #file <- paste0('./MrBayes/',dataSet, '_', perturbMove, '.nex.',i,'.nex.lstat')  ##needed if looping through using i in c(1:100)
  
  treeFile <- str_replace(file, pattern = '.nex.lstat', '_tree.nex')  ##deletes the .nex.lstat from the file name,    need SCO_NNI_chain.nex.4_tree.nex
  tree <- read.nexus(paste0('SingleStartTrees/', treeFile))
  resTreeFile <- str_replace(file, pattern = '.lstat', '.con.tre')  ##delete .nex.lstat from file name, need SCO_NNI_chain.nex.4.nex.con.tre
  resTree <- read.nexus(paste0('MrBayes/',resTreeFile))
  
  outFile <- read_tsv(paste0('MrBayes/',file), comment = '[')   ## reads the data into a tibble of dim nrow=9, ncol=4
  
  #copy whole column of MLs into row of amMat
  bfMat[which(grepl(file, bfMat$file)),c(6:13)] <- t(outFile[c(1:8),'harmonic_mean'])  # harmonic mean column of tibble is extracted as a vector and transposed
  
  ##calculate quartet divergence of start tree to optimal tree
  # need start tree and optimal tree, e.g.  ./SingleSTartTrees/CEA_NNI_chain.nex.1.nex and    ./CEA_optimal_tree.nex
  startQS <- QuartetStatus(tree ,cf = optTree)
  bfMat[which(grepl(file, bfMat$file)),2] <-  QuartetDivergence(startQS)
  resultQS <- QuartetStatus(resTree, cf = optTree)
  bfMat[which(grepl(file, bfMat$file)),3] <- QuartetDivergence(resultQS)
  startresQS <- QuartetStatus(resTree, cf = tree)
  bfMat[which(grepl(file, bfMat$file)),4] <- QuartetDivergence(startresQS)
}

#plot Bayes Factors - difference of ML of each Bayesian result tree against ML of result tree generated from unpartitioned analysis
#this will tell us if model fit is impacted even if the tree that is used for weighting is not ideal!
#we already know that the ability of MC3 to find a good tree is not impacted negatively 

#calculate standard error in Bayes Factors /MLs!
bfMat$se <- apply(bfMat[,c(6:13)],1,sd)
bfMat$se <- bfMat$se/sqrt(8)

#calculate arithmetic mean of Bayes Factors /MLs!
bfMat$BayesFactor <- apply(bfMat[,c(6:13)],1,mean)


#plot means and error bars representing standard error

#bfMat in long format
bfMatLong <- bfMat[,c(1,2,3,4,5,14,15)]
#bfMatLong <- gather(bfMatLong,key = 'start or result?', value = 'QS', c(2,3),factor_key = TRUE)


##########################################################################################################
############################################ Linear Regression ###########################################
##########################################################################################################
# reconstructed/preferred vs partitioning/preferred
scatter.smooth(x=bfMatLong$QS_start_to_published, y=bfMatLong$QS_result_to_published, main="start/result ~ BF?")
cor(x=bfMatLong$QS_start_to_published, y=bfMatLong$QS_result_to_published)


#reconstructed/partitioning vs BF
scatter.smooth(x=bfMatLong$BayesFactor, y=bfMatLong$QS_start_to_result, main="start ~ result?")
cor(x=bfMatLong$BayesFactor, y=bfMatLong$QS_start_to_result)

nlcor(x=bfMatLong$BayesFactor, y=bfMatLong$QS_start_to_result, plt=T)

