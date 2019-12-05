#calculate Goloboff's f for each character over a set of random trees
require(TreeSearch)
require(ape)
require(phangorn)
require(ggplot2)
require(tidyr)
require(readxl)

# Select dataset
datasetName <- "OZLclock"       #HYO, SYL, SCO, CEA, DINO, THER, OZL
rootDir <- "/Volumes/MasterData/dxsb43/Partitioning_Strategies/mutations/" #path
setwd(rootDir)

# Define constants
nPart <- 3 #number of partitions
prop <- 1/nPart # proportion of characters per partition
insertionComment <- "INSERT PARTITIONS HERE"

mrBayesTemplateFile <- paste0(rootDir, '/',datasetName,'/', datasetName, '_SS_TEMPLATE.nex')
dataset <- ReadAsPhyDat(mrBayesTemplateFile)
# Mr Bayes template
mrBayesTemplate <- readLines(mrBayesTemplateFile)
insertLine <- grep(insertionComment, mrBayesTemplate)

###################################################################################################
############################# read in .nex file of 300 random trees ###############################
###################################################################################################
#how many trees to write to file
rep <- seq_len(300)

#functions
File <- function (suffix) paste0('CharStructure/', datasetName, suffix)  ##path and name of the target file

#read optimal tree
inputTree <- read.nexus(paste0('optimal_trees/',datasetName, '_optimal_tree.nex'))   ###optimal tree is in dir 'mutations', not in dataset subdir
inputTree$edge.length <- NULL
inputLabels <- inputTree$tip.label
inputTree <- multi2di(inputTree)
plot(inputTree)

#which trees to calculate CI for
trees <- read.nexus(paste0('CharStructure/', datasetName, '_random_',length(rep),'.nex')) # reads all 300 random trees 


###################################################################################################
##################### calculate homoplasy for each character, 300 times ###########################
###################################################################################################

powerOf2 <- 2^(0:ncol(attr(dataset, "contrast"))) #contrast shows the possible permutations of 
#the character states, i.e. 0, 1, 2, 3, {01}, {02} etc.
decode <- apply(attr(dataset, "contrast"), 1, function(r) 
  sum(powerOf2[as.logical(r)])
)
tab <- t(vapply(dataset, I, dataset[[1]])) # translates lists of taxa and character data into matrix
tab <- tab[, attr(dataset, 'index')]
nChar <- ncol(tab)
chars <- seq_len(nChar)

minSteps <- apply(tab, 2, function(char) 
  TreeSearch:::MinimumSteps(decode[char])
)


############## calculate max steps       g_max(t,n) = t - [t/n]   for even distribution of states
# where t = number of taxa 
# where n = number of character states
ntax <- nrow(tab)  # number of rows in object tab corresponds to number of taxa, as ncol() corresponds
# to number of characters
stats <- apply(tab, 2, unique)  # unique() tells us how many states a char can take (it gives the elements of the vector with duplicates removed)
#Think about whether
# - inapplicable coding should be a considered a state
nstat <- lapply(stats, length) #length() gives the number of elements in the vector
nstat <- unlist(nstat)      #coerce to numeric vector - nstat now is a num vec of the number of states per character
maxSteps <- ntax - ceiling(ntax / nstat)

# g = t - F       for uneven distribution of states
# where F = number of taxa with the most frequent state


k <- 5  ##3 is bit low
fMax <- (k+1)/(maxSteps+k+1+minSteps)



CImat <- matrix(data=NA,nrow=length(rep), ncol=nChar) #preallocate a matrix of NAs, with as many rows as trees and as many columns as characters
colnames(CImat) <- chars


#loop through trees, calculating f for all characters for each tree i
for (i in seq_along(trees)) {
  tree <- trees[[i]]
  parsScore <- Fitch(tree, dataset)
  
  #calculate number to go in exp() for branch lengths prior = inverse of number of changes along tree per character
  expVal <- ncol(tab)/parsScore
  
  obsSteps <- FitchSteps(tree, dataset)
  obsSteps <- obsSteps[attr(dataset, 'index')] #if two characters have the same profile, they are now not collapsed into one
  
  #calculate Goloboff's unbiased measure of homoplasy for a given k (concavity constant) and data set
  k <- 5
  fObs <- (k+1)/(obsSteps+k+1+minSteps)
  f <- fObs/fMax                          ###normalised f using the maximum f for each char as a measure of scale 
  
  #fill ith row with the vector of CIs
  CImat[i, ] <- f     
}


## with CImat in wide format: take average f for each character (column)
meanF <- colMeans(CImat, dims =1)
meanFdata <- data.frame(c(1:nChar), meanF)
colnames(meanFdata) <- c('character','meanF')       ##x$name[order(x$val)]


#rank characters by information content IC and then divide equally into a number of partitions
mat <- rbind(chars,meanF) #combine char no. and homoplasy value into one matrix
orderedMat <- mat[, order(meanF)] #sort columns by homoplasy (f) in ascending order
chunk <- round(prop * nChar) # number of characters per partition

partA <- orderedMat[1, 1:chunk]
partB <- orderedMat[1, (chunk+1) : (2*chunk)]
partC <- orderedMat[1, (2*chunk+1):nChar]

mrBayesOutput <- c(mrBayesTemplate[seq_len(insertLine - 1)], 
                   paste("prset brlenspr = unconstrained: exp(",expVal, ");"),
                   paste("charset partA =", paste(partA, collapse=' '), ";"),
                   paste("charset partB =", paste(partB, collapse=' '), ";"),
                   paste("charset partC =", paste(partC, collapse=' '), ";"),
                   "",
                   "partition chartype=3: partA, partB, partC;",
                   "set partition=chartype;",
                   "",
                   mrBayesTemplate[(insertLine + 1L):length(mrBayesTemplate)])

outputFile <- paste0(datasetName,'/',datasetName, '_F_ss.nex')
writeLines(mrBayesOutput, outputFile)



