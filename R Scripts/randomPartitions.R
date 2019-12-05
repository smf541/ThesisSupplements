#random partitions
#randomly perturb the list of character numbers
#split into 3 equal partitions

require(TreeSearch)
require(ape)
require(phangorn)
require(ggplot2)
require(tidyr)
require(readxl)


datasets <- c('HYO', 'SYL', 'SCO', 'CEA', 'DINO', 'THER', 'OZL')
for (datasetName in datasets) {
  # Select dataset
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
  
  
  powerOf2 <- 2^(0:ncol(attr(dataset, "contrast"))) #contrast shows the possible permutations of 
  #the character states, i.e. 0, 1, 2, 3, {01}, {02} etc.
  decode <- apply(attr(dataset, "contrast"), 1, function(r) 
    sum(powerOf2[as.logical(r)])
  )
  tab <- t(vapply(dataset, I, dataset[[1]])) # translates lists of taxa and character data into matrix
  tab <- tab[, attr(dataset, 'index')]
  
  
  #make vector of character numbers (1 up to number of columns in tab)
  nChar <- ncol(tab)
  chars <- seq_len(nChar)
  
  #order the df by IC (ascending)
  
  #calculate number to go in exp() for branch lengths prior
  optTree <- read.nexus(paste0('optimal_trees/', datasetName,'_optimal_tree.nex'))
  optTree <- multi2di(optTree)  ###################pick tree, should be 'optimal' one
  parsScore <- Fitch(optTree, dataset)
  
  expVal <- ncol(tab)/parsScore
  #if two characters have the same profile, they are now not collapsed into one
  obsSteps <- FitchSteps(optTree, dataset)
  obsSteps <- obsSteps[attr(dataset, 'index')] 
  
  
  #randomise character order
  chars <- sample(chars)
  
  #divide characters into 3 equal partitions
  chunk <- round(prop * nChar) # number of characters per partition
  
  partA <- chars[1:chunk]
  partB <- chars[(chunk+1) : (2*chunk)]
  partC <- chars[(2*chunk+1):nChar]
  
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
  
  outputFile <- paste0(datasetName,'/',datasetName, '_random_ss.nex')
  writeLines(mrBayesOutput, outputFile)
}


