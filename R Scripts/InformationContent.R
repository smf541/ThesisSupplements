# generates .nex files for a partitioned mcmc analysis with character partitions based on Information Content scores
require(TreeSearch)
require(ape)
require(phangorn)
require(ggplot2)
require(tidyr)
require(readxl)

# Select dataset
datasetName <- "OZLclock"
rootDir <- "/Volumes/MasterData/dxsb43/Partitioning_Strategies/mutations/"
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

CharacterInformation <- function (tokens) {
  tokenCounts <- table(tokens)
  # Our character splits our taxa into groups with the same token
  # ?s and -s are best ignored
  splits <- tokenCounts[!(names(tokenCounts) %in% c('?', '-'))]
  
  # Information content = -log2(probability)
  # Probability of a tree being consistent with our character is
  # n trees consistent with character / n trees with that many tips
  # NUnrootedMult(splits) / NUnrooted(splits)
  # As we are working with large numbers we can use logarithms
  lnP <- LnUnrootedMult(splits) - LnUnrooted(sum(splits))
  log2P <- lnP / log(2)
  information <- -log2P
  
  # Return: 
  information
}

tab[!tab %in% powerOf2] <- '?'
IC <- apply(tab, 2, CharacterInformation)  ##vector of IC values - information content of each character

##bung that into a matrix along with char.no, then sort by IC.
## partition by 

#make vector of character numbers (1 up to number of columns in tab)
nChar <- ncol(tab)
chars <- seq_len(nChar)

#order the df by IC (ascending)
    #calculate number to go in exp() for branch lengths prior
    optTree <- read.nexus(paste0('optimal_trees/',datasetName,'_optimal_tree.nex'))
    optTree <- multi2di(optTree)  ###################pick tree, should be 'optimal' one
    parsScore <- Fitch(optTree, dataset)
    
    expVal <- ncol(tab)/parsScore
    #if two characters have the same profile, they are now not collapsed into one
    obsSteps <- FitchSteps(optTree, dataset)
    obsSteps <- obsSteps[attr(dataset, 'index')] 
    
    #rank characters by information content IC and then divide equally into a number of partitions
    mat <- rbind(chars, IC) #combine char no. and homoplasy value into one matrix
    orderedMat <- mat[, order(IC)] #sort columns by homoplasy (f) in ascending order
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
    
    outputFile <- paste0(datasetName,'/',datasetName, '_IC_ss.nex')
    writeLines(mrBayesOutput, outputFile)



