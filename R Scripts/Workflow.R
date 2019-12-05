#calculate homoplasy indices from starting trees
require(TreeSearch)
require(ape)
require(phangorn)

# Define constants
nPart <- 3 #number of partitions
prop <- 1/nPart # proportion of characters per partition
insertionComment <- "INSERT PARTITIONS HERE"
bayesFilesDir <- 'MrBayes_ss'


# Select dataset
datasetName <- "SCO"    #HYO, SYL, SCO, CEA, DINO, THER, OZL, OZLunc, OZLclock
rootDir <- paste0("E:/dxsb43/Partitioning_Strategies/mutations/", datasetName) #path
setwd(rootDir)
mrBayesTemplateFile <- paste0(rootDir, '/', datasetName, '_SS_TEMPLATE.nex')
dataset <- ReadAsPhyDat(mrBayesTemplateFile)

# Mr Bayes template
mrBayesTemplate <- readLines(mrBayesTemplateFile)
insertLine <- grep(insertionComment, mrBayesTemplate)
if (!dir.exists(bayesFilesDir)) dir.create(bayesFilesDir)

powerOf2 <- 2^(0:ncol(attr(dataset, "contrast"))) #contrast shows the possible permutations of 
                                                  #   the character states, i.e. 0, 1, 2, 3, 4, 5, {01}, {02} etc.
decode <- apply(attr(dataset, "contrast"), 1, function(r) #
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


k <- 5  ##concavity constant
fMax <- (k+1)/(maxSteps+k+1+minSteps)
expVal <- 0.435483870967742

#make MrBayes files from perturbed trees
for (ourFile in list.files('StartingTrees', pattern='*.nex')) {
  trees <- read.nexus(paste0('StartingTrees/', ourFile))
  for (i in seq_along(trees)) {
    tree <- trees[[i]]
    
    obsSteps <- FitchSteps(tree, dataset)
    obsSteps <- obsSteps[attr(dataset, 'index')] #if two characters have the same profile, they are now not collapsed into one
    
    #calculate Goloboff's unbiased measure of homoplasy for a given k (concavity constant) and data set

    fObs <- (k+1)/(obsSteps+k+1+minSteps)
    f <- fObs/fMax   

    
    #rank characters by homoplasy values and then divide equally into a number of partitions
    mat <- rbind(chars, f) #combine char no. and homoplasy value into one matrix
    sortedMat <- mat[, order(f)] #sort columns by homoplasy (f) in ascending order
    chunk <- round(prop * nChar) # number of characters per partition
    
    partA <- sortedMat[1, 1:chunk]
    partB <- sortedMat[1, (chunk+1) : (2*chunk)]
    partC <- sortedMat[1, (2*chunk +1) : nChar]
    
    mrBayesOutput <- c(mrBayesTemplate[seq_len(insertLine - 1)], 
                       paste("prset brlenspr = unconstrained: exp(",expVal, "); [using parsScore of THER_optimal_tree]"), #not for OZLclock
                       paste("charset partA =", paste(partA, collapse=' '), ";"),
                       paste("charset partB =", paste(partB, collapse=' '), ";"),
                       paste("charset partC =", paste(partC, collapse=' '), ";"),
                       "",
                       "partition chartype=3: partA, partB, partC;",
                       "set partition=chartype;",
                       "",
                       mrBayesTemplate[(insertLine + 1L):length(mrBayesTemplate)])
    
    outputFile <- paste0(bayesFilesDir, '/', ourFile, '_ss.', i, '.nex')
    writeLines(mrBayesOutput, outputFile)
  }
}

