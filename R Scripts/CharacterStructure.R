#is there a correlation between mean homoplasy scores and information content of a character?


#calculate Goloboff's f for each character over a set of random trees --- mod from fPartitioning.R
require(TreeSearch)
require(ape)
require(phangorn)
require(ggplot2)
require(tidyr)
require(readxl)

# Select dataset
datasetName <- "THER"
rootDir <- "/Volumes/MasterData/dxsb43/Partitioning_Strategies/mutations/" #path
setwd(rootDir)


mrBayesTemplateFile <- paste0(rootDir, '/',datasetName,'/', datasetName, '_SS_TEMPLATE.nex')
dataset <- ReadAsPhyDat(mrBayesTemplateFile)


#how many trees to write to file
rep <- seq_len(300)

#functions
File <- function (suffix) paste0('Rest of files/CharStructure/', datasetName, suffix)  ##path and name of the target file
#######################################################################################################
############# make a .nex file with 300 random trees - only needed once per dataset ###################
#######################################################################################################


 #read optimal tree
 inputTree <- read.nexus(paste0('Rest of files/optimal_trees/',datasetName, '_optimal_tree.nex'))   ###optimal tree is in dir 'mutations', not in dataset subdir
 inputTree$edge.length <- NULL
 inputLabels <- inputTree$tip.label
 inputTree <- multi2di(inputTree)
 plot(inputTree)
 
            #write 300 trees to file, located in dir "mutations/CharStructure"
            write.nexus(structure(lapply(rep, function (i)
              ape::rtree(n = length(inputLabels), br=NULL, tip.label = inputLabels)),
              class='multiPhylo'), file=File(paste0('_random_',length(rep),'.nex')))
 



###################################################################################
###### calculate homoplasy for each character, 300 times ##########################
###################################################################################

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

############## CI function
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

## calculate vector of IC values - 1 for each character
tab[!tab %in% powerOf2] <- '?'
IC <- apply(tab, 2, CharacterInformation)  ##vector of IC values - information content of each character

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



CImat <- matrix(data=NA,nrow=length(rep), ncol=nChar) #preallocate a matrix of NAs, with as many rows as trees and as many columns as characters
colnames(CImat) <- chars
#which trees to calculate CI for
trees <- read.nexus(paste0('Rest of files/CharStructure/', datasetName, '_random_',length(rep),'.nex')) # reads all 300 random trees 

#loop through trees, calculating f for all characters for each tree i
for (i in seq_along(trees)) {
  tree <- trees[[i]]
  parsScore <- Fitch(tree, dataset)
  
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
se <- function(CharCol) {
  sd(CharCol)/sqrt(length(CharCol))
}
meanFse <- apply(CImat,MARGIN=2, FUN= se)

meanFdata <- data.frame(c(1:nChar), meanF, nstat)
colnames(meanFdata) <- c('character','meanF', 'nstates')       ##x$name[order(x$val)]
meanFdata$IC <- IC # add in column of IC values calculated above (under CharInfo function def)
meanFdata$se <- meanFse

meanFdata$character <- factor(meanFdata$character, levels = meanFdata$character[order(meanFdata$meanF)]) #in which order should the characters appear?

          #copy ordered character numbers to the clipboard
          #  writeClipboard(as.character(meanFdata$character[order(meanFdata$meanF)]))
                                                                                        # order() tells me how to get the numbers into ascending order

## this bit can plot all the f values against the character number. 
## it shows overplotting through transparency and the mean as a large Karo shape
          
          #reshape CImat into long format
          fMat <- as.data.frame(CImat)
          fMat <- gather(fMat, key="character", value="f", 1:nChar)
          fMat$character <- factor(fMat$character, levels = meanFdata$character[order(meanFdata$meanF)]) #in which order should the characters appear?
          #CImat$facet <- rep(c(1:4),times=50)

################## plot meanF for each character (size of marker prop to standard error of the mean) 
          ########    vs IC of the character
                  

          
# lm_eqn <- function(df, y, x){
#         formula = as.formula(sprintf('%s ~ %s', y, x))
#         m <- lm(formula, data=df);
#         # formating the values into a summary string to print out
#         # ~ give some space, but equal size and comma need to be quoted
#         #eq <- atop(substitute(italic(target) == a + b %.% italic(input)*,~~italic(r)^2~"="~r2*","~~p~"="~italic(pvalue), 
#         eq <- substitute(italic(r)^2~"="~r2*~"", 
#                          list(target = y,
#                               input = x,
#                               r2 = format(summary(m)$r.squared, digits = 3),
#                               # getting the pvalue is painful
#                               pvalue = format(summary(m)$coefficients[2,'Pr(>|t|)'], digits=1)
#                          )
#         )
#             as.character(as.expression(eq));                 
#           }
# lab <- lm_eqn(meanFdata, 'meanF','IC')
#  
m <-  lm(meanF ~ IC, meanFdata)
rsq <- round(summary(m)$r.squared,2)
pvalue <- format(summary(m)$coefficients[2,'Pr(>|t|)'], 1)

ggplot(meanFdata, aes(x=IC, y=meanF)) +
  geom_point(aes(size=se #nstates
                 ),
             alpha=0.3) +
  geom_smooth(method = 'lm') +
#  geom_text(x=47,y=3.2,label=paste0('r^2='rsq),color='red', size=6, parse=TRUE) +
  theme_light() +
  theme(legend.position='bottom',
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14)) +
  scale_y_continuous(name='Homoplasy (mean F)') +
  scale_x_continuous(name = 'Information Content [% of trees compatible]')


summary(lm(meanF ~ IC, meanFdata))
                 



###specify the size of the plot and path to it (should end up in dir CharStructure)
ggsave(filename = paste0('FvIC_', datasetName,'_se.pdf'),
       device = cairo_pdf, 
       path = paste0(rootDir, 'Rest of files/CharStructure'),
       width = 6,
       height = 5,
       units = 'in'
)

