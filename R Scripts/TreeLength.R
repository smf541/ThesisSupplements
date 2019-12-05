#plots lnML of result tree of a homoplasy-partitioned Bayesian analysis
#   vs longest branch length of result tree 
#
#  to see if worse models infer longer branch lengths
################################################################################

#packages
require(ggplot2)
require(tidyr)
require(stringr)
require(readr)
# require(tibble)
# require(utils)
require(ape)
require(Quartet)
# require(gtools)
require(tidyverse)

#set dataset and perturbation move
dataSet <- 'HYO'         #SCO, THER, HYO

#set working directory 
rootDir <-"/Volumes/MasterData/dxsb43/Partitioning_Strategies/mutations" #path
setwd(paste0(rootDir, '/', dataSet))


#initialise tibble for logMLs
bfMat <- as_tibble(matrix(data=NA,nrow=length(list.files('StartingTrees',pattern='*.nex'))*100,ncol=13))
colnames(bfMat) <- c('file','longBranch', "lengthVar", 'BF1','BF2','BF3','BF4','BF5','BF6','BF7','BF8','BayesFactor','se')
bfMat$file <- str_sort(list.files('MrBayes', pattern = '.lstat'), numeric=TRUE)

#read in lstat file and populate amMat with ML values, indexing by the filename

for (resTreeFile in str_sort(list.files('MrBayes', pattern = '.con.tre'), numeric=TRUE)) {
  
  #resTreeFile <- 'SCO_SPR_chain.nex.4.nex.con.tre'       ## for testing
  resTree <- read.nexus(paste0('MrBayes/',resTreeFile))
  file <- gsub(pattern = ".con.tre", replacement = ".lstat", resTreeFile)
  
  outFile <- read_tsv(paste0('MrBayes/',file), comment = '[')   ## reads the data into a tibble of dim nrow=9, ncol=4
  
  #copy whole column of MLs into row of amMat
  bfMat[which(grepl(file, bfMat$file)),c(4:11)] <- t(outFile[c(1:8),'harmonic_mean'])  # harmonic mean column of tibble is extracted as a vector and transposed

  # grab length of longest branch and number of that branch for each result tree
  bfMat[which(grepl(file, bfMat$file)),2] <- max(unlist(resTree$edge.length))
  bfMat[which(grepl(file, bfMat$file)),3] <- var(unlist(resTree$edge.length))
}


#calculate arithmetic mean of MLs!
bfMat$BayesFactor <- apply(bfMat[,c(4:11)],1,mean)



############### Graph of BF vs tree length

cor(bfMat$BayesFactor, bfMat$longBranch)
cor(bfMat$BayesFactor, bfMat$lengthVar)
linmod <- lm(longBranch ~ BayesFactor, data=bfMat) # is the longest branch length corr. with BF?
# rsq <- summary(linmod)$r.squared
# Fval <- summary(linmod)$fstatistic[1]
# pval <- summary(linmod)$coefficients[2,4]
# str(summary(linmod))


ggplot(data = bfMat, aes(x=bfMat$BayesFactor, y=bfMat$longBranch)) +
  geom_point(alpha=0.5
            , aes(colour=lengthVar)
             ) +
  geom_smooth(method = 'lm') +
  theme_light()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.position='bottom',
        legend.key.width=unit(1.2,'cm'),
        legend.text=element_text(size=12),
        legend.title=element_text(size=15)) + 
  scale_y_continuous(name = 'Length of longest branch') +
  scale_x_continuous(name = 'log Marginal Likelihood') +
  scale_colour_continuous(name = 'Branch length\nvariance')




ggsave(filename = paste0('TreePerturb_longestBranch_', dataSet,'.pdf'),
        device = cairo_pdf, 
        path = paste0('~/Dropbox/MScR Thesis/Results/SimPlots'),
        width = 5,
        height = 5,
        units = 'in'
)


csvMat <- bfMat[,c(1,2,3,4,13)] # don't need all those columns
write_csv(csvMat, 
          path = paste0('~/Dropbox/MScR Thesis/Results/SimPlots/TreeLength_',dataSet , '.csv')
)


