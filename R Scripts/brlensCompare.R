#calculate MLs of unpart analyses for all datasets under all branch length priors

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

#set working directory 
rootDir <-"/Volumes/MasterData/dxsb43/Partitioning_Strategies/mutations" #path
setwd(rootDir)#  dir is e.g. .../mutations/SCO
partStrats <- c('unpart', 'optTree', 'F', 'IC', 'anatomy', 'neotrans', 'neotrans2','random')
for (PS in partStrats) {
MLs <- matrix(data=NA,nrow=4, ncol=6)
colnames(MLs) <- c('brlenspr', 'CEA','OZL','SCO','HYO','THER')
rownames <- c('', 'exp', 'gamma', 'clock')
MLs <- as_tibble(MLs)
MLs$brlenspr <- c('', 'exp', 'gamma', 'clock')

#list of .ss files from viable analyses (15 in total, 3 per dataset)
filesPrep <- list.files(pattern=paste0('_', PS, '_ss.nex.ss'), recursive=TRUE, include.dirs=FALSE)
files <- basename(filesPrep[ !grepl("Rest of files", filesPrep)])

rows <- MLs$brlenspr
cols <- colnames(MLs[,2:6])
for (column in cols) {
  for (row in rows) {
    if (!file.exists(paste0(column,row))) {MLs[which(grepl(row, MLs$brlenspr)),which(grepl(column,colnames(MLs)))] <- NA} else
        #read in unpartitioned tree .ss file - use this for calculating Bayes factor
        {mat <- read_tsv(paste0(column,row,'/',column, row, '_', PS, '_ss.nex.ss'),comment = '[')
        mat <- mat[,c(3:10)]
        #sum run columns over all steps, DISCARD POSITIVE RUNS, then take the mean
        ML <- mean(colSums(mat)[which(colSums(mat) < 0)])

        #place value into MLs tibble
        MLs[which(grepl(row, MLs$brlenspr)),which(grepl(column,colnames(MLs)))] <- ML}
    }
  }
BFs <- MLs
BFs[,2:6] <- MLs[,2:6] - c(MLs[1,2:6])
BFs[1,1] <- 'default'
#restructure data into long format
BFs <- gather(BFs, key='dataset', value='BF', c(2:6), factor_key=TRUE)
BFs$brlenspr <- as.factor(BFs$brlenspr)


## bumps chart:
require(ggrepel)
ggplot(data = BFs, aes(x=dataset, y=BF, group=brlenspr)) +
  geom_line(aes(colour = brlenspr),alpha=0.5, size=1.5) +
  geom_point(aes(colour = brlenspr), 
             size=2, 
             alpha=0.6)+
  scale_size(guide =FALSE) +
  scale_alpha(guide=FALSE)+
  geom_text_repel(data = BFs %>% filter(dataset == "CEA"),
                  aes(label=brlenspr, colour=brlenspr),
                  direction="both",
                  nudge_x=-0.7
                  )+
  geom_text_repel(data = BFs %>% filter(dataset == "THER"),
                  aes(label=brlenspr, colour=brlenspr),
                  direction="y",
                  nudge_x=2)+
  scale_y_continuous(name = "Bayes Factor"
                     #, breaks = seq(floor_any(min(allMLs$BF),10),max(allMLs$BF),10)) + ####limits = c(min(0,min(bfMat$BayesFactor)),max(bfMat$BayesFactor)
                    )+
  scale_x_discrete(name = "Dataset",
                   limits = c("","CEA", "HYO","OZL","SCO","THER","")) +
  theme_light() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()) +
  scale_colour_discrete(guide=FALSE
                        # name='Branch length prior',
                        # breaks = c('default','gamma','exp','clock'),
                        # labels=c('Default','Gamma','Exponential','Clock')
                        )


ggsave(filename = paste0('branch_lengths_withlines_',PS,'.pdf'),
       device = cairo_pdf,
       path = paste0(rootDir, '/Rest of files/BayesFactorPlots'),
       width = 4,
       height = 5,
       units = 'in'
       )
}

