#heatmap for average ranks of different partitioning strategies


require(readr)
require(tidyr)
require(plotly)
require(tidyverse)
setwd("~/Dropbox/MScR Thesis/Results/MLs_BFs/") #path to summary file
dataSets <- c('CEA', 'OZL', 'SCO','HYO','THER')


# initialise table of ranks of ps. columns are ps, rows are dataset-prior combos
ranks <- data.frame(
  dataset = factor(15),
  prior = factor(15),
  unpart = numeric(15),
  opt = numeric(15),
  neo = numeric(15),
  neo2 = numeric(15),
  ana = numeric(15),
  random = numeric(15),
  F = numeric(15),
  IC = numeric(15)
)

  set <- 'THER'
  dat <- read_delim(paste0(set,"_MLsBFs.csv"),';')
   dat_def <- dat[dat$brlenspr == 'default',]
   dat_exp <- dat[dat$brlenspr == 'exponential',]
   dat_gam <- dat[dat$brlenspr == 'gamma',]
   dat_clk <- dat[dat$brlenspr == 'clock',]
  
  dat <- transform(dat,  ### gives me dataframe with mean BFs etc for a dataset
            rank = ave(BF, brlenspr,  ### with ranks 1-8 within each brlenspr
                       FUN = function(x) rank(-x, ties.method = 'first')))
  
  


#enter rank columns by hand
ranks$dataset <- c('CEA','CEA','CEA','OZL','OZL','OZL','SCO','SCO','SCO','HYO','HYO','HYO','THER','THER','THER')
ranks$prior <- c('default','exp','clock',
                 'default','gamma','clock',
                 'default','exp','gamma',
                 'default','exp','clock',
                 'default','gamma','clock')
ranks$unpart <- c(4,3,4,2,2,2,5,1,4,3,1,3,4,4,4)
ranks$opt <- c(1,1,1,1,1,1,1,2,1,4,5,4,1,1,2)
ranks$neo <- c(5,8,5,5,5,6,8,3,8,1,8,2,5,5,5)
ranks$neo2 <- c(2,7,2,4,4,4,6,5,5,2,4,1,6,6,6)
ranks$ana <- c(8,6,8,7,8,8,7,6,6,6,2,6,7,7,7)
ranks$random <- c(7,5,7,8,7,7,2,4,2,7,7,7,8,8,8)
ranks$F <- c(6,4,6,3,3,3,3,7,3,5,3,5,2,2,1)
ranks$IC <- c(3,2,3,6,6,5,4,8,7,8,6,8,3,3,3)

#to only look at data from default priored analyses:
    #ranks <- ranks[which(ranks$prior == 'default'), ]

# initialise df of counts (in how many cases ps in column outperforms ps in row)
counts <- data.frame(
  unpart = numeric(8),
  opt = numeric(8),
  neo = numeric(8),
  neo2 = numeric(8),
  ana = numeric(8),
  random = numeric(8),
  F = numeric(8),
  IC = numeric(8)
)

rownames(counts) <- c('unpart','opt','neo','neo2','ana','random','F','IC')

#count number of times partstrat a is better than partstrat b
ps.rows <- rownames(counts)
ps.cols <- colnames(counts)

#initialising the vector at full length and assigning each element
is.higher <- logical(15)
for (a in ps.cols) {
  for (b in ps.rows) {
    for (row in 1:nrow(ranks)) {
      is.higher[row] <- ranks[row,a] < ranks[row,b] # TRUE if ps a is better than ps b
      counts[b,a] <- sum(is.higher)
    } 
  }  
}





counts <- rownames_to_column(counts, var = 'b')
# gather counts into format  b - a - #a<b
counts.long <- gather(counts, key = 'a', value = 'a.better.than.b', -1)
counts.long$a <- as.factor(counts.long$a)
counts.long$b <- as.factor(counts.long$b)

#heatmap
ggplot(counts.long, aes(b,a)) +
  geom_tile(aes(fill = a.better.than.b), 
            colour = 'white') +
  scale_fill_gradient2('a > b',midpoint = 7.5# low = 'red', mid = 'white', 
    #high = 'purple'
  ) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


##save plot
ggsave(
  filename = paste0('heatmap_partstrats_compared_allpriors.pdf'),
  device = cairo_pdf,
  path = '~/Dropbox/MScR Thesis/Results/Stats/',
  width = 6,
  height = 5,
  units = 'in'
)
