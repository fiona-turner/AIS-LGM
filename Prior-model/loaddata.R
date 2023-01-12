loaddata <- function(dir){
  setwd(dir)
  filenames <- list.files(".", '*.csv')
  ##reorder to help with comparison of analysis
  filenames <- c( 'height5g.csv', 'height6g.csv', 'heightpip.csv', 'heighttar.csv', 'topopollardhadcm3.csv', 'usurfgolledge2012hadcm3.csv',
                  'usurfgolledge2013hadcm3.csv', 'usurfgolledge2014hadcm3.csv', 'topoboer1hadcm3.csv', 'topoboer2hadcm3.csv', 'topoboer3hadcm3.csv', 
                  'topoboer4hadcm3.csv', 'topoboer5hadcm3.csv', 'topoboer6hadcm3.csv', 'topoboer7hadcm3.csv', 'topoboer8hadcm3.csv', 'topoboer9hadcm3.csv',
                  'topoboer10hadcm3.csv', 'topoboer11hadcm3.csv', 'topoboer12hadcm3.csv', 'topoboer13hadcm3.csv', 'topoboer14hadcm3.csv',
                  'topoboer15hadcm3.csv', 'topoboer16hadcm3.csv', 'topoboer17hadcm3.csv', 'topoboer18hadcm3.csv', 'topoboer19hadcm3.csv', 
                  'topoboer20hadcm3.csv', 'topoboer21hadcm3.csv', 'topoboer22hadcm3.csv', 'topoboer23hadcm3.csv', 'topoboer24hadcm3.csv', 
                  'topoboer25hadcm3.csv', 'topoboer26hadcm3.csv', 'topoboer27hadcm3.csv', 'topoboer28hadcm3.csv', 'topoboer29hadcm3.csv', 
                  'topoboer30hadcm3.csv', 'topoboer31hadcm3.csv', 'topoboer32hadcm3.csv')
  ##load the datasets, removing the first column as it is just an index
  X <- sapply(filenames, function(x) {tmp <-read.csv(x)
  return (unlist(tmp[,-1]))
  })
  #change to the transpose, so PCA is applied to the right dimension
  X <- t(X)
  #find the column means, and remove it so the data is centred
  mu = colMeans(X)
  Xc = sweep(X,2, mu)
  #split the data in to the two sets
  Xfourc <- Xc[1:4,]
  Xotherc <- Xc[5:40,]
  return(list(X,mu,Xc,Xfourc,Xotherc))
}