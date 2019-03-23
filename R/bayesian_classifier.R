library(readr)

##############
# Prepare data
##############
# read in data from CSV
assaysDF <- read_csv("/Users/jcasaletto/Desktop/GRAD_SCHOOL/UCSC/WINTER_2019/ROTATION/DATA/data-with-fisher.csv")
summary(assaysDF)

# filter out unnecessary columns
allFeaturesPlusLabel <- c( "in_vivo_ws", "in_vivo_wu", "in_vivo_ms", "in_vivo_mu",  "wt5score", "mu5score", "wt3score", "mu3score",  "lor", "exonlen", "deltaSS3","deltaSS5", "fivePrimeSSscore", "threePrimeSSscore", "n_ess", "n_ese", "esm", "fisher")
allData <- assaysDF[allFeaturesPlusLabel]

# remove NAs
allScores <- allData[complete.cases(allData),]
labeledAllFiniteScores <- allScores[is.finite(rowSums(allScores)),]


nrow(labeledAllFiniteScores[labeledAllFiniteScores$esm == 1,])
nrow(labeledAllFiniteScores[labeledAllFiniteScores$esm == 0,])


all <- c("fivePrimeSSscore", "threePrimeSSscore", "wt5score", "mu5score", "wt3score", "mu3score", "lor", "exonlen", "deltaSS5", "deltaSS3", "n_ese", "n_ess",  "esm")
summary(labeledAllFiniteScores)
#=============================
# Ok, here is how it goes: 
# 1. You sample, and in the two groups of Y=1 and Y=0, if X is continuous, you probably assume that distribution is 
# Gaussian. So for that, you can calculate the mean and standard deviation strictly from data points. 

# 2. P(Y=1 | X) = P(Y)N(mu, sigma) / P(X), and P(Y=0 | X) = P(Y)N(mu, sigma) / P(X), 

# 3. To classify, for the above two formulas, since they both have P(X) as the denominator, you can take them out, and just compare the numerators, and whichever has a larger probability is whichever class (Y value) you classify a new point into. 

# ============================================
# calculate mean and variance of each feature
# ============================================
all_features <- labeledAllFiniteScores[all]
all_featuresPos <- all_features[all_features$esm == 1,]
all_featuresNeg <- all_features[all_features$esm == 0,]

# calculate means for each feature for positive class
musPos <- as.data.frame(sapply(all_featuresPos, FUN=mean), index=all)
# calculate means for each feature for negative class
musNeg <- as.data.frame(sapply(all_featuresNeg, FUN=mean), index=all)
# calculate means for each feature for all
musAll <- as.data.frame(sapply(all_features, FUN=mean), index=all)

# calculate sd for each feature for positive class
sdPos  <- as.data.frame(sapply(all_featuresPos, FUN=sd), index=all)
# calculate sd for each feature for negative class
sdNeg  <- as.data.frame(sapply(all_featuresNeg, FUN=sd), index=all)
# calculate sd for each feature for all
sdAll <- as.data.frame(sapply(all_features, FUN=sd), index=all)
musAll
sdAll
sdPos
sdNeg
#===============================================
# build gaussian distribution for each feature
#===============================================

wt5scoreDistPos <- rnorm(1000, musPos["wt5score",], sdPos["wt5score",])
mu5scoreDistPos <- rnorm(1000, musPos["mu5score",], sdPos["mu5score",])
wt3scoreDistPos <- rnorm(1000, musPos["wt3score",], sdPos["wt3score",])
mu3scoreDistPos <- rnorm(1000, musPos["mu3score",], sdPos["mu3score",])
threeSSscoreDistPos <- rnorm(1000, musPos["threePrimeSSscore",], sdPos["threePrimeSSscore",])
fiveSSscoreDistPos <- rnorm(1000, musPos["fivePrimeSSscore",], sdPos["fivePrimeSSscore",])
lorDistPos <- rnorm(1000, musPos["lor",], sdPos["lor",])
exonDistPos <- rnorm(1000, musPos["exonlen",], sdPos["exonlen",])
deltaSS5DistPos <- rnorm(1000, musPos["deltaSS5",], sdPos["deltaSS5",])
deltaSS3DistPos <- rnorm(1000, musPos["deltaSS3",], sdPos["deltaSS3",])
ese_DistPos <- rnorm(1000, musPos["n_ese",], sdPos["n_ese",])
ess_DistPos <- rnorm(1000, musPos["n_ess",], sdPos["n_ess",])


wt5scoreDistNeg <- rnorm(1000, musNeg["wt5score",], sdNeg["wt5score",])
mu5scoreDistNeg <- rnorm(1000, musNeg["mu5score",], sdNeg["mu5score",])
wt3scoreDistNeg <- rnorm(1000, musNeg["wt3score",], sdNeg["wt3score",])
mu3scoreDistNeg <- rnorm(1000, musNeg["mu3score",], sdNeg["mu3score",])
threeSSscoreDistNeg <- rnorm(1000, musNeg["threePrimeSSscore",], sdNeg["threePrimeSSscore",])
fiveSSscoreDistNeg <- rnorm(1000, musNeg["fivePrimeSSscore",], sdNeg["fivePrimeSSscore",])
lorDistNeg <- rnorm(1000, musNeg["lor",], sdNeg["lor",])
exonDistNeg <- rnorm(1000, musNeg["exonlen",], sdNeg["exonlen",])
deltaSS5DistNeg <- rnorm(1000, musNeg["deltaSS5",], sdNeg["deltaSS5",])
deltaSS3DistNeg <- rnorm(1000, musNeg["deltaSS3",], sdNeg["deltaSS3",])
ese_DistNeg <- rnorm(1000, musNeg["n_ese",], sdNeg["n_ese",])
ess_DistNeg <- rnorm(1000, musNeg["n_ess",], sdNeg["n_ess",])

wt5scoreDistAll <- rnorm(1000, musAll["wt5score",], sdAll["wt5score",])
mu5scoreDistAll <- rnorm(1000, musAll["mu5score",], sdAll["mu5score",])
wt3scoreDistAll <- rnorm(1000, musAll["wt3score",], sdAll["wt3score",])
mu3scoreDistAll <- rnorm(1000, musAll["mu3score",], sdAll["mu3score",])
threeSSscoreDistAll <- rnorm(1000, musAll["threePrimeSSscore",], sdAll["threePrimeSSscore",])
fiveSSscoreDistAll <- rnorm(1000, musAll["fivePrimeSSscore",], sdAll["fivePrimeSSscore",])
lorDistAll <- rnorm(1000, musAll["lor",], sdAll["lor",])
exonDistAll <- rnorm(1000, musAll["exonlen",], sdAll["exonlen",])
deltaSS5DistAll <- rnorm(1000, musAll["deltaSS5",], sdAll["deltaSS5",])
deltaSS3DistAll <- rnorm(1000, musAll["deltaSS3",], sdAll["deltaSS3",])
ese_DistAll <- rnorm(1000, musAll["n_ese",], sdAll["n_ese",])
ess_DistAll <- rnorm(1000, musAll["n_ess",], sdAll["n_ess",])

#############################
# calculate probabilities for bayes
###############################
# priors
p_false = nrow(labeledAllFiniteScores[labeledAllFiniteScores$esm == 0,])/nrow(labeledAllFiniteScores)
p_false
p_true =  nrow(labeledAllFiniteScores[labeledAllFiniteScores$esm == 1,])/nrow(labeledAllFiniteScores)
p_true


ll_x <- function(w5,m5,w3,m3,lo,le,d5,d3,s3,s5,nese,ness,mus,sd) {

  wt5_ <-pnorm(w5,mus["wt5score",], sd["wt5score",])
  mu5_ <-pnorm(m5,mus["mu5score",], sd["mu5score",])
  wt3_ <-pnorm(w3,mus["wt3score",], sd["wt3score",])
  mu3_ <-pnorm(m3,mus["mu3score",], sd["mu3score",])
  ss3_ <-pnorm(s3,mus["threePrimeSSscore",], sd["threePrimeSSscore",])
  ss5_ <-pnorm(s5,mus["fivePrimeSSscore",], sd["fivePrimeSSscore",])
  lor_ <-pnorm(lo,mus["lor",],sd["lor",])
  exonlen_ <-pnorm(le,mus["exonlen",],sd["exonlen",])
  deltaSS5_ <-pnorm(d5,mus["deltaSS5",],sd["deltaSS5",])
  deltaSS3_ <-pnorm(d3,mus["deltaSS3",],sd["deltaSS3",])
  nese_ <- pnorm(nese,mus["n_ese",],sd["n_ese",])
  ness_ <- pnorm(ness,mus["n_ess",],sd["n_ess",])
  
  #llx <- log2(wt5_) +  log2(mu5_) + log2(wt3_)  + log2(mu3_)  + log2(lor_) + log2(exonlen_) + log2(deltaSS5_) + log2(deltaSS3_)  + log2(nese_) + log2(ness_)
  llx <- log2(deltaSS5_) + log2(deltaSS3_)
  return (llx)
}

bayes <- function(df, p_class, mus_class, sd_class, mus_all, sd_all) {
  
  w5 <- df["wt5score"]
  m5 <- df["mu5score"]
  w3 <- df["wt3score"]
  m3 <- df["mu3score"]
  s3 <- df["threePrimeSSscore"]
  s5 <- df["fivePrimeSSscore"]
  lo <- df["lor"]
  le <- df["exonlen"]
  d5 <- df["deltaSS5"]
  d3 <- df["deltaSS3"]
  nese <- df["n_ese"]
  ness <- df["n_ess"]
  #w5,m5,w3,m3,lo,le,d5,d3,diff5,diff3,s3,s5,nese,ness,mus,sd
  ll_data_given_class <- ll_x(w5,m5,w3,m3,lo,le,d5,d3,s3,s5,nese,ness,mus_class,sd_class)

  #ll_data <- ll_x(w5,m5,w3,m3,lo,le,d5,d3,mus_all,sd_all)
  #return (ll_data_given_class + log2(p_class) - ll_data)
  #return (ll_data_given_class - ll_data)
  #return (ll_data_given_class + log2(p_class))
  return (ll_data_given_class)
}


labeledAllFiniteScores$bayes_pos <- apply(labeledAllFiniteScores, 1, bayes, p_true, musPos, sdPos, musAll, sdAll)
labeledAllFiniteScores$bayes_neg <- apply(labeledAllFiniteScores, 1, bayes, p_false, musNeg, sdNeg, musAll, sdAll)


positiveBayes <- labeledAllFiniteScores[labeledAllFiniteScores$bayes_pos > labeledAllFiniteScores$bayes_neg,]
nrow(positiveBayes)
nrow(positiveBayes[positiveBayes$esm == 1,])


negativeBayes <- labeledAllFiniteScores[labeledAllFiniteScores$bayes_pos < labeledAllFiniteScores$bayes_neg,]
nrow(negativeBayes)
nrow(negativeBayes[negativeBayes$esm == 0,])

compare2 <- c("bayes_pos", "bayes_neg", "esm")
labeledAllFiniteScores[compare2]

# ================================
# plot dist of each feature
# ================================
labels <- labeledAllFiniteScores$label_vivo
plot(density(labels))

threePrimeScores <- labeledAllFiniteScores$threePrimeSSscore
plot(density(threePrimeScores))

fivePrimeScores <- labeledAllFiniteScores$fivePrimeSSscore
plot(density(fivePrimeScores))
summary(fivePrimeScores)

wt5scores <- labeledAllFiniteScores$wt5score
plot(density(wt5scores))

mu5scores <- labeledAllFiniteScores$mu5score
plot(density(mu5scores))

wt3scores <- labeledAllFiniteScores$wt3score
plot(density(wt3scores))

wt5scores <- labeledAllFiniteScores$wt5score
plot(density(wt5scores))

lorScores <- labeledAllFiniteScores$lor
plot(density(lorScores))

exonlenScores <- labeledAllFiniteScores$exonlen
plot(density(exonlenScores))

deltaSS3scores <- labeledAllFiniteScores$deltaSS3
plot(density(deltaSS3scores))

deltaSS5scores <- labeledAllFiniteScores$deltaSS5
plot(density(deltaSS5scores))

ese <- labeledAllFiniteScores$n_ese
plot(density(ese))

ess <- labeledAllFiniteScores$n_ess
plot(density(ess))




####################
# r doc example for log likelhood
####################

install.packages("likelihood")
library("likelihood")


