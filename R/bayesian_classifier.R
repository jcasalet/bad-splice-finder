library(readr)

##############
# Prepare data
##############
# parse command-line args
#args <- commandArgs(TRUE)
#inputFile <- args[1]
# read in data from CSV
assaysDF <- read_csv("/Users/jcasaletto/Desktop/GRAD_SCHOOL/UCSC/WINTER_2019/ROTATION/DATA/data-with-skippy.csv")
# calculate ratio of in vivo mutated spliced/total (efficiency)
assaysDF$spliceRatio_mut_vivo <- assaysDF$in_vivo_ms / (assaysDF$in_vivo_ms + assaysDF$in_vivo_mu)
assaysDF$spliceRatio_wt_vivo <- assaysDF$in_vivo_ws / (assaysDF$in_vivo_ws + assaysDF$in_vivo_wu)
# calculate log fold change in vivo of ratio of spliced/total mutant:wild type
assaysDF$l2fcSpliceRatio_vivo <- log(assaysDF$spliceRatio_mut_vivo/assaysDF$spliceRatio_wt_vivo)
# calculate ratio of in vitro mutated spliced/total (efficiency)
assaysDF$spliceRatio_mut_vitro <- assaysDF$vit_ms / (assaysDF$vit_ms + assaysDF$vit_mu)
assaysDF$spliceRatio_wt_vitro <- assaysDF$vit_ws / (assaysDF$vit_ws + assaysDF$vit_wu)
# calculate log fold change in vitro of ratio of spliced/total mutant:wild type
assaysDF$l2fcSpliceRatio_vitro <- log(assaysDF$spliceRatio_mut_vitro/assaysDF$spliceRatio_wt_vitro)
# filter out unnecessary columns
#varsOfInterest <- c("wt5score", "mu5score", "wt3score", "mu3score", "l2fcSpliceRatio_vivo", "l2fcSpliceRatio_vitro")
varsOfInterest <- c("wt5score", "mu5score", "wt3score", "mu3score", "l2fcSpliceRatio_vivo", "l2fcSpliceRatio_vitro", "lor", "exonlen", "deltaSS3","deltaSS5")
spliceScoresWithL2FC <- assaysDF[varsOfInterest]
# filter out rows with NA, Inf, or -Inf in any field
allScores <- spliceScoresWithL2FC[complete.cases(spliceScoresWithL2FC),]
allFiniteScores <- allScores[is.finite(rowSums(allScores)),]


allFiniteScores_5diff <- transform(allFiniteScores, diff_5 = abs(wt5score - mu5score))
allFiniteScores <- allFiniteScores_5diff <- transform(allFiniteScores_5diff, diff_3 = abs(wt3score - mu3score))



# normalize data
normalize <- function(df) {
  for(col in colnames(df)) {
    if (deparse(substitute(col)) == "\"l2fcSpliceRatio_vitro\"" | deparse(substitute(col)) == "\"l2fcSpliceRatio_vivo\"") {
      next
    }
    numerator <- df[[col]] - min(df[[col]])
    denominator <- max(df[[col]] - min(df[[col]]))
    df[[col]] <- numerator / denominator
  }
  return (df)
}

#allFiniteScoresNormalized <- normalize(allFiniteScores)
# label each data point as 0 or 1 
labeledAllFiniteScores = within(allFiniteScores, {
  label_vivo = ifelse(l2fcSpliceRatio_vivo > 0.1 | l2fcSpliceRatio_vivo < -0.1, 1, 0)
  label_vitro = ifelse(l2fcSpliceRatio_vitro > 0.1 | l2fcSpliceRatio_vitro < -0.1, 1, 0)
})
nrow(labeledAllFiniteScores[labeledAllFiniteScores$label_vivo == 1,])
nrow(labeledAllFiniteScores[labeledAllFiniteScores$label_vivo == 0,])


all <- c("label_vivo", "wt5score", "mu5score", "wt3score", "mu3score","lor", "exonlen", "deltaSS3", "deltaSS5", "diff_3", "diff_5")


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

all_featuresPos <- all_features[all_features$label_vivo == 1,]
all_featuresPos
all_featuresNeg <- all_features[all_features$label_vivo == 0,]
all_featuresNeg

# calculate means for each feature for positive class
musPos <- as.data.frame(sapply(all_featuresPos, FUN=mean), index=all)
# calculate means for each feature for negative class
musNeg <- as.data.frame(sapply(all_featuresNeg, FUN=mean), index=all)

musPos
# calculate sd for each feature for positive class
sdPos  <- as.data.frame(sapply(all_featuresPos, FUN=sd), index=all)
# calculate sd for each feature for negative class
sdNeg  <- as.data.frame(sapply(all_featuresNeg, FUN=sd), index=all)

#===============================================
# build gaussian distribution for each feature
#===============================================

wt5scoreDistPos <- rnorm(1000, musPos["wt5score",], sdPos["wt5score",])
mu5scoreDistPos <- rnorm(1000, musPos["mu5score",], sdPos["mu5score",])
wt3scoreDistPos <- rnorm(1000, musPos["wt3score",], sdPos["wt3score",])
mu3scoreDistPos <- rnorm(1000, musPos["mu3score",], sdPos["mu3score",])
lorDistPos <- rnorm(1000, musPos["lor",], sdPos["lor",])
exonDistPos <- rnorm(1000, musPos["exonlen",], sdPos["exonlen",])
deltaSS5DistPos <- rnorm(1000, musPos["deltaSS5",], sdPos["deltaSS5",])
deltaSS3DistPos <- rnorm(1000, musPos["deltaSS3",], sdPos["deltaSS3",])
diff_5DistPos <- rnorm(1000, musPos["diff_5",], sdPos["diff_5",])
diff_3DistPos <- rnorm(1000, musPos["diff_3",], sdPos["diff_3",])

wt5scoreDistNeg <- rnorm(1000, musNeg["wt5score",], sdNeg["wt5score",])
mu5scoreDistNeg <- rnorm(1000, musNeg["mu5score",], sdNeg["mu5score",])
wt3scoreDistNeg <- rnorm(1000, musNeg["wt3score",], sdNeg["wt3score",])
mu3scoreDistNeg <- rnorm(1000, musNeg["mu3score",], sdNeg["mu3score",])
lorDistNeg <- rnorm(1000, musNeg["lor",], sdNeg["lor",])
exonDistNeg <- rnorm(1000, musNeg["exonlen",], sdNeg["exonlen",])
deltaSS5DistNeg <- rnorm(1000, musNeg["deltaSS5",], sdNeg["deltaSS5",])
deltaSS3DistNeg <- rnorm(1000, musNeg["deltaSS3",], sdNeg["deltaSS3",])
diff_5DistNeg <- rnorm(1000, musNeg["diff_5",], sdNeg["diff_5",])
diff_3DistNeg <- rnorm(1000, musNeg["diff_3",], sdNeg["diff_3",])

#############################
# calculate probability of a few data points
###############################
# priors
p_false = nrow(labeledAllFiniteScores[labeledAllFiniteScores$label_vivo == 0,])/nrow(labeledAllFiniteScores)
p_false
p_true =  nrow(labeledAllFiniteScores[labeledAllFiniteScores$label_vivo == 1,])/nrow(labeledAllFiniteScores)
p_true



p_x <- function(p, w5,m5,w3,m3,lo,le,d5,d3,diff5,diff3,mus,sd) {

  wt5_ <-pnorm(w5,mus["wt5score",], sd["wt5score",])
  mu5_ <-pnorm(m5,mus["mu5score",], sd["mu5score",])
  wt3_ <-pnorm(w3,mus["wt3score",], sd["wt3score",])
  mu3_ <-pnorm(m3,mus["mu3score",], sd["mu3score",])
  lor_ <-pnorm(lo,mus["lor",],sd["lor",])
  exonlen_ <-pnorm(le,mus["exonlen",],sd["exonlen",])
  deltaSS5_ <-pnorm(d5,mus["deltaSS5",],sd["deltaSS5",])
  deltaSS3_ <-pnorm(d3,mus["deltaSS3",],sd["deltaSS3",])
  diff_5_ <- pnorm(diff5,mus["diff_5",],sd["diff_5",])
  diff_3_ <- pnorm(diff3,mus["diff_3",],sd["diff_3",])
  
  px <- p * wt5_ * mu5_ * wt3_ * mu3_ * lor_ * exonlen_ * deltaSS5_ * deltaSS3_ * diff_5_ * diff_3_

  return (px)
}

all_px <- function(df, p, mus, sd) {
  
  w5 <- df["wt5score"]
  m5 <- df["mu5score"]
  w3 <- df["wt3score"]
  m3 <- df["mu3score"]
  lo <- df["lor"]
  le <- df["exonlen"]
  d5 <- df["deltaSS5"]
  d3 <- df["deltaSS3"]
  diff5 <- df["diff_5"]
  diff3 <- df["diff_3"]
    
  x <- p_x(p, w5,m5,w3,m3,lo,le,d5,d3,diff5,diff3,mus,sd)
  return (x)
}

for(i in 1:length(labeledAllFiniteScores)) {
  myrow <- labeledAllFiniteScores[i,]
  
}


labeledAllFiniteScores$bayes_pos <- apply(labeledAllFiniteScores, 1, all_px, p_true, musPos, sdPos)
labeledAllFiniteScores$bayes_neg <- apply(labeledAllFiniteScores, 1, all_px, p_false, musNeg, sdNeg)

# now add another label for ratio of pos/neg
labeledAllFiniteScoresBayes = within(labeledAllFiniteScores, {
  label_bayes = ifelse(bayes_pos/bayes_neg > 1, 1, 0)
})


# calculate accuracy
ncorrect <- nrow(labeledAllFiniteScoresBayes[labeledAllFiniteScoresBayes$label_bayes == labeledAllFiniteScoresBayes$label_vivo,])
nincorrect <- nrow(labeledAllFiniteScoresBayes[labeledAllFiniteScoresBayes$label_bayes != labeledAllFiniteScoresBayes$label_vivo,])
ntotal <- nrow(labeledAllFiniteScoresBayes)
accuracy <- ncorrect / ntotal
accuracy


# calculate TPR, TNR, FPR, and FNR
TP <- length(which(labeledAllFiniteScoresBayes$label_bayes == 1 & labeledAllFiniteScoresBayes$label_vivo == 1))
TN <- length(which(labeledAllFiniteScoresBayes$label_bayes == 0 & labeledAllFiniteScoresBayes$label_vivo == 0))
FP<- length(which(labeledAllFiniteScoresBayes$label_bayes == 1 & labeledAllFiniteScoresBayes$label_vivo == 0))
FN<- length(which(labeledAllFiniteScoresBayes$label_bayes == 0 & labeledAllFiniteScoresBayes$label_vivo == 1))

TPR <- TP / (FN + TP)
TPR

FNR <- FN / (FN + TP)
FNR

TNR <- TN / (FP + TN)
TNR

FPR <- FP / (FP + TN)
FPR

#####################
# TLDR;
######################



# ================================
# plot dist of each feature
# ================================

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
l2fcvivoScores <- labeledAllFiniteScores$l2fcSpliceRatio_vivo
plot(density(l2fcvivoScores))
l2fcvitroScores <- labeledAllFiniteScores$l2fcSpliceRatio_vitro
plot(density(l2fcvitroScores))



