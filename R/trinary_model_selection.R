library(readr)
##################
# Prepare data
##################
# read in data from CSV
assaysDF <- read_csv("/Users/jcasaletto/Desktop/GRAD_SCHOOL/UCSC/WINTER_2019/ROTATION/DATA/data-with-skippy.csv")
# determine if any 0's in denominators (no splicing)
noSplicingAtAll <- subset(assaysDF, (in_vivo_ms + in_vivo_ws) ==0)
noMutatedSplicing <- subset(assaysDF, (in_vivo_ms == 0) & (in_vivo_ws !=0))
noWTSplicing <- subset(assaysDF, (in_vivo_ms !=0) & (in_vivo_ws ==0))
nrow(noSplicingAtAll)
nrow(noMutatedSplicing)
noMutatedSplicing
nrow(noWTSplicing)
summary(assaysDF)
assaysDF

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
#varsOfInterest <- c("wt5score", "mu5score", "wt3score", "mu3score", "l2fcSpliceRatio_vivo")
varsOfInterest <- c("wt5score", "mu5score", "wt3score", "mu3score", "l2fcSpliceRatio_vivo", "l2fcSpliceRatio_vitro", "lor", "exonlen", "deltaSS3","deltaSS5")
spliceScoresWithL2FC <- assaysDF[varsOfInterest]
# filter out rows with NA, Inf, or -Inf in any field
allScores <- spliceScoresWithL2FC[is.finite(rowSums(spliceScoresWithL2FC)),]
allFiniteScores <- allScores[complete.cases(allScores),]


# filter out scores with identical wt=mu 5' and 3'
allUnique5Scores <- allFiniteScores[allFiniteScores$wt5score != allFiniteScores$mu5score,]
allUniqueScores <- allUnique5Scores[allUnique5Scores$wt3score != allUnique5Scores$mu3score,]
allUniqueScores


allUniqueScores_5diff <- transform(allUniqueScores, diff_5 = abs(wt5score - mu5score))
allUniqueScores <- allUniqueScores_5diff <- transform(allUniqueScores_5diff, diff_3 = abs(wt3score - mu3score))

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
allFiniteScoresNormalized <- normalize(allUniqueScores)
# label each data point as 0 or 1 
labeledAllFiniteScores = within(allFiniteScoresNormalized, {
  label_vivo = ifelse(l2fcSpliceRatio_vivo > 0.1 , 1, ifelse(l2fcSpliceRatio_vivo < -0.1, -1, 0))
})

enhancedSubset <- subset(labeledAllFiniteScores, l2fcSpliceRatio_vivo > 0.1)
suppressedSubset <- subset(labeledAllFiniteScores, l2fcSpliceRatio_vivo < -0.1)
neutralSubset <- subset(labeledAllFiniteScores, l2fcSpliceRatio_vivo >= -0.1 & l2fcSpliceRatio_vivo <= 0.1)
nrow(enhancedSubset)
nrow(suppressedSubset)
nrow(neutralSubset)

# create train and test data sets from labeled data
require(caTools)
set.seed(123)
sample = sample.split(labeledAllFiniteScores$label_vivo, SplitRatio = 0.80)
trainingData <- subset(labeledAllFiniteScores, sample == 1)
testingData <- subset(labeledAllFiniteScores, sample == 0)


#####################
# Set y and x values for model
#####################
# build logistic regression model
y <- "label_vivo"

allScores <- c("wt5score", "mu5score", "wt3score", "mu3score",  "lor", "exonlen", "deltaSS3", "deltaSS5")

onlyMaxent <- c("wt5score", "mu5score", "wt3score", "mu3score")

onlySkippy <- c( "lor", "exonlen", "deltaSS3", "deltaSS5")

all3 <-  c("wt3score", "mu3score", "deltaSS3")

all5 <-  c("wt5score", "mu5score", "deltaSS5")

best <- c("deltaSS3", "exonlen")
all <- c("wt5score", "mu5score", "wt3score", "mu3score","lor", "exonlen", "deltaSS3", "deltaSS5", "diff_3", "diff_5")
bestCombo <- c("mu3score", "wt3score", "exonlen", "diff_5", "diff_3")





###############
# random forest
###############
library(randomForest)
set.seed(5123512)
x <- labeledAllFiniteScores[bestCombo]
y <- as.factor(labeledAllFiniteScores$label_vivo)
RFmodel <- randomForest(x=x, y=y, ntree=100, nodesize=7, importance=TRUE)
print(RFmodel)
varImp <- importance(RFmodel)
varImp
varImp[1:8,]
varImpPlot(RFmodel, type=1)
















