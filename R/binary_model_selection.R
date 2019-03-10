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
  label_vivo = ifelse(l2fcSpliceRatio_vivo > 0.1 | l2fcSpliceRatio_vivo < -0.1, TRUE, FALSE)
  label_vitro = ifelse(l2fcSpliceRatio_vitro > 0.1 | l2fcSpliceRatio_vitro < -0.1, TRUE, FALSE)
})
nrow(labeledAllFiniteScores[labeledAllFiniteScores$label_vivo == 1,])
nrow(labeledAllFiniteScores[labeledAllFiniteScores$label_vivo == 0,])
# create train and test data sets from labeled data
require(caTools)
set.seed(123)
sample = sample.split(labeledAllFiniteScores$label_vivo, SplitRatio = 0.80)
trainingData <- subset(labeledAllFiniteScores, sample == 1)
testingData <- subset(labeledAllFiniteScores, sample == 0)
y <- "label_vivo"
all <- c("wt5score", "mu5score", "wt3score", "mu3score","lor", "exonlen", "deltaSS3", "deltaSS5", "diff_3", "diff_5")
bestCombo <- c("mu3score", "wt3score", "exonlen", "diff_5", "diff_3")
onlyMaxent <-  c("wt5score", "mu5score", "wt3score", "mu3score")
onlySkippy <-  c("lor", "exonlen", "deltaSS3", "deltaSS5")



#####################
# logistic regression
#####################
formulaString <- paste(y, paste(bestCombo, collapse="+"), sep="~")
LRmodel <- glm(formulaString, data=trainingData, family=binomial(link="logit"))
# make predictions using model
trainingData$prediction <- predict(LRmodel, newdata=trainingData, type="response")
testingData$prediction <- predict(LRmodel, newdata=testingData, type="response")
labeledAllFiniteScores$prediction <- predict(LRmodel, newdata=labeledAllFiniteScores, type="response")
# plot distribution of prediction score grouped by known outcome
library(ggplot2)
ggplot(trainingData, aes(x=prediction, color=label_vivo, linetype=label_vivo)) + geom_density()
# build confusion matrix
confusion.matrix <- table(pred=labeledAllFiniteScores$prediction>0.15, label=labeledAllFiniteScores$label_vivo)
confusion.matrix
#ctab.test <- table(pred=labeledAllFiniteScores$prediction>0.5, label=labeledAllFiniteScores$label_vivo)
#ctab.test
# calculate TPR, FPR, TNR, and FNAR and enrichment
TPR <- confusion.matrix[2,2]/sum(confusion.matrix[2,])
TNR <- confusion.matrix[1,1]/sum(confusion.matrix[1,])
FPR <- 1 - TNR
FNR <- 1 - TPR
TPR
FPR
TNR
FNR
AIC(LRmodel)
# examine coefficients of model
coefficients(LRmodel)
# examine LRmodel summary
summary(LRmodel)


###############
# random forest
###############
library(randomForest)
set.seed(5123512)
x <- labeledAllFiniteScores[bestCombo]
y <- as.factor(labeledAllFiniteScores$label_vivo)
# RF balanced?
RFmodel <- randomForest(x=x, y=y, ntree=100, nodesize=7, importance=TRUE)
print(RFmodel)
varImp <- importance(RFmodel)
varImp
varImp[1:8,]
varImpPlot(RFmodel, type=1)
RFmodel$confusion
RFmodel$classes









