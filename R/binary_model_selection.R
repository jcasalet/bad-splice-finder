library(readr)
##############
# Prepare data
##############
# read in data from CSV
assaysDF <- read_csv("/Users/jcasaletto/Desktop/GRAD_SCHOOL/UCSC/WINTER_2019/ROTATION/DATA/data-with-fisher.csv")
nrow(assaysDF)

# try to calculate soemedi ratio from figure 1
assaysDF$a_over_c <- assaysDF$in_vivo_wu / assaysDF$in_vivo_mu
assaysDF$b_over_d <- assaysDF$in_vivo_ws / assaysDF$in_vivo_ms
assaysDF$soemedi_ratio <- assaysDF$a_over_c / assaysDF$b_over_d
assaysDF$l2_soemedi_ratio <- log(assaysDF$soemedi_ratio)

# try to calculate soemedi ratio from allelic imbalance analses (this should be it, but what is "input" and "output"?)
assaysDF$mo <- assaysDF$in_vivo_ms
assaysDF$mi <- assaysDF$in_vivo_mu
assaysDF$wo <- assaysDF$in_vivo_ws
assaysDF$wi <- assaysDF$in_vivo_wu
assaysDF$allelic_ratio <- ((assaysDF$mo / assaysDF$mi) / (assaysDF$wo / assaysDF$wi))
assaysDF$l2_allelic_ratio <- log(assaysDF$allelic_ratio)
assaysDF$l2_allelic_ratio

# calculate ratio of in vivo mutated spliced/total (efficiency)
assaysDF$spliceRatio_mut_vivo <- assaysDF$in_vivo_ms / (assaysDF$in_vivo_ms + assaysDF$in_vivo_mu)
assaysDF$spliceRatio_wt_vivo <- assaysDF$in_vivo_ws / (assaysDF$in_vivo_ws + assaysDF$in_vivo_wu)
# calculate log fold change in vivo of ratio of spliced/total mutant:wild type
assaysDF$l2fcSpliceRatio_vivo <- log(assaysDF$spliceRatio_mut_vivo/assaysDF$spliceRatio_wt_vivo)


# filter out unnecessary columns
allFeaturesPlusLabel <- c("l2_allelic_ratio", "l2_soemedi_ratio", "in_vivo_ws", "in_vivo_wu", "in_vivo_ms", "in_vivo_mu", "spliceRatio_mut_vivo", "spliceRatio_wt_vivo", "l2fcSpliceRatio_vivo", "wt5score", "mu5score", "wt3score", "mu3score",  "lor", "exonlen", "deltaSS3","deltaSS5", "fivePrimeSSscore", "threePrimeSSscore", "n_ess", "n_ese", "esm", "fisher")
allData <- assaysDF[allFeaturesPlusLabel]

# adjust p-value for FDR
allData$pAdj <- p.adjust(allData$fisher, method="fdr")

# remove NAs
allScores <- allData[complete.cases(allData),]
allFiniteScores <- allScores[is.finite(rowSums(allScores)),]

nrow(allFiniteScores)

# normalize data
#normalize <- function(df) {
#  for(col in colnames(df)) {
#    if (deparse(substitute(col)) == "\"l2fcSpliceRatio_vitro\"" | deparse(substitute(col)) == "\"l2fcSpliceRatio_vivo\"") {
#      next
#    }
#    numerator <- df[[col]] - min(df[[col]])
#    denominator <- max(df[[col]] - min(df[[col]]))
#    df[[col]] <- numerator / denominator
#  }
#  return (df)
#}


#allFiniteScoresNormalized <- normalize(allUniqueScores)
# label each data point as 0 or 1 
#labeledAllFiniteScores = within(allFiniteScoresNormalized, {
#  label_vivo = ifelse(l2fcSpliceRatio_vivo > 0.5 | l2fcSpliceRatio_vivo < -0.5, TRUE, FALSE)
#  label_vitro = ifelse(l2fcSpliceRatio_vitro > 0.5 | l2fcSpliceRatio_vitro < -0.5, TRUE, FALSE)
#})
labeledAllFiniteScores = within(allFiniteScores, {
  label_allele = ifelse(pAdj < 0.05 & (l2_allelic_ratio > 1.5 | l2_allelic_ratio < -1.5), TRUE, FALSE)
  label_soemedi = ifelse(pAdj < 0.05 & (l2_soemedi_ratio > 1.5 | l2_soemedi_ratio < -1.5), TRUE, FALSE)
  label_fisher = ifelse(pAdj < 0.05 & (l2fcSpliceRatio_vivo > 1.5 | l2fcSpliceRatio_vivo < -1.5), TRUE, FALSE )
  label_esm = ifelse(pAdj < 0.05 & esm == 1, TRUE, FALSE)
  
  #rfv = fisher.test(rbind(c(in_vivo_wu, in_vivo_ws), c(in_vivo_mu, in_vivo_ms)))
})

head(allData)
fisher.test(rbind(c(5515, 385), c(780, 10)))

compare3 <- c("esm", "l2_allelic_ratio", "pAdj")
compareData <- allData[compare3]
compareData[compareData$pAdj > 0.05 & compareData$esm == 1,]


nrow(labeledAllFiniteScores[labeledAllFiniteScores$label_esm == T,])
nrow(labeledAllFiniteScores[labeledAllFiniteScores$label_fisher == T,])
nrow(labeledAllFiniteScores[labeledAllFiniteScores$label_soemedi == T,])
nrow(labeledAllFiniteScores[labeledAllFiniteScores$label_allele == T,])

nrow(labeledAllFiniteScores[labeledAllFiniteScores$label_esm == F,])
nrow(labeledAllFiniteScores[labeledAllFiniteScores$label_fisher == F,])
nrow(labeledAllFiniteScores[labeledAllFiniteScores$label_soemedi == F,])
nrow(labeledAllFiniteScores[labeledAllFiniteScores$label_allele == F,])

compareLabels <- c("esm", "label_esm", "label_nolog", "label_allele", "label_fisher", "label_soemedi", "l2fcSpliceRatio_vivo", "l2_soemedi_ratio", "pAdj")
comp <- labeledAllFiniteScores[compareLabels]
comp[comp$label_esm != comp$label_allele,]

# create train and test data sets from labeled data
require(caTools)
set.seed(123)
allFeatures <- c("wt5score", "mu5score", "wt3score", "mu3score",  "lor", "exonlen", "deltaSS3","deltaSS5", "fivePrimeSSscore", "threePrimeSSscore", "n_ese", "n_ess")

sample = sample.split(labeledAllFiniteScores$label_esm, SplitRatio = 0.80)
trainingData <- subset(labeledAllFiniteScores, sample == TRUE)
testingData <- subset(labeledAllFiniteScores, sample == FALSE)


y <- "label_fisher"
all <- c("wt5score", "mu5score", "wt3score", "mu3score", "exonlen", "deltaSS3", "deltaSS5","n_ese", "n_ess")
justESX <- c("n_ese", "n_ess")
justDeltas <- c("deltaSS5", "deltaSS3")
deltaPlusESX <- c(justDeltas, justESX)
maxentScan <- c("wt5score", "wt3score", "mu5score", "mu3score")
bestLR <- c("exonlen", "deltaSS5", "deltaSS3")
bestRF <- c("exonlen", "deltaSS5", "deltaSS3", "n_ess", "n_ese")



#####################
# logistic regression
#####################
formulaString <- paste(y, paste(bestLR, collapse="+"), sep="~")
LRmodel <- glm(formulaString, data=trainingData, family=binomial(link="logit"))
# make predictions using model
trainingData$prediction <- predict(LRmodel, newdata=trainingData, type="response")
testingData$prediction <- predict(LRmodel, newdata=testingData, type="response")
labeledAllFiniteScores$prediction <- predict(LRmodel, newdata=labeledAllFiniteScores, type="response")
# plot distribution of prediction score grouped by known outcome
library(ggplot2)
ggplot(labeledAllFiniteScores, aes(x=prediction, color=label_esm, linetype=label_esm)) + geom_density()
# build confusion matrix
confusion.matrix <- table(pred=labeledAllFiniteScores$prediction>0.20, label=labeledAllFiniteScores$label_esm)
confusion.matrix
# calculate TPR, FPR, TNR, and FNAR and enrichment
TP <- confusion.matrix[2,2]
TN <- confusion.matrix[1,1]
P <- sum(confusion.matrix[2,])
N <- sum(confusion.matrix[1,])
TPR <- TP / P
TNR <- TN / N
FPR <- 1 - TNR
FNR <- 1 - TPR
TPR
FPR
TNR
FNR
accuracy <- (TP + TN) / (P + N)
accuracy
AIC(LRmodel)
# examine coefficients of model
coefficients(LRmodel)
# examine LRmodel summary
summary(LRmodel)
# plot ROC curve
library(pROC)
library(ROCR)
pROC::roc(trainingData$label_esm, trainingData$prediction)
plot(pROC::roc(trainingData$label_esm, trainingData$prediction))
ROCR::plot

###############
# random forest
###############
library(randomForest)
#set.seed(5123512)
x <- trainingData[bestRF]
y <- as.factor(trainingData$label_esm)
xTest <- testingData[bestRF]
yTest <- as.factor(testingData$label_esm)

# RF balanced?
#RFmodel <- randomForest(x=x, y=y, ntree=1000, importance=TRUE, xtest=xTest, ytest=yTest)
RFmodel <- randomForest(x=x, y=y, ntree=1000, importance=TRUE)

print(RFmodel)
varImp <- importance(RFmodel)
varImp
varImp[1:8,]
varImpPlot(RFmodel, type=1)
RFmodel$confusion
RFmodel$test$confusion
mean(RFmodel$test$err.rate)
# get AUC
class(RFmodel)
rf_predict_train <- predict(RFmodel, type="prob")[,2]
rf_prediction_train <- prediction(rf_predict_train, trainingData$label_esm)
rf_auc_train <- performance(rf_prediction_train, measure="auc")@y.values[[1]]
rf_auc_train

rf_predict_test <- predict(RFmodel, newdata=xTest, type="prob")
rf_prediction_test <- prediction(rf_predict_test[,2], testingData$label_esm)
rf_auc_test <- performance(rf_prediction_test, measure="auc")@y.values[[1]]
rf_auc_test





