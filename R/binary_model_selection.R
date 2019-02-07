library(readr)

# read in data from CSV
#assaysDF <- read_csv("/Users/jcasaletto/Desktop/GRAD_SCHOOL/UCSC/WINTER_2019/ROTATION/DATA/final-data-with-best-skippy-features.csv")
assaysDF <- read_csv("/Users/jcasaletto/Desktop/GRAD_SCHOOL/UCSC/WINTER_2019/ROTATION/DATA/final-data.csv")

names(assaysDF)
summary(assaysDF)

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
assaysDF$l2fcSpliceRatio_vitro
# filter out unnecessary columns
varsOfInterest <- c("wt5score", "mu5score", "wt3score", "mu3score", "l2fcSpliceRatio_vivo", "l2fcSpliceRatio_vitro")
#varsOfInterest <- c("wt5score", "mu5score", "wt3score", "mu3score", "l2fcSpliceRatio_vivo", "l2fcSpliceRatio_vitro", "logoddsratio", "exonlen", "deltaSS3","deltaSS5")

spliceScoresWithL2FC <- assaysDF[varsOfInterest]

summary(spliceScoresWithL2FC)

# filter out rows with NA, Inf, or -Inf in any field
allScores <- spliceScoresWithL2FC[complete.cases(spliceScoresWithL2FC),]
allFiniteScores <- allScores[is.finite(rowSums(allScores)),]
summary(allFiniteScores)

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

allFiniteScoresNormalized <- normalize(allFiniteScores)


summary(allFiniteScoresNormalized)

# label each data point as 0 or 1 
labeledAllFiniteScores = within(allFiniteScoresNormalized, {
  label_vivo = ifelse(l2fcSpliceRatio_vivo > 1.5 | l2fcSpliceRatio_vivo < -1.5, TRUE, FALSE)
  label_vitro = ifelse(l2fcSpliceRatio_vitro > 1.5 | l2fcSpliceRatio_vitro < -1.5, TRUE, FALSE)
})

nrow(labeledAllFiniteScores[labeledAllFiniteScores$label_vivo == 1,])
nrow(labeledAllFiniteScores[labeledAllFiniteScores$label_vivo == 0,])


# create train and test data sets from labeled data
require(caTools)
set.seed(123)
sample = sample.split(labeledAllFiniteScores$label_vivo, SplitRatio = 0.80)
labeledAllFiniteScores
trainingData <- subset(labeledAllFiniteScores, sample == 1)
testingData <- subset(labeledAllFiniteScores, sample == 0)




# build logistic regression model
y <- "label_vivo"

allMaxent_allSkippy <- c("wt5score", "mu5score", "wt3score", "mu3score", "l2fcSpliceRatio_vivo", "l2fcSpliceRatio_vitro", "logoddsratio", "exonlen", "deltaSS3",
  "deltaSS5")

bestCombo <- c("mu5score", "mu3score", "logoddsratio", "exonlen", "deltaSS3", "deltaSS5")

onlyMaxent <-  c("wt5score", "mu5score", "wt3score", "mu3score")

onlySkippy <-  c("logoddsratio", "exonlen", "deltaSS3", "deltaSS5")




formulaString <- paste(y, paste(onlyMaxent, collapse="+"), sep="~")
LRmodel <- glm(formulaString, data=trainingData, family=binomial(link="logit"))

# make predictions using model
trainingData$prediction <- predict(LRmodel, newdata=trainingData, type="response")
testingData$prediction <- predict(LRmodel, newdata=testingData, type="response")

# plot distribution of prediction score grouped by known outcome
library(ggplot2)
ggplot(trainingData, aes(x=prediction, color=label_vivo, linetype=as.factor(label_vivo))) + geom_density()

# plot precision and recall

library(ROCR)
library(grid)
predObj <- prediction(trainingData$prediction, trainingData$label_vivo)
precObj <- performance(predObj, measure="prec")
recObj <- performance(predObj, measure="rec")
precision <- (precObj@y.values)[[1]]
prec.x <- (precObj@x.values)[[1]]
recall <- (recObj@y.values)[[1]]
rocFrame <- data.frame(threshold=prec.x, precision=precision, recall=recall)
nplot <- function(plist) {
  n <- length(plist)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(n,1)))
  vplayout=function(x,y) {viewport(layout.pos.row=x, layout.pos.col=y)}
  for(i in 1:n) {
    print(plist[[i]], vp=vplayout(i,1))
  }
}

pnull <- mean(as.numeric(trainingData$label_vivo))

p1 <- ggplot(rocFrame, aes(x=threshold)) +
  geom_line(aes(y=precision/pnull)) +
  coord_cartesian(xlim=c(0,1), ylim=c(0,10))
 
p2 <- ggplot(rocFrame, aes(x=threshold)) +
  geom_line(aes(y=recall)) +
  coord_cartesian(xlim = c(0,1))

nplot(list(p1, p2))


# build confusion matrix
ctab.test <- table(pred=testingData$prediction<0.10, label=testingData$label_vivo)
ctab.test

ctab.test <- table(pred=trainingData$prediction<0.10, label=trainingData$label_vivo)
ctab.test

# calculate TPR, FPR, TNR, and FNAR and enrichment
TPR <- ctab.test[2,2]/sum(ctab.test[2,])
TPR
FPR <- ctab.test[2,1]/sum(ctab.test[2,])
FPR
TNR <- ctab.test[1,1]/sum(ctab.test[1,])
TNR
FNR <- ctab.test[1,2]/sum(ctab.test[1,])
FNR

enrichment <- precision/mean(as.numeric(testingData$label_vivo))
testingData$label_vivo
enrichment

AIC(LRmodel)

# examine coefficients of model
coefficients(LRmodel)

# examine LRmodel summary
summary(LRmodel)


# compare pseudo R^2
ll.null <- LRmodel$null.deviance/-2
ll.proposed <- LRmodel$deviance/-2
(ll.null-ll.proposed)/ll.null

# plot the data
predicted.data <- data.frame(
  probability.of.hd=LRmodel$fitted.values,
  label=trainingData$label_vivo)

predicted.data <- data.frame(probability.of.hd=LRmodel$fitted.values, label=trainingData$label_vivo)

predicted.data <- predicted.data[
  order(predicted.data$probability.of.hd, decreasing=FALSE),]

predicted.data$rank <- 1:nrow(predicted.data)

ggplot(data=predicted.data, aes(x=rank, y=probability.of.hd)) +
  geom_point(aes(color=label), alpha=1, shape=4, stroke=2) +
  xlab("Index") +
  ylab("Predicted probability of altered splicing")

# compute deviance
loglikelihood <- function(y, py) {
  sum(y * log(py) + (1-y)*log(1-py))
}

pnull <- mean(as.numeric(trainingData$label_vivo))
null.dev <- 2*loglikelihood(as.numeric(trainingData$label_vivo), pnull)
pnull
null.dev
LRmodel$null.deviance

pred <- predict(LRmodel, newdata=trainingData, type="response")
resid.dev <- -2*loglikelihood(as.numeric(trainingData$label_vivo), pred)

resid.dev

LRmodel$deviance

testy <- as.numeric(testingData$label_vivo)
testpred <- predict(LRmodel, newdata=testingData, type="response")

pnull.test <- mean(testy)
null.dev.test <- -2*loglikelihood(testy, pnull.test)
resid.dev.test <- -2*loglikelihood(testy, testpred)

pnull.test

null.dev.test

resid.dev.test

















