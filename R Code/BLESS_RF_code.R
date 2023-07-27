# BLESS: Bagged Logistic Regression Algorithm
# Kyle Gardiner
# July 27, 2023

# load in libraries
library(caret)
library(randomForest)


# load in data
load("data.rda")


# making empty dataframe to store selected SNPs
predictor.df<-data.frame(Predictor=character(0))


############################ BLESS #############################

# Loop for bagged logistic regression
for(ii in 1:5000){
  
  # subsampling observations and predictors
  set.seed(ii)# setting seed each iteration for reproducibility
  
  obs.id<-sample(nrow(BLESS_data),round(nrow(BLESS_data)*0.9)) # subsampling 90% of observations
  sample.SNPs<-sample(ncol(BLESS_data[,-1001]),30) # subsampling 30 predictors
  sub.data<-BLESS_data[obs.id,c(sample.SNPs,1011)]
  
  # fitting logistic model with subset data
  myfit<-glm(Y~., data=sub.data, family="binomial")
  
  # determining which predictor have p-values < 0.05 
  selected.predictor.id<-which(coef(summary(myfit))[-1,4]<0.05)
  selected.predictor<-row.names(coef(summary(myfit))[-1,])[selected.predictor.id]
  
  # adding selected predictors from iteration to iteration
  temp.predictor.df<-data.frame(Predictor=selected.predictor)
  predictor.df<-rbind(predictor.df,temp.predictor.df)
  
}

# sorting frequency of important predictors and converting to data.frame
char.counts<-sort(table(predictor.df),decreasing=TRUE)
predictor.freq_5000<-as.data.frame(char.counts)



###################### Random Forest ########################
set.seed(1) # setting seed for reproducibility

# fitting random forest model
rf.model_5000<-randomForest(Y~.,data=BLESS_data,
                                ntree=5000,
                                importance=TRUE)

# plotting predictor importance
varImpPlot(rf.model_5000,main="RF for n=5000 trees")

