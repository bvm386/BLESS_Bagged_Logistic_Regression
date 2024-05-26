# BLESS: Bagged Logistic Regression Algorithm
# Kyle Gardiner, Xuekui Zhang, Li Xing

 # load in the necessary libraries
library(dplyr)


# load in the data
load("data.rda")
BLESS_data[1:5, 1:5]
BLESS_data[1:5, c(1001:1004, 1011)]


# Create a dataframe to store the feature names and pvalues
feature.df<-data.frame(Feature = character(0), pvalues=numeric(0))



############# BLESS Algorithm #####################

for(ii in 1:1000){ # 1000 is the number of iterations performed. 
  ## This number can be changed depending on the data being analyzed
  
  ## Setting seed each iteration for reproducibility
  set.seed(ii) 
  
  # Step 1) Subsampling Observations and Selecting Subset of Features
  ## Randomly subsampling 90% of observations
  obs.id<-sample(nrow(BLESS_data),round(nrow(BLESS_data)*0.9), replace=FALSE)
  
  ## Randomly subsampling 30 predictors
  ### 30 predictors is roughly square root the numbers of predictors
  #### Only samples the SNPs (columns 1:1000)
  sample.SNPs<-sample(ncol(BLESS_data[,1:1000]),30, replace=FALSE)
  
  ## Combining selected observations and predictors with the outcome Y
  ### Include covariates of interest (for now, include Cov.1-Cov.4)
  ### Column 1011 is the outcome variable
  sub.data<-BLESS_data[obs.id,c(sample.SNPs,1001:1004,1011)]
  
  
  # Step 2) Fitting a logistic model with each subset data
  myfit<-glm(Y~., data=sub.data, family="binomial")
  
  
  ## Extracting the feature names and p-values for this logistic regression model
  selected.feature<-row.names(coef(summary(myfit))[2:30,])
  selected.pvalue<-coef(summary(myfit))[2:30,4]
  
  
  ## Storing the extracted features and p-values in a temporary dataframe
  temp.feature.df<-data.frame(Feature=selected.feature, pvalue=selected.pvalue)
  
  ## Binding all the outputs from the logistic regression models over all iterations
  feature.df<- rbind(feature.df, temp.feature.df)
}


# Step 3) Aggregating P-values
## Grouping the output by feature name
grouped <- split(feature.df, feature.df$Feature)

## Create an empty dataframe to store the adjusted p-values
adj.feature<-data.frame(Feature=character(0), adj.pvalue=numeric(0))


## Apply FDR to each grouped feature
for(ss in 1:length(grouped)){ # looping for how many features are in the `grouped` object
  
  ## applying FDR testing correction to the grouped feature p-values
  adj.p.value<-p.adjust(grouped[[ss]]$pvalue, method="fdr")
  
  
  grouped[[ss]]$adj.pvalue<-adj.p.value # creating a new column for the adjusted p-values
  
  
  ## Checking if any of the newly created adjusted p-values <0.05
  ### If one of the adjusted p.values<0.05, then the feature is added to result
  if(any(grouped[[ss]]$adj.pvalue<0.05)){
    
    ## If the feature passes the threshold, extract the minimum adjust p-value for that feature
    smallest.id<-which(grouped[[ss]]$adj.pvalue==min(grouped[[ss]]$adj.pvalue))[1]
    
    ## Storing the feature name and smallest adjusted p-value in a temporary dataframe
    temp<-grouped[[ss]][smallest.id,c(1,3)]
    row.names(temp)<-NULL
    
    
    ## Storing the selected feature and adjusted p-value across all features in the `grouped` object
    adj.feature<-rbind(adj.feature, temp)
  }
}

# Removing any duplicate values
adj.feature<-adj.feature[!duplicated(adj.feature$Feature),][-1,]

# Step 4) Obtaining Final Ranked Output Based on Adjusted P-values
# Ranking the output based on FDR adjusted p-values
ranked.adj.feature<-adj.feature[order(adj.feature$adj.pvalue, decreasing=FALSE),]
row.names(ranked.adj.feature)<-NULL

head(ranked.adj.feature) # Displaying the top identified features























