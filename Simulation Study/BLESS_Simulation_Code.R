# Kyle Gardiner
# BLESS Simulation Code
# Feb. 22, 2024

library(tidyverse)

# Load in the data
load("data_for_sim.rda")

# Extract the top 30 SNPs for later comparison
top.30.SNP<-data.frame(SNPs=names(simulation.data)[1:30])
metrics<-data.frame(recall=rep(NA,100),precision=rep(NA,100), accuracy=rep(NA,100))

# creating simulation parameters
prop_observation <- c(0.8, 0.9, 1)  # values for prop_observation
num_feat <- c(10, 30, 50)       # values for num_snps
parameter_combinations <- expand.grid(prop_observation = prop_observation, num_feat = num_feat)


# Cycling through each combinations of simulation parameters
for (param_index in 1:nrow(parameter_combinations)) {
  prop_observation <- parameter_combinations$prop_observation[param_index]
  num_feat <- parameter_combinations$num_feat[param_index]
  
  
  
  # Looping for each of the 100 simulated datasets
  for(jj in 1:100){
    
    # creating file path for the simulated data
    file_path<- paste0("sim_data_",jj,".rda")
    
    # reading all the simulated data files
    assign(paste0("simulated_data_", jj), load(file_path), envir = .GlobalEnv)
    
    
    #################### BLESS #################################
    
    # creating empty dataframe to store selected SNPs
    assign(paste0("selected_SNP_df_sim_", jj), data.frame(SNPs = character(0), pvalues=numeric(0)), envir = .GlobalEnv)
    
    
    
    # Loop for bagged logistic regression (1000 iterations)
    for(ii in 1:1000){
      
      # selecting 90% of the observations and 30 SNPs
      set.seed(ii)
      
      obs.id<-sample(nrow(get(paste0("sim.data.", jj))),round(nrow(get(paste0("sim.data.", jj)))*prop_observation), replace=FALSE)
      
      sample.SNPs<-sample(ncol(get(paste0("sim.data.", jj))[,-772]),num_feat, replace=FALSE) # column 772 is the outcome
      
      sub.data<-get(paste0("sim.data.", jj))[obs.id,c(sample.SNPs,772)]
      
      # now run logisic regression on sub.data
      myfit<-glm(simulated.outcome~., data=sub.data, family="binomial")
      
      # selecting all the names and SNPs
      selected.SNP<-row.names(coef(summary(myfit))[-1,])
      
      selected.pvalue<-coef(summary(myfit))[-1,4]
      
      # storing selected SNPs in temporary dataframe
      temp.snp.df<-data.frame(SNP=selected.SNP,pvalue=selected.pvalue)
      # saving all selected SNPs from each iteration to a dataframe
      assign(paste0("selected_SNP_df_sim_", jj), rbind(get(paste0("selected_SNP_df_sim_", jj)), temp.snp.df), envir = .GlobalEnv)
      
      
      print(ii)
    }
    
    
    
    
    # sorting frequency of important snps and converting to data.frame 
    final.output.df<-get(paste0("selected_SNP_df_sim_", jj))
    row.names(final.output.df)<-NULL
    
    # grouping the data by SNP name
    grouped <- split(final.output.df, final.output.df$SNP)
    
    # applying FDR correction
    adj.selected<-character(0)
    
    for(ss in 1:length(grouped)){
      adj.p.value<-p.adjust(grouped[[ss]]$pvalue, method="fdr")
      
      # If one of the adjusted p.values<0.05, then its added to result
      if(any(adj.p.value<0.05)){
        adj.selected<-rbind(adj.selected,grouped[[ss]]$SNP[1])
      }
    }
    
    # Calculating Recall, Precision, and Accuracy
    real.positive<-length(top.30.SNP$SNPs)
    real.negative<-length(which(!colnames(simulation.data[,-1]) %in% top.30.SNP$SNPs))
    
    pred.positive<-length(adj.selected)
    pred.negative<-length(which(!colnames(simulation.data[,-1]) %in% adj.selected))
    
    true.positive<-length(which(adj.selected %in% top.30.SNP$SNPs))
    false.negative<-length(top.30.SNP$SNPs)-true.positive
    true.negative<-pred.negative-false.negative
    false.positive<-real.negative-true.negative
    
    
    
    metrics$recall[jj]<-(true.positive/real.positive)
    metrics$precision[jj]<-(true.positive/pred.positive)
    metrics$accuracy[jj]<- ((true.positive+true.negative)/dim(simulation.data[,-1])[2])
    
    # sorting the frequency, thinking I don't need this anymore
    char.counts<-sort(table(final.output.df$SNP[which(adj.selected %in% final.output.df$SNP)]), decreasing=TRUE)
    assign(paste0("snp_freq_sim_", jj), as.data.frame(char.counts), envir = .GlobalEnv)
    
    
    # saving each ordered frequency table as rda files
    save(list= paste0("snp_freq_sim_", jj), file = paste0("SNP_freq_sim_", jj, ".rda"))
    
    print(jj)
  }
  
  # saving the simulation results for each paramter combination
  save(metrics, file= paste0("sim_metrics_",num_feat,"feat_",prop_observation,"obs.rda" ))
}
