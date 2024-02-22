# Kyle Gardiner
# BLESS Data Simulation Code
# Feb. 22, 2024


# Load in the data used to simulate
load("data_for_sim.rda")

# creating sub data with top 30 SNPs and outcome
sub.SNP<-simulation.data[,1:30]
Y<-simulation.data[,772]
sub.data<-cbind(sub.SNP, Y)


# fitting logistic regression model
fit<-glm(Y~.,data=sub.data, family="binomial")


# creating 0.5 threshold for randomly generating binary outcome
threshold<-0.5

############# Generating 100 simulated outcome data sets ###################

for(jj in 1:100){
  # setting seed for reproducibility
  set.seed(jj)
  
  # adding noise to fixed effects
  noisy_coeff<-fit$coefficients + rnorm(length(fit$coefficients), mean=0.3, sd=0.1)
  # setting the intercept beta to 0
  noisy_coeff[1]<-0
  
  # constructing linear model with noisy coefficients
  linear_predictor<- as.numeric(noisy_coeff[1] + as.matrix(sub.SNP) %*% noisy_coeff[-1])
  
  # calculating probability from logit function
  pred.prob<-(exp(linear_predictor)/(1+exp(linear_predictor)))
  
  # simulating outcome data
  simulated.outcome<-rbinom(n=length(pred.prob), size=1, prob=pred.prob >= threshold)
  
  # combining simulated outcome with original SNP data
  sim.data<-cbind(simulation.data[-772], simulated.outcome)
  
  # assigning each data frame its own name
  assign(paste0("sim.data.",jj), sim.data)
  
  # constructing the file name
  file_name<-paste0("sim_data_", jj, ".rda")
  
  # save simulated data to .rds file
  save(list = paste0("sim.data.", jj), file = file_name)
}
