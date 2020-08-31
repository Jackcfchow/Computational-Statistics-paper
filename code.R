rm(list = ls())
set.seed(1)
#setwd("Set Your Own Directory")
#setwd("Compstat_Chow")


library(MASS)
library(grf)
library(sandwich)
library(lmtest)
library(Hmisc)
library(ggplot2)
library(dplyr)
library(readr)
library(haven)
library(rddtools)
library(rddtools)
library(data.table)

#section 1 deterministic and probabilistic assignment of treatment 
#


############################ section  5.1 Deterministic vs probabilistic ########################################
#import datset 
ori_data <- read_csv("20140162_Data_updated.csv")

#simulate data 
N= 10000
gen_dataset<- function(N, ori_data){
  deno.exp = model.matrix(~ factor(ori_data$denomination)) #factorise categorical variable: denomination
  colnames(deno.exp)<-paste0("denomination.", sort(unique(ori_data$denomination)))
  deno.exp=deno.exp[,2:5]
  ori_data <- subset(ori_data, select = -c(city, area,group ,attendance,denomination,bins, n_bins, binnr,winnersvsalwayslosers, losersvsalwayswinners))
  #exclude unusable columns 
  
  ori_data <- cbind(ori_data,deno.exp)
  ori_data <- na.omit(ori_data) # dropping missing value column 
  
  dataset_sim <- as.data.frame(mvrnorm(n = N, mu = colMeans(ori_data), Sigma = cov(ori_data)))
  # simulating the  dataset using the mean and covariance of the original dataset 
  for (i in {1:27}) { #rounding the continous columns into categorical columns, capinf by the maximum and minimum of the values in the original dataset
    maxi = max(ori_data[,i])
    mini = min(ori_data[,i])
    dataset_sim[i][dataset_sim[i]>maxi] <- maxi
    dataset_sim[i][dataset_sim[i]<mini] <- mini
    dataset_sim[i] = round(dataset_sim[i],digits =0)
  }
  
  for (i in {1:nrow(dataset_sim)}) {
    # for denomination, one can only have one religion. 
    #So we choose the the religion with the highest value to be the denomination of an observation. Other denomination columns are assgined to zero
    maxi_reli <- max(dataset_sim[i,28:31])
    for (j in {28:31}) {
      dataset_sim[i,j][dataset_sim[i,j]==maxi_reli] <- 1
      dataset_sim[i,j][dataset_sim[i,j]<maxi_reli] <- 0
    }
  }
  return(dataset_sim)
}
dataset <- gen_dataset(N,ori_data)

#simulation the genertaion of outcome Y
gen_treatment_outcome_sharp_fuzzy <- function(dataset_sim){#we simulate 2 situations, treatment are determinisitic vs probabilistic 
  mu <- rnorm(nrow(dataset_sim),0,35) #error for treatment assignment 
  dataset_sim$expenditures_temp <- dataset_sim$expenditures+ mu #add randomness to expenditure when deciding treatment
  dataset_sim$W <- ifelse(dataset_sim$expenditures <=200, 1, 0)
  dataset_sim$W_sharp <- ifelse(dataset_sim$expenditures <=200, 1, 0)     #treatment is not affected by mu 
  dataset_sim$W_fuzzy <- ifelse(dataset_sim$expenditures_temp <=200, 1, 0)  #treatment is affected by mu 
  
  mu <- rnorm(nrow(dataset_sim),0,5) #error for Y
  beta_w <<- 7 
  
  dataset_sim$Y_sharp<- (dataset_sim$communal_meal_parti +dataset_sim$common_purchase_food_parti+dataset_sim$community_workshop_parti+dataset_sim$work_exchange_parti
                         +dataset_sim$fundraising_parti + dataset_sim$communal_child_care_parti + dataset_sim$prepare_for_fund_parti + dataset_sim$communal_construction_parti
                         +dataset_sim$property_invasion_parti + dataset_sim$security_committee_parti + dataset_sim$election_champaign_parti+
                           # participation dummy variable have a beta =1
                           + 1.4*dataset_sim$schooling_resp  + 0.4*dataset_sim$householdsize + 0.2*dataset_sim$ageresponder +
                           +  beta_w*dataset_sim$W_sharp + #true beta is set to be 7
                           +  0.05*dataset_sim$expenditures
                         +  mu)
  
  dataset_sim$Y_fuzzy<- (dataset_sim$communal_meal_parti +dataset_sim$common_purchase_food_parti+dataset_sim$community_workshop_parti+dataset_sim$work_exchange_parti
                         +dataset_sim$fundraising_parti + dataset_sim$communal_child_care_parti + dataset_sim$prepare_for_fund_parti + dataset_sim$communal_construction_parti
                         +dataset_sim$property_invasion_parti + dataset_sim$security_committee_parti + dataset_sim$election_champaign_parti+
                           # participation dummy variable have a beta =1
                           + 1.4*dataset_sim$schooling_resp  + 0.4*dataset_sim$householdsize + 0.2*dataset_sim$ageresponder +
                           +  beta_w*dataset_sim$W_fuzzy + #true beta is set to be 7
                           +  0.05*dataset_sim$expenditures
                         +  mu)
  
  True_ATE<<- beta_w
  return(dataset_sim)
}

dataset_sharp_fuzzy<- gen_treatment_outcome_sharp_fuzzy(dataset)

#visualise treatment recived in sharp and fuzzy case
X <-as.matrix (cbind(dataset_sharp_fuzzy[,2:12], dataset_sharp_fuzzy$householdsize, dataset_sharp_fuzzy$ageresponder, dataset_sharp_fuzzy$schooling_resp)  )

rdd_w_sharp <- rdd_data(x=dataset_sharp_fuzzy$expenditures, y=dataset_sharp_fuzzy$W ,z=ifelse(dataset_sharp_fuzzy$expenditures<200,1,0), covar =  X,   data=dataset_sharp_fuzzy, cutpoint=200)
png(filename="rddgraph_sharp.png")
rddgraph_sharp <- plot(rdd_w_sharp, nbins=N/10,xlab="Expenditure", ylab = "Treatment",main = "Figure 1a. Deterministic treatment assignment" , ylim=c(0,1))
dev.off()

rdd_w_fuzzy <- rdd_data(x=dataset_sharp_fuzzy$expenditures, y=dataset_sharp_fuzzy$W_fuzzy ,z=ifelse(dataset_sharp_fuzzy$expenditures<200,1,0), covar =  X,   data=dataset_sharp_fuzzy, cutpoint=200)
png(filename="rddgraph_fuzzy.png")
rddgraph_fuzzy <- plot(rdd_w_fuzzy, nbins=N/10,xlab="Expenditure", ylab = "Treatment",main = "Figure 1b. Probabilistic treatment assignment" , ylim=c(0,1))
dev.off()



#############RDD estimation ###################
RDD_sharp_fuzzy <- function(dataset_sim) {
  
  fit_Y_sharp <- lm(Y_sharp ~ W_sharp   + expenditures + householdsize + ageresponder + 
                      communal_meal_parti +common_purchase_food_parti+community_workshop_parti+work_exchange_parti
                    +fundraising_parti + communal_child_care_parti + prepare_for_fund_parti + communal_construction_parti
                    +property_invasion_parti + security_committee_parti + election_champaign_parti              
                    , data=dataset_sim)
  
  
  fit_Y_fuzzy <- lm(Y_fuzzy ~ W_fuzzy  + expenditures + householdsize + ageresponder + 
                      communal_meal_parti +common_purchase_food_parti+community_workshop_parti+work_exchange_parti
                    +fundraising_parti + communal_child_care_parti + prepare_for_fund_parti + communal_construction_parti
                    +property_invasion_parti + security_committee_parti + election_champaign_parti              
                    , data=dataset_sim)
  
  ATE_RDD_sharp <- as.numeric (fit_Y_sharp$coefficients[2] )
  ATE_RDD_fuzzy <- as.numeric (fit_Y_fuzzy$coefficients[2] )
  Bias_ATE_RDD_sharp = ATE_RDD_sharp - True_ATE
  Bias_ATE_RDD_fuzzy = ATE_RDD_fuzzy - True_ATE
  
  return (cbind( True_ATE= True_ATE, ATE_RDD_sharp,Bias_ATE_RDD_sharp,ATE_RDD_fuzzy,Bias_ATE_RDD_fuzzy))
}


#############CF estimation ###################
forest_sharp_fuzzy<-function(forest_data){
  Y_sharp <- forest_data$Y_sharp
  Y_fuzzy <- forest_data$Y_fuzzy 
  W<- forest_data$W
  W_sharp<-forest_data$W_sharp 
  W_fuzzy<- forest_data$W_fuzzy 
  
  
  X <-as.matrix (cbind(forest_data[,1:12], forest_data$householdsize, forest_data$ageresponder, forest_data$schooling_resp)  )
  Y_sharp.forest = regression_forest(X, Y_sharp)
  Y_sharp.hat = predict(Y_sharp.forest)$predictions
  
  Y_fuzzy.forest = regression_forest(X, Y_fuzzy)
  Y_fuzzy.hat = predict(Y_fuzzy.forest)$predictions
  cf_sharp = causal_forest(X , Y_sharp, W_sharp,
                           Y.hat = Y_sharp.hat,
                           # W.hat = W_sharp.hat,
                           tune.parameters = "all" )
  
  ATE_CF_sharp = as.numeric (average_treatment_effect(cf_sharp, target.sample = "all")[1])
  ATE_overlap_CF_sharp = as.numeric (average_treatment_effect(cf_sharp, target.sample = "overlap")[1])
  
  Bias_ATE_CF_sharp = True_ATE- ATE_CF_sharp
  Bias_ATE_overlap_CF_sharp = ATE_overlap_CF_sharp - True_ATE
  
  
  cf_fuzzy = causal_forest(X , Y_fuzzy, W_fuzzy,
                           Y.hat = Y_fuzzy.hat,
                           tune.parameters = "all")
  
  tau_fuzzy.hat = predict(cf_fuzzy)$predictions
  ATE_CF_fuzzy = as.numeric (average_treatment_effect(cf_fuzzy, target.sample = "all")[1])
  ATE_overlap_CF_fuzzy = as.numeric (average_treatment_effect(cf_fuzzy, target.sample = "overlap")[1]) 
  
  Bias_ATE_CF_fuzzy = ATE_CF_fuzzy- True_ATE 
  Bias_ATE_overlap_CF_fuzzy = ATE_overlap_CF_fuzzy - True_ATE
  
  
  return(cbind(ATE_CF_sharp,Bias_ATE_CF_sharp, 
               ATE_overlap_CF_sharp, Bias_ATE_overlap_CF_sharp,
               ATE_CF_fuzzy,Bias_ATE_CF_fuzzy,
               ATE_overlap_CF_fuzzy,Bias_ATE_overlap_CF_fuzzy ))
}



for (n in c(500,1000,5000,10000)) {
  
  result_sharp_fuzzy <- as.data.frame (NULL)
  
  for (i in 1:100) {
    print(i)
    set.seed(i)
    N=n
    dataset <- gen_dataset(N,ori_data)
    dataset_sharp_fuzzy<- gen_treatment_outcome_sharp_fuzzy(dataset)
    RDD_result_sharp_fuzzy <- RDD_sharp_fuzzy(dataset_sharp_fuzzy)
    cf_result_sharp_fuzzy<- forest_sharp_fuzzy(dataset_sharp_fuzzy)
    result_sharp_fuzzy_temp <- cbind(True_ATE= True_ATE, sample_size=n  , RDD_result_sharp_fuzzy ,cf_result_sharp_fuzzy )
    result_sharp_fuzzy <- rbind(result_sharp_fuzzy, result_sharp_fuzzy_temp)
    
  }
  write.csv(result_sharp_fuzzy, file=paste("result_sharp_fuzzy_", n ,"sample.csv" , sep=""))
  
}


######merge all RDD anf CF results into one statistics file

#result for section  5.1 
result_sharp_fuzzy_500sample <- read_csv("result_sharp_fuzzy_500sample.csv")
result_sharp_fuzzy_1000sample <- read_csv("result_sharp_fuzzy_1000sample.csv")
result_sharp_fuzzy_5000sample <- read_csv("result_sharp_fuzzy_5000sample.csv")
result_sharp_fuzzy_10000sample <- read_csv("result_sharp_fuzzy_10000sample.csv")
# merge resultant file into one file 
result_sharp_fuzzy <- rbind(result_sharp_fuzzy_500sample,result_sharp_fuzzy_1000sample,result_sharp_fuzzy_5000sample ,result_sharp_fuzzy_10000sample)
result_sharp_fuzzy <- na.omit(result_sharp_fuzzy) 

result_sharp_fuzzy$Bias_ATE_RDD_sharp_per <- abs(result_sharp_fuzzy$Bias_ATE_RDD_sharp)/ result_sharp_fuzzy$True_ATE 
result_sharp_fuzzy$Bias_ATE_RDD_fuzzy_per <- abs(result_sharp_fuzzy$Bias_ATE_RDD_fuzzy)/ result_sharp_fuzzy$True_ATE 

result_sharp_fuzzy$Bias_ATE_CF_sharp_per <- abs(result_sharp_fuzzy$Bias_ATE_CF_sharp)/ result_sharp_fuzzy$True_ATE
result_sharp_fuzzy$Bias_ATE_overlap_CF_sharp_per <- abs(result_sharp_fuzzy$Bias_ATE_overlap_CF_sharp)/ result_sharp_fuzzy$True_ATE
result_sharp_fuzzy$Bias_ATE_CF_fuzzy_per <- abs(result_sharp_fuzzy$Bias_ATE_CF_fuzzy)/ result_sharp_fuzzy$True_ATE
result_sharp_fuzzy$Bias_ATE_overlap_CF_fuzzy_per <- abs(result_sharp_fuzzy$Bias_ATE_overlap_CF_fuzzy)/ result_sharp_fuzzy$True_ATE

aggregate_result_sharp_fuzzy <- as_data_frame(NULL)
stat_result_sharp_fuzzy <- as_data_frame(NULL)

for (n in c(500,1000,5000,10000)) {
  print(n)
  result_subset = subset(result_sharp_fuzzy,sample_size==n )
  result_subset$True_ATE_average <- mean(result_subset$True_ATE)
  result_subset$RMSE_ATE_RDD_sharp <- (result_subset$ATE_RDD_sharp-result_subset$True_ATE_average)^2/nrow(result_subset)
  result_subset$RMSE_ATE_RDD_fuzzy <- (result_subset$ATE_RDD_fuzzy-result_subset$True_ATE_average)^2/nrow(result_subset)
  
  result_subset$RMSE_ATE_CF_sharp <- (result_subset$ATE_CF_sharp-result_subset$True_ATE_average)^2/nrow(result_subset)
  result_subset$RMSE_ATE_overlap_CF_sharp <- (result_subset$ATE_overlap_CF_sharp-result_subset$True_ATE_average)^2/nrow(result_subset)
  
  result_subset$RMSE_ATE_CF_fuzzy <- (result_subset$ATE_CF_fuzzy-result_subset$True_ATE_average)^2/nrow(result_subset)
  result_subset$RMSE_ATE_overlap_CF_fuzzy <- (result_subset$ATE_overlap_CF_fuzzy-result_subset$True_ATE_average)^2/nrow(result_subset)
  
  aggregate_result_sharp_fuzzy <- rbind(aggregate_result_sharp_fuzzy,result_subset)
  
  stat_result_sharp_fuzzy_temp <- cbind(True_ATE = mean(result_subset$True_ATE),sample_size=n
                                        ,Bias_ATE_RDD_sharp_per_mean = mean(result_subset$Bias_ATE_RDD_sharp_per)
                                        ,RMSE_ATE_RDD_sharp_sum = sum(result_subset$RMSE_ATE_RDD_sharp)
                                        ,Bias_ATE_RDD_fuzzy_per_mean = mean(result_subset$Bias_ATE_RDD_fuzzy_per)
                                        ,RMSE_ATE_RDD_fuzzy_sum = sum(result_subset$RMSE_ATE_RDD_fuzzy)
                                        
                                        ,Bias_ATE_CF_sharp_per_mean = mean(result_subset$Bias_ATE_CF_sharp_per)
                                        ,RMSE_ATE_CF_sharp_sum = sum(result_subset$RMSE_ATE_CF_sharp)
                                        ,Bias_ATE_overlap_CF_sharp_per_mean = mean(result_subset$Bias_ATE_overlap_CF_sharp_per)
                                        ,RMSE_ATE_overlap_CF_sharp_sum = sum(result_subset$RMSE_ATE_overlap_CF_sharp)
                                        
                                        ,Bias_ATE_CF_fuzzy_per_mean = mean(result_subset$Bias_ATE_CF_fuzzy_per)
                                        ,RMSE_ATE_CF_fuzzy_sum = sum(result_subset$RMSE_ATE_CF_fuzzy)
                                        ,Bias_ATE_overlap_CF_fuzzy_per_mean = mean(result_subset$Bias_ATE_overlap_CF_fuzzy_per)
                                        ,RMSE_ATE_overlap_CF_fuzzy_sum = sum(result_subset$RMSE_ATE_overlap_CF_fuzzy)
                                        
  )
  
  stat_result_sharp_fuzzy <- rbind(stat_result_sharp_fuzzy,stat_result_sharp_fuzzy_temp)
  
}
t_stat_result_sharp_fuzzy <- transpose(stat_result_sharp_fuzzy)
colnames(t_stat_result_sharp_fuzzy) <- c("Sample size =500","Sample size =1000","Sample size =5000","Sample size =10000")
rownames(t_stat_result_sharp_fuzzy) <- colnames(stat_result_sharp_fuzzy)
write.csv(t_stat_result_sharp_fuzzy, file=paste("stat_result_sharp_fuzzy.csv" ))



############################ section  5.2 Hetergeneity and nonlinearity ########################################

#import datset 
ori_data <- read_csv("20140162_Data_updated.csv")
#ori_data <- read_csv("The Effect of Income on Religiousness/Data/20140162_INEC.csv")

#simulate data 

gen_dataset<- function(N, ori_data){
  deno.exp = model.matrix(~ factor(ori_data$denomination)) #factorise categorical variable: denomination
  colnames(deno.exp)<-paste0("denomination.", sort(unique(ori_data$denomination)))
  deno.exp=deno.exp[,2:5]
  ori_data <- subset(ori_data, select = -c(city, area,group ,attendance,denomination,bins, n_bins, binnr,winnersvsalwayslosers, losersvsalwayswinners))
  #exclude unusable columns 
  ori_data <- cbind(ori_data,deno.exp)
  ori_data <- na.omit(ori_data) # dropping missing value column 
  
  dataset_sim <- as.data.frame(mvrnorm(n = N, mu = colMeans(ori_data), Sigma = cov(ori_data)))
  # simulating the  dataset using the mean and covariance of the original dataset 
  for (i in {1:27}) { #rounding the continous columns into categorical columns, capinf by the maximum and minimum of the values in the original dataset
    maxi = max(ori_data[,i])
    mini = min(ori_data[,i])
    dataset_sim[i][dataset_sim[i]>maxi] <- maxi
    dataset_sim[i][dataset_sim[i]<mini] <- mini
    dataset_sim[i] = round(dataset_sim[i],digits =0)
  }
  
  for (i in {1:nrow(dataset_sim)}) {
    # for denomination, one can only have one religion. 
    #So we choose the the religion with the highest value to be the denomination of an observation. Other denomination columns are assgined to zero
    maxi_reli <- max(dataset_sim[i,28:31])
    for (j in {28:31}) {
      dataset_sim[i,j][dataset_sim[i,j]==maxi_reli] <- 1
      dataset_sim[i,j][dataset_sim[i,j]<maxi_reli] <- 0
    }
  }
  return(dataset_sim)
}

N=10000

dataset <- gen_dataset(N,ori_data)


#simulation 
#treatment effect is correlated with expenditure 
gen_treatment_outcome_heter <- function(dataset_sim){
  mu <- rnorm(nrow(dataset_sim),0,35)#error for treatment assignment 
  dataset_sim$expenditures_temp <- dataset_sim$expenditures+ mu #add randomness to expenditure when deciding treatment
  dataset_sim$W <- ifelse(dataset_sim$expenditures_temp <=200, 1, 0) #treatment is affected by mu 
  
  mu <- rnorm(nrow(dataset_sim),0,5)#error for Y
  beta_w <<- 7 
  beta_w_exp <- 0.05
  beta_w_exp_sq <- -0.00007
  
  dataset_sim$Y<- (dataset_sim$communal_meal_parti +dataset_sim$common_purchase_food_parti+dataset_sim$community_workshop_parti+dataset_sim$work_exchange_parti
                   +dataset_sim$fundraising_parti + dataset_sim$communal_child_care_parti + dataset_sim$prepare_for_fund_parti + dataset_sim$communal_construction_parti
                   +dataset_sim$property_invasion_parti + dataset_sim$security_committee_parti + dataset_sim$election_champaign_parti+
                     # participation dummy variable have a beta =1
                     + 1.4*dataset_sim$schooling_resp  + 0.4*dataset_sim$householdsize + 0.2*dataset_sim$ageresponder +
                     +  beta_w*dataset_sim$W +
                     -0.0045*dataset_sim$expenditures 
                   + beta_w_exp*dataset_sim$expenditures*dataset_sim$W + beta_w_exp_sq*(dataset_sim$expenditures^2)*dataset_sim$W +
                     +  mu)
  True_ATE<<- beta_w + beta_w_exp*mean(dataset_sim$expenditures)  + beta_w_exp_sq*mean(dataset_sim$expenditures^2) #interaction between treatment and expenditure
  
  return(dataset_sim)
}


#w =treatment 
#visualise running variable and Y 
dataset_heter<- gen_treatment_outcome_heter(dataset)
X <-as.matrix (cbind(dataset_heter[,2:12], dataset_heter$householdsize, dataset_heter$ageresponder, dataset_heter$schooling_resp)  )

rdd <- rdd_data(x=dataset_heter$expenditures, y=dataset_heter$Y ,z=ifelse(dataset_heter$expenditures<200,1,0), covar =  X,   data=dataset_heter, cutpoint=200)
png(filename="rddgraph1.png")
rddgraph1 <- plot(rdd, nbins=N/20,xlab="Expenditure", ylab = "Church Attendence",main = "Figure 2. Graphical presentation of RDD", ylim=c(12,40))
dev.off()





reg_para_order1 <- rdd_reg_lm(rdd_object=rdd,order=1)
png(filename="rddgraph_order1.png")
rddgraph_order1 <- plot(reg_para_order1,xlab="Expenditure", ylab = "Church Attendence",main = "Figure 3a. RDD estimation with a 1st order polynomial", ylim=c(12,35)) 
dev.off()

# 
reg_para_order2 <- rdd_reg_lm(rdd_object=rdd,order=2)
# print(reg_para_order2)
png(filename="rddgraph_order2.png")
rddgraph_order2 <- plot(reg_para_order2,xlab="Expenditure", ylab = "Church Attendence",main = "Figure 3b. RDD estimation with a 2nd order polynomial", ylim=c(12,35))
dev.off()

reg_para_order3 <- rdd_reg_lm(rdd_object=rdd,order=3)
# print(reg_para_order2)
png(filename="rddgraph_order3.png")
rddgraph_order3 <- plot(reg_para_order3,xlab="Expenditure", ylab = "Church Attendence",main = "Figure 3c. RDD estimation with a 3rd order polynomial", ylim=c(12,35))
dev.off()



#estimate the ATE using RDD 
RDD_heter <- function(dataset_sim) {
  dataset_sim$W_times_expenditures <- dataset_sim$W * dataset_sim$expenditures
  dataset_sim$W_times_expenditures_sq <- dataset_sim$W * (dataset_sim$expenditures^2)
  dataset_sim$W_times_expenditures_cube <- dataset_sim$W * (dataset_sim$expenditures^3)
  
  fit_order1 <- lm(Y ~ W + W_times_expenditures + expenditures + schooling_resp + householdsize + ageresponder + 
                     communal_meal_parti +common_purchase_food_parti+community_workshop_parti+work_exchange_parti
                   +fundraising_parti + communal_child_care_parti + prepare_for_fund_parti + communal_construction_parti
                   +property_invasion_parti + security_committee_parti + election_champaign_parti       
                   , data=dataset_sim) #1st order polynomial 
  
  
  fit_order2 <- lm(Y ~ W + W_times_expenditures + W_times_expenditures_sq + expenditures + schooling_resp + householdsize + ageresponder + 
                     communal_meal_parti +common_purchase_food_parti+community_workshop_parti+work_exchange_parti
                   +fundraising_parti + communal_child_care_parti + prepare_for_fund_parti + communal_construction_parti
                   +property_invasion_parti + security_committee_parti + election_champaign_parti              
                   , data=dataset_sim) #2nd order polynomial 
  
  fit_order3 <- lm(Y ~ W + W_times_expenditures + W_times_expenditures_sq + W_times_expenditures_cube + expenditures + schooling_resp + householdsize + ageresponder + 
                     communal_meal_parti +common_purchase_food_parti+community_workshop_parti+work_exchange_parti
                   +fundraising_parti + communal_child_care_parti + prepare_for_fund_parti + communal_construction_parti
                   +property_invasion_parti + security_committee_parti + election_champaign_parti              
                   , data=dataset_sim) #3rd order polynomial 
  
  ATE_order1 <- as.numeric (fit_order1$coefficients[2] + fit_order1$coefficients[3]*mean(dataset_sim$expenditures) )
  ATE_order2 <- as.numeric ( fit_order2$coefficients[2] + fit_order2$coefficients[3]*mean(dataset_sim$expenditures) 
                             + fit_order2$coefficients[4]*mean(dataset_sim$expenditures^2) )
  ATE_order3 <- as.numeric ( fit_order3$coefficients[2] + fit_order3$coefficients[3]*mean(dataset_sim$expenditures) 
                             + fit_order3$coefficients[4]*mean(dataset_sim$expenditures^2) 
                             + fit_order3$coefficients[5]*mean(dataset_sim$expenditures^3) )
  Bias_ATE_order1 <- ATE_order1 - True_ATE
  Bias_ATE_order2 <- ATE_order2 - True_ATE
  Bias_ATE_order3 <- ATE_order3 - True_ATE
  Y_hat_order1 = predict(fit_order1)
  Y_hat_order2 = predict(fit_order2)
  Y_hat_order3 = predict(fit_order3)
  
  sample_size <- nrow(dataset_sim)
  return (cbind(ATE_order1 , Bias_ATE_order1
                , ATE_order2,Bias_ATE_order2
                , ATE_order3, Bias_ATE_order3
  ))
}


forest_heter<-function(forest_data){
  Y <- forest_data$Y 
  W<- forest_data$W
  X <-as.matrix (cbind(forest_data[,1:12], forest_data$householdsize, forest_data$ageresponder, forest_data$schooling_resp)  )
  Y.forest = regression_forest(X, Y)
  Y.hat = predict(Y.forest)$predictions
  
  cf = causal_forest(X, Y, W
                     ,Y.hat = Y.hat
  )
  tau.hat = predict(cf)$predictions
  ATE_CF_heter <<- as.numeric (average_treatment_effect(cf, target.sample = "overlap")[1]) # a lot of computational power is needed! 
  # ATT_CF_heter <<- as.numeric (average_treatment_effect(cf_fuzzy, target.sample = "treated")[1]) # a lot of computational power is needed! 
  Bias_ATE_CF_heter = ATE_CF_heter - True_ATE
  return ( cbind(ATE_CF_heter, Bias_ATE_CF_heter ) )
}




for (n in c(500,1000,5000,10000)) {
  result_heter <- as.data.frame (NULL)
  for (i in 1:100) {
    print(i)
    set.seed(i)
    N=n
    dataset <- gen_dataset(N,ori_data)
    dataset_heter<- gen_treatment_outcome_heter(dataset)
    RDD_rseult_heter <- RDD_heter(dataset_heter)
    cf_result_heter = forest_heter(dataset_heter)
    
    result_heter_temp <- cbind(True_ATE= True_ATE , sample_size=n , RDD_rseult_heter,cf_result_heter)
    result_heter <- rbind(result_heter, result_heter_temp)
    
  }
  write.csv(result_heter, file=paste("result_heter_", n ,"sample.csv" , sep=""))
}

######merge all RDD anf CF results into one statistics file
#result for section  5.2 
result_heter_500sample <- read_csv("result_heter_500sample.csv")
result_heter_1000sample <- read_csv("result_heter_1000sample.csv")
result_heter_5000sample <- read_csv("result_heter_5000sample.csv")
result_heter_10000sample <- read_csv("result_heter_10000sample.csv")

result_heter <- rbind(result_heter_500sample,result_heter_1000sample,result_heter_5000sample ,result_heter_10000sample)

result_heter$Bias_ATE_order1_per <- abs(result_heter$Bias_ATE_order1)/ result_heter$True_ATE
result_heter$Bias_ATE_order2_per <- abs(result_heter$Bias_ATE_order2)/ result_heter$True_ATE
result_heter$Bias_ATE_order3_per <- abs(result_heter$Bias_ATE_order3)/ result_heter$True_ATE
result_heter$Bias_ATE_CF_heter_per <- abs(result_heter$Bias_ATE_CF_heter)/ result_heter$True_ATE

aggregate_result_heter <- as_data_frame(NULL)
stat_result_heter <- as_data_frame(NULL)
for (n in c(500,1000,5000,10000)) {
  print(n)
  result_subset = subset(result_heter,sample_size==n )
  result_subset$True_ATE_average <- mean(result_subset$True_ATE)
  result_subset$RMSE_ATE_order1 <- (result_subset$ATE_order1-result_subset$True_ATE_average)^2/nrow(result_subset)
  result_subset$RMSE_ATE_order2 <- (result_subset$ATE_order2-result_subset$True_ATE_average)^2/nrow(result_subset)
  result_subset$RMSE_ATE_order3 <- (result_subset$ATE_order3-result_subset$True_ATE_average)^2/nrow(result_subset)
  result_subset$RMSE_ATE_CF_heter <- (result_subset$ATE_CF_heter-result_subset$True_ATE_average)^2/nrow(result_subset)
  
  aggregate_result_heter <- rbind(aggregate_result_heter,result_subset)
  stat_result_heter_temp <- cbind(True_ATE = mean(result_subset$True_ATE),sample_size=n
                                  ,Bias_ATE_order1_per_mean = mean(result_subset$Bias_ATE_order1_per)
                                  ,RMSE_ATE_order1_sum = sum(result_subset$RMSE_ATE_order1)
                                  
                                  ,Bias_ATE_order2_per_mean = mean(result_subset$Bias_ATE_order2_per)
                                  ,RMSE_ATE_order2_sum = sum(result_subset$RMSE_ATE_order2)
                                  
                                  ,Bias_ATE_order3_per_mean = mean(result_subset$Bias_ATE_order3_per)
                                  ,RMSE_ATE_order3_sum = sum(result_subset$RMSE_ATE_order3)
                                  
                                  ,Bias_ATE_CF_heter_per_mean = mean(result_subset$Bias_ATE_CF_heter_per)
                                  ,RMSE_ATE_CF_heter_sum = sum(result_subset$RMSE_ATE_CF_heter)
  )
  
  stat_result_heter <- rbind(stat_result_heter,stat_result_heter_temp)
  
}
t_stat_result_heter <- transpose(stat_result_heter)

# get row and colnames in order
colnames(t_stat_result_heter) <- c("Sample size =500","Sample size =1000","Sample size =5000","Sample size =10000")
rownames(t_stat_result_heter) <- colnames(stat_result_heter)
write.csv(t_stat_result_heter, file="stat_result_heter.csv")


############################ section  5.3 omitted variable ########################################

rm(list = ls())
set.seed(1)

library(MASS)
library(grf)
library(sandwich)
library(lmtest)
library(Hmisc)
library(ggplot2)
library(dplyr)
library(readr)
library(haven)
library(rddtools)

#import datset 
ori_data <- read_csv("20140162_Data_updated.csv")
#ori_data <- read_csv("The Effect of Income on Religiousness/Data/20140162_INEC.csv")

#simulate data 
N= 1000

gen_dataset<- function(N, ori_data){
  deno.exp = model.matrix(~ factor(ori_data$denomination)) #factorise categorical variable: denomination
  colnames(deno.exp)<-paste0("denomination.", sort(unique(ori_data$denomination)))
  deno.exp=deno.exp[,2:5]
  ori_data <- subset(ori_data, select = -c(city, area,group ,attendance,denomination,bins, n_bins, binnr,winnersvsalwayslosers, losersvsalwayswinners))
  #exclude unusable columns 
  
  ori_data <- cbind(ori_data,deno.exp)
  ori_data <- na.omit(ori_data) # dropping missing value column 
  
  dataset_sim <- as.data.frame(mvrnorm(n = N, mu = colMeans(ori_data), Sigma = cov(ori_data)))
  # simulating the  dataset using the mean and covariance of the original dataset 
  for (i in {1:27}) { #rounding the continous columns into categorical columns, capinf by the maximum and minimum of the values in the original dataset
    maxi = max(ori_data[,i])
    mini = min(ori_data[,i])
    dataset_sim[i][dataset_sim[i]>maxi] <- maxi
    dataset_sim[i][dataset_sim[i]<mini] <- mini
    dataset_sim[i] = round(dataset_sim[i],digits =0)
  }
  
  for (i in {1:nrow(dataset_sim)}) {
    # for denomination, one can only have one religion. 
    #So we choose the the religion with the highest value to be the denomination of an observation. Other denomination columns are assgined to zero
    maxi_reli <- max(dataset_sim[i,28:31])
    for (j in {28:31}) {
      dataset_sim[i,j][dataset_sim[i,j]==maxi_reli] <- 1
      dataset_sim[i,j][dataset_sim[i,j]<maxi_reli] <- 0
    }
  }
  return(dataset_sim)
}



#simulation 
gen_treatment_outcome_omit <- function(dataset_sim){
  mu <- rnorm(nrow(dataset_sim),0,35)#error for treatment assignment 
  dataset_sim$expenditures_temp <- dataset_sim$expenditures+ mu #add randomness to expenditure when deciding treatment
  dataset_sim$W <- ifelse(dataset_sim$expenditures_temp <=200, 1, 0)#treatment is affected by mu 
  mu <- rnorm(nrow(dataset_sim),0,5) #error for Y
  beta_w <- 7 
  
  
  dataset_sim$Y<- (dataset_sim$communal_meal_parti +dataset_sim$common_purchase_food_parti+dataset_sim$community_workshop_parti+dataset_sim$work_exchange_parti
                   +dataset_sim$fundraising_parti + dataset_sim$communal_child_care_parti + dataset_sim$prepare_for_fund_parti + dataset_sim$communal_construction_parti
                   +dataset_sim$property_invasion_parti + dataset_sim$security_committee_parti + dataset_sim$election_champaign_parti+
                     # participation dummy variable have a beta =1
                     + 1.4*dataset_sim$schooling_resp #
                   + 0.4*dataset_sim$householdsize + 0.2*dataset_sim$ageresponder +
                     +  beta_w*dataset_sim$W + #true beta is set to be 7
                     +  rnorm(nrow(dataset_sim),0,5))
  
  
  True_ATE<<- beta_w 
  return(dataset_sim)
}



#############RDD estimation ###################
RDD_omitted <- function(dataset_sim) {
  #schooling is excluded from the model
  fit_Y <- lm(Y~ W   + expenditures + householdsize + ageresponder + #schooling is excluded
                communal_meal_parti +common_purchase_food_parti+community_workshop_parti+work_exchange_parti
              +fundraising_parti + communal_child_care_parti + prepare_for_fund_parti + communal_construction_parti
              +property_invasion_parti + security_committee_parti + election_champaign_parti              
              , data=dataset_sim)
  
  ATE_RDD <- as.numeric (fit_Y$coefficients[2] )
  
  Bias_ATE_RDD = ATE_RDD - True_ATE
  
  return (cbind(ATE_RDD,Bias_ATE_RDD))
}


forest_omitted<-function(forest_data){
  Y <- forest_data$Y
  W<- forest_data$W
  X <-as.matrix (cbind(forest_data[,1:12], forest_data$householdsize, forest_data$ageresponder
                       #  , forest_data$schooling_rep    #schooling is excluded  
  )  )
  Y.forest = regression_forest(X, Y)
  Y.hat = predict(Y.forest)$predictions
  
  cf = causal_forest(X , Y, W,
                     Y.hat = Y.hat,
                     # W.hat = W.hat,
                     tune.parameters = "all"
  )
  
  
  tau.hat = predict(cf)$predictions
  ATE_CF = as.numeric (average_treatment_effect(cf,target.sample = "overlap")[1])
  Bias_ATE_CF = (ATE_CF - True_ATE)
  
  return(cbind(ATE_CF,Bias_ATE_CF ))
}



ori_data <- read_csv("20140162_Data_updated.csv")

for (n in c(500,1000,5000,10000)) {
  result_omit  <- as.data.frame (NULL)
  
  for (i in 1:100) {
    print(i)
    set.seed(i)
    N=n
    
    data_sim_omit<- gen_dataset(N,ori_data)
    data_sim_omit<- gen_treatment_outcome_omit(data_sim_omit)
    RDD_result_omit <- RDD_omitted(data_sim_omit)
    cf_result_omit<- forest_omitted(data_sim_omit)   
    result_omit_temp <- cbind(True_ATE= True_ATE ,sample_size=N , RDD_result_omit,cf_result_omit)
    result_omit<- rbind(result_omit, result_omit_temp)
  }
  
  write.csv(result_omit, file=paste("result_omit_",n ,"sample.csv" , sep=""))
  
}


######merge all RDD anf CF results into one statistics file
#result for section  5.3
rm(list = ls())
# merge resultant file into one file 
result_omit_500sample <- read_csv("result_omit_500sample.csv")
result_omit_1000sample <- read_csv("result_omit_1000sample.csv")
result_omit_5000sample <- read_csv("result_omit_5000sample.csv")
result_omit_10000sample <- read_csv("result_omit_10000sample.csv")

result_omit <- rbind(result_omit_500sample,result_omit_1000sample,result_omit_5000sample ,result_omit_10000sample)
result_omit <- na.omit(result_omit) 


result_omit$Bias_ATE_RDD_per <- abs(result_omit$Bias_ATE_RDD)/ result_omit$True_ATE
result_omit$Bias_ATE_CF_per <- abs(result_omit$Bias_ATE_CF)/ result_omit$True_ATE


aggregate_result_omit <- as.data.frame(NULL)
stat_result_omit <- as.data.frame(NULL)
for (n in c(500,1000,5000,10000)) {
  print(n)
  
  result_subset = subset(result_omit,(sample_size==n))
  
  print(result_subset)
  result_subset$True_ATE_average <- mean(result_subset$True_ATE)
  result_subset$RMSE_ATE_RDD <- (result_subset$ATE_RDD-result_subset$True_ATE_average)^2/nrow(result_subset)
  result_subset$RMSE_ATE_CF <- (result_subset$ATE_CF-result_subset$True_ATE_average)^2/nrow(result_subset)
  # 
  aggregate_result_omit <- rbind(aggregate_result_omit,result_subset)
  stat_result_omitted_temp <- cbind(True_ATE = mean(result_subset$True_ATE),sample_size=n
                                    ,Bias_ATE_RDD_per_mean = mean(result_subset$Bias_ATE_RDD_per)
                                    ,Bias_ATE_CF_per_mean = mean(result_subset$Bias_ATE_CF_per)
                                    ,RMSE_ATE_RDD_sum = sum(result_subset$RMSE_ATE_RDD)
                                    ,RMSE_ATE_CF_sum = sum(result_subset$RMSE_ATE_CF)
  )
  
  stat_result_omit <- rbind(stat_result_omit,stat_result_omitted_temp)
}

t_stat_result_omit <- transpose(stat_result_omit)

# get row and colnames in order
colnames(t_stat_result_omit) <- c("Sample size =500","Sample size =1000","Sample size =5000","Sample size =10000")
rownames(t_stat_result_omit) <- colnames(stat_result_omit)
write.csv(t_stat_result_omit, file="stat_result_omit.csv")

