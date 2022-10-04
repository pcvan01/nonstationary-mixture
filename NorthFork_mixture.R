######################################################################
# Load Packages, functions, set working directory
######################################################################
setwd("./Site Analysis/09405500_NorthFork")
source('mixture_lib.R')

#
library('dplyr')
library('data.table')
library('rstan') 
library('evd')
library('extRemes')
library('parallel')
library('stats4')
library('DEoptim')
library('sm')
library('fBasics')
library("tibble")

#
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

######################################################################
# Choose Station and Load Data 
stn="09405500" # Input a USGS station for analysis to load data 
######################################################################

#Load Daily Mean Peak Flow
file <- paste("./Inputs/",stn,"_obs_mxflwmean.txt",sep="")
Mxfl <- read.table(file,header=T) 

#Load Future Predictors (CMIP 1950-2099)
file <- paste("./Inputs/",stn,"_cmip_pr.txt",sep="")
cm_Pr <- read.table(file,header=T)
file <- paste("./Inputs/",stn,"_cmip_tasmin.txt",sep="")
cm_MinT <- read.table(file,header=T)
file <- paste("./Inputs/",stn,"_cmip_tasmax.txt",sep="")
cm_MaxT <- read.table(file,header=T)
file <- paste("./Inputs/",stn,"_cmip_flow.txt",sep="")
cm_Flow <- read.table(file,header=T)

#Load Historical Covariates (Livneh)
file <- paste("./Inputs/",stn,"_obs_pr.txt",sep="")
Pr <- read.table(file,header=T) 
file <- paste("./Inputs/",stn,"_obs_tasmin.txt",sep="")
MinT <- read.table(file,header=T)
file <- paste("./Inputs/",stn,"_obs_tasmax.txt",sep="")
MaxT <- read.table(file,header=T)
file <- paste("./Inputs/",stn,"_obs_liv_flow.txt",sep="")
livFlow <- read.table(file,header=T)

######################################################################
# Organize Historical Data For Fitting
######################################################################

#Create Uniform Data Set
df_list <- list(Pr, MaxT,MinT,Mxfl,livFlow) #Load any covariates that have Year Month and Day columns into this list
gg <- Reduce(function(x, y) merge(x, y,by=c("Year","Month","Day")), df_list)
df_list <- list(gg) #Load any covariates that have only Year and Month into this list, along with matrix gg from above
Daily_Data <- Reduce(function(x, y) merge(x, y,by=c("Year","Month")), df_list)
colnames(Daily_Data) <- c("Year","Month","Day","Pr","MaxT","MinT","Mxfl","Liv_Flow")

#Water Year
Daily_Data <- WaterYear(Daily_Data)

######################################################################
# Select Season and Organize Data 
season=c(4,5,6) # select season for max flow/avg covariates
pre_season=c(11,12,1,2,3) # select previous season for avg covariates
#Select the months of interest for seasonal block maxima
DATA <- getSeasonal_from_daily(Daily_Data,season,pre_season)

NiNondj12 <- read.table("./SST_Data/NDJFM_Nino12.txt",header = T)
NiN0ann12 <- read.table("./SST_Data/AMJ_Nino12.txt",header = T)
NiNondj34 <- read.table("./SST_Data/NDJFM_Nino34.txt",header = T)
NiNoann34 <- read.table("./SST_Data/AMJ_Nino34.txt",header = T)
NAOndjfm <- read.table("./SST_Data/NDJFM_NAO.txt",header = T)
NAOamj <- read.table("./SST_Data/AMJ_NAO.txt",header = T)
PDOndjfm <- read.table("./SST_Data/NDJFM_PDO.txt",header = T)
PDOamj <- read.table("./SST_Data/AMJ_PDO.txt",header = T)
AMOndjfm <- read.table("./SST_Data/NDJFM_AMO.txt",header = T)
AMOamj <- read.table("./SST_Data/AMJ_AMO.txt",header = T)

df_list <- list(DATA,NiNondj12,NiN0ann12,NiNondj34,NiNoann34,NAOndjfm,NAOamj,PDOndjfm,PDOamj,AMOndjfm,AMOamj) 
DATA <- Reduce(function(x, y) merge(x, y,by=c("Year")), df_list)

cor.test(DATA$Mxfl,DATA$Pre_Pr)
cor.test(DATA$Mxfl,DATA$Pr)
cor.test(DATA$Mxfl,DATA$NDJFM_NINO34)
cor.test(DATA$Mxfl,DATA$NDJFM_PDO)

###################################
# Stationary  GEV
###################################
Y <- DATA$Mxfl
lower <- c(0,0,-2)
upper <- c(3000,1000,2)
fit.gevS <- test_stationary_GEV(Y,lower,upper)
saveRDS(fit.gevS,"NFK_mixture_stationary_results_loc_scale_shape_final.rds")
#fit.gevS <- readRDS(file="NFK_mixture_stationary_results_loc_scale_shape_final.rds")
my_QQ1 <- my_QQ_plot(fit.gevS,DATA,stationary=TRUE)
my_dens_plot1 <- marg_distr_plot(fit.gevS,DATA,stationary = TRUE,yplotmax=.0025)
my_plots1 <- my_stationary_sim(fit.gevS,DATA)

###################################
# Climate 
##################################
# Direct Maximization
Y <- DATA$Mxfl
covs <- covs <- c("Pre_Pr","Pr","NDJFM_NINO34","NDJFM_PDO")
lower <- c(0,30,0,0,-.5,-.5,-65,-65,-65,-65,-65)
upper <- c(3000,3100,750,750,.5,.5,65,65,65,65,65)
Index_results <- test_Mixture_Models(DATA,Y,covs,lower,upper)
# 
# 
# 
# 
# ###################################
# # MLE Covariance Matrix
# ###################################
 LL_Hes <- function(par){

          # Constants
          Y <- as.vector(DATA$Mxfl)
          X <- as.vector(DATA$Pre_Pr)

          # Set up prelim variables
          likehood <- c()

          # Pull out variables
          location1 <- par[1]
          location2 <- par[2]
          scale1 <- par[3]
          scale2 <- par[4]
          shape1 <- par[5]
          shape2 <- par[6]
          beta1 <- par[7]
          beta2 <- par[8]

          for (tt in 1:length(Y)) {
             # logistic (a in state two; b in state 1)
            this_beta <- beta1 + beta2*X[tt]
            a <- exp(this_beta)/(1+exp(this_beta))
            b <- 1-a

            # multinomial transform to get TPM proabilities
            state1_prob <- b * devd(Y[tt],loc=location1, scale=scale1, shape=shape1)
            state2_prob <- a * devd(Y[tt],loc=location2, scale=scale2, shape=shape2)

            # log(sumphi) is the log-liklihood for the current point
            sumphi <- state1_prob + state2_prob
            likehood[tt] <- sumphi
            # if(pevd(-.01,loc=location1, scale=scale1, shape=shape1)>=.0001){likehood[tt] <- 1E-2}
            # if(pevd(-.01,loc=location2, scale=scale2, shape=shape2)>=.0001){likehood[tt] <- 1E-2}
          }
          totalnegative_LL <- -log(prod(likehood))

          return(totalnegative_LL)
 }

init_best_model <- as.numeric(unlist(strsplit(Index_results$Params[Index_results$AIC == min(Index_results$AIC)],split=",")))

par <- init_best_model#Approx DEoptim Solution
final_best_model <- optim(par=par,fn=LL_Hes,hessian=T,control=list(reltol=1e-20,maxit=20000))
Index_results[1,2] <- toString(final_best_model$par)
Index_results[1,3] <- final_best_model$value
AIC <- 2*final_best_model$value + 2*length(final_best_model$par)
Index_results[1,4] <- AIC

fisher_info <- solve(final_best_model$hessian)
prop_sigma <- sqrt(diag(fisher_info))
params <- data.frame(value=final_best_model$par, sigma=prop_sigma)

saveRDS(Index_results,"NFK_mixture_Index__climate_results_final_winter.rds")
#Index_results <- readRDS(file="NFK_mixture_Index__climate_results_final_winter.rds")


saveRDS(params,"NFK_mixture_Index__climate_results_final_winter_params.rds")
#params <- readRDS(file="NFK_mixture_Index__climate_results_final_winter_params.rds")


#Plot
my_plots3 <- Best_Model_Plots(Index_results,DATA,fit.gevS)
my_QQ3 <- my_QQ_plot(Index_results,DATA,stationary=FALSE)
my_dens_plot3 <- marg_distr_plot(Index_results,DATA,yplotmax=.0025)


######################################
# Comparison Plot
######################################
my_comparison(Index_results,fit.gevS,DATA)
multiplot(my_QQ1,my_QQ3,cols=2)


#####################################
# Cross Valiidate
#####################################
best_params <- as.numeric(unlist(strsplit(Index_results$Params[Index_results$AIC == min(Index_results$AIC)],split=",")))
parameter_results <- as.data.frame(matrix(nrow=length(DATA$Mxfl),ncol=length(best_params)))
samp_data <- DATA[,c(2,7)]


for(i in 4:length(DATA$Mxfl)){

    samp <- samp_data[-(i),]
    my_lower <- c(0,30,10,10,-.5,-.5,-50,-50)
    my_upper <- c(3000,3100,750,750,.5,.5,50,50)

    Y <- as.vector(samp[,1])
    XX <- as.matrix(samp[,2])
    m <- length(my_lower)-7
    boot_LL <- function(optimize){
                    # Set up prelim variables
                    likehood <- c()
                    # Pull out variables
                    location1 <- optimize[1]
                    location2 <- optimize[2]
                    scale1 <- optimize[3]
                    scale2 <- optimize[4]
                    shape1 <- optimize[5]
                    shape2 <- optimize[6]
                    beta1 <- optimize[7]
                    if(m>=1){beta2 <- optimize[8]}
                    if(m>=2){beta3 <- optimize[9]}
                    if(m>=3){beta4 <- optimize[10]}
                    if(m>=4){beta5 <- optimize[11]}
                    if(m>=5){beta6 <- optimize[12]}
                    if(m>=6){beta7 <- optimize[13]}
                    if(m>=7){beta8 <- optimize[14]}
                    if(m>=8){beta9 <- optimize[15]}


                    for (tt in 1:length(Y)) {
                       # logistic (a in state two; b in state 1)
                      this_beta <- beta1
                      for (ggg in 1:m){
                        this_var <- paste("beta",ggg+1,sep="")
                        this_beta <- this_beta+(get(this_var)*XX[tt,ggg])
                      }
                      a <- exp(this_beta)/(1+exp(this_beta))
                      b <- 1-a

                      # multinomial transform to get TPM proabilities
                      state1_prob <- b * devd(Y[tt],loc=location1, scale=scale1, shape=shape1)
                      state2_prob <- a * devd(Y[tt],loc=location2, scale=scale2, shape=shape2)

                      # log(sumphi) is the log-liklihood for the current point
                      sumphi <- state1_prob + state2_prob
                      likehood[tt] <- sumphi
                      if(pevd(-.01,loc=location1, scale=scale1, shape=shape1)>=.0001){likehood[tt] <- 1E-1000}
                      if(pevd(-.01,loc=location2, scale=scale2, shape=shape2)>=.0001){likehood[tt] <- 1E-1000}
                    }
                    totalnegative_LL <- -log(prod(likehood))
                    return(totalnegative_LL)
                }
    itermax=900
    my_optim <- DEoptim(boot_LL,lower=my_lower,upper=my_upper,
                              control=list(itermax=itermax,CR=.8,steptol=50,reltol=.00001, parallelType=2,
                                           packages=c("extRemes")))
    bestmod <- my_optim$optim$bestmem

    parameter_results[i,] <- bestmod
}

saveRDS(parameter_results,"NFK_mixture_index_climate_final_wintercrosval.rds")
# parameter_results <- readRDS(file="NFK_mixture_index_climate_final_wintercrosval.rds")

for(i in 1:nrow(parameter_results)){
  this_row <- parameter_results[i,]
  if(parameter_results[i,1] >= parameter_results[i,2]){
    parameter_results[i,1] <- this_row[1,2]
    parameter_results[i,2] <- this_row[1,1]
    parameter_results[i,3] <- this_row[1,4]
    parameter_results[i,4] <- this_row[1,3]
    parameter_results[i,5] <- this_row[1,6]
    parameter_results[i,6] <- this_row[1,5]
    parameter_results[i,7] <- this_row[1,7]*(-1)
    parameter_results[i,8] <- this_row[1,8]*(-1)
  }
}

##########################################
# Simulate from Cross Validation
##########################################
m <- 1
y <- DATA$Mxfl
X <- subset(DATA,select=Pre_Pr)
nsims=10000
predictions <- as.data.frame(matrix(0,nrow=length(DATA$Year),ncol=nsims+1))
sim_preds <- as.data.frame(matrix(nrow=nrow(DATA),ncol=2))
sim_preds$V1 <- DATA$Year
###loop

for(i in 1:nrow(parameter_results)){
    this_best_model <- as.vector(parameter_results[i,])
    location1 <- this_best_model[1]
    location2 <- this_best_model[2]
    scale1 <- this_best_model[3]
    scale2 <- this_best_model[4]
    shape1 <- this_best_model[5]
    shape2 <- this_best_model[6]
    beta1 <- this_best_model[7]
    if(m>=1){beta2 <- this_best_model[8]}
    if(m>=2){beta3 <- this_best_model[9]}
    if(m>=3){beta4 <- this_best_model[10]}
    if(m>=4){beta5 <- this_best_model[11]}
    if(m>=5){beta6 <- this_best_model[12]}
    if(m>=6){beta7 <- this_best_model[13]}
    if(m>=7){beta8 <- this_best_model[14]}
    if(m>=8){beta9 <- this_best_model[15]}
    
    this_beta <- beta1
    for (ggg in 1:m){
      this_var <- paste("beta",ggg+1,sep="")
      this_beta <- this_beta+(get(this_var)*X[i,ggg])
    }
    state_probs <- exp(this_beta)/(1+exp(this_beta))
    sim_preds[i,2] <- state_probs 
    
    predictions[i,1] <- DATA$Year[i] 
    
    for (gg in 1:nsims){
        rr <- runif(1,min=0,max=1)
        if(rr >= state_probs){
          sim_state=1
          predictions[i,gg+1]=revd(1,loc=as.numeric(location1),scale=as.numeric(scale1),shape=as.numeric(shape1))
        }
        if(rr < state_probs){
          sim_state=2
          predictions[i,gg+1]=revd(1,loc=as.numeric(location2),scale=as.numeric(scale2),shape=as.numeric(shape2))
        }
    }
}

names(predictions)[1] <- c("Year")
predictions_plot <- melt(predictions,id.vars="Year")
rownames(predictions) <- predictions[,1]
predictions[,1] <- NULL
  
plotdataframe <- DATA
years <- as.data.frame(c(min(plotdataframe$Year):max(plotdataframe$Year)))
names(years)[1] <- c("Year")
plotdataframe2 <- merge(years,plotdataframe,by="Year",all=T)



ggplot() +
    stat_summary(fun.data = f, geom="boxplot",aes(x=predictions_plot$Year,y=predictions_plot$value,color="Simulated")) +
    geom_line(aes(x=plotdataframe2$Year,y=plotdataframe2$Mxfl,color="Observed"),size=1.2) +
    xlab("Year") +
    ylab("Flow (CFS)") +
    theme(plot.title = element_text(hjust = 0.5),legend.position = "bottom")+
    scale_color_manual(name="",values = c('Simulated' = 'Black', 'Observed' = 'Red')) 


sim_preds2 <- arrange(as.data.frame(sim_preds),V2)
sim_preds3 <- add_column(sim_preds2,xx=1:nrow(sim_preds2),.after="V1")
sim_preds4 <- add_column(sim_preds3,yy=nrow(sim_preds3):1,.after="V1")
names(sim_preds4)[1] <- "Year"
plot_flow <- merge(DATA,sim_preds4)
predictions_plot2 <- merge(predictions_plot,sim_preds4)

ggplot() +
      stat_summary(fun.data = f, geom="boxplot",aes(x=predictions_plot2$xx,y=predictions_plot2$value,color="Simulated")) +
      geom_point(aes(x=plot_flow$xx,y=plot_flow$Mxfl,color="Observed"),size=2.0) +
      xlab("Risk Year") +
      ylab("Flow (CFS)") +
      theme(plot.title = element_text(hjust = 0.5),legend.position = "bottom")+
      scale_color_manual(name="",values = c('Simulated' = 'Black', 'Observed' = 'Red')) +
      scale_x_discrete(labels=as.character(sim_preds3$Year))






