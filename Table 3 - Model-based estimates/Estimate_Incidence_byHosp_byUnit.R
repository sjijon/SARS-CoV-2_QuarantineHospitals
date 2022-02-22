####
#### ANRS Project - NOSOCOVID 
#### MESuRS Lab - Cnam
####
#### Paris, 2022
####
#### Incidence estimation among healthcare workers 
#### by hospital unit


#### 0. SETUP ################################################################

## CLEAR ALL
rm(list = ls())

## PACKAGES 
library(tidyverse)
library(dplyr)
library(FME)


## USER-DETERMINED OPTIONS

# Hospitals
HOSP = "Hosp1"
# HOSP = "Hosp2"
# HOSP = "Hosp3"

## Affectation
UNIT = "ICU"
# UNIT = "NoICU"

## Save results?
# SAVE_RES = "YES";
SAVE_RES = "NO";

## View objects?
# VIEW_RES = "YES";
VIEW_RES = "NO";

## Run mode
RUN_MODE = "Run"
# RUN_MODE = "Load" # Load results presented in the Article


#### 1. READ DATA ############################################################

# Recruitments and infections by shift
Recruitments = readRDS(file = paste0("Data/Recruitments_",HOSP,".rds")) 
Infections = readRDS(file = paste0("Data/Infections_",HOSP,".rds")) 

# Shifts
if(HOSP=="Hosp1"){
    Shifts = readRDS(file = "Data/Shifts_Hosp1.rds") %>% 
        mutate(Hospital = HOSP) %>%
        select(Hospital,Duration,Susceptible,Infected)
    
    # Shift length
    d = 14
}else if(HOSP=="Hosp3"){
    Shifts = readRDS(file = "Data/Shifts_Hosp3.rds") %>% 
        mutate(Hospital = HOSP) %>%
        select(Hospital,Duration,Susceptible,Infected)
    
    # Shift length
    d = 7
}else if(HOSP=="Hosp2"){
    Shifts = readRDS(file = "Data/Shifts_Hosp2.rds") %>% 
        mutate(Hospital = HOSP) %>%
        select(Hospital,Duration,Susceptible,Infected)
    
    # Shift length
    d = 14
}

# Redefine susceptible and infected in Shifts
if(UNIT == "ICU"){
    Shifts$Susceptible = Recruitments$ICU
    Shifts$Infected = Infections$ICU
} else if(UNIT == "NoICU"){
    Shifts$Susceptible = Recruitments$NoICU
    Shifts$Infected = Infections$NoICU
}

NumShifts = dim(Shifts)[1] # Number of Shifts

NumShifts
sum(Shifts$Infected)

## Observed infections, susceptible, days until test
Obs = NULL
j=1 # Auxiliary counter
for(k in seq(1,NumShifts,1)){
  for (i in 1:Shifts$Duration[k]){
    Obs$Shift[j] = k
    Obs$Hospital[j] = Shifts$Hospital[k]
    Obs$Day[j] = i
    Obs$DaysTillTest[j] = Shifts$Duration[k]-i
    Obs$Susceptible[j] = Shifts$Susceptible[k]
    if (i<Shifts$Duration[k]){
      Obs$Infected[j] = 0
    }else if (i==Shifts$Duration[k]){
      Obs$Infected[j] = Shifts$Infected[k]
    }
    j=j+1
  }
}
Obs = as_tibble(Obs)
# view(Obs)


## Total number of observed infections
Total_Obs_Inf = sum(Shifts$Infected)
Total_Obs_Inf


## Time-of-infection dependent sensitivity PCR
Sensi = readRDS(file = "Data/DailyPCRSensitivity.rds")
Sensi

#### 2. MODEL ###############################################################


#### INFECTION MODEL ##############################


## FUNCTION: Simulated number of infections
siminci_fun = function(lambda){
  # for a lambda value, siminci_fun returns a list
  # for each element of the list we have a simulation by shift
  # one simulation = vector of daily infections which takes into account the daily number of susceptible HCW
  SimInci = list()
  for(k in 1:NumShifts){
    Susceptible = Shifts$Susceptible[k]
    DailyInci = NULL
    for (i in 1:Shifts$Duration[k]){
      NewInf = rbinom(1, Susceptible, p=lambda)
      Susceptible = Susceptible - NewInf # update nb of Susceptible
      DailyInci = c(DailyInci, NewInf)
    }
    SimInci[[k]] = DailyInci
  }
  return(SimInci)
}

#### OBSERVATION MODEL ############################

sim_captinf = function(Obs, SimInci){
  Obs = mutate(Obs, 
               SimInci = unlist(SimInci),
               Sensi = Sensi$Sensi[Obs$DaysTillTest+1],
               CapturedInf = NA)
  for (i in 1:nrow(Obs)){
    Obs$CapturedInf[i] = rbinom(1, Obs$SimInci[i], prob = Obs$Sensi[i])
  }
  return(Obs)
}

## WRAP RESULTS ################################
## As in the Shifts table
wrapmodel = function(model){
  ModelShifts = Shifts
  ModelShifts$Infected = NA
  for (k in 1:NumShifts){
    sub = filter(model,Shift==k)
    ModelShifts$Infected[k] = sum(sub$CapturedInf)
  }
  ModelShifts$PropPos = ModelShifts$Infected / Shifts$Susceptible
  return(ModelShifts)
}


## TOTAL MODEL ################################
model = function(lambda){
  # registerDoParallel(cores=nbcores)
  SimInci = siminci_fun(lambda)
  model = sim_captinf(Obs, SimInci)
  sim_shift = wrapmodel(model)
  return(sim_shift)
}


#### TEST model functions ----------------------
lambda = 0.01 # Daily probability of getting infected
SimInci = siminci_fun(lambda); SimInci
Test_ObsModel = sim_captinf(Obs, SimInci); Test_ObsModel
Test_Model = model(lambda); Test_Model
Test_TotalInf = sum(Test_Model$Infected); Test_TotalInf
sum(Test_ObsModel$Infected)

#### 3. MCMC SIMULATIONS ####################################################


## Using modMCMC of the FME package
## modMCMC(f, p, ..., jump = NULL,  lower = -Inf, upper = +Inf, 
# prior = NULL, var0 = NULL, wvar0 = NULL, n0 = NULL, 
# niter = 1000, outputlength = niter, burninlength = 0, 
# updatecov = niter, covscale = 2.4^2/length(p),
# ntrydr = 1, drscale = NULL, verbose = TRUE)

## f = -2 * log-likelihood
## [required by the FEM package]
loglike = function(lambda){
  mo = model(lambda)
  loglik = dbinom(sum(Shifts$Infected), sum(Shifts$Susceptible), prob=sum(mo$Infected)/sum(mo$Susceptible), log=TRUE)
  return(-2*loglik)
}

## -2 * log(prior)
## [required by the FEM package]
prior = function(lambda){
  ln_prior = dunif(lambda, min=0, max=1, log=TRUE)
  return(-2*ln_prior)
}

N_mcmc = 11000
p = 0.01


if (RUN_MODE=="Run"){
  Sim_MCMC = modMCMC(loglike,
                     p = p,
                     jump = p/10,
                     prior = prior,
                     niter = N_mcmc,
                     burninlength = 1000,
                     updatecov = 1000) # For now, default value
  
  lambda = Sim_MCMC$pars
}else if (RUN_MODE=="Load"){
  ## Load results
  Sim_MCMC = readRDS(file = paste0("Table 3 - Model-based estimates/Output/Incidence by unit/SimMCMC_",UNIT,"_",HOSP,".rds"))
  lambda = Sim_MCMC$pars[,1]
}

#### EXAMINE MCMC CHAINS #####################

summary(Sim_MCMC)

pdf(file = paste0("Table 3 - Model-based estimates/Output/Incidence by unit/lambda_",UNIT,"_",HOSP,".pdf"))
par(mfrow = c(1,2))
## Plot trace
plot(lambda,
     type = "l",
     xlab = "iteration",
     ylab = "",
     main = "lambda (trace)")
## Plot histogram
hist(lambda, freq=TRUE,
     xlab = "lambda",
     main = "Histogram of lambda")
dev.off()

#### PROB OF INFECTION DURING A SHIFT ##############################################

summary(lambda)
lambda_95CI = c(quantile(lambda, probs = 0.025), quantile(lambda, probs = 0.975))


Prob_Inf_7dShift = 1-(1-lambda)^7
Prob_Inf_7dShift_95CI = c(quantile(Prob_Inf_7dShift, probs = 0.025), quantile(Prob_Inf_7dShift, probs = 0.975))
summary(Prob_Inf_7dShift)


#### 4. PRINT AND SAVE RESULTS #####################################################

## View tables and objects
if (VIEW_RES=="YES"){
  view(Obs)
  view(Sensi)
}

## Print results
print(paste0("Results for ",UNIT," at ",HOSP," Hospital(s)"))

print(paste0("Total observed infected: ",Total_Obs_Inf))

print(paste0("Accepted ",Sim_MCMC$count["num_accepted"]," of ",N_mcmc,
             "(~",round(100*Sim_MCMC$count["num_accepted"]/N_mcmc,0),"%)"))

print(paste0('lambda: ',100*round(median(lambda),4),
             ' (95% CrI: [',100*round(lambda_95CI [1],4),',',100*round(lambda_95CI [2],4),'])'))

print(paste0('Prob infection for a 7-day shift: ',100*round(median(Prob_Inf_7dShift),4),
             '% (95% CrI: [',100*round(Prob_Inf_7dShift_95CI [1],4),'%-',100*round(Prob_Inf_7dShift_95CI [2],4),'%])'))


## Save results
if (SAVE_RES=="YES"){
  dir = "Table 3 - Model-based estimates/Output/Incidence by unit/"
  saveRDS(lambda, paste0(dir, "lambda_",UNIT,"_",HOSP,".rds"))
  saveRDS(Sim_MCMC, paste0(dir, "SimMCMC_",UNIT,"_",HOSP,".rds"))
}

print(paste0("Results for ",UNIT," at ",HOSP))