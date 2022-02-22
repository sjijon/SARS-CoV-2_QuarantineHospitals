####
#### ANRS Project - NOSOCOVID
#### MESuRS Lab - Cnam
####
#### Paris, 2022
####
#### Estimating the rate of Patient-HCW transmission in Hosp1

#### 0. SETUP ################################################################
## CLEAR ALL
rm(list = ls())

## PACKAGES
library(tidyverse)
library(FME)

library(latex2exp)
library(ggplot2)
library(gridExtra)

library(tictoc)

## Hospital and categories
HOSP = "Hosp1"
# HOSP = "Hosp2"
# HOSP = "Hosp3"
CatsLabels = list(Unit = c("ICU","NoICU"))


## USER-DETERMINED OPTIONS
## Run simulations or load saved results?
RunMode = "Run"
# RunMode = "Load" # Load results presented in the Article

## Save results?
# SAVE_RES = "YES";
SAVE_RES = "NO";

## View objects?
# VIEW_RES = "YES";
VIEW_RES = "NO";

#### 2. DATA ############################################################

#### 2.1. Read data #############
d = 14

## Shifts
Shifts = readRDS(file = paste0("Data/Shifts_",HOSP,".rds"))

# Recruitments and infections by shift
Recruitments = readRDS(file = paste0("Data/Recruitments_",HOSP,".rds"))
Recruitments$Hospital =  Recruitments$Hospital %>% as_factor()

Infections = readRDS(file = paste0("Data/Infections_",HOSP,".rds"))
Infections$Hospital =  Infections$Hospital %>% as_factor()

## Number of shifts
NumShifts = nrow(Shifts)

## PD
ShiftDays = Shifts %>%
  select(Duration) %>%
  sum()

## Time-of-infection dependent sensibility PCR
Sensi = readRDS(file = "Data/DailyPCRSensitivity.rds")


if(VIEW_RES == "YES"){
  # sum(Shifts$Duration)

  view(Sensi)
  view(Shifts)
  view(Recruitments)
  view(Infections)
}

#### 3. MODEL ###############################################################

## Load model functions
source("Table 3 - Model-based estimates/model_functions.R")

## Generate the observations table, for each CAT (ICU, NoICU, Total)
for (CAT in CatsLabels$Unit){
  Obs = CreateObs(CAT)
  assign(paste0("Obs_", CAT), Obs)
}

Obs = tibble(Shift = Obs_ICU$Shift,
             ShiftDay = Obs_ICU$ShiftDay,
             DaysTillTest = Obs_ICU$DaysTillTest,
             S_ICU = Obs_ICU$Susceptible,
             S_NoICU = Obs_NoICU$Susceptible,
             S_Tot = S_ICU + S_NoICU,
             I_ICU = Obs_ICU$Infected,
             I_NoICU = Obs_NoICU$Infected,
             I_Tot = I_ICU + I_NoICU)
Obs


#### 3.1. Infection model ##############################


## Function: Simulated number of infections
siminci_fun = function(p){
  lambda = p["lambda"]

  SimInci = list()
  for (CAT in CatsLabels$Unit){
    Aux = list()
    for(k in 1:NumShifts){
      Susceptible = Recruitments[[k,CAT]]
      DailyInci = NULL
      for (i in 1:Shifts$Duration[k]){
        j = Shifts$Start[k] + i - 1
        NewInf = rbinom(1, Susceptible, p = p["lambda"])
        Susceptible = Susceptible - NewInf # update nb of Susceptible
        DailyInci = c(DailyInci, NewInf)
      }
      Aux[[k]] = DailyInci
    }
    assign(paste0("SimInci_", CAT), Aux)
  }
  SimInci = list(SimInci_ICU,
                 SimInci_NoICU)
  names(SimInci) = CatsLabels$Unit

  return(SimInci)
}

#### 3.2. Observation model ############################

sim_captinf <- function(Obs, SimInci){
  # Add simulated infections and sensitivity to table Obs
  Obs = Obs %>%
    mutate(SimInci_ICU = unlist(SimInci$ICU),
           SimInci_NoICU = unlist(SimInci$NoICU),
           SimInci_Tot = SimInci_ICU + SimInci_NoICU,
           Sensi = Sensi$Sensi[Obs$DaysTillTest + 1],
           CaptInf_ICU = NA,
           CaptInf_NoICU = NA,
           CaptInf_Tot = NA)

  # Diagnosed infections
  for (i in 1:nrow(Obs)){
    Obs$CaptInf_ICU[i] = rbinom(1, Obs$SimInci_ICU[i], prob = Obs$Sensi[i])
    Obs$CaptInf_NoICU[i] = rbinom(1, Obs$SimInci_NoICU[i], prob = Obs$Sensi[i])
  }
  Obs$CaptInf_Tot = Obs$CaptInf_ICU + Obs$CaptInf_NoICU
  return(Obs)
}

## 3.3. WRAP RESULTS ################################
## As in the Shifts table
wrapmodel <- function(model){
  ModelShifts = Shifts %>%
    select(Start, End, Duration) %>%
    mutate(S_ICU = Recruitments$ICU,
           S_NoICU = Recruitments$NoICU,
           S_Tot = Recruitments$Total)
  ModelShifts$I_ICU = NA
  ModelShifts$I_NoICU = NA
  ModelShifts$I_Tot = NA
  for (k in 1:NumShifts){
    sub = filter(model,Shift==k)
    ModelShifts$I_ICU[k] = sum(sub$CaptInf_ICU)
    ModelShifts$I_NoICU[k] = sum(sub$CaptInf_NoICU)
    ModelShifts$I_Tot[k] = sum(sub$CaptInf_Tot)
  }
  return(ModelShifts)
}

## 3.4. TOTAL MODEL ################################
model <- function(p){
  SimInci = siminci_fun(p)
  model = sim_captinf(Obs, SimInci)
  sim_shift = wrapmodel(model)
  return(sim_shift)
}


#### TEST model function ######
p = c(lambda = 0.01)

Test_SimInci = siminci_fun(p); Test_SimInci
Test_ObsModel = sim_captinf(Obs, Test_SimInci); Test_ObsModel
Test_Model = model(p); Test_Model


## Total observed infections
sum(Infections[CatsLabels$Unit])

## Number of observed (diagnosed) infections in simulation
Test_Model %>% select(starts_with("I_")) %>% sum()

## By affectation
## ICU
Infections[,"ICU"] %>% sum()
Test_Model %>% select(I_ICU) %>% sum()
## NoICU
Infections[,"NoICU"] %>% sum()
Test_Model %>% select(I_NoICU) %>% sum()



#### 4. MCMC SIMULATIONS ####################################################

## Using modMCMC of the FME package

## f = -2 * log-likelihood
## [required by the FEM package]
loglike = function(p){

  mo = model(p)

  ## log-likelihood
  loglik = dbinom(sum(Infections$ICU),sum(Recruitments$ICU),prob=sum(mo$I_ICU)/sum(mo$S_ICU),log=TRUE) + dbinom(sum(Infections$NoICU),sum(Recruitments$NoICU),prob=sum(mo$I_NoICU)/sum(mo$S_NoICU),log=TRUE)

  return(-2*loglik)
}

## -2 * log(prior)
## [required by the FEM package]
prior = function(p){
  ln_prior = dunif(p, min=0, max=1, log=TRUE)
  return(-2*ln_prior)
}

N_mcmc = 11000
N_burn = 1000
p_init = c(lambda = 0.02)

## Run or load simMCMC simulations
if (RunMode == "Run"){
  tic("Running MCMC")
  Sim_MCMC = modMCMC(loglike,
                     p = p_init,
                     lower = 0,
                     upper = 1,
                     jump = p_init/10,
                     updatecov = N_mcmc/50,
                     prior = prior,
                     niter = N_mcmc,
                     burninlength = N_burn)
  RunTime = toc()
}else if (RunMode == "Load"){
  Sim_MCMC = readRDS(file = paste0("Table 3 - Model-based estimates/Output/Incidence overall/SimMCMC_",HOSP,".rds"))
}

head(Sim_MCMC$pars)
summary(Sim_MCMC$pars)

lambda = Sim_MCMC$pars[,1]
dev_params = Sim_MCMC$SS

lambda_95CI = c(quantile(lambda, probs = 0.025), quantile(lambda, probs = 0.975))

Prob_Inf_Shift = 1-(1-lambda)^d
Prob_Inf_Shift_95CI = c(quantile(Prob_Inf_Shift, probs = 0.025), quantile(Prob_Inf_Shift, probs = 0.975))
summary(Prob_Inf_Shift)

Prob_Inf_7dShift = 1-(1-lambda)^7
Prob_Inf_7dShift_95CI = c(quantile(Prob_Inf_7dShift, probs = 0.025), quantile(Prob_Inf_7dShift, probs = 0.975))
# summary(Prob_Inf_7dShift)

Prob_Inf_14dShift = 1-(1-lambda)^14
Prob_Inf_14dShift_95CI = c(quantile(Prob_Inf_14dShift, probs = 0.025), quantile(Prob_Inf_14dShift, probs = 0.975))
# summary(Prob_Inf_14dShift)


#### EXAMINE MCMC CHAINS #####################

100*summary(Sim_MCMC)


params = tibble(lambda = as.vector(lambda))

## Plot traces

p1 = ggplot(data = params, aes(x = seq(1,nrow(params)))) +
  geom_line(aes(y = lambda),size=0.5) +
  labs(title = TeX("Trace of $\\lambda$"),
       x = "\n Iteration",
       y = TeX("$\\lambda$\n")) +
  theme_minimal()


## Plot histograms
p2 = ggplot(params, aes(x=lambda)) +
  geom_histogram(aes(y=..count..), fill="white", color="black")+
  geom_vline(aes(xintercept=mean(lambda)), color="blue",
             linetype="dashed")+
  labs(title = TeX("histogram of $\\lambda$"),
       x = TeX("$\\lambda$"),
       y = "Count")+
  theme_classic()

p = grid.arrange(p1, p2,
                 nrow = 1,
                 top = paste("Incidence rate among HCW in", HOSP))


#### 4. PRINT AND SAVE RESULTS #####################################################

## View tables and objects
if (VIEW_RES=="YES"){
  view(Obs)
  view(Sensi)
}

## Print results
print(paste0("Results for Hospital",HOSP))

print(paste0("Accepted ",Sim_MCMC$count["n/um_accepted"]," of ",N_mcmc,
             "(~",round(100*Sim_MCMC$count["num_accepted"]/N_mcmc,0),"%)"))

print(paste0('lambda median (95% CrI): ',100*round(median(lambda),4),
             ' (',100*round(lambda_95CI [1],4),'-',100*round(lambda_95CI [2],4),')'))

print(paste0('Prob infection by shift: ',100*round(median(Prob_Inf_Shift),4),
             '% (95% CrI: [',100*round(Prob_Inf_Shift_95CI [1],4),'%-',100*round(Prob_Inf_Shift_95CI [2],4),'%])'))

print(paste0('Prob infection for a 7-day shift: ',100*round(median(Prob_Inf_7dShift),4),
             '% (95% CrI: [',100*round(Prob_Inf_7dShift_95CI [1],4),'%-',100*round(Prob_Inf_7dShift_95CI [2],4),'%])'))
print(paste0('Prob infection for a 14-day shift: ',100*round(median(Prob_Inf_14dShift),4),
             '% (95% CrI: [',100*round(Prob_Inf_14dShift_95CI [1],4),'%-',100*round(Prob_Inf_14dShift_95CI [2],4),'%])'))


if (RunMode == "Run"){
  print(paste0("Time running MCMC: ",round((RunTime$toc-RunTime$tic)/60,3), " minutes elapsed"))
}

## Save results
if (SAVE_RES=="YES"){
  ## Results directory
  dir = paste0("Table 3 - Model-based estimates/Output/Incidence overall/")

  ## Save objects
  saveRDS(Obs, paste0(dir,"/Obs_",HOSP,".rds"))
  saveRDS(Sim_MCMC, paste0(dir,"/SimMCMC_",HOSP,".rds"))
  saveRDS(lambda, paste0(dir,"/lambda_",HOSP,".rds"))
  
  ## Save plots
  ggsave(paste0(dir,"/lambda_",HOSP,".png"),
         plot=p, height=10, width=14, units=c("cm"), dpi=600)
}
