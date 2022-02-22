####
#### ANRS Project - NOSOCOVID 
#### MESuRS Lab - Cnam
####
#### Paris, 2022
####
#### Number of infected HCW (captured and non-captured by routine testing)
#### Running the model for all the estimated parameters


#### 0. SETUP ################################################################
## CLEAR ALL
rm(list = ls())

library(tidyverse)
library(latex2exp)
library(ggplot2)
library(gridExtra)
library(tictoc)

# HOSP ="Hosp1"
HOSP ="Hosp2"
# HOSP ="Hosp3"

CatsLabels = list(Unit = c("ICU","NoICU"))

## USER-DETERMINED OPTIONS
## Run simulations or load saved results?
# RunMode = "Run"
RunMode = "Load" # Load results presented in the Article

## Save results?
# SAVE_RES = "YES";
SAVE_RES = "NO";

## View objects?
# VIEW_RES = "YES";
VIEW_RES = "NO";

dir = paste0("Table 3 - Model-based estimates/Output/Number of infections")

#### 1. MODEL ###############################################################

## Load model functions
source("Table 3 - Model-based estimates/model_functions.R")

#### 1.1. Infection model ##############################
## Function: Simulated number of infections
siminci_fun = function(p){
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

#### 1.2. Observation model ############################

sim_captinf <- function(Obs, SimInci){
    # From a vector of daily infections SimInci,
    # the function sim_captinf returns a vector of diagnosed infections
    # taking into account the time-dependent test sensitivity
    
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

## 1.3. WRAP RESULTS ################################
## As in the Shifts table
wrapmodel <- function(model, SimInci){
    ModelShifts = Shifts %>%
        select(Start, End, Duration) %>%
        mutate(S_ICU = Recruitments$ICU,
               S_NoICU = Recruitments$NoICU,
               S_Tot = Recruitments$Total)
    ## Infected 
    ModelShifts$X_ICU = NA
    ModelShifts$X_NoICU = NA
    for (k in 1:NumShifts){
        sub = filter(model,Shift==k)
        ModelShifts$X_ICU[k] = sum(unlist(SimInci$ICU[k]))
        ModelShifts$X_NoICU[k] = sum(unlist(SimInci$NoICU[k]))
    }
    ModelShifts$X_Tot = ModelShifts$X_ICU + ModelShifts$X_NoICU
    ## Diagnosed
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

## 1.4. TOTAL MODEL ################################
model <- function(p){
    SimInci = siminci_fun(p)
    model = sim_captinf(Obs, SimInci)
    sim_shift = wrapmodel(model,SimInci)
    return(sim_shift)
}

## Time-of-infection dependent sensitivity PCR
Sensi = readRDS(file = "Data/DailyPCRSensitivity.rds")

#### 2. RUN MODEL  #####################################################

tic("Running the model")

N_CaptInf = NULL
TotalInf = NULL

#### 2.1. READ DATA #####################

## Shifts
Shifts = readRDS(file = paste0("Data/Shifts_",HOSP,".rds"))

# Recruitments and infections by shift
Recruitments = readRDS(file = paste0("Data/Recruitments_",HOSP,".rds"))
Infections = readRDS(file = paste0("Data/Infections_",HOSP,".rds"))

Recruitments$Hospital =  Recruitments$Hospital %>% as_factor()
Infections$Hospital =  Infections$Hospital %>% as_factor()

NumShifts = nrow(Shifts)

## PD
ShiftDays = Shifts %>%
    select(Duration) %>%
    sum()


## Generate the observations table, for each CAT (ICU, NoICU, NoICU, Total)
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
head(Obs)

## Parameter estimation
SimMCMC = readRDS(file = paste0("Table 3 - Model-based estimates/Output/Incidence overall/SimMCMC_",HOSP,".rds"))
Params = SimMCMC$pars %>%
    as_tibble() %>%
    mutate(Hospital = HOSP)

summary(Params)

lambda_P50 = median(Params$lambda)
lambda_P025 = quantile(Params$lambda, probs = 0.025)
lambda_P975 = quantile(Params$lambda, probs = 0.975)
print(paste0('HOSP: ',HOSP,';  lambda, median (95% CrI): ',100*round(lambda_P50,4),
             ' (',100*round(lambda_P025,4),'-',100*round(lambda_P975,4),')'))

#### TEST model ######
p = c(lambda = 0.01)

Test_SimInci = siminci_fun(p); Test_SimInci
Test_ObsModel = sim_captinf(Obs, Test_SimInci); Test_ObsModel
Test_Model = model(p); Test_Model

## Total observed infections
sum(Infections[CatsLabels$Unit])

## Number of observed (diagnosed) infections in simulation
Test_Model %>% select(starts_with("I_")) %>% sum()

## Number of infections in simulation
Test_Model %>% select(starts_with("X_")) %>% sum()

#### 2.2. RUN MODEL #####################

print(paste('Running the model for Hosp:',HOSP))

if (RunMode == "Run"){
    
    N_Inf_HOSP = NULL
    N_CaptInf_HOSP = NULL
    
    ## Run for the estimated values
    for(j in seq(1,nrow(Params),1)){
        
        # Infections
        AUX = model(c(lambda = Params$lambda[j])) %>%
            select(starts_with("X_"))
        colSums(AUX)
        N_Inf_HOSP = rbind(N_Inf_HOSP, colSums(AUX))
        
        # Diagosed infections
        AUX = model(c(lambda = Params$lambda[j])) %>%
            select(starts_with("I_"))
        N_CaptInf_HOSP = rbind(N_CaptInf_HOSP, colSums(AUX))
    }
    
    ## Build table of number of infections by hospital unit ##############
    ## All infections
    N_Inf = N_Inf_HOSP %>% 
        as_tibble() %>% 
        mutate(Hospital = HOSP) %>%
        relocate(Hospital,.before="X_ICU") 
    
    colnames(N_Inf) = c("Hospital", CatsLabels$Unit, "Total")
    
    N_Inf = N_Inf %>%
        pivot_longer(cols = c(CatsLabels$Unit,"Total"),
                     names_to = "Unit",
                     values_to = "Infected")
    N_Inf
    
    ## Captured
    N_CaptInf = N_CaptInf_HOSP %>% 
        as_tibble() %>% 
        mutate(Hospital = HOSP) %>%
        relocate(Hospital,.before="I_ICU") 
    
    colnames(N_CaptInf) = c("Hospital", CatsLabels$Unit, "Total")
    
    N_CaptInf = N_CaptInf %>%
        pivot_longer(cols = c(CatsLabels$Unit,"Total"),
                     names_to = "Unit",
                     values_to = "Infected")
    N_CaptInf
} else if (RunMode == "Load"){
    N_Inf = readRDS(paste0(dir,"/N_Inf_",HOSP,".rds"))
    N_CaptInf = readRDS(paste0(dir,"/N_CaptInf_",HOSP,".rds"))
}

## Observed
Obs_Inf = Obs %>% 
    select(starts_with("I_")) %>%
    colSums() %>%
    as_tibble_row() 
colnames(Obs_Inf) = c(CatsLabels$Unit,"Total")
Obs_Inf = Obs_Inf %>%
    pivot_longer(cols = c(CatsLabels$Unit,"Total"),
                 names_to = "Unit",
                 values_to = "Infected")
Obs_Inf

RunTime = toc()
print(paste0("Time running MCMC: ",round((RunTime$toc-RunTime$tic)/60,3), " minutes elapsed"))

#### 3. Display results ##########################################

#### Number of infections
I_ICU = N_CaptInf %>% 
    filter(Unit == "ICU")

I_NonICU = N_CaptInf %>% 
    filter(Unit == "NoICU")

I_Total = N_CaptInf %>% 
    filter(Unit == "Total")

# print(paste0("Median number of infections in ICU: ", 
#              median(I_ICU$Infected), " (95% CrI: ",
#              quantile(I_ICU$Infected, probs = 0.025),"-",
#              quantile(I_ICU$Infected, probs = 0.975),")"))
# 
# print(paste0("Median number of infections in Non-ICU: ", 
#              median(I_NonICU$Infected), " (95% CrI: ",
#              quantile(I_NonICU$Infected, probs = 0.025),"-",
#              quantile(I_NonICU$Infected, probs = 0.975),")"))

print(paste0("Median number of diagnosed infections in ", HOSP, " :", 
             median(I_Total$Infected), " (95% CrI: ",
             quantile(I_Total$Infected, probs = 0.025),"-",
             quantile(I_Total$Infected, probs = 0.975),")"))


# X_ICU = N_Inf %>% 
#     filter(Unit == "ICU")
# 
# X_NonICU = N_Inf %>% 
#     filter(Unit == "NoICU")

X_Total = N_Inf %>% 
    filter(Unit == "Total")

print(paste0("Median number of infections in ", HOSP, ": ", 
             median(X_Total$Infected), " (95% CrI: ",
             quantile(X_Total$Infected, probs = 0.025),"-",
             quantile(X_Total$Infected, probs = 0.975),")"))


#### PLOT ##########################################

ColorOne = "#00648C" # Blue
ColorTwo = "#A0A0A0" # Gray
ColorThree = "#FFAA00" # Yellow


p = ggplot(N_CaptInf, aes(x=Unit, y=Infected)) + 
    geom_boxplot(aes(fill=Unit),show.legend = FALSE) + 
    scale_fill_manual(values=c(ColorOne, ColorTwo, ColorThree)) +
    xlab("Hospital unit") +
    ylim(0, 60) +
    ylab("") +
    labs(title = "Number of diagnosed infections",
         subtitle = paste("Model ran using the estimated parameters in",HOSP)) +
    geom_point(data=Obs_Inf, aes(x=Unit,y=Infected),shape=4, size=4, show.legend =TRUE) + 
    annotate("text", x=3.3, y=55, label= "X  Observed") +
    theme_minimal() +
    theme(legend.position = "none")
p



#### SAVE RESULTS
if (SAVE_RES=="YES"){
    
    saveRDS(N_Inf, paste0(dir,"NumberInfected/",HOSP,"_Inf.rds"))
    saveRDS(N_CaptInf, paste0(dir,"NumberInfected/",HOSP,"_CaptInf.rds"))
    
    ggsave(paste0(dir,"NumberInfected/",HOSP,"_",".png"),
           plot=p, height=10, width=12, units=c("cm"), dpi=600)
}