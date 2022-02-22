#### 2. Build daily observations table #############

CreateObs <- function(CAT){
## Daily observed infections, susceptible, days until test
Obs = NULL
j=1 # Auxiliary counter
for(k in seq(1,NumShifts,1)){
  for (i in 1:Shifts$Duration[k]){
    Obs$Shift[j] = k
    Obs$ShiftDay[j] = i
    Obs$DaysTillTest[j] = Shifts$Duration[k]-i
    Obs$Susceptible[j] = Recruitments[[k,CAT]] # Susceptible at the beginning of the shift
    if (i<Shifts$Duration[k]){
      Obs$Infected[j] = 0
    }else if (i==Shifts$Duration[k]){
      Obs$Infected[j] = Infections[[k,CAT]] # Infected at the end of the shift
    }
    j=j+1
  }
}

Obs = Obs %>% 
  as_tibble() 

return(Obs)
}

#### 3. MODEL ###############################################################


#### 3.2. Observation model ############################

sim_captinf <- function(Obs, SimInci){
  # From a vector of daily infections SimInci,
  # the function sim_captinf returns a vector of diagnosed infections
  # taking into account the time-dependent test sensitivity
  
  # Add simulated infections and sensitivity to table Obs
  Obs = Obs %>%
    mutate(SimInci = unlist(SimInci),
           Sensi = Sensi$Sensi[Obs$DaysTillTest + 1],
           CapturedInf = NA)
  
  # Diagnosed infections
  for (i in 1:nrow(Obs)){
    Obs$CapturedInf[i] = rbinom(1, Obs$SimInci[i], prob = Obs$Sensi[i])
  } ## On peut pas vectorizer car rbinom n'est pas fonction vectoris?e
  return(Obs)
}

## 3.3. WRAP RESULTS ################################
## As in the Shifts table
wrapmodel <- function(model){
  ModelShifts = Shifts
  ModelShifts$Susceptible = Recruitments[,CAT]
  ModelShifts$Infected = NA
  for (k in 1:NumShifts){
    sub = filter(model,Shift==k)
    ModelShifts$Infected[k] = sum(sub$CapturedInf)
  }
  # ModelShifts$PropPos = ModelShifts$Infected / ModelShifts$Susceptible
  return(ModelShifts)
}


## 3.4. TOTAL MODEL ################################
model <- function(beta_P){
  SimInci = siminci_fun(beta_P)
  model = sim_captinf(Obs, SimInci)
  sim_shift = wrapmodel(model)
  return(sim_shift)
}

#### 4. functions required by modMCMC  #############
## f = -2 * log-likelihood 
## -> deviance
loglike = function(beta_P){
  mo = model(beta_P)
  loglik = dbinom(sum(Infections[,CAT]), sum(Recruitments[,CAT]), prob=sum(mo$Infected)/sum(mo$Susceptible), log=TRUE)
  return(-2*loglik)
}

## -2 * log(prior)
prior = function(beta_P){
  ln_prior = dunif(beta_P, min=0, max=1, log=TRUE)
  return(-2*ln_prior)
}


