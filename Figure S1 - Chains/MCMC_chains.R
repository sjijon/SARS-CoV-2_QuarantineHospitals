####
#### ANRS Project - NOSOCOVID 
#### MESuRS Lab - Cnam
####
#### Paris, 2022
####
#### Supplementary figures

####
#### 0. SETUP ################################################################
####

## CLEAR ALL
rm(list = ls())
library(tidyverse)
library(latex2exp)
library(ggplot2)
library(gridExtra)


#### 1. READRESULTS ############################################################

params = NULL
for(HOSP in c("Hosp1","Hosp2","Hosp3")){
    Sim_MCMC = readRDS(paste0("Table 3 - Model-based estimates/Output/Incidence overall/SimMCMC_",HOSP,".rds"))
    # Sim_MCMC = readRDS(paste0("Table 3 - Model-based estimates/Output/Incidence overall/SimMCMC_",HOSP,".rds"))
    params = cbind(params,Sim_MCMC$pars[,1])
}

colnames(params) = c("Hosp1","Hosp2","Hosp3")
params = params %>% as_tibble()
params

## Plot traces

p1 = ggplot(data = params, aes(x = seq(1,nrow(params)))) +
    geom_line(aes(y = Hosp1),size=0.5) +
    labs(title = TeX("Hosp1"),
         x = "\n Iteration",
         y = TeX("Force of infection ($\\lambda$)\n")) +
    theme_minimal()

p2 = ggplot(data = params, aes(x = seq(1,nrow(params)))) +
    geom_line(aes(y = Hosp2),size=0.5) +
    labs(title = TeX("Hosp2"),
         x = "\n Iteration",
         y = TeX("Force of infection ($\\lambda$)\n")) +
    theme_minimal()

p3 = ggplot(data = params, aes(x = seq(1,nrow(params)))) +
    geom_line(aes(y = Hosp3),size=0.5) +
    labs(title = TeX("Hosp3"),
         x = "\n Iteration",
         y = TeX("Force of infection ($\\lambda$)\n")) +
    theme_minimal()


p = grid.arrange(p1, p2, p3,
                 nrow = 1,
                 top = "")

## Save plots
ggsave("Figure S1 - Chains/SupplementaryFig1_MCMCchains.png",
       plot=p, height=10, width=28, units=c("cm"), dpi=600)
