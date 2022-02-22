####
#### ANRS Project - NOSOCOVID 
#### MESuRS Lab - Cnam
####
#### Paris, 2022
####
#### Comparing incidence rate distributions


#### 0. SETUP ################################################################

## CLEAR ALL
rm(list=ls())

## PACKAGES 
library(tidyverse)
library(dplyr)
library(FME)
library(ggplot2)
library(gridExtra)



#### 1. By Hospital ############################################################

lambda_Hosp1 <- readRDS(file = "Table 3 - Model-based estimates/Output/Incidence overall/lambda_Hosp1.rds") %>%
  as_tibble() %>%
  rename(Risk = "value") %>%
  mutate(Hospital = "Hosp1")

lambda_Hosp2 <- readRDS(file = "Table 3 - Model-based estimates/Output/Incidence overall/lambda_Hosp2.rds") %>%
  as_tibble() %>%
  rename(Risk = "value") %>%
  mutate(Hospital = "Hosp2")

lambda_Hosp3 <- readRDS(file = "Table 3 - Model-based estimates/Output/Incidence overall/lambda_Hosp3.rds") %>%
    as_tibble() %>%
    rename(Risk = "value") %>%
    mutate(Hospital = "Hosp3")

LAMBDA = rbind(lambda_Hosp1, lambda_Hosp2, lambda_Hosp3)

mean_risk_hosp = LAMBDA %>%
  group_by(Hospital) %>%
  summarize(mean = mean(Risk))

p1 = ggplot(LAMBDA, aes(Risk, color = Hospital, fill = Hospital)) +
  geom_density(alpha = 0.2, size = 0.5) +
  geom_vline(data = mean_risk_hosp, aes(xintercept = mean, 
                                       color = Hospital), size=0.5) +
  labs(title = "By hospital") +
  theme_minimal()
p1


## Test Hosp1 vs Hosp3
ks.test(lambda_Hosp1$Risk,lambda_Hosp3$Risk)

## Test Hosp1 vs Hosp2
ks.test(lambda_Hosp1$Risk,lambda_Hosp2$Risk)

# wilcox.test(lambda_Hosp1,lambda_Hosp3)