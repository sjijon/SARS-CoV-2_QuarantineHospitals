####
#### ANRS Project - NOSOCOVID 
#### MESuRS Lab - Cnam
#### Paris, 2021
####
#### Generate a summary table of the infections in both hospitals

####
#### 0. SETUP ################################################################
####

## CLEAR ALL
rm(list = ls())

## PACKAGES 
library(readxl)
library(tidyverse)
library(geepack)

## USER-DETERMINED OPTIONS

## Uncomment which data to read 
## Save results?
SAVE_RES = "YES";
# SAVE_RES = "NO";

####
#### 1. BUILD TABLE ############################################################
####

SummaryGEEGLM_v2 <- function(fit, alpha=.05, dig=2, p.dig=4){
    # output a summary data frame of GEE results
    # fit is a fitted geese() object. dig is number of digits to report.
    zq       <- qnorm(1-alpha/2)
    estimate <- fit$coefficients[-1]
    lower    <- fit$coefficients[-1] - zq*coef(summary(fit))[,"Std.err"][-1]
    upper    <- fit$coefficients[-1] + zq*coef(summary(fit))[,"Std.err"][-1]
    
    p <- coef(summary(fit))[,"Pr(>|W|)"][-1]
    p <- round(p, digits=p.dig)
    #Sig <- ifelse(p<0.05, '*', ifelse(p<0.01,'**', ifelse(p<0.001, '***', '')))
    RR       <- round(exp(estimate), dig) #incidence rate ratio
    RR.lower <- round(exp(lower), dig)
    RR.upper <- round(exp(upper), dig)
    
    return(data.frame(cbind(RR, paste0(RR.lower, "-", RR.upper) , p)))	
}

####
#### 2. READ DATA ############################################################
####

Data = readRDS(file = "Data/LineList_InfectedByUnit.rds")

colnames(Data) = c("Hospital","Unit","Days","Infected")

head(Data)


####
#### 3. CRUDE RATES ############################################################
####

## Total number of infections
sum(Data$Infected[Data$Hospital=="Hosp1"]) + sum(Data$Infected[Data$Hospital=="Hosp3"]) + sum(Data$Infected[Data$Hospital=="Hosp2"])
sum(Data$Infected[Data$Hospital=="Hosp1"])
sum(Data$Infected[Data$Hospital=="Hosp3"])
sum(Data$Infected[Data$Hospital=="Hosp2"])


## Total number of days of study
sum(Data$Days[Data$Hospital=="Hosp1"])
sum(Data$Days[Data$Hospital=="Hosp3"])
sum(Data$Days[Data$Hospital=="Hosp2"])


## Estimates rates
Rate_Hosp3 = sum(Data$Infected[Data$Hospital=="Hosp3"]) / sum(Data$Days[Data$Hospital=="Hosp3"]);100*Rate_Hosp3
sd_Hosp3 = sqrt( sum(Data$Infected[Data$Hospital=="Hosp3"])/(sum(Data$Days[Data$Hospital=="Hosp3"])^2))
100*(Rate_Hosp3 - 1.96*sd_Hosp3)
100*(Rate_Hosp3 + 1.96*sd_Hosp3)

Rate_Hosp1 = sum(Data$Infected[Data$Hospital=="Hosp1"]) / sum(Data$Days[Data$Hospital=="Hosp1"]);100*Rate_Hosp1
sd_Hosp1 = sqrt( sum(Data$Infected[Data$Hospital=="Hosp1"])/(sum(Data$Days[Data$Hospital=="Hosp1"])^2))
100*(Rate_Hosp1 - 1.96*sd_Hosp1)
100*(Rate_Hosp1 + 1.96*sd_Hosp1)

Rate_Hosp2 = sum(Data$Infected[Data$Hospital=="Hosp2"]) / sum(Data$Days[Data$Hospital=="Hosp2"]);100*Rate_Hosp2
sd_Hosp2 = sqrt( sum(Data$Infected[Data$Hospital=="Hosp2"])/(sum(Data$Days[Data$Hospital=="Hosp2"])^2))
100*(Rate_Hosp1 - 1.96*sd_Hosp2)
100*(Rate_Hosp1 + 1.96*sd_Hosp2)

IRR_1 = Rate_Hosp1/Rate_Hosp3;IRR_1
IRR_2 = Rate_Hosp2/Rate_Hosp3;IRR_2


### write results as table
res_rec_rate_Hosp3 = c(sum(Data$Infected[Data$Hospital=="Hosp3"]),
                       round(sum(Data$Days[Data$Hospital=="Hosp3"])),
                       signif(100*Rate_Hosp3, 3),
                       paste(signif(100*(Rate_Hosp3 - 1.96*sd_Hosp3), 2), signif(100*(Rate_Hosp3 + 1.96*sd_Hosp3), 3), sep="-"))
res_rec_rate_Hosp1 = c(sum(Data$Infected[Data$Hospital=="Hosp1"]),
                       round(sum(Data$Days[Data$Hospital=="Hosp1"])),
                       signif(100*Rate_Hosp1, 3),
                       paste(signif(100*(Rate_Hosp1 - 1.96*sd_Hosp1), 2), signif(100*(Rate_Hosp1 + 1.96*sd_Hosp1), 3), sep="-"))

res_rec_rate_Hosp2 = c(sum(Data$Infected[Data$Hospital=="Hosp2"]),
                       round(sum(Data$Days[Data$Hospital=="Hosp2"])),
                       signif(100*Rate_Hosp2, 3),
                       paste(signif(100*(Rate_Hosp2 - 1.96*sd_Hosp2), 2), signif(100*(Rate_Hosp2 + 1.96*sd_Hosp2), 3), sep="-"))

res_rec_rate = rbind(res_rec_rate_Hosp3, res_rec_rate_Hosp1, res_rec_rate_Hosp2)
row.names(res_rec_rate)= c("Hosp3", "Hosp1", "Hosp2")
res_rec_rate


var_vec = c("Hospital", "Unit")
res_table = NULL
for (var in var_vec){
    levs = unique(Data[,var])
    for(cla in levs){
        Events = sum(Data[Data[,var]==cla, "Infected"])
        PD = sum(Data[Data[,var]==cla, "Days"])
        Rate = Events/PD
        sd_rate = sqrt( sum(Data[Data[,var]==cla, "Infected"])/(sum(Data[Data[,var]==cla, "Days"]) )^2 )
        IC_low = (Rate - 1.96*sd_rate)
        IC_sup = (Rate + 1.96*sd_rate)
        res = c(Events, round(PD), signif(100*Rate, 3),
                paste(signif(100*IC_low, 2), signif(100*IC_sup, 3), sep="-"))
        res_table = rbind(res_table, res)
        rownames(res_table)[nrow(res_table)] = paste0(var, "_",cla)
    }
}

colnames(res_table) = c("Events","PD","Rate","CI")
res_table
write.csv(res_table, file = "Table 2 - Observed rates/Output/crude_rates.csv")



####
#### 4. ADJUSTED POISSON ####################################################
####

adj_gee2 = geeglm(Infected~Hospital+Unit, 
                  offset = log(Days), 
                  data = Data, 
                  family = poisson(link='log'),
                  id = 1:nrow(Data), corstr = "exchangeable")
SummaryGEEGLM_v2(adj_gee2)
lnL_mod_tot = sum( adj_gee2$y*log(adj_gee2$fitted.values)- adj_gee2$fitted.values)
ddl_mod_tot = adj_gee2$df.residual

res_table = SummaryGEEGLM_v2(adj_gee2)

write.csv(res_table, file = "Table 2 - Observed rates/Output/adjusted_poisson.csv")


####
#### 5. GLOBAL p-VALUES ####################################################
####
# pval DAA
adj_gee_w_DAA = geeglm(Infected~Days+Unit, offset = log(Days), data=Data, family=poisson(link='log'),
                       id=1:nrow(Data), corstr = "exchangeable")

1-pchisq(2*(lnL_mod_tot - sum( adj_gee_w_DAA$y*log(adj_gee_w_DAA$fitted.values)- adj_gee_w_DAA$fitted.values)), 
         adj_gee_w_DAA$df.residual-ddl_mod_tot)

