####
#### ANRS Project - NOSOCOVID 
#### MESuRS Lab - Cnam
#### Paris, 2021
####
#### Generating the Epidemic curve

####
#### 0. SETUP ################################################################
####

## CLEAR ALL
rm(list = ls())

## PACKAGES 
library(readxl)
library(ggplot2)
library(tidyverse)
library(EpiCurve)
library(gridExtra)
library(lubridate)


## USER-DETERMINED OPTIONS

## Colors
ColorOne = 	"#00648C" # Blue
ColorTwo = 	"#A0A0A0" # Gray
ColorThree = "#FFAA00" # Yellow
# ColorThree = "#A01E18" # Red

## Save results?
SAVE_RES = "YES";
# SAVE_RES = "NO";

## View objects?
# VIEW_RES = "YES";
VIEW_RES = "NO";

####
#### 1. READ DATA ############################################################
####

####
#### 1.1 Hosp3 #########################################
####

Data_Hosp3 = readRDS("Data/Data_HCW_Hosp3.rds") %>%
  as_tibble()
Data_Hosp3


Inf_Hosp3 = readRDS("Data/PP_Data_Hosp3.rds") %>%
  select(-c(Hospital,Shift))
Inf_Hosp3$Unit = Inf_Hosp3$Unit %>% as.factor()
Inf_Hosp3$Date = Data_Hosp3$StartShift
Inf_Hosp3

Inf_Hosp3$Date[58] = Inf_Hosp3$Date[57] -4;

## Per-shift attack rate

Infected_Hosp3 = NULL
Susceptible_Hosp3 = NULL
k = 1
for (day in unique(Inf_Hosp3$Date[-58])){
  Shift = Inf_Hosp3 %>% filter(Date == day)
  Susceptible_Hosp3[k] = Shift %>% nrow()
  Infected_Hosp3[k] = sum(Shift$Infected)
  if (k==3){Infected_Hosp3[k] = Infected_Hosp3[k] +1}
  k=k+1
}
PerShiftAttackRate_Hosp3 = Infected_Hosp3/Susceptible_Hosp3

print("Per-shift attack rate")
print(format(PerShiftAttackRate_Hosp3,digits=2))


####
#### 1.2 Hosp 1  ######################################
####

Data_Hosp1 = readRDS("Data/HospShifts_ReInfCat_Hosp1.rds")
StartShift = readRDS("Data/Dates_Hosp1.rds")

Inf_Hosp1 = Data_Hosp1 %>% filter(Status == "Infected")
Rec_Hosp1 = Data_Hosp1 %>% filter(Status == "Recruited")
Uninf_Hosp1 = tibble(Shift = Inf_Hosp1$Shift,
               Status = rep("Uninfected",10),
               ICU = Rec_Hosp1$ICU - Inf_Hosp1$ICU,
               NoICU = Rec_Hosp1$NoICU - Inf_Hosp1$NoICU)

Data_Hosp1_ = rbind(Inf_Hosp1,Uninf_Hosp1) %>%
  mutate(Date = c(StartShift$Date,StartShift$Date),
         Hospital = "Hosp1",
         Obs_All = ICU + NoICU)
Data_Hosp1_$Status = Data_Hosp1_$Status %>% as.factor()
Data_Hosp1_

Susceptible_Hosp1 = Data_Hosp1_ %>%
  filter(Status == "Uninfected") %>%
  select(c(ICU,NoICU)) %>%
  rowSums()
Infected_Hosp1 = Data_Hosp1_ %>%
  filter(Status == "Infected") %>%
  select(c(ICU,NoICU)) %>%
  rowSums()
PerShiftAttackRate_Hosp1 = Infected_Hosp1/Susceptible_Hosp1


Inf_Hosp1_Unit = Data_Hosp1_ %>%
  select(Date,Status,ICU,NoICU) %>%
  filter(Status=="Infected") %>%
  pivot_longer(cols = c("ICU","NoICU"),
               names_to = "Unit",
               values_to = "Obs")

Inf_Hosp1_Unit$Unit = Inf_Hosp1_Unit$Unit %>% as.factor()
Inf_Hosp1_Unit

####
#### 1.3 Hosp 2  ######################################
####

Data_Hosp2 = readRDS("Data/HospShifts_ReInfCat_Hosp2.rds")
StartShift = readRDS("Data/Dates_Hosp2.rds")
NumShifts = 9

Inf_Hosp2 = Data_Hosp2 %>% filter(Status == "Infected")
Rec_Hosp2 = Data_Hosp2 %>% filter(Status == "Recruited")
Uninf_Hosp2 = tibble(Shift = Inf_Hosp2$Shift,
                   Status = rep("Uninfected",NumShifts),
                   ICU = Rec_Hosp2$ICU - Inf_Hosp2$ICU,
                   NoICU = Rec_Hosp2$NoICU - Inf_Hosp2$NoICU)

Data_Hosp2_ = rbind(Inf_Hosp2,Uninf_Hosp2) %>%
  mutate(Date = c(StartShift$Date,StartShift$Date),
         Hospital = "Hosp2")
Data_Hosp2_ = Data_Hosp2_ %>%  
  mutate(Obs_All = Data_Hosp2_$ICU + Data_Hosp2_$NoICU)
Data_Hosp2_$Status = Data_Hosp2_$Status %>% as.factor()
Data_Hosp2_

Susceptible_Hosp2 = Data_Hosp2_ %>% 
  filter(Status == "Uninfected") %>%
  select(c(ICU,NoICU)) %>%
  rowSums()
Infected_Hosp2 = Data_Hosp2_ %>% 
  filter(Status == "Infected") %>%
  select(c(ICU,NoICU)) %>%
  rowSums()
PerShiftAttackRate_Hosp2 = Infected_Hosp2/Susceptible_Hosp2


Inf_Hosp2_Unit = Data_Hosp2_ %>%
  select(Date,Status,ICU,NoICU) %>%
  filter(Status=="Infected") %>%
  pivot_longer(cols = c("ICU","NoICU"),
               names_to = "Unit",
               values_to = "Obs")

Inf_Hosp2_Unit$Unit = Inf_Hosp2_Unit$Unit %>% as.factor()
Inf_Hosp2_Unit


## 3. Figs article manuscript ###############################
## Epidemic curves for the 3 hosps
## By unit

# Reorganize tables and 
Inf_Hosp1 = Inf_Hosp1_Unit %>% 
  mutate(Hospital = as.factor("Hosp1")) %>%
  relocate(Hospital, .after = Date)

Inf_Hosp2 = Inf_Hosp2_Unit %>% 
  mutate(Hospital = as.factor("Hosp2"))%>%
  relocate(Hospital, .after = Date)

Inf_Hosp3 = Inf_Hosp3 %>%
  mutate(Status = as.factor("Infected")) %>%
  mutate(Hospital = as.factor("Hosp3")) %>%
  select(Date,Hospital,Status,Unit,Infected) %>% 
  rename(Obs = Infected)

#Complete tables 
# with all dates and NA beyond study period
Aux_DatesHosp1 = tibble(Date = as_date(setdiff(Inf_Hosp2$Date,Inf_Hosp1$Date))) %>%
  mutate(Hospital = as.factor("Hosp1"),
         Status  = as.factor("Infected"),
         Unit = NA,
         Obs = 0)

Aux_DatesHosp2 = tibble(Date = as_date(setdiff(Inf_Hosp1$Date,Inf_Hosp2$Date))) %>%
  mutate(Hospital = as.factor("Hosp2"),
         Status  = as.factor("Infected"),
         Unit = NA,
         Obs = 0)

Aux_DatesHosp31 = tibble(Date = as_date(setdiff(Inf_Hosp1$Date,Inf_Hosp3$Date))) %>%
  mutate(Hospital = as.factor("Hosp3"),
         Status  = as.factor("Infected"),
         Unit = NA,
         Obs = 0)
Aux_DatesHosp32 = tibble(Date = as_date(setdiff(Inf_Hosp2$Date,Inf_Hosp3$Date))) %>%
  mutate(Hospital = as.factor("Hosp3"),
         Status  = as.factor("Infected"),
         Unit = NA,
         Obs = 0)


INF = rbind(Inf_Hosp1,
            Aux_DatesHosp1,
            Inf_Hosp2,
            Aux_DatesHosp2,
            Inf_Hosp3,
            Aux_DatesHosp31,
            Aux_DatesHosp32) %>%
  arrange(Date)

# Replace labels
INF = INF %>%
    mutate(Unit = str_replace(Unit, "NoICU", "Non-ICU"))

## Plot
datebreaks= c("2020-W11","2020-W13","2020-W15","2020-W17","2020-W19","2020-W21","2020-W23","2020-W25","2020-W27","2020-W29")

p1 = EpiCurve(filter(INF,Hospital=="Hosp1"),
             date ="Date",
             period="day",
             to.period = "week",
             cutvar = "Unit",
             colors = c(ColorThree,ColorOne,ColorTwo),
             freq="Obs",
             # ylabel = "Per-shift number of infections",
             title = "Hosp1") +
    scale_x_discrete(breaks = datebreaks) +
    scale_y_continuous(limits = c(0,18), 
                       expand = c(0, 0)) +
    theme(legend.position = "none") +
    labs(
        earliest_date = format(min(INF$Date, na.rm=T), format = '%W'),
        latest_date = format(max(INF$Date, na.rm=T), format = '%W')
    ) +
    annotate("text",
           x = seq(1,19,2),
           y = Infected_Hosp1[c(1,seq(3,10,1),2)] + 1.2,
           size = 3,
           angle = 90,
           label = format(PerShiftAttackRate_Hosp1[c(1,seq(3,10,1),2)],digits=1))

p2 =  EpiCurve(filter(INF,Hospital=="Hosp2"),
               date ="Date",
               period="day",
               to.period = "week",
               cutvar = "Unit",
               colors = c(ColorThree,ColorOne,ColorTwo),
               freq="Obs",
               # ylabel = "Per-shift number of infections",
               title = "Hosp2") +
    scale_x_discrete(breaks = datebreaks) +
    scale_y_continuous(limits = c(0,18), 
                       expand = c(0, 0)) +
    theme(legend.position = "none") +
    annotate("text",
           x = seq(1,2*NumShifts-1,2) + 3,
           y = Infected_Hosp2 + 1.2,
           size = 3,
           angle = 90,
           label = format(PerShiftAttackRate_Hosp2,digits=1))

p3 = EpiCurve(filter(INF,Hospital=="Hosp3"),
              date ="Date",
              period="day",
              to.period = "week",
              cutvar = "Unit",
              colors = c(ColorThree,ColorOne,ColorTwo),
              freq="Obs",
              # ylabel = "Per-shift number of infections",
              title = "Hosp3") +
    scale_x_discrete(breaks = datebreaks) +
    scale_y_continuous(limits = c(0,18), 
                       expand = c(0, 0)) +
    labs(fill="Hospital unit") +
    # scale_fill_discrete(name = "Hospital unit", 
                        # labels = c("A", "B")) +
    annotate("text", 
           x = seq(1,5,1) + 13 ,
           y = Infected_Hosp3 + 1.2,
           size = 3,
           angle = 90,
           label = format(PerShiftAttackRate_Hosp3,digits=1))

# p1
# p2
# p3

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend <- get_legend(p3)

p3 = p3 + theme(legend.position = "none")

p = grid.arrange(p1,p2,p3,legend,
                 ncol=4,
                 widths = c(6,6,6,5), 
                 heights=unit(8, "cm"))

p

#### save figure

ggsave("Figure 1 - Epidemic curves/Output/EpiCurves.png",
       plot=p, height=8, width=28, units=c("cm"), dpi=600)

ggsave("Figure 1 - Epidemic curves/Output/EpiCurves.pdf",
       plot=p, height=8, width=28, units=c("cm"), dpi=800)

## 2. Outbreaks table ###############################

weeks_Hosp1 = format(unique(Inf_Hosp1$Date),"W%W")
weeks_Hosp1 = weeks_Hosp1[order(weeks_Hosp1)]
weeks_Hosp2 = format(unique(Inf_Hosp2$Date),"W%W")
weeks_Hosp3 = unique(format(unique(Inf_Hosp3$Date),"W%W"))

AttRat = tibble(
  Hospital = c(rep("Hosp1",length(Infected_Hosp1)),rep("Hosp2",length(Infected_Hosp2)),rep("Hosp3",length(Infected_Hosp3))),
  Shift = c(seq(1,length(Infected_Hosp1),1),seq(1,length(Infected_Hosp2),1),seq(1,length(Infected_Hosp3),1)),
  Week = c(weeks_Hosp1, weeks_Hosp2,weeks_Hosp3),
  Susceptible = c(Susceptible_Hosp1,Susceptible_Hosp2,Susceptible_Hosp3),
  Infected = c(Infected_Hosp1,Infected_Hosp2,Infected_Hosp3),
  AttackRate = c(PerShiftAttackRate_Hosp1,PerShiftAttackRate_Hosp2,PerShiftAttackRate_Hosp3),
  AR = round(Infected/Susceptible,2))
AttRat

saveRDS(AttRat,"Figure 1 - Epidemic curves/Output/AttackRates")

## Show outbreaks
Outbreaks = AttRat %>% filter(AttackRate >= 0.20)
# Outbreaks = AttRat %>% filter(AttackRate >= 0.30)
Outbreaks

print(paste0("The outbreak-related infections represented the ", round(100*sum(Outbreaks$Infected)/sum(AttRat$Infected),0), 
      "% of all infections (", sum(Outbreaks$Infected),"/", sum(AttRat$Infected),")"))

