setwd(getwd())

library(igraph)
library(readr)
library(readxl)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(scales)
options(digits = 4)
options(max.print = 100)

load("RData/00-simulation-envir.RData")

w <- 0
tre <- 0
n_10 <- c("ER_10", "SB_10", "RF_10", "RGG_10")
n_40 <- c("ER_40", "SB_40", "RF_40", "RGG_40")
n_10_partial <- c("ER_partial", "SB_partial", "RF_partial", "RGG_partial")
n_10_all <- c(n_10,  n_10_partial)

# Dataframe total construction (df_tot)  ---------------------------

df_tot <- data.frame()
for (type in n_10_all) {
  # import data of n=10
  df <- as.data.frame(read.csv(
    paste0("../input/raw_data/data", type, ".csv"),
    header = TRUE,
    sep = ";"
  ))
  
  # remove clicks of "history checking" and period 0
  df <- df[df$choice != 0 & df$time != 0, ]
  # remove duplicates
  df <- distinct(df, id_subj, id_grp, 
                 grp, round, period, 
                 sig, choice, .keep_all = TRUE)
  # insert type column
  df$type <- type
  
  df_tot <- rbind(df_tot, df)
}
for (i in c(1:4)) {
  for (type in n_40) {
    # import data of n=40 across 4 rounds
    df <-
      as.data.frame(read.csv(
        paste0("../input/raw_data/data", type, "_", i, ".csv"),
        header = TRUE,
        sep = ";"
      ))
    
    # remove clicks of "history checking" and period 0
    df <- df[df$choice != 0 & df$time != 0, ]
    # remove duplicates
    df <- distinct(df, id_subj, id_grp,
                   grp, round, period,
                   sig, choice, .keep_all = TRUE)
    # insert type column
    df$type <- type
    df$grp <- i
    
    df_tot <- rbind(df_tot, df)
  }
}
df_tot <- rename(df_tot, group = grp)
# recode empty signal as NA (for partial signals)
df_tot$sig[df_tot$sig == 0] <-NA 

# Recode neighbours of agents from string to vector
df_tot$links <- gsub("\\[", "", df_tot$links)
df_tot$links <- gsub("\\]", "", df_tot$links)
df_tot$links <-
  lapply(strsplit(as.character(df_tot$links), ","), as.numeric)

# max number of group, rounds, periods
max_grp <- max(df_tot$group)
max_rnd <- max(df_tot$round)
max_prd <- max(df_tot$period)

# Import simulation to compare with data ---------------------------

# Import DeGroot simulation results
df_simulation <-
  read_csv(
    paste0("dataframe/df_simulation.csv"),
    col_types = cols(X1 = col_skip())
  )
df_simulation <- as.data.frame(df_simulation)

# Import Bayesian simulation results
df_simulation2 <-
  read_csv(
    paste0(
      "C:/Users/tyyto/Desktop/Social Learning/Simulation_DeGroot/Simulation_data_Bayesian,n_10.csv"
    ),
    col_types = cols(X1 = col_skip())
  )
df_simulation2 <- as.data.frame(df_simulation2)


# % of those who guesses round majority signal (df_summary)  ---------------------------

df_summary <- data.frame(
    group = NA, round = NA, period = NA,
    mean = NA, type = NA, good = NA,
    initsig = NA
  )
sum_sig <- 0
num_sig <- 0

for (ntype in c(n_10, n_10_partial, n_40)) {
  # Remember that Bayesian simulation only available on n_10
  for (grp in c(1:max_grp)) {
    for (rnd in c(1:max_rnd)) {
      df_temp <- filter(df_tot, type == ntype, group == grp, round == rnd)
      
      # treat majority of signal as truth
      sum_sig_temp <- #add the sum of signals (1 being wrong signal, 2 being correct, 0 being no signal)
        sum(as.numeric(lapply(df_temp$sig, sum)), na.rm = TRUE)
      num_sig_temp <- #count the number of signals
        sum(as.numeric(df_temp$sig != 0), na.rm = TRUE)
      
      num_sig <- num_sig + num_sig_temp #total number of signals observed
      
      truth <- round(sum_sig_temp / num_sig_temp) #majority of signals as truth
      
      if (truth <= 1.5) {
        sum_sig <- sum_sig + 3 * num_sig_temp - sum_sig_temp
        Ns <-
          3 - sum_sig_temp / num_sig_temp - 1 #Set number of signals matching majority signal
      }
      else {
        sum_sig <- sum_sig + sum_sig_temp
        Ns <-
          sum_sig_temp / num_sig_temp - 1 #Set number of signals matching majority signal
      }
      
      # check if the set of signal used in group-round is bad
      y <- (grp - 1) * max_rnd + rnd #signal y used in group grp, round rnd
      z <- #whether signal y is reached correct consensus, is "good"
        filter(df_simulation, type == ntype, signal == y, period == 12)$good
      
      for (prd in c(1:max_prd)) {
        df_temp2 <- df_temp[df_temp$period == prd, ]
        x <- mean(df_temp2$choice == truth) #record the mean of guess across the network at period prd
        df_summary <-
          rbind(df_summary, list(grp, rnd, prd, x, ntype, z, Ns))
      }
    }
  }
}
df_summary <- df_summary[-1, ]
df_summary$type <- as.factor(as.character(df_summary$type))

# Export df_tot
write.csv(df_summary,"dataframe/df_summary.csv", row.names = FALSE)

# # % of those who guesses their initial signal (df_summary2)  ---------------------------
# df_summary2=data.frame(group = NA, round = NA, period = NA, mean = NA, type =NA)
# for (ntype in n_10){
#   for (grp in c(1:max_grp)){
#     for (rnd in c(1:max_rnd)){
#       df_temp <- filter(df_tot, type==ntype, group==grp, round==rnd)
#
#       #treat initial signal as truth
#       sum_sig_temp <- as.numeric(lapply(df_temp$sig,sum))
#       num_sig_temp <- as.numeric(lapply(df_temp$sig,length))
#
#       df_temp$init_signal <- round(sum_sig_temp/num_sig_temp)
#       df_temp<- df_temp[complete.cases(df_temp$init_signal), ]
#
#       for (prd in c(1:max_prd)){
#         x <- mean(df_temp[df_temp$period==prd,]$choice == df_temp[df_temp$period==prd,]$init_signal)
#         df_summary2<-rbind(df_summary2,list(grp,rnd,prd,x,ntype))
#       }
#     }
#   }
# }
# df_summary2 <- df_summary2[-1,]
# df_summary2$type <- as.factor(as.character(df_summary2$type))

# For partial info, what percentage overall received a signal? (expected 2/3)
(nrow(df_tot) - num_sig) / (max_grp * max_prd * max_rnd * 10 * 2)

# Total sum of true input/signals/ number of signals (expected 0.7)
signals <- 
  as.matrix(read_csv(paste0("../input/signals/Signals_n_10.csv"), col_names = FALSE))
mean(signals, na.rm = TRUE) - 1

signals <- 
  as.matrix(read_csv(paste0("../input/signals/Signals_n_10_partial.csv"), col_names = FALSE))
signals[signals == 0] <- NA
mean(signals, na.rm = TRUE) - 1

signals <- 
  as.matrix(read_csv(paste0("../input/signals/Signals_n_40.csv"), col_names = FALSE))
mean(signals, na.rm = TRUE) - 1

# Saving environment for later use
save.image("RData/01-process_data-envir.RData")
