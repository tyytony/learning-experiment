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

load("RData/01-process_data-envir.RData")

# import signals (df_signal) ---------------------------

df_signal <-   data.frame(
  "typesize" = NA, "group" = NA,
  "round" = NA, "Ns" = NA
)
for (ntypelist in c("n_10", "n_10_partial", "n_40")) {
  balls <-
    as.matrix(read_csv(paste0("input/signals/Signals_", ntypelist, ".csv"), 
                       col_names = FALSE)) - 1
  balls[balls == -1] <- 0.5
  for (ntype in get(ntypelist)) {
    for (grp in c(1:max_grp)) {
      for (rnd in c(1:max_rnd)) {
        i <- (grp - 1) * max_rnd + rnd
        Ns <- ifelse(mean(balls[i, ]) < 0.5, 
                     1 - mean(balls[i, ]), mean(balls[i, ]))
        df_signal <- rbind(df_signal, c(ntype, grp, rnd, Ns))
      }
    }
  }
}
df_signal$Ns <- as.numeric(df_signal$Ns)
df_signal <- df_signal[-1, ]

# Good and Bad signals (goods, bads) ---------------------------

# Create classification of networks
# subset simulation into those that fail to reach correct consensus
df_class <- filter(df_simulation, period==12, mean<0.9) 

# Find out set of signals that would give correct consensus for all networks
# Collectively "goods" signals for all networks
good10 <-
  c(1:24)[!(c(1:24) %in% c(filter(df_class, type == "ER_10")$signal)) &
            !(c(1:24) %in% c(filter(df_class, type == "SB_10")$signal)) &
            !(c(1:24) %in% c(filter(df_class, type == "RF_10")$signal))]
good40 <-
  c(1:24)[!(c(1:24) %in% c(filter(df_class, type == "ER_40")$signal)) &
            !(c(1:24) %in% c(filter(df_class, type == "SB_40")$signal)) &
            !(c(1:24) %in% c(filter(df_class, type == "RF_40")$signal))]

goods <- data.frame(type = NA, goodbad = NA)
for (ntype in c(n_10, n_40)) {
  temp <- if (ntype %in% n_10) {list(good10)} else {list(good40)}
  goods <-
    rbind(goods, list(ntype, temp))
}
goods <- goods[-1, ]

df_goods <- df_summary[1,]
df_goods$signal <- df_goods$round+(df_goods$group-1)*6
for (gb in c("goods")){
  for (ntypelist in c("n_10","n_40")){
    df_temp <- filter(df_summary, type %in% get(ntypelist), period ==10)
    df_temp$signal <- df_temp$round+(df_temp$group-1)*6
    for (ntype in get(ntypelist)){
      # Set of group-rounds that are in the "goods"
      x <- filter(df_temp, type==ntype, signal %in% filter(get(gb), type==ntype)$goodbad[[1]])
      df_goods <- rbind(df_goods, x)
    }
  }
}
df_goods <- df_goods[-1,]
df_goods$mean <- as.numeric(df_goods$mean)

# Find out set of signals that would not give correct consensus for all networks
# Pairwise "bads" signals between ntype and ER
# ie. ER would reach correct consensus, but ntype would not
bads  <- data.frame(type = NA, goodbad = NA)
for (ntype in c(n_10, n_40)) {
  
  # get set of signals that fail to reach correct consensus in ER10 and ER40
  bad_sim_ER_sig <- c(filter(df_class, type == (
    if (ntype %in% n_10) { "ER_10" } 
    else{ "ER_40" }
  ))$signal)
  # get set of signals that fail to reach correct consensus in network ntype
  ntype_sim_sig <- c(filter(df_class, type == ntype)$signal)
  # remove signals that fail in both ER and ntype 
  # from list of bad signals in ntype
  # to form "bads" signals for ntype
  x <- unique(ntype_sim_sig[!(ntype_sim_sig %in% bad_sim_ER_sig)])
  bads <- rbind(bads, list(ntype, list(x)))
}
bads <- bads[-1, ]

# Convergence time of data (df_converge) ---------------------------

width <- c(0.1,0.2)
df_converge <- data.frame(width = NA, from_period = NA, mean = NA, type =NA, good =NA)
for (wd in width){
  for (ntype in c(n_10,n_40)){
    for (grp in c(1:max_grp)){
      for (rnd in c(1:max_rnd)){
        for (prd in c(1:12)){
          df_temp <- filter(df_summary, type==ntype, group==grp, round==rnd, period >= prd)
          # check if from period prd onwards, 
          # the difference in guesses are less than width.
          # if so, then record prd
          if (max(df_temp$mean)-min(df_temp$mean)<=wd){
            df_converge <- rbind(df_converge, c(wd, prd, mean(df_temp$mean), ntype, unique(df_temp$good)))
            break
          }
        }
      }
    }
  }
}
df_converge <- df_converge[-1,]
df_converge$from_period <- as.numeric(df_converge$from_period)

# Stable Consensus, truth, stability (df_"style") ---------------------------
# reaches "style" from period prd onwards
width <- c(0.1,0.15,0.2)

for (style in c("consensus","truth","stability")){
  assign(paste0("df_",style), data.frame(width = NA, from_period = NA, mean = NA, type =NA))
  for (wd in width){
    for (prd in c(1:12)){
      for (ntype in c(n_10_all,n_40)){
        y <- c()
        for (grp in c(1:max_grp)){
          for (rnd in c(2:max_rnd)){
            df_temp <- filter(df_summary, type==ntype, group==grp, round==rnd, period >= prd)
            
            if (style == "consensus"){
              # check if consensus to 1 or 0 with error of width
              x <- as.numeric(max(df_temp$mean)<=wd | min(df_temp$mean)>=1-wd)
            }
            if (style == "truth"){
              # check if obtain close to truth
              x <- as.numeric(min(df_temp$mean)>=1-wd)
            }
            if (style == "stability"){
              # check if the difference in guesses are less than width
              x <- as.numeric(max(df_temp$mean)-min(df_temp$mean)<=wd)
            }
            y <- c(y,x)
          }
        }
        assign(paste0("df_",style), rbind(get(paste0("df_",style)), c(wd, prd, mean(y), ntype)))
      }
    }
  }
  assign(paste0("df_",style), get(paste0("df_",style))[-1,])
  assign(paste0("df_",style), transform(get(paste0("df_",style)), type = as.factor(type)))
  assign(paste0("df_",style), transform(get(paste0("df_",style)), mean = as.numeric(mean)))
  assign(paste0("df_",style), transform(get(paste0("df_",style)), from_period = as.numeric(from_period)))
}

#Individual level consensus (df_indivcon) ---------------------------
df_indivcon=data.frame(type=NA, group=NA, round=NA, id_grp=NA, from_period = NA)
for (ntype in c(n_10,n_40)){
  for (grp in c(1:max_grp)){
    for (rnd in c(1:max_rnd)){
      for (i in c(1:if(grepl("40",ntype)){40} else{10})){
        for (prd in c(1:12)){
          df_temp <- filter(df_tot, type==ntype, group==grp, round==rnd, id_grp==i, period >= prd)
          
          # check if all guesses after period prd is the same
          if(max(df_temp$choice)-min(df_temp$choice)==0){
            df_indivcon <- rbind(df_indivcon, c(ntype, grp, rnd, i, prd))
            break
          }
        }
      }
    }
  }
}
df_indivcon <- df_indivcon[-1,]
df_indivcon$type <- as.factor(as.character(df_indivcon$type))
df_indivcon$from_period <- as.numeric(df_indivcon$from_period)


# Saving environment for later use
save.image("RData/02-analysis-envir.RData")
