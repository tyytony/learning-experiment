# To print all output into "output" folder
# Use ctrl-F to:
# replace "#pdf(" with "pdf("
# replace "#dev.off()" with "dev.off()"

# To print all dataframes into "dataframe" folder
# Use ctrl-F to:
# replace "#write.csv" with "write.csv"

setwd(getwd())

library(igraph)
library(readr)
library(readxl)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(scales)
library(MASS)
options(digits=4)
options(max.print=100)
source("function/Function-influence_ratio.R")
source("function/Function-Spectral_homophily.R")
source("function/Function-r_cohesiveness.R")
source("function/Function-DeGroot_once,DeGroot,Simulation.R")
source("function/Function-Bayesian.R")

w <- 0
tre <- 0
k <- 1
c <- 0

n_10 <- c("ER_10","SB_10","RF_10","RGG_10")
n_40 <- c("ER_40","SB_40","RF_40","RGG_40")
n_10_partial <- c("ER_partial","SB_partial","RF_partial","RGG_partial")
n_10_all <- c(n_10, n_10_partial)


# # Bayesian simulation for RF, 3 periods only
# periods <- 10
# g <- glist$ER_40
# s <- c(1,0,0,0,1,1,1,1,1,1,1,0,0,0,1,1,1,1,1,1,1,0,0,0,1,1,1,1,1,1,1,0,0,0,1,1,1,1,1,1)
# 
# x <- Bayesian(glist$ER_40, rbinom(40,1,0.7), periods =10)
# y <- x[[1]]
# x[[2]]

#*** Part 1 - Creates dataframes from simulations ---------------------------

# Import network (glist) ---------------------------

glist <- list()
# import n=40 network
for (ntype in c("ER", "SB", "RF", "RGG"))
{
  network = paste0("n=", 40, ",d=4_", ntype)
  g <- as.matrix(read_csv(paste0("../input/networks/", network, ".csv"), col_names = FALSE))
  colnames(g) <- 1:ncol(g)
  glist[[paste0(ntype, "_40")]] <- g
}

# import n=10 network
for (ntype in c("ER", "SB", "RF", "RGG", "ER2", "SB2", "RF2", "RGG2"))
{
  network = paste0("n=", 10, ",d=4_", ntype)
  g <- as.matrix(read_csv(paste0("../input/networks/", network, ".csv"), col_names = FALSE))
  colnames(g) <- 1:ncol(g)
  glist[[ntype]] <- g
}
names(glist) <- c(n_40, n_10_all)

# print networks
for (ntype in c(n_10,n_40)) {
  g <- as.matrix(glist[ntype][[1]])
  g <- graph_from_adjacency_matrix(g-diag(ncol(g)))
  #pdf(paste0("../output/", ntype, ".pdf"), width=9, height=9)
  plot(g, main = paste0(ntype),
       edge.arrow.size= (if(ntype %in% c("RF_10","RF_40")){0.4} else{0}),
       vertex.size=7,
       vertex.color="#4f85db",
       vertex.label.dist=3, vertex.label.cex=2,
       vertex.label.color= "black",
       vertex.label.degree=0,
       #layout= layout_in_circle,
       layout= layout_with_kk)
  #dev.off()
}

# Network statistics (network_stats) ---------------------------
statistics <- c("type","degree", "diameter", "density", "pathlength", 
                "clustering", "connected", "laplacian", 
                "cohesion", "homophily", "influence")
m <- data.frame(matrix(nrow=length(c(n_10_all,n_40)), 
                       ncol=length(statistics), 
                       dimnames=list(c(n_10_all,n_40), statistics)))

#df_clan <- read_delim("../output/Simulation_clan.csv", ";", escape_double = FALSE, col_types = cols(X1 = col_skip()), trim_ws = TRUE)

for (ntype in c(n_10_all,n_40)){
  g <- graph_from_adjacency_matrix(as.matrix(glist[ntype][[1]])-diag(ncol(glist[ntype][[1]])))
  m[ntype, "degree"] <- mean(degree(g, mode="out", loops=FALSE))
  m[ntype, "diameter"] <- diameter(g,directed= if(grepl("^RF",ntype)){TRUE} else{FALSE})
  m[ntype, "density"] <- edge_density(g)
  m[ntype, "pathlength"] <- average.path.length(g,directed=if(grepl("^RF",ntype)){TRUE} else{FALSE})
  m[ntype, "clustering"] <- transitivity(g, type="global")
  m[ntype, "connected"] <- as.numeric(is.connected(g, mode = "strong"))
  m[ntype, "laplacian"] <- as.numeric(sort(eigen(laplacian_matrix(g),symmetric=TRUE,only.values=TRUE)[1][["values"]])[2] >= sqrt(2)/2)
  m[ntype, "cohesion"] <- 0
  m[ntype, "homophily"] <- assortativity_degree(g,directed = if(grepl("^RF",ntype)){TRUE} else{FALSE})
  m[ntype, "influence"] <- Influence_ratio(g)
  m[ntype, "type"] <- ntype
}

# 0.5-cohesiveness calculated separatedly due to high computation time

#m[,"cohesion"] <- c(0.36,0.31,0.50,0.32,0.36,0.31,0.50,0.32,0.02,0.13,0.07,0.00)
m[,"cohesion"] <- c(0.029,0.125,0.115,0.039,0.029,0.125,0.115,0.039,0,0.03,0.015,0)
m$laplacian<- as.character(m$laplacian)
m$connected<- as.character(m$connected)
m[sapply(m, is.infinite)] <- 1000000
network_stats <- m
remove(m)
#write.csv(network_stats ,file = paste0("../output/Network_statistics.csv"))

# Randomly generate 1000 signals and simulate DeGroot (df_1000) ---------------------------
# for guess distribution and convergence speed histogram

set.seed(1)
df_1000 <- data.frame(con_time = NA, mean = NA, type = NA)
for (ntype in c(n_10,n_40)){
  g <- glist[ntype][[1]]
  df_temp <- Simulation(rep=1000, periods=12, g, p=0.7, w, tre, k, c)
  df_temp$type <- ntype
  print(c(ntype,
          signif(mean(df_temp$con_time, na.rm=TRUE),4),
          signif(max(df_temp$con_time, na.rm=TRUE),4)))
  df_1000  <- rbind(df_1000 ,df_temp)
}
df_1000 <- df_1000[-1,]
df_1000 <- 
  df_1000 %>% 
  group_by(type) %>% 
  mutate(avg.mean = mean(mean),
         avg.con_time = mean(con_time, na.rm=TRUE),
         type = as.factor(type))

for (ntypelist in c("n_10","n_40")){
  # print guess distribution of 1000 signals
  final <- ggplot(filter(df_1000, type %in% get(ntypelist)), aes(x=mean)) + 
    geom_histogram(breaks = seq(0,1,0.025))+
    stat_bin(breaks = seq(0,1,0.025), 
             geom="text", 
             colour="black", 
             size=3, 
             aes(label= ifelse(..count.. > 50, ..count.., ""), 
                 group=type, y=30+(..count..))) +
    #ylim(NA, 830) + 
    scale_x_continuous(breaks=seq(0,1,0.1), limits = c(-0.05, 1.05))+
    facet_wrap(~ factor(type, levels=get(ntypelist)))+
    ggtitle(paste0("DeGroot simulation 1000 signals, ",ntypelist,
                   if (w!=0|tre!=0|k!=1|c!=0){
                     paste0(",weight=",w,",
                            tremble=",tre,",
                            logistic growth rate=",k,",
                            clustering=",c)
                     })) +
    
    xlab("Average guess at period 12") + 
    ylab("Frequency of signals (out of 1000)") +
    theme_bw()+
    theme(legend.title = element_blank(), 
          strip.text.x = element_text(size = 13))
  tname <- paste0("../output/Simulation_distribution_1000sig,", ntypelist, 
                  if (w!=0|tre!=0|k!=1|c!=0){
                    paste0(",weight=",w,
                           ",tremble=",tre,
                           ",logistic growth rate=",k,
                           ",clustering=",c)
                    }, ".pdf")
  #pdf(tname, width=7, height=5)
  print(final)
  #dev.off()
  
  # print converge time of 1000 sig simulation
  final <- ggplot(filter(df_1000, type %in% get(ntypelist)), aes(x=con_time)) +
    geom_histogram(binwidth = 1, color="white", )+
    stat_bin(binwidth = 1, 
             geom="text", 
             colour="black", 
             size=3, 
             aes(label= ifelse(..count.. > 50, ..count.., ""), y=60+(..count..))) +
    geom_vline(aes(xintercept=avg.con_time, group=type), colour="red") +

    #ylim(NA, 600) +
    scale_x_continuous(breaks=seq(1,12,1), limits = c(0.5, 12.5))+
    facet_wrap(~ factor(type, levels=get(ntypelist)))+
    ggtitle(paste0("DeGroot simulation 1000 signals, ",ntypelist,
                   if (w!=0|tre!=0|k!=1|c!=0){
                     paste0(",weight=",w,
                            ",tremble=",tre,
                            ",logistic growth rate=",k,
                            ",clustering=",c)
                     })) +

    xlab("Convergence time") + 
    ylab("Frequency of signals (out of 1000)")+
    theme_bw()+
    theme(legend.title = element_blank(), strip.text.x = element_text(size = 13))
  tname <- paste0("../output/Simulation_distribution_1000sig_contime,", ntypelist, 
                  if (w!=0|tre!=0|k!=1|c!=0){
                    paste0(",weight=",w,
                           ",tremble=",tre,
                           ",logistic growth rate=",k,
                           ",clustering=",c)
                    }, ".pdf")
  #pdf(tname, width=7, height=5)
  print(final)
  #dev.off()
}

# Create DeGroot simulation of 24 signals (df_simulation) ---------------------------
nseed <- 7
set.seed(nseed)

df_simulation=data.frame(period = NA, mean = NA, converge = NA, 
                         signal = NA, type =NA, good =NA)
for (ntypelist in c("n_10","n_10_partial","n_40")){
  balls <- as.matrix(read_csv(paste0("../input/signals/Signals_",ntypelist,".csv"), 
                              col_names = FALSE))-1
  balls[balls == -1] <- 0.5

  # plot distribution of initial signal
  x <- c(rowMeans(balls,na.rm=TRUE))
  final <- ggplot(mapping = aes(x)) +
    geom_histogram(breaks = seq(0,1,0.025), colour="white")+
    stat_bin(breaks = seq(0,1,0.025), 
             geom="text", 
             colour="black", 
             size=2.5, 
             aes(label= ifelse(..count.. > 0, ..count.., ""), y=0.3+(..count..))) +
    geom_vline(aes(xintercept=mean(x)), colour="red") +

    ylim(NA, 8) +
    scale_x_continuous(breaks=seq(0,1,0.1), limits = c(-0.05, 1.05))+
    ggtitle(paste0("Signal distribution, ",ntypelist)) +

    xlab("Average signal") + 
    ylab("Frequency of networks (out of 24)")+
    theme_bw()+
    theme(legend.title = element_blank())
  tname <- paste0("../output/Signal_distribution,", ntypelist, ".pdf")
  #pdf(tname, width=7, height=5)
  print(final)
  #dev.off()
  
  # Create df_simulation
  for (ntype in get(ntypelist)){
    g <- glist[ntype][[1]]
    for (i in c(1:nrow(balls))){
      s<-balls[i,]
      
      if (mean(s)<0.5){s<-1-s} #Set majority signal as truth
      
      record<-DeGroot(g,s,12,w,tre,k,c) #Perform DeGroot for 12 periods
      
      for (prd in c(1:12)){
        if(mean(record[prd,]) >1){print(c(ntype,i,prd,mean(record[prd,])))}
        x <- mean(record[prd,])
        z <- if(prd==1){NA} else{as.numeric(all(record[prd,]==record[prd-1,]))}
        df_simulation <- rbind(df_simulation, list(prd,x,z,i,ntype,NA))
      }
    }
  }
}
df_simulation <- df_simulation[-1,] %>%
  mutate(signal = as.numeric(signal),
         period = as.numeric(period),
         type = as.factor(as.character(type)),
         )

# all details of simulation saved in csv
df_temp <- filter(df_simulation, period==12, mean>0.90)
df_temp2 <- filter(df_simulation, period==12, mean<0.10)
for (ntype in c(n_10_all,n_40)){
  for (i in c(1:24)){
    #check if the set of signal used in group-round is bad
    if (i %in% filter(df_temp, type==ntype)$signal)       {y <- "Correct consensus"} 
    else if (i %in% filter(df_temp2, type==ntype)$signal) {y <- "Incorrect consensus"} 
    else                                                  {y <- "Breakdown"}
    df_simulation$good[df_simulation$type==ntype & df_simulation$signal==i] <- y
  }
}
df_simulation$good <- factor(df_simulation$good,
                              levels = c("Correct consensus", "Incorrect consensus", "Breakdown"))

#write.csv(df_simulation,file = paste0("../output/dataframe/df_simulation.csv"))

# Create Bayesian simulation of 24 signals (df_bay_sim) ---------------------------
# nseed <- 7
# set.seed(nseed)
# 
# df_bay_sim=data.frame(period = NA, mean = NA, converge_time = NA, 
#                       signal = NA, type =NA, good =NA)
# for (ntypelist in c("n_10","n_40")){
#   balls <- as.matrix(read_csv(paste0("Signals_",ntypelist,".csv"), col_names = FALSE))-1
#   balls[balls == -1] <- 0.5
# 
#   for (ntype in get(ntypelist)){
#     g <- glist[ntype][[1]]
#     for (i in c(1:nrow(balls))){
#       s<-balls[i,]
#       
#       if (mean(s)<0.5){s<-1-s} #Set majority signal as truth
#       
#       record<-Bayesian(g, s, periods =12) #Perform Bayesian for 12 periods
#       
#       for (prd in c(1:12)){
#         #if(mean(record[prd,]) >1){print(c(ntype,i,prd,mean(record[prd,])))}
#         x <- mean(record[[1]][prd,])
#         z <- record[[2]]
#         df_bay_sim <- rbind(df_bay_sim,list(prd,x,z,i,ntype,NA))
#       }
#     }
#   }
# }
# df_bay_sim <- df_bay_sim[-1,]
# df_bay_sim$signal<- as.numeric(df_bay_sim$signal)
# df_bay_sim$type<- as.factor(as.character(df_bay_sim$type))
# df_bay_sim$period<- as.numeric(df_bay_sim$period)
# 
# # all details of simulation saved in csv
# 
# df_temp <- filter(df_bay_sim, period==12, mean>0.90)
# df_temp2 <- filter(df_bay_sim, period==12, mean<0.10)
# for (ntype in c(n_10)){
#   for (i in c(1:24)){
#     #check if the set of signal used in group-round is bad
#     if (i %in% filter(df_temp, type==ntype)$signal)         {y <- "Correct consensus"} 
#     else if (i %in% filter(df_temp2, type==ntype)$signal)   {y <- "Incorrect consensus"} 
#     else                                                    {y <- "Breakdown"}
#     df_bay_sim$good[df_bay_sim$type==ntype & df_bay_sim$signal==i] <- y
#   }
# }
# #write.csv(df_bay_sim,file = paste0("../output/dataframe/df_bay_sim.csv"))


#*** Part 2 - Summarises and prints output ---------------------------

# Plot 24sig guess distribution in histogram for ntype ---------------------------

# for (ntypelist in c("n_10","n_40")){
#   for (ntype in get(ntypelist)){
#     df_temp <- filter(df_simulation, period==12, type ==ntype)
#     assign(paste0("final_",ntype), ggplot(df_temp, aes(x=mean, fill=good)) +
#       geom_histogram(breaks = seq(0,1,0.025)) +
#       stat_bin(breaks = seq(0,1,0.025), geom="text", colour="black", size=3.5, aes(label= ifelse(..count.. > 0, ..count.., ""), group=type, y=1+(..count..))) +
#       ylim(NA, 24) +
#       scale_x_continuous(breaks=seq(0,1,0.1), limits = c(-0.05, 1.05))+
#       ggtitle(paste0("DeGroot simulation, majority signal as truth, ",ntype,if (w!=0|tre!=0|k!=1|c!=0){paste0(",weight=",w,",tremble=",tre,",logistic growth rate=",k,",clustering=",c)})) +
#       xlab("Average guess at period 12") + ylab("Frequency of networks (out of 24)")+
#       theme_bw() +
#       theme(legend.title = element_blank(), legend.justification=c(1,1), legend.position=c(0.3,0.9))
#     )
#     tname <- paste0("../output/Simulation_distribution,", ntype, if (w!=0|tre!=0|k!=1|c!=0){paste0(",weight=",w,",tremble=",tre,",logistic growth rate=",k,",clustering=",c)}, ".pdf")
#     #pdf(tname, width=7, height=5)
#     print(get(paste0("final_",ntype)))
#     #dev.off()
#   }
# }

# Plot 24sig guess distribution in histogram for ntypelist, by 4 ntype ---------------------------

# DeGroot simulation
for (ntypelist in c("n_10","n_40")){
  df_temp <- filter(df_simulation, period==12, type %in% get(ntypelist))
  final <- ggplot(df_temp, aes(x=mean)) + #add "fill=good" for classification in colour
    geom_histogram(breaks = seq(0,1,0.025))+
    stat_bin(breaks = seq(0,1,0.025), 
             geom="text", 
             colour="black", 
             size=2.5, 
             aes(label= ifelse(..count.. > 0, ..count.., ""), 
                 group=type, y=1+(..count..))) +
    #ylim(NA, 24) +
    scale_x_continuous(breaks=seq(0,1,0.1), limits = c(-0.05, 1.05))+
    facet_wrap(~ factor(type, levels=get(ntypelist)))+
    ggtitle(paste0("DeGroot simulation, majority signal as truth, ",ntypelist,
                   if (w!=0|tre!=0|k!=1|c!=0){
                     paste0(",weight=",w,
                            ",tremble=",tre,
                            ",logistic growth rate=",k,
                            ",clustering=",c)
                   }))+
    xlab("Average guess at period 12") + 
    ylab("Frequency of networks (out of 24)")+
    
      theme_bw()+
    theme(legend.title = element_blank(), 
          strip.text.x = element_text(size = 13))
  tname <- paste0("../output/Simulation_distribution,", ntypelist, 
                  if (w!=0|tre!=0|k!=1|c!=0){
                    paste0(",weight=",w,
                           ",tremble=",tre,
                           ",logistic growth rate=",k,
                           ",clustering=",c)
                  }, ".pdf")
  #pdf(tname, width=7, height=5)
  print(final)
  #dev.off()
}

# Bayesian simulation
# for (ntypelist in c("n_10")){
#   df_temp <- filter(df_bay_sim, period==12, type %in% get(ntypelist))
#   final <- ggplot(df_temp, aes(x=mean)) + 
#     geom_histogram(breaks = seq(0,1,0.025))+
#     stat_bin(breaks = seq(0,1,0.025), geom="text", colour="black", size=2.5, aes(label= ifelse(..count.. > 0, ..count.., ""), group=type, y=1+(..count..))) +
#     ylim(NA, 25) +
#     scale_x_continuous(breaks=seq(0,1,0.1), limits = c(-0.05, 1.05))+
#     facet_wrap(~ factor(type, levels=get(ntypelist)))+
#     ggtitle(paste0("Bayesian simulation, majority signal as truth, ",ntypelist,if (w!=0|tre!=0|k!=1|c!=0){paste0(",weight=",w,",tremble=",tre,",logistic growth rate=",k,",clustering=",c)})) +
#     
#     xlab("Average guess at period 12") + ylab("Frequency of networks (out of 24)")+
#     theme_bw()+
#     theme(legend.title = element_blank(), strip.text.x = element_text(size = 13))
#   tname <- paste0("../output/Simulation_distribution,", ntypelist, ",Bayesian", if (w!=0|tre!=0|k!=1|c!=0){paste0(",weight=",w,",tremble=",tre,",logistic growth rate=",k,",clustering=",c)}, ".pdf")
#   #pdf(tname, width=7, height=5)
#   print(final)
#   #dev.off()
# }

# Plot extended form of guess across periods ---------------------------
for (ntypelist in c("n_10","n_40")){
  for (ntype in get(ntypelist)){
    df_plot <- filter(df_simulation,type==ntype)
    df_plot$signal<- as.factor(as.character(df_plot$signal))
    
    {final <- ggplot(df_plot, aes(x=period, y=mean, group=signal)) +
        geom_line(aes(color=good)) +
        geom_point(aes(shape=signal, color=good), size = 1, stroke = 1) +
        scale_shape_manual(values=c(1:24))+
        
        scale_x_continuous(name="Period", 
                           limits=c(1, 12), 
                           breaks = c(0:12), 
                           expand = c(0.02, 0) ) +
        scale_y_continuous(name="fraction of agents", 
                           limits=c(0, 1), 
                           expand = c(0.02, 0)) +
        
        coord_cartesian(clip = 'off') +
        theme_bw() +
        theme(legend.justification=c(1,1), legend.position="right")+
        ggtitle(paste0("Detailed simulation of ", ntype, 
                       if (w!=0|tre!=0|k!=1|c!=0){
                         paste0(",weight=",w,
                                ",tremble=",tre,
                                ",logistic growth rate=",k,
                                ",clustering=",c)
                         }))
      
      #final <- final + geom_hline(yintercept=avg_sig,linetype="dashed", color = "red")
    }
    tname <- paste0("../output/Simulation_extended,", ntype, 
                    if (w!=0|tre!=0|k!=1|c!=0){
                      paste0(",weight=",w,
                             ",tremble=",tre,
                             ",logistic growth rate=",k,
                             ",clustering=",c)
                      }, ".pdf")
    #pdf(tname, width=10, height=7)
    print(final)
    #dev.off()
  }
}

# Generate convergent time series- Pool all rounds ---------------------------

#Pool across rounds per network treatment
agg <- aggregate(df_simulation$mean, 
                 list(df_simulation$type, df_simulation$period), 
                 FUN=function(x) c(mean = mean(x), sd = sd(x)))
df_plot<- do.call(data.frame, agg)
colnames(df_plot) <- c("type","period","mean","sd")

for (ntypelist in c("n_10","n_10_partial","n_40")){
  {
    final <- ggplot(subset(df_plot, type %in% get(ntypelist)), 
                    aes(x=period, y=mean, group=type)) +
      geom_line(aes(color=type)) +
      geom_point(aes(shape=type, color=type), size = 3, stroke = 1) +
      scale_shape_manual(values=c(15,16,17,18))+
      
      scale_x_continuous(name="Period", 
                         limits=c(1, 12), 
                         breaks = c(0:12), 
                         expand = c(0.02, 0) ) +
      scale_y_continuous(name="fraction of agents", 
                         limits=c(0.5,1), 
                         expand = c(0.02, 0)) +
      
      coord_cartesian(clip = 'off') +
      theme_bw() +
      ggtitle(paste0("DeGroot simulation, mean guess pooling across group-rounds, ",ntypelist)) +
      theme(legend.title = element_blank(), 
            legend.justification=c(1,1), 
            legend.position=c(0.9,1))
    
    #final <- final + geom_hline(yintercept=avg_sig,linetype="dashed", color = "red")
  }
  tname <- paste0("../output/Simulation-pooled_mean_",ntypelist,".pdf")
  #pdf(tname, width=10, height=7)
  print(final)
  #dev.off()
}

for (ntypelist in c("n_10","n_10_partial","n_40")){
  {
    final <- ggplot(subset(df_plot, type %in% get(ntypelist)), 
                    aes(x=period, y=sd, group=type)) +
      geom_line(aes(color=type)) +
      geom_point(aes(shape=type, color=type), size = 3, stroke = 1) +
      scale_shape_manual(values=c(15,16,17,18))+
      
      scale_x_continuous(name="Period", 
                         limits=c(1, 12), 
                         breaks = c(0:12), 
                         expand = c(0.02, 0) ) +
      scale_y_continuous(name="fraction of agents", 
                         limits=c(0,0.5), 
                         expand = c(0.02, 0)) +
      
      coord_cartesian(clip = 'off') +
      theme_bw() +
      ggtitle(paste0("DeGroot simulation, sd of guess pooling across group-rounds, ",ntypelist)) +
      theme(legend.title = element_blank(), 
            legend.justification=c(1,1), 
            legend.position=c(0.9,1))
    
    #final <- final + geom_hline(yintercept=avg_sig,linetype="dashed", color = "red")
  }
  tname <- paste0("../output/Simulation_pooled_sd_",ntypelist,".pdf")
  #pdf(tname, width=10, height=7)
  print(final)
  #dev.off()
}

#Converge? Check period and period-1 guess consistent
df_converge <- aggregate(df_simulation$converge, 
                         list(df_simulation$type, df_simulation$period), 
                         FUN=mean, na.rm=TRUE)
colnames(df_converge) <- c("type","period","fraction_converge")

for (ntypelist in c("n_10","n_10_partial","n_40")){
  {
    final <- ggplot(subset(df_converge, type %in% get(ntypelist)), 
                    aes(x=period, y=fraction_converge, group=type)) +
      geom_line(aes(color=type)) +
      geom_point(aes(shape=type, color=type), size = 3, stroke = 1) +
      scale_shape_manual(values=c(15,16,17,18))+
      
      scale_x_continuous(name="Between Period t-1 and t", 
                         limits=c(1, 12), 
                         breaks = c(0:12), 
                         expand = c(0.02, 0) ) +
      scale_y_continuous(name="Fraction of networks", 
                         limits=c(0,1), 
                         expand = c(0.02, 0)) +
      
      coord_cartesian(clip = 'off') +
      theme_bw() +
      ggtitle(paste0("DeGroot simulation, Network convergence, ",ntypelist)) +
      theme(legend.title = element_blank(), 
            legend.justification=c(1,1), 
            legend.position=c(0.9,0.3))
    
    #final <- final + geom_hline(yintercept=avg_sig,linetype="dashed", color = "red")
  }
  tname <- paste0("../output/Simulation_Convergence_",ntypelist,".pdf")
  #pdf(tname, width=10, height=7)
  print(final)
  #dev.off()
}

# 1000 signals simulated, main stats (consensus, truth, falsehood, breakdown) ---------------------------

width <- c(0.1)
for (style in c("consensus","truth")){
  for (wd in width){
    if (style == "consensus"){ #check if consensus to 1 or 0 with error of width
      df_plot <- filter(df_1000, (mean>=(1-wd)|mean<=wd))
    }
    if (style == "truth"){ #check if obtain close to correct consensus
      df_plot <- filter(df_1000, mean>=(1-wd))
    }
    df_plot <- df_plot %>% 
      group_by(type, .drop = FALSE) %>%
      summarise(counts = n()/1000) %>%
      mutate(label = paste0(signif(counts*100,3),"%"))
    
    for (ntypelist in c("n_10","n_40")){
      final <- ggplot(filter(df_plot, type %in% get(ntypelist)), 
                      aes(x=type, y=counts)) +
        geom_bar(position="dodge", stat="identity", width=0.6) +
        geom_text(aes(label = label), 
                  vjust = -0.3, 
                  position = position_dodge(width = 0.9)) + 
        scale_y_continuous(labels = scales::percent_format(accuracy = 1L), 
                           limits = c(0,1.05)) +
        
        ggtitle(paste0("../output/Simulation of 1000 signals, ", 
                       str_to_sentence(style),
                       ", width=",wd,", ", ntypelist))+
        xlab("Network type") + 
        ylab("Frequency of networks")+
        theme_bw()
      tname <- paste0("../output/Simulation_no.",style,",", 
                      ntypelist,",width=",wd, 
                      ".pdf")
      #pdf(tname, width=7, height=5)
      print(final)
      #dev.off()
    }
  }
}

width <- c(0.1,0.4)
for (style in c("falsehood")){
  for (wd in width){
    if (style == "falsehood"){ #check if obtain close to incorrect consensus
      df_plot <- filter(df_1000, mean<=wd)
    }
    df_plot <- df_plot %>% 
      group_by(type, .drop = FALSE) %>%
      summarise(counts = n()/1000) %>%
      mutate(label = paste0(signif(counts*100,3),"%"))
    
    for (ntypelist in c("n_10","n_40")){
      final <- ggplot(filter(df_plot, type %in% get(ntypelist)), 
                      aes(x=type, y=counts)) +
        geom_bar(position="dodge", stat="identity", width=0.6) +
        geom_text(aes(label = label), 
                  vjust = -0.3,
                  position = position_dodge(width = 0.9)) + 
        scale_y_continuous(labels = scales::percent_format(accuracy = 1L), 
                           limits = c(0,1.05)) +
        
        ggtitle(paste0("../output/Simulation of 1000 signals, ",
                       str_to_sentence(style),
                       ", width=",wd,", ", ntypelist))+
        xlab("Network type") + 
        ylab("Frequency of networks")+
        theme_bw()
      tname <- paste0("../output/Simulation_no.",style,",", 
                      ntypelist,",width=",wd, ".pdf")
      #pdf(tname, width=7, height=5)
      print(final)
      #dev.off()
    }
  }
}

#*** Part 3 - Others ---------------------------

# # Stable Consensus,truth,stability from period end onwards ---------------------------
# 
# width <- c(0,0.05,0.1,0.15,0.2)
# end <- c(1:12)
# 
# for (style in c("consensus","truth","stability")){
#   assign(paste0("df_",style), data.frame(width = NA, from_period = NA, mean = NA, type =NA))
#   for (wd in width){
#     for (prd in end){
#       for (ntypelist in c("n_10","n_10_partial","n_40")){
#         for (ntype in get(ntypelist)){
#           y <- c()
#           for (sig in c(1:24)){
#             df_temp <- filter(df_simulation, type==ntype, signal==sig, period >= prd)
#             
#             if (style == "consensus"){
#               #check if consensus to 1 or 0 with error of width
#               x <- as.numeric(max(df_temp$mean)<=wd | min(df_temp$mean)>=1-wd)
#             }
#             if (style == "truth"){
#               #check if obtain close to truth
#               x <- as.numeric(min(df_temp$mean)>=1-wd)
#             }
#             if (style == "stability"){
#               #check if the difference in guesses are less than width
#               x <- as.numeric(max(df_temp$mean)-min(df_temp$mean)<=wd)
#             }
#             y <- c(y,x)
#           }
#           assign(paste0("df_",style), rbind(get(paste0("df_",style)), c(wd, prd, mean(y), ntype)))
#         }
#       }
#     }
#   }
#   assign(paste0("df_",style), get(paste0("df_",style))[-1,])
#   assign(paste0("df_",style), transform(get(paste0("df_",style)), type = as.factor(type)))
#   assign(paste0("df_",style), transform(get(paste0("df_",style)), mean = as.numeric(mean)))
#   assign(paste0("df_",style), transform(get(paste0("df_",style)), from_period = as.numeric(from_period)))
# }
# 
# # Print truth, consensus, stability of different width
# for (ntypelist in c("n_10","n_10_partial","n_40")){
#   for (style in c("consensus","truth")){
#     for (wd in c(0,0.05,0.1)){
#       final <- ggplot(filter(subset(get(paste0("df_",style)), type %in% get(ntypelist)), width ==wd)[,-1], 
#                       aes(x=from_period, y=mean, group=type)) +
#         geom_line(aes(color=type)) +
#         geom_point(aes(shape=type, color=type), size = 3, stroke = 1) +
#         scale_shape_manual(values=c(15,16,17,18))+
#         
#         scale_x_continuous(name="From period onwards",
#                            limits=c(1, 12),
#                            breaks = c(0:12),
#                            expand = c(0.02, 0) ) +
#         scale_y_continuous(name="Percentage of network",
#                            limits=c(0,1),
#                            expand = c(0.02, 0),
#                            labels=scales::percent) +
#         
#         coord_cartesian(clip = 'off') +
#         theme_bw() +
#         theme(legend.title = element_blank(), 
#               legend.justification=c(1,1), 
#               legend.position=c(0.9,0.5)) +
#         ggtitle(paste0("../output/Simulation ", style,", width = ",wd))
#       
#       #final <- final + geom_hline(yintercept=avg_sig,linetype="dashed", color = "red")
#       
#       tname <- paste0("../output/Simulation_", style,"_", 
#                       ntypelist,"_width=", wd,".pdf")
#       #pdf(tname, width=10, height=7)
#       print(final)
#       #dev.off()
#     }
#   }
# }


# #Find the clans ---------------------------
# df_clan <- data.frame(clan = NA, type = NA)
# for (ntype in c(n_10_complete,n_40))
# {
#   g <- as.matrix(glist[ntype][[1]])
#   g <- graph_from_adjacency_matrix(g-diag(ncol(g)))
#   xx <- r_cohesiveness(g,5)
#   print(c(ntype, round(xx[[1]],3)))
#   if(xx[[1]] != 0){
#     df_temp <- as.data.frame(as.matrix(xx[[2]]))
#     df_temp$type <- ntype
#     colnames(df_temp) <- c("clan","type")
#     df_clan <- rbind(df_clan,df_temp)
#   }
# }
# df_clan <- df_clan[-1,]
# 
# #print list of clans
# write.csv(as.matrix(df_clan),file = "../output/dataframe/df_clan_dirty.csv")

# Saving environment for later use
save.image("RData/00-simulation-envir.RData")