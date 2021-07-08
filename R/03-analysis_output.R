# To print all output into "output" folder
# Use ctrl-F to:
# replace "#pdf(" with "pdf("
# replace "#dev.off()" with "dev.off()"

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

load("RData/02-analysis-envir.RData")

# Experimental time-series data plotted with simulation theory (data_theory) ---------------------------

# aggregate data across rounds i to 6, by network type and signal class (good)
i <- 1
df_temp <- filter(df_summary,round %in% seq(i,6))
agg <- aggregate(df_temp$mean, list(df_temp$type, df_temp$period, df_temp$good), 
                 FUN=function(x) c(mean = mean(x), sd = sd(x)))
df_plot<- do.call(data.frame, agg)
colnames(df_plot) <- c("type","period","good","mean","sd")

# aggregate simulation rounds i to 6, by network type and signal class (good)
sig_interest <- c(seq(i,6),seq(i+6,6*2),seq(i+6*2,6*3),seq(i+6*3,6*4))
df_temp <- filter(df_simulation, signal %in% sig_interest)
agg <- aggregate(df_temp$mean, list(df_temp$type, df_temp$period, df_temp$good), 
                 FUN=function(x) c(mean = mean(x), sd = sd(x)))
df_plot2<- do.call(data.frame, agg)
colnames(df_plot2) <- c("type","period","good","mean","sd")

# plot time-series data against simulation, by network type and signal class (good)
for (ntype in c(n_10,n_40)){
  final <- 
    ggplot(subset(df_plot, type == ntype), 
           aes(x=period, y=mean, group=good, color=good)) +
    geom_line() +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.3,position=position_dodge(0.05))+
    geom_point(aes(color=good), size = 3, stroke = 1) +
    scale_shape_manual(values=c(15,16,17,18))+
    
    geom_line(data=subset(df_plot2, type == ntype), 
              aes(x=period, y=mean, group=good), 
              linetype = "dashed") +
    
    # geom_errorbar(data=subset(df_plot2, type == ntype), 
    #               aes(ymin=mean-sd, ymax=mean+sd), 
    #               width=.2,position=position_dodge(0.1), 
    #               linetype = "dashed")+
    
    scale_x_continuous(name="Period", 
                       limits=c(0.5, 12.5), 
                       breaks = c(0:12), 
                       expand = c(0.02, 0) ) +
    scale_y_continuous(name="fraction of players", 
                       limits=c(-0.05,1.05), 
                       expand = c(0.02, 0)) +
    
    coord_cartesian(clip = 'off') +
    theme_bw() +
    theme(legend.title = element_blank(), 
          legend.justification=c(1,1), 
          legend.position=c(0.2,0.2))+
    
    ggtitle(paste0("Fraction of players guessing truth (majority signals),",ntype, 
                   if(w>0){paste0(",weight=",w)},
                   if(tre>0){paste0(",tremble=",tre)}))
  
  #final <- final + geom_hline(yintercept=avg_sig,linetype="dashed", color = "red")
  tname <- paste0("../output/data_theory,",ntype,
                  if(w>0){paste0(",weight=",w)},
                  if(tre>0){paste0(",tremble=",tre)},
                  ".pdf")
  #pdf(tname, width=10, height=7)
  print(final)
  #dev.off()
}


# Average guesses of each network treatment over periods (truth) ---------------------------
# pooled across all groups, rounds 

# Pool across rounds per network treatment
df_temp <- filter(df_summary, round %in% seq(1,6))
df_temp <- df_temp[!(df_temp$good=="Incorrect consensus" & grepl("RF",df_temp$type)),] #remove all RF signals with incorrect consensus

agg <- aggregate(df_temp$mean, list(df_temp$type, df_temp$period),
                 FUN=function(x) c(mean = mean(x), sd = sd(x)))
df_plot<- do.call(data.frame, agg)
colnames(df_plot) <- c("type","period","mean","sd")
df_plot$type <- as.factor(as.character(df_plot$type))

# Generate convergent time series- Pool all rounds
for (ntypelist in c("n_10","n_40")){
  final <- 
    ggplot(subset(df_plot, type %in% get(ntypelist)), 
           aes(x=period, y=mean, group=type, color=type)) +
    geom_line() +
    geom_point(aes(shape=type, color=type), size = 3, stroke = 1) +
    scale_shape_manual(values=c(15,16,17,18))+
    #geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.3,position=position_dodge(0.05))+
    geom_hline(yintercept = 0.7, color = "grey", linetype="dashed")+
    
    scale_x_continuous(name="Period", 
                       limits=c(1, 12), 
                       breaks = c(0:12), 
                       expand = c(0.02, 0) ) +
    scale_y_continuous(name="Fraction of players", 
                       limits=c(0.6,1), 
                       expand = c(0.02, 0)) +
    
    coord_cartesian(clip = 'off') +
    theme_bw() +
    theme(legend.title = element_blank(), 
          legend.justification=c(1,1), 
          legend.position=c(0.2,0.95))+
    
    ggtitle(paste0("Fraction of players guessing truth (majority signals),", ntypelist)) + 
    labs(caption = "*Bad signals in RF removed")
  
  #final <- final + geom_hline(yintercept=avg_sig,linetype="dashed", color = "red")
  
  tname <- paste0("../output/Truth,", ntypelist, ".pdf")
  #pdf(tname, width=7, height=5)
  print(final)
  #dev.off()
}

#Generate convergent time series- Pool all rounds - standard deviation
for (ntypelist in c("n_10", "n_40")){
  final <- ggplot(subset(df_plot, type %in% get(ntypelist)), 
                  aes(x=period, y=sd, group=type, color=type)) +
    geom_line() +
    geom_point(aes(shape=type, color=type), size = 3, stroke = 1) +
    scale_shape_manual(values=c(15,16,17,18))+
    #geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.3,position=position_dodge(0.05))+
    
    scale_x_continuous(name="Period", 
                       limits=c(1, 12), 
                       breaks = c(0:12), 
                       expand = c(0.02, 0) ) +
    scale_y_continuous(name="Standard deviation of mean guess", 
                       limits=c(0,0.5), 
                       expand = c(0.02, 0)) +
    
    coord_cartesian(clip = 'off') +
    theme_bw() +
    theme(legend.title = element_blank(), 
          legend.justification=c(1,1), 
          legend.position=c(0.9,1.05))+
    ggtitle(paste0("Standard deviation of network guesses across all group-rounds,", ntypelist))
  
  #final <- final + geom_hline(yintercept=avg_sig,linetype="dashed", color = "red")
  
  tname <- paste0("../output/Truth-sd,", ntypelist, ".pdf")
  #pdf(tname, width=10, height=7)
  print(final)
  #dev.off()
}

# Period 10/11 distribution of average guess across groups and rounds (guess_distribution) ---------------------------
prd <- 10
for (cl in c("class-","")){ #with classification or without
  for (ntypelist in c("n_10","n_40")){
    df_plot<-filter(df_summary, type %in% get(ntypelist), period ==prd)
    if (cl == "class-"){
      final <- ggplot(df_plot, aes(x=mean, group=good, fill=good))
    }
    else {
      final <- ggplot(df_plot, aes(x=mean))
    }
    final <- final +
      geom_histogram(breaks = seq(0,1,0.1), color="white", closed = "left")+
      stat_bin(breaks = seq(0,1,0.1), closed = "left", geom="text", 
               colour="black", size=2.5, 
               aes(label= ifelse(..count.. > 0, ..count.., ""), 
                   group=type, y=1+(..count..))) +
      geom_vline(aes(xintercept = 1, color = "Bayesian"), linetype="dashed")+
      geom_vline(aes(xintercept = 0.7, color = "No learning"), linetype="dashed")+
      
      facet_wrap(~ factor(type, levels=get(ntypelist)))+
      #ylim(NA, 15) +
      scale_colour_manual(values = c("blue","red"))+
      scale_x_continuous(breaks=seq(0,1,0.1), limits = c(-0.05, 1.05))+
      ggtitle(paste0("Period ",prd," distribution of guessing truth, ", ntypelist))+
      xlab("Average guess of network") + 
      ylab("Frequency of networks (out of 24)")+
      theme_bw() +
      theme(legend.title = element_blank(), 
            strip.text.x = element_text(size = 13), 
            legend.justification=c(1,1), 
            legend.text = element_text(size=10), 
            legend.position=c(0.2,0.92))
    
    tname <- paste0("../output/guess_distribution-",cl, ntypelist, ".pdf")
    #pdf(tname, width=10, height=7)
    print(final)
    #dev.off()
  }
}

# Compare 4 main learning outcomes across networks (no."style") ---------------------------
# consensus, truth, falsehood, breakdown

width <- c(0.1, 0.2) #margin of error allowed
for (style in c("consensus", "truth")) {
  for (wd in width) {
    if (style == "consensus") {
      # check if consensus to 1 or 0 with error of width
      df_plot <-
        filter(df_summary, type %in% c(n_10,n_40), 
               period >= 7, (mean >= (1 - wd) | mean <= wd))
      df_plot2 <-
        filter(df_simulation, type %in% c(n_10,n_40), 
               period >= 7, (mean >= (1 - wd) | mean <= wd))
    }
    if (style == "truth") {
      # check if obtain close to correct consensus
      df_plot <- filter(df_summary, type %in% c(n_10,n_40), 
                        period >= 7, mean >= (1 - wd))
      df_plot2 <- filter(df_simulation, type %in% c(n_10,n_40), 
                         period >= 7, mean >= (1 - wd))
    }
    df_plot <- 
      df_plot %>% group_by(type, .drop = FALSE) %>% summarise(mean = mean(mean))
    df_plot2 <-
      df_plot2 %>% group_by(type, .drop = FALSE) %>% summarise(mean = mean(mean))
    
    df_plot$id <- "Data"
    df_plot2$id <- "Simulation"
    df_temp <- df_plot
    #df_temp <- rbind(df_plot2, df_plot)
    
    # print results bar charts
    for (ntypelist in c("n_10", "n_40")) {
      final <-
        ggplot(filter(df_temp, type %in% get(ntypelist)),
               aes( x = type, y = mean, group = id, fill = id )) +
        geom_bar(
          position = "dodge",
          stat = "identity",
          color = "white",
          width = 0.6
        ) +
        geom_text(
          aes(label = mean, group = id),
          vjust = -0.3,
          position = position_dodge(width = 0.6)
        ) +
        
        #ylim(NA, 25) +
        ggtitle(paste0(str_to_sentence(style), ", width=", wd, ", ", ntypelist)) +
        xlab("Network type") + ylab("Frequency of networks (out of 24)") +
        theme_bw()
      tname <- paste0("../output/no.",style,",",ntypelist,",width=",wd,".pdf")
      #pdf(tname, width=7, height=5)
      print(final)
      #dev.off()
    }
  }
}

width <- c(0.2, 0.3)
for (style in c("falsehood", "breakdown")) {
  for (wd in width) {
    if (style == "falsehood") {
      # check if obtain close to incorrect consensus
      df_plot <- filter(df_summary, period == 10, mean <= wd)
      df_plot2 <- filter(df_simulation, period == 10, mean <= wd)
    }
    if (style == "breakdown") {
      # check if consensus to 1 or 0 with error of width
      df_plot <-
        filter(df_summary, period == 10, (mean < (1 - wd) & mean > wd))
      df_plot2 <-
        filter(df_simulation, period == 10, (mean < (1 - wd) & mean > wd))
    }
    df_plot <-
      df_plot %>% group_by(type, .drop = FALSE) %>% summarise(counts = n())
    df_plot2 <-
      df_plot2 %>% group_by(type, .drop = FALSE) %>% summarise(counts = n())
    df_plot$id <- "Data"
    df_plot2$id <- "Simulation"
    df_temp <- df_plot
    #df_temp <- rbind(df_plot2, df_plot)
    
    # print result bar charts
    for (ntypelist in c("n_10", "n_40")) {
      final <-
        ggplot(filter(df_temp, type %in% get(ntypelist)),
               aes(x = type, y = counts, group = id, fill = id )) +
        geom_bar(
          position = "dodge",
          stat = "identity",
          color = "white",
          width = 0.6
        ) +
        geom_text(aes(label = counts),
                  vjust = -0.3,
                  position = position_dodge(width = 0.9)) +
        #ylim(NA, 25) +
        ggtitle(paste0(str_to_sentence(style), ", width=", wd, ", ", ntypelist)) +
        xlab("Network type") + ylab("Frequency of networks (out of 24)") +
        theme_bw()
      tname <-
        paste0("../output/no.",style,",",ntypelist,",width=",wd,".pdf")
      #pdf(tname, width=7, height=5)
      print(final)
      #dev.off()
    }
  }
}

# Compare 4 main learning outcomes across networks (p."style") ---------------------------
# consensus, truth, falsehood, breakdown

width <- c(0.15, 0.2) #margin of error allowed
for (style in c("consensus", "truth", "falsehood", "breakdown")) {
  for (wd in width) {
    if (style == "consensus") {
      # check if consensus to 1 or 0 with error of width
      df_plot <- df_summary %>%
        mutate(consensus = ifelse(mean >= 1 - wd | mean <= wd, 1,0))
      df_plot2 <- df_simulation %>%
        mutate(consensus = ifelse(mean >= 1 - wd | mean <= wd, 1,0))
    }
    if (style == "truth") {
      # check if obtain close to correct consensus
      df_plot <- df_summary %>%
        mutate(consensus = ifelse(mean >= 1 - wd, 1,0))
      df_plot2 <- df_simulation %>%
        mutate(consensus = ifelse(mean >= 1 - wd, 1,0))
    }
    if (style == "falsehood") {
      # check if consensus to 1 or 0 with error of width
      df_plot <- df_summary %>%
        mutate(consensus = ifelse(mean <= wd, 1,0))
      df_plot2 <- df_simulation %>%
        mutate(consensus = ifelse(mean <= wd, 1,0))
    }
    if (style == "breakdown") {
      # check if obtain close to correct consensus
      df_plot <- df_summary %>%
        mutate(consensus = ifelse(mean < (1 - wd) & mean > wd, 1,0))
      df_plot2 <- df_simulation %>%
        mutate(consensus = ifelse(mean < (1 - wd) & mean > wd, 1,0))
    }
    
    df_plot <- df_plot %>% 
      filter(type %in% c(n_10,n_40), period >= 7) %>%
      group_by(type, .drop = FALSE) %>% 
      summarise(mean = mean(consensus),
                count = n()) %>% 
      mutate(sd = (mean * (1 - mean) / count)^(0.5))
    
    df_plot2 <- df_plot2 %>% 
      filter(type %in% c(n_10,n_40), period >= 7) %>%
      group_by(type, .drop = FALSE) %>% 
      summarise(mean = mean(consensus),
                count = n()) %>% 
      mutate(sd = (mean * (1 - mean) / count)^(0.5))
    
    df_plot$id <- "Data"
    df_plot2$id <- "Simulation"
    df_temp <- df_plot #only data bar plot
    #df_temp <- rbind(df_plot2, df_plot) #combine data and simulation
    
    # print results bar charts
    for (ntypelist in c("n_10", "n_40")) {
      final <-
        ggplot(filter(df_temp, type %in% get(ntypelist)),
               aes( x = type, y = mean, group = id, fill = id )) +
        geom_bar(
          position = "dodge",
          stat = "identity",
          color = "white",
          width = 0.6
        ) +
        geom_text(
          aes(label = signif(mean,3), group = id),
          vjust = -1.7,
          position = position_dodge(width = 0.6)
        ) +
        geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), 
                      width = .2,
                      position = position_dodge(0.6)) +
        
        ylim(NA, 1) +
        ggtitle(paste0(str_to_sentence(style), ", width=", wd, ", ", ntypelist)) +
        xlab("Network type") + ylab("Fraction of networks") +
        theme_bw()
      tname <- paste0("../output/p.",style,",",ntypelist,",width=",wd,".pdf")
      #pdf(tname, width=7, height=5)
      print(final)
      #dev.off()
    }
  }
}

# Convergence time of data, histogram using df_converge (contime) ---------------------------
width <- c(0.1,0.2)

for (wd in width){
  for (ntypelist in c("n_10","n_40")){
    final <- ggplot(filter(df_converge, type %in% get(ntypelist), width==wd), 
                    aes(x=from_period)) + 
      geom_histogram(binwidth = 1, color="white")+
      stat_bin(binwidth = 1, geom="text", colour="black", size=3.5, 
               aes(label= ifelse(..count.. > 5, ..count.., ""), y=2+(..count..))) +
      #geom_vline(aes(xintercept=avg.con_time, group=type), colour="red") +
      
      #ylim(NA, 24) +
      scale_x_continuous(breaks=seq(1,12,1), limits = c(0.5, 12.5))+
      facet_wrap(~ factor(type, levels=get(ntypelist)))+
      ggtitle(paste0("Convergence time, ",ntypelist,", width=",wd,
                     if (w!=0|tre!=0){paste0(",weight=",w,",tremble=",tre)})) +
      
      xlab("Convergence time") + ylab("Frequency of networks (out of 24)")+
      theme_bw()+
      theme(legend.title = element_blank(), 
            strip.text.x = element_text(size = 13))
    
    tname <- paste0("../output/contime,", ntypelist,",width=",wd, 
                    if (w!=0|tre!=0){paste0(",weight=",w,",tremble=",tre)}, 
                    ".pdf")
    #pdf(tname, width=10, height=7)
    print(final)
    #dev.off()
  }
}

# Stable Consensus, truth, stability (stable_"style") ---------------------------
# reaches "style" from period prd onwards

width <- c(0.1,0.15,0.2)
for (ntypelist in c("n_10","n_40")){
  for (style in c("truth","consensus")){
    for (wd in width){
      final <- ggplot(filter(subset(get(paste0("df_",style)), type %in% get(ntypelist)), width ==wd)[,-1], 
                      aes(x=from_period, y=mean, group=type)) +
        geom_line(aes(linetype=type)) +
        geom_point(aes(shape=type, color=type), size = 3, stroke = 1) +
        scale_shape_manual(values=c(15,16,17,18))+
        
        scale_x_continuous(name="From period onwards", 
                           limits=c(1, 12), 
                           breaks = c(0:12), 
                           expand = c(0.02, 0) ) +
        scale_y_continuous(name="Percentage of network", 
                           limits=c(0,1), 
                           expand = c(0.02, 0)) +
        
        coord_cartesian(clip = 'off') +
        theme_bw() +
        theme(legend.title = element_blank(), 
              legend.justification=c(1,1), 
              legend.position=c(0.2,.8))+
        ggtitle(paste0("Stable ", style,", width = ",wd))
      
      #final <- final + geom_hline(yintercept=avg_sig,linetype="dashed", color = "red")
      
      tname <- paste0("../output/stable_", style,",", ntype,",width=", wd,".pdf")
      #pdf(tname, width=10, height=7)
      print(final)
      #dev.off()
    }
  }
}

# Raw time-series plot with no pooling (rawaverage) ---------------------------
# by  network treatment
agg <- aggregate( df_summary$mean,
                  list( df_summary$type, df_summary$group, df_summary$round,
                        df_summary$period, df_summary$good ),
                  FUN = function(x) c(mean = mean(x))
)

df_plot<- do.call(data.frame, agg)
colnames(df_plot) <- c("type","group","round","period","good","mean")
df_plot$round<- as.factor(as.character(df_plot$round))

for (ntype in c(n_10,n_40)){
  for (grp in c(1:max_grp)){
    {final <- ggplot(filter(df_plot, group ==grp, type==ntype), 
                     aes(x=period, y=mean, group=round)) +
      geom_line(aes(linetype=round)) +
      geom_point(aes(shape = round, color = good), size = 3, stroke = 1) +
      scale_shape_manual( values = c(15,16,17,18,3,4))+
      
      scale_x_continuous(name="Period", 
                         limits=c(1, 12), 
                         breaks = c(0:12), 
                         expand = c(0.02, 0) ) +
      scale_y_continuous(name="fraction of players",
                         limits=c(0, 1), 
                         expand = c(0.02, 0)) +
      
      coord_cartesian(clip = 'off') +
      theme_bw() +
      theme(legend.justification=c(1,1), legend.position="right")+
      ggtitle(paste0(ntype," Group",grp))
    
    final <- final + 
      scale_colour_discrete(name  ="",labels=c("Correct\nConsensus", "Other"))
    }
    tname <- paste0("../output/rawaverage",ntype,"-Group",grp,".pdf")
    #pdf(tname, width=10, height=7)
    print(final)
    #dev.off()
  }
}

# Pooling individual consensus by network type (converge) ---------------------------

df_plot=data.frame(type=NA, from_period = NA, mean = NA)
for (ntype in c(n_10,n_40)){
  for (prd in c(1:12)){
    df_temp <- filter(df_indivcon, type==ntype, from_period<=prd)
    
    # x is fraction of players stop changing guesses from period prd onward
    x<- nrow(df_temp)/((if(grepl("40",ntype)){40} else{10})*max_grp*max_rnd)
    df_plot <- rbind(df_plot, c(ntype, prd, x))
  }
}
df_plot <- df_plot[-1,]
df_plot$type <- as.factor(as.character(df_plot$type))
df_plot$mean <- as.numeric(df_plot$mean)
df_plot$from_period <- as.numeric(df_plot$from_period)

for (ntypelist in c("n_10","n_40")){
  final <- ggplot(filter(df_plot, type %in% get(ntypelist)), 
                  aes(x=from_period, y=mean, group=type)) +
    geom_line(aes(linetype=type)) +
    geom_point(aes(shape=type, color=type), size = 3, stroke = 1) +
    scale_shape_manual(values=c(15,16,17,18))+
    
    scale_x_continuous(name="From period onwards", 
                       limits=c(1, 12), 
                       breaks = c(0:12), 
                       expand = c(0.02, 0) ) +
    scale_y_continuous(name="Fraction of players", 
                       limits=c(0.4,1), 
                       expand = c(0.02, 0)) +
    
    coord_cartesian(clip = 'off') +
    theme_bw() +
    theme(legend.title = element_blank(), 
          legend.justification=c(1,1), 
          legend.position=c(0.2,.8))+
    ggtitle(paste0("Fraction of players converged,", ntypelist))
  
  #final <- final + geom_hline(yintercept=avg_sig,linetype="dashed", color = "red")
  
  tname <- paste0("../output/converge,",ntypelist,".pdf")
  #pdf(tname, width=10, height=7)
  print(final)
  #dev.off()
}

# #individual consensus per ntype per group
# df_plot=data.frame(type=NA, group=NA, round=NA, from_period = NA, mean = NA)
# for (ntype in c(n_10_all,n_40)){
#   for (grp in c(1:max_grp)){
#     for (rnd in c(1:max_rnd)){
#       for (prd in c(1:12)){
#         df_temp <- filter(df_indivcon, type==ntype, group==grp, round==rnd, from_period<=prd)
# 
#         #x is fraction of players stop changing guesses from period prd onward
#         x<- nrow(df_temp)/(if(grepl("40",ntype)){40} else{10})
#         df_plot <- rbind(df_plot, c(ntype, grp, rnd, prd, x))
#       }
#     }
#   }
# }
# df_plot <- df_plot[-1,]
# df_plot$type <- as.factor(as.character(df_plot$type))
# df_plot$mean <- as.numeric(df_plot$mean)
# df_plot$from_period <- as.numeric(df_plot$from_period)
# 
# for (ntype in c(n_10,n_40)){
#   for (grp in c(1:max_grp)){
#     final <- ggplot(filter(df_plot, type==ntype, group==grp)[,-2], aes(x=from_period, y=mean, group=round)) +
#       geom_line() +
#       geom_point(aes(shape=round, color=round), size = 3, stroke = 1) +
#       scale_shape_manual(values=c(15:21))+
# 
#       scale_x_continuous(name="From period onwards", limits=c(1, 12), breaks = c(0:12), expand = c(0.02, 0) ) +
#       scale_y_continuous(name="Fraction of guesses", limits=c(0,1), expand = c(0.02, 0)) +
# 
#       coord_cartesian(clip = 'off') +
#       theme_bw() +
#       theme(legend.title = element_blank(), legend.justification=c(1,1), legend.position=c(0.2,.8))+
#       ggtitle(paste0("Individual consensus, ",ntype, " Group",grp))
# 
#     #final <- final + geom_hline(yintercept=avg_sig,linetype="dashed", color = "red")
# 
#     tname <- paste0("../output/Individual consensus, ",ntype, " Group",grp,".pdf")
#     #pdf(tname, width=10, height=7)
#     print(final)
#     #dev.off()
#   }
# }

# Good and Bad signals (goods, bads) ---------------------------
# Using df_goods
for (gb in c("goods")){
  for (ntypelist in c("n_10","n_40")){
    final <- ggplot(filter(df_goods, type %in% get(ntypelist)), 
                    aes(x=mean, fill=type)) +
      geom_histogram(breaks = seq(0,1,0.2), color="white", closed = "left")+
      stat_bin(breaks = seq(0,1,0.2), 
               geom="text", 
               closed = "left", 
               aes(label= ifelse(..count.. > 0, ..count.., ""), y=0.5+(..count..))) +
      facet_wrap(~ factor(type, levels=get(ntypelist)))+
      
      ylim(NA, max(unlist(lapply(filter(get(gb), type %in% get(ntypelist))$goodbad, length)))) +
      scale_x_continuous(breaks=seq(0,1,0.1), limits = c(-0.05, 1.05))+
      xlab("Average guess of network") + ylab(paste0("Frequency of networks"))+
      
      ggtitle(paste0("Comparing signals - ",gb,",", ntypelist))+
      theme_bw()+
      theme(legend.title = element_blank(), legend.position="none")
    
    tname <- paste0("../output/",gb,",", ntypelist, ".pdf")
    #pdf(tname, width=7, height=5)
    print(final)
    #dev.off()
  }
}

for (gb in c("bads")){
  for (ntypelist in c("n_10","n_40")){
    df_plot <- filter(df_summary, type %in% get(ntypelist), period ==10)
    df_plot$signal <- df_plot$round+(df_plot$group-1)*6
    
    for (ntype in get(ntypelist)[-1]){
      df_temp <- df_plot[1,]
      df_temp <- rbind(df_temp, filter(df_plot, type==ntype, signal %in% filter(get(gb), type==ntype)$goodbad[[1]]))
      df_temp <- rbind(df_temp, filter(df_plot, type==(if(ntypelist=="n_10"){"ER_10"} else{"ER_40"}), signal %in% filter(get(gb), type==ntype)$goodbad[[1]]))
      df_temp <- df_temp[-1,]
      df_temp$mean <- as.numeric(df_temp$mean)
      
      final <- ggplot(df_temp, aes(x=mean, fill=type)) +
        geom_histogram(breaks = seq(0,1,0.2), color="white", closed = "left")+
        stat_bin(breaks = seq(0,1,0.2), geom="text", closed = "left", 
                 aes(label= ifelse(..count.. > 0, ..count.., ""), y=0.5+(..count..))) +
        facet_wrap(~ factor(type))+
        
        ylim(NA, max(unlist(lapply(filter(get(gb), type %in% get(ntypelist))$goodbad, length)))) +
        scale_x_continuous(breaks=seq(0,1,0.1), limits = c(-0.05, 1.05))+
        xlab("Average guess of network") + 
        ylab(paste0("Frequency of networks"))+
        
        ggtitle(paste0("Comparing signals - ",gb,",", ntypelist))+
        theme_bw()+
        theme(legend.title = element_blank(), legend.position="none")
      
      tname <- paste0("../output/",gb,",", ntype, ".pdf")
      #pdf(tname, width=7, height=5)
      print(final)
      #dev.off()
    }
  }
}

# Influencers of RF (influencer) ---------------------------

df_influencer <- data.frame(typelist = NA, signal = NA)
for (ntypelist in c("n_10","n_40")){
  balls <- as.matrix(read_csv(paste0("../input/signals/Signals_",ntypelist,".csv"), col_names = FALSE))-1
  balls <- balls[,c(1,2,3)]
  r <- c()
  # record list of signals where influencers get wrong signals
  for (i in 1:nrow(balls)){
    if(mean(balls[i,c(1,2,3)]) <= 0.5) {r <- c(r,i)}
  }
  df_influencer <- rbind(df_influencer, list(ntypelist,list(r)))
}
df_influencer <- df_influencer[-1,]

for (ntypelist in c("n_10","n_40")){
  df_plot <- filter(df_summary, type %in% get(ntypelist), period >= 7)
  df_plot$signal <- df_plot$round+(df_plot$group-1)*6
  
  for (ntype in get(ntypelist)[grepl("RF",get(ntypelist))]){
    df_temp <- df_plot[1,]
    df_temp <- rbind(df_temp, 
                     filter(df_plot, 
                            type==ntype, 
                            signal %in% filter(df_influencer, typelist==ntypelist)$signal[[1]]))
    df_temp <- rbind(df_temp, 
                     filter(df_plot, 
                            type==get(ntypelist)[grepl("ER",get(ntypelist))], 
                            signal %in% filter(df_influencer, typelist==ntypelist)$signal[[1]]))
    df_temp <- df_temp[-1,]
    df_temp$mean <- as.numeric(df_temp$mean)
    
    final <- ggplot(df_temp, aes(x=mean, fill=type)) +
      geom_histogram(breaks = seq(0,1,0.1), 
                     color="white", 
                     closed = "left")+
      stat_bin(breaks = seq(0,1,0.1), 
               geom="text", 
               closed = "left", 
               aes(label= ifelse(..count.. > 0, ..count.., ""), 
                   y=0.5+(..count..))) +
      facet_wrap(~ factor(type))+
      
      #ylim(NA, length(filter(df_influencer, typelist==ntypelist)$signal[[1]])) +
      scale_x_continuous(breaks=seq(0,1,0.1), 
                         limits = c(-0.05, 1.05))+
      xlab("Average guess of network") + 
      ylab(paste0("Frequency of networks"))+
      
      ggtitle(paste0("Effects of bad influencers, ", ntypelist))+
      theme_bw()+
      theme(legend.title = element_blank(), 
            legend.position="none")
    
    tname <- paste0("../output/influencer,", ntypelist, ".pdf")
    #pdf(tname, width=7, height=5)
    print(final)
    #dev.off()
  }
}
