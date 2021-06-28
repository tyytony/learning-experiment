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
source("function/Function-DeGroot_once,DeGroot.R")

# *** Part 1 - Creates dataframes from df_tot ---------------------

load("RData/02-analysis-envir.RData")

# forming period by period guess vector from data (df_detail, df_truth)  ---------------------------
df_detail <- data.frame(group = NA, round = NA, period = NA, type =NA, 
                        guess =NA, id_subj = NA, id_grp = NA)
df_truth <- data.frame(group = NA, round = NA, type =NA, truth = NA)

for (ntype in c(n_10_all,n_40)){
  for (grp in c(1:max_grp)){
    for (rnd in c(1:max_rnd)){
      df_temp <- filter(df_tot, type==ntype & group == grp & round == rnd)
      tru <- round(mean(df_temp$sig),0)-1 #Majority of signal from group-round: 1 if mean signal >=0.5
      df_truth <- rbind(df_truth, c(grp, rnd, ntype, tru))
      
      for (prd in c(1:max_prd)){
        df_temp2<- filter(df_temp, period==prd)
        s<- df_temp2[order(df_temp2$id_grp),]$choice #the vector of guesses from...
        subj <- df_temp2[order(df_temp2$id_grp),]$id_subj #the subjects, in order
        location <- df_temp2[order(df_temp2$id_grp),]$id_grp
        df_detail <- rbind(df_detail,list(grp, rnd, prd, ntype, list(s-1),list(subj),list(location)))
      }
    }
  }
}
df_detail <- df_detail[-1,]
df_truth <- df_truth[-1,]

# performing DeGroot of next period on guess vector using network glist (df_DeGroot, df_DeGrootmore)  ---------------------------
df_DeGroot=data.frame(group = NA, round = NA, period = NA, 
                      type =NA, DeGroot =NA, id_subj = NA)
df_DeGrootmore=data.frame(group = NA, round = NA, period = NA, 
                          type =NA, DeGroot =NA, id_subj = NA)

for (ntype in c(n_10_all,n_40)){
  g <- glist[ntype][[1]]
  for (grp in c(1:max_grp)){
    for (rnd in c(1:max_rnd)){
      for (prd in c(1:(max_prd-1))){
        s <- as.numeric(filter(df_detail, type==ntype & group == grp & round == rnd & period == prd)$guess[[1]]) #Guess vector from data at period t
        s_next <- as.numeric(filter(df_detail, type==ntype & group == grp & round == rnd & period == (prd+1))$guess[[1]]) #Guess vector from data at period t+1
        predict <- DeGroot_once(g,s,w,tre) #Simulation from guess vector at period t for next period
        subj <- filter(df_detail, type==ntype & group == grp & round == rnd & period == prd)$id_subj[[1]]
        df_DeGroot <- rbind(df_DeGroot,list(grp, rnd, prd+1, ntype, list(as.numeric(s_next==predict)),list(subj))) #check if the simulation matches with guess at period t+1
        for (j in 1:length(subj)){
          df_DeGrootmore <- rbind(df_DeGrootmore,c(grp, rnd, prd+1, ntype, as.numeric(s_next==predict)[j],subj[j]))
        }
      }
    }
  }
}
df_DeGroot <- df_DeGroot[-1,]
df_DeGrootmore <- df_DeGrootmore[-1,]
df_DeGrootmore$DeGroot <- as.numeric(df_DeGrootmore$DeGroot)

# Repeat for just guessing signal, how well that matches with DeGroot (df_sigdetail, df_sigDeGroot)  ---------------------------
df_sigdetail=data.frame(group = NA, round = NA, period = NA, type =NA, guess =NA, id_subj = NA)
for (ntype in c(n_10_all,n_40)){
  for (grp in c(1:max_grp)){
    for (rnd in c(1:max_rnd)){
      df_temp <- filter(df_tot, type==ntype & group == grp & round == rnd)
      
      for (prd in c(1:max_prd)){
        df_temp2<- filter(df_temp, period==prd)
        s<- df_temp2[order(df_temp2$id_grp),]$sig
        subj <- df_temp2[order(df_temp2$id_grp),]$id_subj
        df_sigdetail <- rbind(df_sigdetail,list(grp, rnd, prd, ntype, list(s-1),list(subj)))
      }
    }
  }
}
df_sigdetail <- df_sigdetail[-1,]

df_sigDeGroot=data.frame(group = NA, round = NA, period = NA, type =NA, DeGroot =NA, id_subj = NA)
for (ntype in c(n_10,n_40)){
  g <- glist[ntype][[1]]
  for (grp in c(1:max_grp)){
    for (rnd in c(1:max_rnd)){
      for (prd in c(1:(max_prd-1))){
        s <- as.numeric(filter(df_sigdetail, type==ntype & group == grp & round == rnd & period == prd)$guess[[1]])
        s_next <- as.numeric(filter(df_sigdetail, type==ntype & group == grp & round == rnd & period == (prd+1))$guess[[1]])
        predict <- DeGroot_once(g,s)
        subj <- filter(df_sigdetail, type==ntype & group == grp & round == rnd & period == prd)$id_subj
        df_sigDeGroot <- rbind(df_sigDeGroot,list(grp, rnd, prd+1, ntype, list(as.numeric(s_next==predict)),subj))
      }
    }
  }
}
df_sigDeGroot <- df_sigDeGroot[-1,]


# Analyze per network type's percentage of DeGroot match, and other matches (df_match)  ---------------------------
df_match <- data.frame(type=NA,match=NA,percentage=NA)
for (ntype in c(n_10,n_40)){
  sum_sig<- sum(as.numeric(lapply(filter(df_DeGroot, type ==ntype)$DeGroot,sum)))
  num_sig<- sum(as.numeric(lapply(filter(df_DeGroot, type ==ntype)$DeGroot,length)))
  df_match <- rbind(df_match,c(ntype, "Data match with DeGroot", signif(sum_sig/num_sig,3)))
  
  sum_sig<- sum(as.numeric(lapply(filter(df_sigDeGroot, type ==ntype)$DeGroot,sum)))
  num_sig<- sum(as.numeric(lapply(filter(df_sigDeGroot, type ==ntype)$DeGroot,length)))
  df_match <- rbind(df_match,c(ntype, "Guessing signal match with DeGroot", signif(sum_sig/num_sig,3)))
}
df_match <- df_match[-1,]
df_match$percentage <- as.numeric(df_match$percentage)

# Percentage of agents switching guess between prd, prd-1 (df_switch) ---------------------------
df_switch <- data.frame(type=NA, group=NA, round=NA, period = NA, switch =NA )
for (ntype in c(n_10,n_40)){
  for (prd in c(2:12)){
    for (grp in c(1:max_grp)){
      for (rnd in c(1:max_rnd)){
        x <- 
          filter(df_detail, type==ntype, group==grp, round==rnd, period ==prd)$guess[[1]] - 
          filter(df_detail, type==ntype, group==grp, round==rnd, period ==(prd-1))$guess[[1]]
        df_switch <- rbind(df_switch, list(ntype, grp, rnd, prd, mean(abs(x))))
      }
    }
  }
}
df_switch <- df_switch[-1,]
df_switch$typelist <- str_split(df_switch$type, "_", simplify = TRUE)[,2]
df_switch$typelist <- as.factor(as.character(df_switch$typelist))

# Guess signal instead of DeGroot, when DeGroot contradict with your signal (df_guess_signal) ---------------------------

df_guess_signal <- 
  data.frame(group = NA, round = NA, period = NA, type =NA, DeGroot =NA, id_subj = NA, id_grp =NA)

for (ntype in c(n_10,n_40)){
  g <- glist[ntype][[1]]
  for (grp in c(1:max_grp)){
    for (rnd in c(1:max_rnd)){
      df_temp <- 
        filter(df_tot, type==ntype, group == grp, round == rnd, period ==1)
      signal <- df_temp[order(df_temp$id_grp),]$sig -1
      for (prd in c(1:(max_prd-1))){
        df_temp <- df_temp
        s <- as.numeric(df_temp$guess[[1]]) #Guess vector from data at period t
        s_next <- as.numeric(filter(df_detail, type==ntype & group == grp & round == rnd & period == (prd+1))$guess[[1]]) #Guess vector from data at period t+1
        predict <- DeGroot_once(g,s,w,tre)
        subj <- df_temp$id_subj[[1]]
        fit <- (1-as.numeric(s_next==predict))*(s_next==signal)
        id <- as.numeric(df_temp$id_grp[[1]])

        for(j in id){
          df_guess_signal <- rbind(df_guess_signal,list(grp, rnd, prd+1, ntype, fit[j],subj[j], j)) #check if the simulation matches with guess at period t+1
        }
      }
    }
  }
}
df_guess_signal <- df_guess_signal[-1,]
df_guess_signal_n_10 <- filter(df_guess_signal, type %in% n_10, id_subj!=0)
df_guess_signal_n_40 <- filter(df_guess_signal, type %in% n_40, id_subj!=0)

# Saving environment for later use
save.image("RData/05-degroot-envir.RData")


# *** Part 2 - Summarises and prints output -----------------------

load("RData/05-degroot-envir.RData")

# Print output in bar chart comparing the match -----------------------
# between how "data" and "guessing signal" matches with DeGroot
for (ntypelist in c("n_10","n_40")){
  final <- ggplot(filter(df_match, type %in% get(ntypelist)), 
                  aes(x=type, y=percentage, group=match, fill=str_wrap(match,15))) +
    geom_bar(position="dodge", stat="identity", color = "white", width=0.6) +
    geom_text(aes(label=scales::percent(percentage,1), group=match), 
              vjust = -0.3, position = position_dodge(width = 0.6)) + 
    
    scale_y_continuous(labels=scales::percent, limits = c(0,1)) + 
    coord_cartesian(ylim=c(0.5,1)) +
    
    ggtitle(paste0("Data matching with decision rules, ",ntypelist))+
    xlab("Network type") + 
    ylab("Percentage")+
    
    theme_bw()+
    theme(legend.title = element_blank(),legend.key.height=unit(1.5, "cm"))
  tname <- paste0("../output/DeGroot_match,", ntypelist,".pdf")
  #pdf(tname, width=7, height=5)
  print(final)
  #dev.off()
}


# Analyze % of data matching DeGroot per network over rounds ---------------------------
df_plot=data.frame(type = NA, round = NA, DeGroot_per =NA)
for (ntype in c(n_10_all,n_40)){
  for (rnd in c(1:max_rnd)){
    sum_sig<- sum(as.numeric(lapply(filter(df_DeGroot, type ==ntype & round ==rnd)$DeGroot,sum)))
    num_sig<- sum(as.numeric(lapply(filter(df_DeGroot, type ==ntype & round ==rnd)$DeGroot,length)))
    df_plot<-rbind(df_plot,list(ntype, rnd, signif(sum_sig/num_sig,3)))
  }
}
df_plot <- df_plot[-1,]
df_plot$type <- as.factor(as.character(df_plot$type))

# Print timeseries: % of data matching DeGroot per network over rounds
for (ntypelist in c("n_10", "n_40")){
  final <- ggplot(subset(df_plot, type %in% get(ntypelist)), 
                  aes(x=round, y=DeGroot_per, group=type)) +
    geom_line(aes(linetype=type)) +
    geom_point(aes(shape=type, color=type), size = 3, stroke = 1) +
    scale_shape_manual(values=c(15,16,17,18))+
    
    scale_x_continuous(name="Round", 
                       limits=c(1, max_rnd), 
                       breaks = c(0:max_rnd), 
                       expand = c(0.02, 0) ) +
    scale_y_continuous(name="Fraction of guesses matching DeGroot", 
                       limits=c(0.70,1), 
                       expand = c(0.02, 0)) +
    
    coord_cartesian(clip = 'off') +
    theme_bw() +
    theme(legend.title = element_blank(), 
          legend.justification=c(1,1), 
          legend.position=c(.9,0.3))+
    ggtitle(paste0("DeGroot percentages per network by round, ", ntypelist))
  
  #final <- final + geom_hline(yintercept=avg_sig,linetype="dashed", color = "red")
  
  tname <- paste0("../output/DeGroot_percentage_by_rounds,", ntypelist,".pdf")
  #pdf(tname, width=7, height=5)
  print(final)
  #dev.off()
}

# Analyze % of data matching DeGroot per network over periods ---------------------------
df_plot=data.frame(type = NA, period = NA, DeGroot_per =NA)
for (ntype in c(n_10_all,n_40)){
  for (prd in c(2:max_prd)){
    sum_sig<- sum(as.numeric(lapply(filter(df_DeGroot, type ==ntype & period ==prd)$DeGroot,sum)))
    num_sig<- sum(as.numeric(lapply(filter(df_DeGroot, type ==ntype & period ==prd)$DeGroot,length)))
    df_plot<-rbind(df_plot,list(ntype, prd, signif(sum_sig/num_sig,3)))
  }
}
df_plot <- df_plot[-1,]
df_plot$type <- as.factor(as.character(df_plot$type))

# Print timeseries: % of data matching DeGroot per network over periods 
for (ntypelist in c("n_10", "n_40")){
  final <- ggplot(subset(df_plot, type %in% get(ntypelist)), aes(x=period, y=DeGroot_per, group=type)) +
    geom_line(aes(linetype=type)) +
    geom_point(aes(shape=type, color=type), size = 3, stroke = 1) +
    scale_shape_manual(values=c(15,16,17,18))+
    
    scale_x_continuous(name="Period", limits=c(2, 12), breaks = c(0:12), expand = c(0.02, 0) ) +
    scale_y_continuous(name="fraction of DeGroot guesses", limits=c(0.7, 1), expand = c(0.02, 0)) +
    
    coord_cartesian(clip = 'off') +
    theme_bw() +
    theme(legend.title = element_blank(), legend.justification=c(1,1), legend.position=c(.9,0.4))+
    ggtitle(paste0("DeGroot percentages by period, ", ntypelist))
  
  #final <- final + geom_hline(yintercept=avg_sig,linetype="dashed", color = "red")
  
  tname <- paste0("../output/DeGroot_percentage_by_periods,", ntypelist,".pdf")
  #pdf(tname, width=7, height=5)
  print(final)
  #dev.off()
}

# Analyze % of data matching DeGroot by network position ---------------------------

for (n in c(10,40)){
  df_netpos = matrix(data=NA, nrow=4, ncol =n) 
  rownames(df_netpos) <- get(paste0("n_",n))
  for (ntype in get(paste0("n_",n))){
    df_temp <- filter(df_DeGroot, type ==ntype)
    df_temp2 <-rep(0,n) 
    for (i in c(1:nrow(df_temp))){
      # sum up all rows: whether each position's guess match with DeGroot
      df_temp2 <- df_temp2 + df_temp$DeGroot[[i]]
    }
    df_netpos[ntype,] <- signif(df_temp2/nrow(df_temp),3) #average all rows
    
    g <- graph_from_adjacency_matrix(as.matrix(glist[ntype][[1]])-diag(n))
    #pdf(paste0("../output/DeGroot_percent,",ntype, ".pdf"), width=9, height=9)
    plot(g, edge.arrow.size= if(grepl("RF",ntype)){0.4} else{0},
         main = ntype,
         vertex.size=7,
         vertex.color="#4f85db",
         vertex.label= df_netpos[ntype,],
         vertex.label.dist=1, vertex.label.cex=1,
         vertex.label.color= "black",
         vertex.label.degree=0,
         #layout= layout_in_circle,
         layout= layout_with_kk)
    #dev.off()
  }
}

# Analyze % of data matching DeGroot per network by person id_subj ---------------------------

{
  df_temp <- filter(df_DeGrootmore, type %in% n_10)
  agg <- aggregate(df_temp$DeGroot, list(df_temp$id_subj, df_temp$type), 
                   drop=FALSE, FUN=function(x) mean = mean(x))
  df_temp <- do.call(data.frame, agg)
  colnames(df_temp) <- c("id_subj","type","mean") 
  df_temp <- filter(df_temp, id_subj !=0 )
  df_temp <- df_temp %>% group_by(type) %>% mutate(avg.mean = mean(mean))
  df_plot <- filter(df_temp, id_subj !=0 )
  
  df_temp <- filter(df_DeGrootmore, type %in% n_40)
  agg <- aggregate(df_temp$DeGroot, list(df_temp$id_subj, df_temp$type, df_temp$group), 
                   drop=FALSE, FUN=function(x) mean = mean(x))
  df_temp <- do.call(data.frame, agg)
  colnames(df_temp) <- c("id_subj","type","group","mean") 
  df_temp <- filter(df_temp, id_subj !=0 )
  df_temp <- df_temp %>% group_by(type) %>% mutate(avg.mean = mean(mean))
  df_plot2 <- filter(df_temp, id_subj !=0 )
  
  final <- ggplot(df_plot, aes(x=mean)) +
    geom_histogram(breaks = seq(0,1,0.05), color="white", closed = "left")+
    stat_bin(breaks = seq(0,1,0.05), 
             closed = "left", 
             geom="text", 
             colour="black", 
             size=2.5, 
             aes(label= ifelse(..count.. > 0, ..count.., ""), group=type, y=1+(..count..))) +
    #geom_vline(aes(xintercept=avg.mean, group=type, color="mean")) +
    #geom_vline(data=filter(df_match, match=="Guessing signal match with DeGroot", type %in% n_10), aes(xintercept=percentage, group=type, color = "Guessing signal match with DeGroot"), linetype="dashed")+
    
    facet_wrap(~ factor(type, levels=n_10))+
    #ylim(NA, 20) +
    scale_colour_manual(values = c("#00bfc4","#f8766d") )+
    scale_x_continuous(breaks=seq(0,1,0.1), 
                       limits = c(-0.05, 1.05), 
                       labels = scales::percent_format(accuracy = 1))+
    
    ggtitle(paste0("Distribution of DeGroot match across networks, n_10"))+
    xlab("Percentage of guesses matching with DeGroot") + 
    ylab("Frequency of subjects (out of 40)")+
    
    theme_bw() +
    theme(legend.title = element_blank(), 
          strip.text.x = element_text(size = 13), 
          legend.justification=c(1,1), 
          legend.text = element_text(size=10), 
          legend.position=c(0.3,0.92))
  
  tname <- paste0("../output/DeGroot_match_distr,n_10.pdf")
  #pdf(tname, width=10, height=7)
  print(final)
  #dev.off()
  
  final <- ggplot(df_plot2, aes(x=mean)) +
    geom_histogram(breaks = seq(0,1,0.05), color="white", closed = "left")+
    stat_bin(breaks = seq(0,1,0.05), 
             closed = "left", 
             geom="text", 
             colour="black", 
             size=2.5, 
             aes(label= ifelse(..count.. > 0, ..count.., ""), group=type, y=1+(..count..))) +
    #geom_vline(aes(xintercept=avg.mean, group=type, color="mean")) +
    geom_vline(data=filter(df_match, 
                           match=="Guessing signal match with DeGroot", 
                           type %in% n_40),
               aes(xintercept=percentage, 
                   group=type, 
                   color = "Guessing signal match with DeGroot"), 
               linetype="dashed")+
    
    facet_wrap(~ factor(type, levels=n_40))+
    #ylim(NA, 60) +
    scale_colour_manual(values = c("#00bfc4","#f8766d"))+
    scale_x_continuous(breaks=seq(0,1,0.1), 
                       limits = c(-0.05, 1.05), 
                       labels = scales::percent_format(accuracy = 1))+
    
    ggtitle(paste0("Distribution of DeGroot match across networks, n_40"))+
    xlab("Percentage of guesses matching with DeGroot") + 
    ylab("Frequency of subjects (out of 160)")+
    theme_bw() +
    theme(legend.title = element_blank(), 
          strip.text.x = element_text(size = 13), 
          legend.justification=c(1,1), 
          legend.text = element_text(size=10), 
          legend.position=c(0.3,0.92))
  
  tname <- paste0("../output/DeGroot_match_distr,n_40.pdf")
  #pdf(tname, width=10, height=7)
  print(final)
  #dev.off()
}

# Analyze % of guessing signal against DeGroot per network by person id_subj ---------------------------
# Guessing signal when DeGroot predict otherwise

{
  agg <- aggregate(df_guess_signal_n_10$DeGroot, list(df_guess_signal_n_10$type, df_guess_signal_n_10$group, df_guess_signal_n_10$id_subj), FUN=function(x) c(mean = mean(x), sd = sd(x)))
  df_plotn_10 <- do.call(data.frame, agg)
  colnames(df_plotn_10) <- c("type","group","id_subj","mean","sd")
  #df_plot <- filter(df_plot, mean!=0)

  agg <- aggregate(df_guess_signal_n_40$DeGroot, list(df_guess_signal_n_40$type, df_guess_signal_n_40$group, df_guess_signal_n_40$round, df_guess_signal_n_40$id_subj), FUN=function(x) c(mean = mean(x), sd = sd(x)))
  df_plotn_40 <- do.call(data.frame, agg)
  colnames(df_plotn_40) <- c("type","group","round","id_subj","mean","sd")
  #df_plot2 <- filter(df_plot2, mean!=0)

  for (ntypelist in c("n_10","n_40")){
    final <- ggplot(get(paste0("df_plot",ntypelist)), aes(x=mean)) +
      geom_histogram(breaks = seq(0,1,0.1), color="white", closed = "left")+
      stat_bin(breaks = seq(0,1,0.1), closed = "left", geom="text", colour="black", size=2.5, aes(label= ifelse(..count.. > 0, ..count.., ""), group=type, y=1+(..count..))) +
      #geom_vline(aes(xintercept=avg.mean, group=type, color="mean")) +
      #geom_vline(data=filter(df_match, match=="Guessing signal match with DeGroot", type %in% n_10), aes(xintercept=percentage, group=type, color = "Guessing signal match with DeGroot"), linetype="dashed")+

      facet_wrap(~ factor(type, levels=get(ntypelist)))+
      #ylim(NA, 25) +
      scale_colour_manual(values = c("#00bfc4","#f8766d"))+
      scale_x_continuous(breaks=seq(0,1,0.1), limits = c(-0.05, 1.05), labels = scales::percent_format(accuracy = 1))+
      ggtitle(paste0("Distribution of Data matching with Guessing signal across networks, ",ntypelist))+
      xlab("Percentage of guesses matching with guessing signal, against DeGroot") + ylab("Frequency of subjects")+
      theme_bw() +
      theme(legend.title = element_blank(), strip.text.x = element_text(size = 13), legend.justification=c(1,1), legend.text = element_text(size=10), legend.position=c(0.3,0.92))

    tname <- paste0("../output/DeGroot_match_guesssig_distr,",ntypelist,".pdf")
    #pdf(tname, width=10, height=7)
    print(final)
    #dev.off()
  }
}

# DeGroot simulation results of signals for 12 periods ---------------------------

# Import signal
balls <- as.matrix(read_csv("Signals_n_10.csv", col_names = FALSE))-1
for (ntype in n_10){
  g <- glist[ntype][[1]]
  xy <- data.frame()
  # For each set of signals, loop
  for (i in c(1:nrow(balls))){
    s<-balls[i,]
    #Set majority signal as truth
    if (mean(s)<0.5){s<-1-s}
    # Perform DeGroot_once for 12 periods, iteratively
    record<-DeGroot(g,s,12,tremble=tre)
    s<-record[nrow(record),]
    print(s)
    xy<- rbind(xy,mean(s))
  }
  colnames(xy) <- "weight"
  final <- ggplot(xy, aes(x=weight)) + geom_histogram(binwidth=0.05)+ylim(NA, nrow(balls)) +
    scale_x_continuous(breaks=seq(0,1,0.1), limits = c(-0.05, 1.05))+
    ggtitle(paste0("DeGroot simulation, majority signal as truth, ",ntype)) +
    xlab("Average guess at period 12") + ylab("Frequency of networks (out of 24)")
  #pdf(paste("Simulation_DeGroot_majorsig", ntype, ".pdf"), width=7, height=5)
  print(final)
  #dev.off()
}

for (tre in c(0)){
  for (ntype in n_10){
    g <- glist[ntype][[1]]
    df_plot=data.frame(period = NA, mean = NA, signal = NA)
    for (i in c(1:nrow(balls))){
      s<-balls[i,]
      # Set majority signal as truth
      if (mean(s)<0.5){s<-1-s}
      # Perform DeGroot_once for 12 periods, iteratively
      record<-DeGroot(g,s,max_prd,tremble=tre)
      
      for (prd in c(1:max_prd)){
        x <- mean(record[prd,])
        df_plot<-rbind(df_plot,list(prd,x,i))
      }
    }
    df_plot <- df_plot[-1,]
    df_plot$signal<- as.factor(as.character(df_plot$signal))
    
    {final <- ggplot(df_plot, aes(x=period, y=mean, group=signal)) +
        geom_line(aes(color=signal)) +
        geom_point(aes(shape=signal, color=signal), size = 1, stroke = 1) +
        scale_shape_manual(values=c(1:24))+
        
        scale_x_continuous(name="Period", 
                           limits=c(1, 12), 
                           breaks = c(0:12), 
                           expand = c(0.02, 0) ) +
        scale_y_continuous(name="fraction of players", 
                           limits=c(0, 1), 
                           expand = c(0.02, 0)) +
        
        coord_cartesian(clip = 'off') +
        theme_bw() +
        theme(legend.title = element_blank(), 
              legend.justification=c(1,1), 
              legend.position=c(0.9,0.4),
              legend.key.size = unit(0.4, 'cm'))+
        
        ggtitle(paste0("Simulation",ntype))
      
      #final <- final + geom_hline(yintercept=avg_sig,linetype="dashed", color = "red")
    }
    tname <- paste0("Simulation_DeGroot_majorsig ",ntype,"_extensive.pdf")
    #pdf(tname, width=10, height=7)
    print(final)
    #dev.off()
  }
}

# Clan DeGroot checking ---------------------------
# # import list of clans of each ntype
# df_clan <- read_delim("C:/Users/tyyto/Desktop/Social Learning/Simulation_DeGroot/Simulation_clan.csv", ";", 
#                       escape_double = FALSE, 
#                       col_types = cols(X1 = col_skip()), 
#                       trim_ws = TRUE)
# 
# df_stuckness <- data.frame(type=NA, group=NA, round=NA, from_period = NA, clan =NA, guess = NA, stuck=NA )
# for (ntype in c(n_10,n_40)){
#   df_temp <- filter(df_clan, type==ntype)
#   
#   if (nrow(df_temp) !=0){ #if there is a non-empty list of clans in ntype
#     for(i in c(1: nrow(df_temp))){
#       cl <- as.numeric(strsplit(gsub('c|(|)', '', df_temp$clan[i]), ',')[[1]]) #split string of clans
#       
#       for (grp in c(1:max_grp)){
#         for (rnd in c(1:max_rnd)){
#           df_temp2 <- filter(df_detail, type==ntype, group==grp, round==rnd)
#           m_temp <- matrix(unlist(df_temp2$guess), ncol = length(df_temp2$guess[[1]]), byrow = TRUE) #matrix of individual guesses over 12 periods
#           #print(m_temp[,cl])
#           for (prd in c(1:12)){
#             m_temp2 <- matrix(m_temp[,cl][prd:12,], ncol = length(cl)) #matrix of guesses of clan from period t onwards
#             tru <- filter(df_truth, type ==ntype, group==grp, round==rnd)$truth #the true state for the group-round
#             
#             #if all guesses at period t is the same and on the wrong state, and if across all periods, all guesses don't change
#             if (all(m_temp2[1,] !=tru) & all(apply(m_temp2, 2, function(x) length(unique(x)) == 1) == TRUE)){ 
#               df_stuckness <- rbind(df_stuckness, list(ntype, grp, rnd, prd, list(cl), list(m_temp2[1,]), "stuck"))
#             }
#             else if (all(m_temp2[1,] !=tru)){ #elseif all guesses at period t is the same and on the wrong state 
#               df_stuckness <- rbind(df_stuckness, list(ntype, grp, rnd, prd, list(cl), list(m_temp2[1,]), "wrong"))
#             }
#             else{
#               df_stuckness <- rbind(df_stuckness, list(ntype, grp, rnd, prd, list(cl), list(m_temp2[1,]), "other"))
#             }
#           }
#         }
#       }
#     }
#   }
# }
# df_stuckness <- df_stuckness[-1,]
# 
# df_temp <- df_stuckness[c(lapply(df_stuckness$clan, length) <= 5),]
# 
# agg <- aggregate(df_temp$stuck,  
#                  list(df_temp$type, df_temp$stuck), 
#                  drop=FALSE,
#                  FUN=function(x) count = length(x))
# df_temp<- do.call(data.frame, agg)
# colnames(df_temp) <- c("type","stuckornot","count")
# 
# for (ntype in c(n_10,n_40)){
#   x <- filter(df_temp, type == ntype, stuckornot == "stuck")$count
#   y <- filter(df_temp, type == ntype, stuckornot == "wrong")$count
#   z <- filter(df_temp, type == ntype, stuckornot == "other")$count
#   
#   print(c(ntype,percent(x/(x+y)),percent(x/(x+y+z)),percent((x+y)/(x+y+z))))
# }
# 

