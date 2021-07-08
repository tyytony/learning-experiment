# To print all output into "output" folder
# Use ctrl-F to:
# replace "#pdf(" with "pdf("
# replace "#dev.off()" with "dev.off()"
# replace "#sink(" with "sink("

setwd(getwd())

library(igraph)
library(readr)
library(readxl)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(scales)
library(miceadds)
library(texreg)
library(Zelig)
options(digits = 4)
options(max.print = 100)
source("function/Function-DeGroot_once,DeGroot,Simulation.R")
source("function/Function-Bayesian.R")

# *** Part 1 - Creates dataframes from df_tot ---------------------
load("RData/01-process_data-envir.RData")
load("RData/02-analysis-envir.RData")

# forming period by period guess vector from data (df_detail, df_truth)  ---------------------------
df_detail <- data.frame(group = NA, round = NA, period = NA, type =NA, 
                        guess =NA, id_subj = NA, id_grp = NA)
df_truth <- data.frame(group = NA, round = NA, type =NA, truth = NA)

for (ntype in c(n_10,n_40)){
  for (grp in c(1:max_grp)){
    for (rnd in c(1:max_rnd)){
      
      df_temp <- filter(df_tot, type==ntype, group == grp, round == rnd)
      df_temp2 <- filter(df_temp, period == 1)
      # Majority of signal from group-round: 1 if mean signal >=0.5
      tru <- round(mean(df_temp2$sig), 0) - 1 
      df_truth <- rbind(df_truth, c(grp, rnd, ntype, tru))
      
      for (prd in c(0:max_prd)){
        df_temp2 <- df_temp %>% filter(period == max(prd, 1))
        df_temp2 <- df_temp2[order(df_temp2$id_grp),] #reorder rows by group id 1-n
        subj <- df_temp2$id_subj #the subjects, in order
        location <- df_temp2$id_grp #the locations of subjects, in order
        
        if (prd == 0){
          if (tru == 0) {s <- 1- (df_temp2$sig - 1)}
          else {s <- df_temp2$sig - 1} #signals 
          
        }
        else{
          if (tru == 0) {s <- 1 - (df_temp2$choice - 1)}
          else {s <- df_temp2$choice - 1} #guesses at prd
           
        }

        df_detail <- df_detail %>%
          rbind(list(grp, rnd, prd, ntype, list(s), list(subj), list(location)))
      }
    }
  }
}
df_detail <- df_detail[-1,]
df_truth <- df_truth[-1,]

# performing DeGroot of next period on guess vector using network glist (df_DeGroot, df_DeGrootmore)  ---------------------------
df_DeGroot <- data.frame(group = NA, round = NA, period = NA, 
                      type =NA, DeGroot =NA, id_subj = NA)

for (ntype in c(n_10,n_40)){
  g <- glist[[ntype]]
  for (grp in c(1:max_grp)){
    for (rnd in c(1:max_rnd)){
      for (prd in c(0:(max_prd-1))){
        df_temp <- filter(df_detail, type == ntype, group == grp, round == rnd)
        # Guess vector from data at period t
        s <- as.numeric(filter(df_temp, period == prd)$guess[[1]])
        # Guess vector from data at period t+1
        s_next <- as.numeric(filter(df_temp, period == (prd+1))$guess[[1]])
        # Simulation from guess vector at period t for next period
        if (prd == 0){
          predict <- s #at period 0, DeGroot should guess signal
        }
        else {
          predict <- DeGroot_once(g, s, w = 0, tre = 0) 
        }
        fit <- 1 - abs(s_next - predict) #check if simulation matches with guess at period t+1
        subj <- filter(df_temp, period == prd)$id_subj[[1]] #subjects
        df_DeGroot <- df_DeGroot %>%
          rbind(list(grp, rnd, prd+1, ntype, 
                     list(fit), 
                     list(subj)))
      }
    }
  }
}
df_DeGroot <- df_DeGroot[-1,]

df_temp <- df_DeGroot %>% 
  separate(type, into = c(NA, "size"), sep = -2, remove = FALSE) %>%
  mutate(size = as.numeric(size)) %>%
  select(group, round, period, type, size)

df_DeGrootmore <- as.data.frame(lapply(df_temp, rep, df_temp$size))
df_DeGrootmore <- df_DeGrootmore %>% 
  cbind(DeGroot = unlist(df_DeGroot$DeGroot), id_subj = unlist(df_DeGroot$id_subj)) %>%
  select(-"size")

# Repeat for just guessing signal, how well that matches with DeGroot (df_sigDeGroot, df_sigDeG_more)  ---------------------------

df_sigDeGroot <- data.frame(group = NA, round = NA, period = NA, 
                            type =NA, DeGroot =NA, id_subj = NA)
for (ntype in c(n_10,n_40)){
  g <- glist[[ntype]]
  for (grp in c(1:max_grp)){
    for (rnd in c(1:max_rnd)){
      for (prd in c(0:(max_prd-1))){
        df_temp <- filter(df_detail, type == ntype, group == grp, round == rnd)
        # Signal vector at period 0
        sig <- as.numeric(filter(df_temp, period == 0)$guess[[1]])
        # Guess vector from data at period t+1
        s_next <- as.numeric(filter(df_temp, period == (prd+1))$guess[[1]])
        # Simulation from guess vector at period t for next period
        predict <- sig #at period prd, agents guess signal
        fit <- 1 - abs(s_next - predict) #check if simulation matches with guess at period t+1
        subj <- filter(df_temp, period == prd)$id_subj[[1]] #subjects
        df_sigDeGroot <- df_sigDeGroot %>%
          rbind(list(grp, rnd, prd+1, ntype, 
                     list(fit),
                     list(subj)))
      }
    }
  }
}
df_sigDeGroot <- df_sigDeGroot[-1,]

df_temp <- df_sigDeGroot %>% 
  separate(type, into = c(NA, "size"), sep = -2, remove = FALSE) %>%
  mutate(size = as.numeric(size)) %>%
  select(group, round, period, type, size)

df_sigDeG_more <- as.data.frame(lapply(df_temp, rep, df_temp$size))
df_sigDeG_more <- df_sigDeG_more %>% 
  cbind(DeGroot = unlist(df_sigDeGroot$DeGroot), id_subj = unlist(df_sigDeGroot$id_subj)) %>%
  select(-"size")

# Analyze per network type's percentage of DeGroot match, and other matches (df_match)  ---------------------------
df_temp <- df_DeGrootmore %>% 
  filter(type %in% c(n_10,n_40)) %>%
  group_by(type, .drop = FALSE) %>% 
  summarise(percentage = mean(DeGroot),
            count = n()) %>% 
  mutate(sd = (percentage * (1 - percentage) / count)^(0.5),
         match = "Data match with DeGroot")

df_temp2 <- df_sigDeG_more %>% 
  filter(type %in% c(n_10,n_40)) %>%
  group_by(type, .drop = FALSE) %>% 
  summarise(percentage = mean(DeGroot),
            count = n()) %>% 
  mutate(sd = (percentage * (1 - percentage) / count)^(0.5),
         match = "Guessing signal match with DeGroot")

df_match <- rbind(df_temp, df_temp2)

# Percentage of agents switching guess between prd, prd-1 (df_switch) ---------------------------
df_switch <- data.frame(type = NA, group = NA, round = NA, 
                        period = NA, switch = NA)
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
  data.frame(group = NA, round = NA, period = NA, type =NA, 
             DeGroot =NA, id_subj = NA, id_grp =NA)

for (ntype in c(n_10,n_40)){
  g <- glist[[ntype]]
  for (grp in c(1:max_grp)){
    for (rnd in c(1:max_rnd)){
      for (prd in c(1:(max_prd-1))){
        df_temp <- filter(df_detail, type == ntype, group == grp, round == rnd)
        # Signal vector at period 0
        sig <- as.numeric(filter(df_temp, period == 0)$guess[[1]])
        
        # Guess vector from data at period t
        s <- as.numeric(filter(df_temp, period == prd)$guess[[1]])
        # Guess vector from data at period t+1
        s_next <- as.numeric(filter(df_temp, period == (prd+1))$guess[[1]])
        # Simulation from guess vector at period t for next period
        predict <- DeGroot_once(g, s, w = 0, tre = 0)
        
        # Only care about the cases where signal against DeGroot
        x <- abs(sig - predict) #indicator of signal and prediction not matching
        x[x == 0 | x == 0.5] <- NA #remove cases where sig and predict match, or could match 50% of the time

        # When DeGroot against signal, follow signal
        fit <- x * (s_next == sig)
        
        subj <- filter(df_temp, period == prd)$id_subj[[1]] #subjects
        id <- as.numeric(filter(df_temp, period == prd)$id_grp[[1]]) #subject location in network
        df_guess_signal <- df_guess_signal %>%
          rbind(list(grp, rnd, prd+1, ntype, list(fit), list(subj), list(id)))
      }
    }
  }
}
df_guess_signal <- df_guess_signal[-1,]

df_temp <- df_guess_signal %>% 
  separate(type, into = c(NA, "size"), sep = -2, remove = FALSE) %>%
  mutate(size = as.numeric(size)) %>%
  select(group, round, period, type, size)

df_gsig_more <- as.data.frame(lapply(df_temp, rep, df_temp$size))
df_gsig_more <- df_gsig_more %>% 
  cbind(DeGroot = unlist(df_guess_signal$DeGroot), id_subj = unlist(df_guess_signal$id_subj)) %>%
  select(-"size")

df_gsig_n_10 <- filter(df_gsig_more, type %in% n_10, id_subj!=0)
df_gsig_n_40 <- filter(df_gsig_more, type %in% n_40, id_subj!=0)

# Share of time info. dominated fail to copy dominant (df_leader) ---------------------------

df_leader <- data.frame(type = NA, leader = NA, follower = NA)
for (ntype in c(n_10,n_40)){
  g <- glist[[ntype]]
  for (i in 1:(ncol(g)-1)){
    for (j in (i+1):ncol(g)){
      if (!is.na(match(-1,g[i,]-g[j,])) & #we can find someone that j observes but i doesn't, and
          is.na(match(-1,g[j,]-g[i,]))){  #we cannot find someone that i observes but j doesn't
        df_leader <- rbind(df_leader, c(ntype, j,i)) #then j is an information leader of i
      }
      if (!is.na(match(-1,g[j,]-g[i,])) & #vice versa
          is.na(match(-1,g[i,]-g[j,]))){  
        df_leader <- rbind(df_leader, c(ntype, i,j)) #then i is an information leader of j
      }
    }
  }
}
df_leader <- df_leader[-1,]

# if i is a leader of j, but a follower of k
# j should remove i as a leader and follow k instead
# note that this is a directed graph
for (ntype in c(n_10,n_40)){

  df_temp <- filter(df_leader, type == ntype)
  # el <- cbind(df_temp$follower,df_temp$leader) #create edgelist of follower to leader
  # g <- graph.edgelist(el)
  # plot(g)
  
  # remove leaders that are followers of others
  for (i in df_temp$follower){
    df_leader <- df_leader[!(df_leader$type == ntype & df_leader$leader == i), ]
  }
  
  df_temp <- filter(df_leader, type == ntype)
  
  # people who follow more than one leader, should choose one leader by degree
  n_occur <- data.frame(table(df_temp$follower))
  multi_follow <- df_temp[df_temp$follower %in% n_occur$Var1[n_occur$Freq > 1],] #dataframe of people following more than one person
  
  if(nrow(multi_follow)==0){next} #If no more agent following multiple people, then move on

  for (j in n_occur[n_occur$Freq > 1,]$Var1){
    df_temp2 <- filter(df_temp, follower==j)
    x <- degree(graph.adjacency(glist[[ntype]]),df_temp2$leader) #named vector of degrees of j's leaders
    y <- names(x)[-(which.max(x))] #name of the neighbour with non-max degree
    df_leader <- df_leader[!(df_leader$type == ntype & #remove other leaders
                             df_leader$follower == j &
                             df_leader$leader %in% y) , ]
  }
  
  df_temp <- filter(df_leader, type == ntype)
  
  n_occur <- data.frame(table(df_temp$follower))
  multi_follow <- df_temp[df_temp$follower %in% n_occur$Var1[n_occur$Freq > 1],] #dataframe of people following more than one person
  if(nrow(multi_follow)==0){next} #If no more agent following multiple people, then move on
  else{ #otherwise, stop and show where it went wrong
    print(paste0("Oh no...", ntype))
    break
  }
}

# Following leader instead of DeGroot, when DeGroot contradict with your leader (df_info_dom) ---------------------------

df_info_dom <-
  data.frame(group = NA, round = NA, period = NA, type =NA, 
             DeGroot =NA, id_subj = NA, id_grp =NA)

for (ntype in c(n_10,n_40)){
  g <- glist[[ntype]]
  
  # leader-follower graph - who i follows in network g
  df_temp2 <- filter(df_leader, type == ntype)
  el <- cbind(df_temp2$follower,df_temp2$leader) #create edgelist of follower to leader
  m <- as.matrix(as_adjacency_matrix(make_graph(el)))
  
  # Zero matrix, replace the row and columns that has data from m
  m0 <- matrix(0, ncol(g), ncol(g))
  location_mat <- cbind(
    as.numeric(rep(rownames(m), ncol(m))), 
    as.numeric(rep(colnames(m), each=nrow(m))) 
  )
  m0[location_mat] <- as.vector(m)
  
  # if follower has no leader, we don't care about them for now
  m0[rowSums(m0)==0,] <- NA
  
  for (grp in c(1:max_grp)){
    for (rnd in c(1:max_rnd)){
      for (prd in c(1:(max_prd-1))){
        df_temp <- filter(df_detail, type == ntype, group == grp, round == rnd)
        # Guess vector from data at period t
        s <- as.numeric(filter(df_temp, period == prd)$guess[[1]])
        # Guess vector from data at period t+1
        s_next <- as.numeric(filter(df_temp, period == (prd+1))$guess[[1]])
        # Simulation from guess vector at period t for next period
        predict <- DeGroot_once(g, s, w = 0, tre = 0)

        # What you would guess at t+1 if you followed your leader's guess at t
        s_leader <- c(m0 %*% s)
          
        # Only care about the cases where following leader is against DeGroot
        x <- abs(s_leader - predict) #indicator of signal and prediction not matching
        x[x == 0 | x == 0.5] <- NA #remove cases where sig and predict match, or could match 50% of the time
        
          
        # When DeGroot against following leader, follow leader
        fit <- x * (s_next == s_leader)
        
        subj <- filter(df_temp, period == prd)$id_subj[[1]] #subjects
        id <- as.numeric(filter(df_temp, period == prd)$id_grp[[1]]) #subject location in network
        df_info_dom <- df_info_dom %>%
          rbind(list(grp, rnd, prd+1, ntype, list(fit), list(subj), list(id)))
      }
    }
  }
}
df_info_dom <- df_info_dom[-1,]

df_temp <- df_info_dom %>% 
  separate(type, into = c(NA, "size"), sep = -2, remove = FALSE) %>%
  mutate(size = as.numeric(size)) %>%
  select(group, round, period, type, size)

df_info_dom_more <- as.data.frame(lapply(df_temp, rep, df_temp$size))
df_info_dom_more <- df_info_dom_more %>% 
  cbind(DeGroot = unlist(df_info_dom$DeGroot), id_subj = unlist(df_info_dom$id_subj)) %>%
  select(-"size") %>%
  filter(id_subj != 0) #remove bots

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
              vjust = -0.8, position = position_dodge(width = 0.6)) +
    geom_errorbar(aes(ymin=percentage-sd, ymax=percentage+sd), 
                  width = .2, position = position_dodge(0.6)) +
    
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

# Analyze % of data matching DeGroot per network over periods (Match_mistake.tex) ---------------------------
df_plot <- data.frame(type = NA, period = NA, DeGroot_per =NA)
for (ntype in c(n_10_all,n_40)){
  for (prd in c(1:max_prd)){
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
    
    scale_x_continuous(name="Period", limits=c(1, 12), breaks = c(0:12), expand = c(0.02, 0) ) +
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

# Binomial probability of period 1 and 2 guess goes against Bayesian and DeGroot (Match_mistake.tex)
for (ntypelist in c("n_10", "n_40")){
  df_temp <- filter(df_DeGrootmore, type %in% get(ntypelist), period <= 2, DeGroot != 0.5)
  df_reg <- df_temp %>%  #Recoding type, size, group, round
    mutate(DeGroot = 1 - DeGroot, #indicator for mistake made
           typesize = df_temp$type,
           type = str_split_fixed(df_temp$type, "_",2)[,1],
           size = str_split_fixed(df_temp$type, "_",2)[,2],
           typesizegroup = paste0(as.character(type), 
                                  as.character(size), 
                                  as.character(group),
                                  sep = ''),
           typesizegroupround = paste0(as.character(type),
                                       as.character(size),
                                       as.character(group),
                                       as.character(round),sep = ''),
    ) %>%
    mutate_at(vars(group, round, type, size, 
                   typesizegroup, typesizegroupround), as.factor) %>%
    filter(!is.na(DeGroot))
  
  lin.1 <- lm.cluster(data = df_reg, formula = DeGroot ~ typesize, cluster = "typesizegroup")
  logit.1 <- glm.cluster(data = df_reg, formula = DeGroot ~ typesize, family = "binomial", cluster="typesizegroup")
  print(
    screenreg(list(lin.1, logit.1),
              ci.test = 0, ci.force.level = 0.95,
              custom.model.names = c("OLS (Bayesian & DeGroot predicts 0)", "Logit"),
              custom.header = list("Guess against majority in period 1,2" = 1:2))
  )
  
  # Print reg into LaTeX
  #sink(paste0("../output/latex/Match_mistake,",ntypelist,".tex"))
  print(
    texreg(list(lin.1, logit.1), 
           ci.test = 0, ci.force.level = 0.95,
           custom.model.names = c("OLS (Bayesian, DeGroot predicts 0)", "Logit"),
           custom.header = list("Guess against majority in period 1,2" = 1:2),
           booktabs = TRUE, dcolumn = FALSE, use.packages = FALSE,
           label = paste0("table:match_mistake,", ntypelist),
           caption = print(paste0("Fraction of guesses against Bayesian and DeGroot prediction")),
    )
  )
  #sink()
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
    
    g <- graph_from_adjacency_matrix(as.matrix(glist[[ntype]])-diag(n))
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
df_plot <- data.frame(id_subj = NA, type = NA, mean = NA, match = NA, avg.mean = NA)
for (ntypelist in c("n_10","n_40")){
  for (mth in c("df_DeGrootmore", "df_sigDeG_more")){
    df_temp <- filter(get(mth), type %in% get(ntypelist))
    df_temp <- aggregate(df_temp$DeGroot, list(df_temp$id_subj, df_temp$type),
                         drop=FALSE, FUN=function(x) mean = mean(x))
    df_temp <- cbind(df_temp, match = mth)
    colnames(df_temp) <- c("id_subj", "type", "mean", "match")
    
    df_temp2 <- df_temp %>% 
      group_by(type) %>% 
      mutate(avg.mean = mean(mean)) %>%
      filter(id_subj != 0)
    
    df_plot <- rbind(df_plot, df_temp2)
  }
}
df_plot <- df_plot[-1,]


for (ntypelist in c("n_10","n_40")){
  df_plot2 <- filter(df_plot, type %in% get(ntypelist))
  final <- ggplot(df_plot2, aes(x = mean, fill = match)) +
    geom_histogram(breaks = seq(0,1,0.05), color="white", closed = "left", 
                   position = "identity", alpha = 0.5)+
    # stat_bin(breaks = seq(0,1,0.05), 
    #          closed = "left", 
    #          geom="text", 
    #          colour="black", 
    #          size=2.5, 
    #          aes(label= ifelse(..count.. > 0, ..count.., ""), group=type, y=1+(..count..))) +
    #geom_vline(aes(xintercept=avg.mean, group=type, color="mean")) +
    # geom_vline(data=filter(df_match, 
    #                        match=="Guessing signal match with DeGroot", 
    #                        type %in% get(ntypelist)),
    #            aes(xintercept=percentage, 
  #                group=type, 
  #                color = "Guessing signal match with DeGroot"), 
  #            linetype="dashed")+
  
  facet_wrap(~ factor(type, levels = get(ntypelist)))+
    #ylim(NA, 60) +
    scale_colour_manual(values = c("#00bfc4","#f8766d"))+
    scale_fill_discrete(labels = c("Data match \nwith DeGroot", "Guessing signal match \nwith DeGroot")) +
    scale_x_continuous(breaks=seq(0,1,0.1), 
                       limits = c(-0.05, 1.05), 
                       labels = scales::percent_format(accuracy = 1)) +
    
    ggtitle(paste0("Distribution of DeGroot match across networks, ", 
                   ntypelist))+
    xlab("Percentage of guesses matching with DeGroot") + 
    ylab("Frequency of subjects")+
    theme_bw() +
    theme(legend.title = element_blank(), 
          strip.text.x = element_text(size = 13), 
          legend.justification=c(1,1), 
          legend.text = element_text(size=10), 
          legend.position=c(0.3,0.92),
          legend.key.size = unit(1, "cm"),
          legend.key.width = unit(0.5,"cm")
          )
  
  tname <- paste0("../output/DeGroot_match_distr,", ntypelist, ".pdf")
  #pdf(tname, width=10, height=7)
  print(final)
  #dev.off()
}

# Analyze % of guessing signal against DeGroot per network by person id_subj (Match_fol_sig.tex) ---------------------------
# Guessing signal when DeGroot predict otherwise

# Initialize file for learning rule matching with data, against DeGroot.

{
  df_plotn_10 <- 
    aggregate(df_gsig_n_10$DeGroot, 
              list(df_gsig_n_10$type, df_gsig_n_10$group, df_gsig_n_10$id_subj), 
              FUN=function(x) c(mean = mean(x, na.rm=TRUE)))
  colnames(df_plotn_10) <- c("type","group","id_subj","mean")
  #df_plot <- filter(df_plot, mean!=0)

  df_plotn_40 <- 
    aggregate(df_gsig_n_40$DeGroot, 
              list(df_gsig_n_40$type, df_gsig_n_40$group, df_gsig_n_40$round, df_gsig_n_40$id_subj), 
              FUN=function(x) c(mean = mean(x, na.rm=TRUE)))
  colnames(df_plotn_40) <- c("type","group","round","id_subj","mean")
  #df_plot2 <- filter(df_plot2, mean!=0)
  
  # Mean of each network of guessing signal against DeGroot
  aggregate(df_gsig_n_10$DeGroot,
            list(df_gsig_n_10$type), 
            FUN=function(x) c(mean = mean(x, na.rm=TRUE)))
  aggregate(df_gsig_n_40$DeGroot, 
            list(df_gsig_n_40$type), 
            FUN=function(x) c(mean = mean(x, na.rm=TRUE)))
  
  # # Histogram of each network of guessing signal against DeGroot
  # for (ntypelist in c("n_10","n_40")){
  #   final <- ggplot(get(paste0("df_plot",ntypelist)), aes(x=mean)) +
  #     geom_histogram(breaks = seq(0,1,0.1), color="white", closed = "left")+
  #     stat_bin(breaks = seq(0,1,0.1), closed = "left", geom="text", colour="black", size=2.5, aes(label= ifelse(..count.. > 0, ..count.., ""), group=type, y=1+(..count..))) +
  #     #geom_vline(aes(xintercept=avg.mean, group=type, color="mean")) +
  #     #geom_vline(data=filter(df_match, match=="Guessing signal match with DeGroot", type %in% n_10), aes(xintercept=percentage, group=type, color = "Guessing signal match with DeGroot"), linetype="dashed")+
  # 
  #     facet_wrap(~ factor(type, levels=get(ntypelist)))+
  #     #ylim(NA, 25) +
  #     scale_colour_manual(values = c("#00bfc4","#f8766d"))+
  #     scale_x_continuous(breaks=seq(0,1,0.1), limits = c(-0.05, 1.05), labels = scales::percent_format(accuracy = 1))+
  #     ggtitle(paste0("Distribution of Data matching with Guessing signal across networks, ",ntypelist))+
  #     xlab("Percentage of guesses matching with guessing signal, against DeGroot") + ylab("Frequency of subjects")+
  #     theme_bw() +
  #     theme(legend.title = element_blank(), strip.text.x = element_text(size = 13), legend.justification=c(1,1), legend.text = element_text(size=10), legend.position=c(0.3,0.92))
  # 
  #   tname <- paste0("../output/DeGroot_match_guesssig_distr,",ntypelist,".pdf")
  #   #pdf(tname, width=10, height=7)
  #   print(final)
  #   #dev.off()
  # }
  
  # Regression of Indicator on type 
  df_fol_sig <- data.frame(type = NA, fol_sig = NA, fol_sig_sd = NA)
  for (ntypelist in c("n_10", "n_40")){
    df_temp <- get(paste0("df_gsig_", ntypelist))
    df_reg <- df_temp %>%  #Recoding type, size, group, round
      mutate(typesize = df_temp$type,
             type = str_split_fixed(df_temp$type, "_",2)[,1],
             size = str_split_fixed(df_temp$type, "_",2)[,2],
             typesizegroup = paste0(as.character(type), 
                                    as.character(size), 
                                    as.character(group),
                                    sep = ''),
             typesizegroupround = paste0(as.character(type),
                                         as.character(size),
                                         as.character(group),
                                         as.character(round),sep = ''),
      ) %>%
      mutate_at(vars(group, round, type, size, 
                     typesizegroup, typesizegroupround), as.factor) %>%
      filter(!is.na(DeGroot))
    
    lin.1 <- lm.cluster(data = df_reg, formula = DeGroot ~ typesize, cluster = "typesizegroup")
    logit.1 <- glm.cluster(data = df_reg, formula = DeGroot ~ typesize, family = "binomial", cluster="typesizegroup")
    print(
      screenreg(list(lin.1, logit.1),
                ci.test = 0, ci.force.level = 0.95,
                custom.model.names = c("OLS (Stubbornness predicts 1)", "Logit"),
                custom.header = list("Always follow signal" = 1:2))
    )
    
    # Print reg into LaTeX (Match_fol_sig.tex)
    #sink(paste0("../output/latex/Match_fol_sig,",ntypelist,".tex"))
    print(
      texreg(list(lin.1, logit.1), 
             ci.test = 0, ci.force.level = 0.95,
             custom.model.names = c("OLS (Stubbornness predicts 1)", "Logit"),
             custom.header = list("Always follow signal" = 1:2),
             booktabs = TRUE, dcolumn = FALSE, use.packages = FALSE,
             label = paste0("table:fol_sig,", ntypelist),
             caption = print(paste0("Fraction of guesses following signal against DeGroot prediction")),
      )
    )
    #sink()
    df_temp <- 
      as.data.frame(unname(cbind(get(ntypelist), 
                                 summary(lin.1)[,1] + c(0, rep(summary(lin.1)[1,1],3)),
                                 summary(lin.1)[,2])))
    df_fol_sig <- rbind(df_fol_sig, setNames(df_temp, names(df_fol_sig)))
  }
  df_fol_sig <- df_fol_sig[-1,]
}

# Analyze % of following leader against DeGroot per network by person id_subj (Match_imi_lead.tex) ---------------------------
# Guessing signal when DeGroot predict otherwise

{
  df_temp <- filter(df_info_dom_more, type %in% n_10)
  df_plotn_10 <- 
    aggregate(df_temp$DeGroot, 
              list(df_temp$type, df_temp$group, df_temp$id_subj), 
              FUN=function(x) c(mean = mean(x, na.rm=TRUE)))
  colnames(df_plotn_10) <- c("type","group","id_subj","mean")
  
  df_temp <- filter(df_info_dom_more, type %in% n_40)
  df_plotn_40 <- 
    aggregate(df_temp$DeGroot, 
              list(df_temp$type, df_temp$group, df_temp$round, df_temp$id_subj), 
              FUN=function(x) c(mean = mean(x, na.rm=TRUE)))
  colnames(df_plotn_40) <- c("type","group","round","id_subj","mean")
  
  # Mean of each network of guessing signal against DeGroot
  df_temp <- filter(df_info_dom_more, type %in% n_10)
  aggregate(df_temp$DeGroot,
            list(df_temp$type), 
            FUN=function(x) c(mean = mean(x, na.rm=TRUE)))
  df_temp <- filter(df_info_dom_more, type %in% n_40)
  aggregate(df_temp$DeGroot, 
            list(df_temp$type), 
            FUN=function(x) c(mean = mean(x, na.rm=TRUE)))
  
  # # Histogram of each network of following leader against DeGroot
  # for (ntypelist in c("n_10","n_40")){
  #   final <- ggplot(get(paste0("df_plot",ntypelist)), aes(x=mean)) +
  #     geom_histogram(breaks = seq(0,1,0.1), color="white", closed = "left")+
  #     stat_bin(breaks = seq(0,1,0.1), closed = "left", geom="text", colour="black", size=2.5, aes(label= ifelse(..count.. > 0, ..count.., ""), group=type, y=1+(..count..))) +
  #     #geom_vline(aes(xintercept=avg.mean, group=type, color="mean")) +
  #     #geom_vline(data=filter(df_match, match=="Guessing signal match with DeGroot", type %in% n_10), aes(xintercept=percentage, group=type, color = "Guessing signal match with DeGroot"), linetype="dashed")+
  # 
  #     facet_wrap(~ factor(type, levels=get(ntypelist)))+
  #     #ylim(NA, 25) +
  #     scale_colour_manual(values = c("#00bfc4","#f8766d"))+
  #     scale_x_continuous(breaks=seq(0,1,0.1), limits = c(-0.05, 1.05), labels = scales::percent_format(accuracy = 1))+
  #     ggtitle(paste0("Distribution of Data matching with Guessing signal across networks, ",ntypelist))+
  #     xlab("Percentage of guesses matching with guessing signal, against DeGroot") + ylab("Frequency of subjects")+
  #     theme_bw() +
  #     theme(legend.title = element_blank(), strip.text.x = element_text(size = 13), legend.justification=c(1,1), legend.text = element_text(size=10), legend.position=c(0.3,0.92))
  # 
  #   tname <- paste0("../output/DeGroot_match_folead_distr,",ntypelist,".pdf")
  #   #pdf(tname, width=10, height=7)
  #   print(final)
  #   #dev.off()
  # }
  

  # Regression of Indicator on type
  df_pre_reg <- df_info_dom_more %>%  #Recoding type, size, group, round
    mutate(typesize = df_info_dom_more$type,
           type = str_split_fixed(df_info_dom_more$type, "_",2)[,1],
           size = str_split_fixed(df_info_dom_more$type, "_",2)[,2],
           typesizegroup = paste0(as.character(type), 
                                  as.character(size), 
                                  as.character(group),
                                  sep = ''),
           typesizegroupround = paste0(as.character(type),
                                       as.character(size),
                                       as.character(group),
                                       as.character(round),sep = ''),
           ) %>%
    mutate_at(vars(group, round, type, size, 
                   typesizegroup, typesizegroupround), as.factor) %>%
    filter(!is.na(DeGroot))
  
  df_imi_lead <- data.frame(type = NA, imi_lead = NA, imi_lead_sd = NA)
  for (ntypelist in c("n_10", "n_40")){
    df_reg <- filter(df_pre_reg, typesize %in% get(ntypelist))
    lin.1 <- lm.cluster(data = df_reg, formula = DeGroot ~ typesize, cluster = "typesizegroup")
    logit.1 <- glm.cluster(data = df_reg, formula = DeGroot ~ typesize, family = "binomial", cluster="typesizegroup")
    print(
      screenreg(list(lin.1, logit.1), 
                ci.test = 0, ci.force.level = 0.95,
                custom.model.names = c("OLS (Bayesian predicts 1)", "Logit"),
                custom.header = list("Correctly Imitate Leader" = 1:2))
    )
    # Print reg into LaTeX (Match_imi_lead.tex)
    #sink(paste0("../output/latex/Match_imi_lead,",ntypelist,".tex"))
    print(
    texreg(list(lin.1, logit.1), 
           ci.test = 0, ci.force.level = 0.95,
           custom.model.names = c("OLS (Bayesian predicts 1)", "Logit"),
           custom.header = list("Correctly follow leader" = 1:2),
           booktabs = TRUE, dcolumn = FALSE, use.packages = FALSE,
           label = paste0("table:match_imi_lead,", ntypelist),
           caption = print(paste0("Fraction of guesses imitate leader against DeGroot prediction")),
    )
    )
    #sink()
    
    df_temp <- 
      as.data.frame(unname(cbind(get(ntypelist), 
                                 summary(lin.1)[,1] + c(0, rep(summary(lin.1)[1,1],3)),
                                 summary(lin.1)[,2])))
    df_imi_lead <- rbind(df_imi_lead, setNames(df_temp, names(df_imi_lead)))
  }
  df_imi_lead <- df_imi_lead[-1,]
}

# Plot scatter plot between df_fol_sig, df_imi_lead ---------------------------
df_plot <- cbind(df_fol_sig, df_imi_lead[,-1])
df_plot$type <- as.factor(df_plot$type)
df_plot[,2:5] <- sapply(df_plot[,2:5], as.numeric)

final <- ggplot(df_plot,aes(x=fol_sig, y=imi_lead, shape=type, color=type )) + 
  geom_point(size = 3) +
  geom_errorbar(aes(ymin=pmax(imi_lead-imi_lead_sd,0), ymax=imi_lead+imi_lead_sd), width=.002,) +
  geom_errorbar(aes(xmin=fol_sig-fol_sig_sd, xmax=fol_sig+fol_sig_sd), width=.002,) +
  scale_shape_manual(values=c(0,15,1,16,2,17,4,18)) +
  
  scale_y_continuous(breaks=seq(0,1,0.1), 
                     limits = c(0, 0.4), 
                     labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(breaks=seq(0,1,0.1), 
                     limits = c(0,0.4), 
                     labels = scales::percent_format(accuracy = 1)) +
  
  ggtitle(paste0("Learning rule - Guesses going against DeGroot"))+
  ylab("Fraction correctly imitate leader") + 
  xlab("Fraction always follow signal")+
  theme_bw()

tname <- paste0("../output/Match_sig_lead.pdf")
#pdf(tname, width=10, height=7)
print(final)
#dev.off()


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

