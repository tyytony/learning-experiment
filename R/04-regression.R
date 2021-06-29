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
# library(aod)
# library(sandwich)
# library(multiwayvcov)
# library(lmtest)
# library(VGAM)
options(digits=4)
options(max.print=100)


load("RData/02-analysis-envir.RData")

# #import network statistics
# m <- read_csv("C:/Users/tyyto/Desktop/Social Learning/Simulation_DeGroot/Network_statistics.csv",
#               col_types = cols(X1 = col_skip(), connected = col_factor(levels = c("1","0")),
#               laplacian = col_factor(levels = c("1","0")) ) )


#Create dataframe for regression
df_pre_reg <- filter(df_summary, period>=7) #Selecting rounds
# df_pre_reg <- aggregate(df_pre_reg$mean, drop=FALSE, list(df_pre_reg$group, df_pre_reg$round, df_pre_reg$type), FUN=function(x) mean = mean(x))
# colnames(df_pre_reg) <- c("group", "round", "type","mean")

df_pre_reg <- df_pre_reg %>%  #Recoding type, size, group, round
  mutate(typesize = df_pre_reg$type,
         type = str_split_fixed(df_pre_reg$type, "_",2)[,1],
         size = str_split_fixed(df_pre_reg$type, "_",2)[,2],
         typesizegroup = paste0(as.character(type), 
                                as.character(size), 
                                as.character(group),
                                sep = ''),
         typesizegroupround = paste0(as.character(type),
                                     as.character(size),
                                     as.character(group),
                                     as.character(round),sep = ''),
         signal = df_pre_reg$round+(df_pre_reg$group-1)*6, #Recoding in signal used by group# and round#
         #Catagorize signals to be goods and non-goods
         good_signal = 
           ifelse(size == "10", 
                  ifelse(signal %in% goods$goodbad[goods$type=="ER_10"][[1]],"Goods", "Non-goods"),
                  ifelse(signal %in% goods$goodbad[goods$type=="ER_40"][[1]],"Goods", "Non-goods")),
  )
df_pre_reg <- merge(df_pre_reg,df_signal)
df_pre_reg <- mutate_at(df_pre_reg, vars(group, round, type, size, 
                                         typesizegroup, typesizegroupround, 
                                         good_signal), as.factor)

# Period 7-12 distribution of average guess across groups and rounds ---------------------------
for (ntypelist in c("n_10","n_40")){
  for (gs in list(c("Goods", "Non-goods"),"Goods", "Non-goods")){
    df_plot<-filter(df_pre_reg, 
                    typesize %in% get(ntypelist), 
                    period >= 7, 
                    good_signal %in% gs)
    
    final <- ggplot(df_plot, aes(x=mean)) +
      geom_histogram(breaks = seq(0,1,0.1), color="white", closed = "left")+
      stat_bin(breaks = seq(0,1,0.1), 
               closed = "left", 
               geom="text", 
               colour="black", 
               size=2.5, 
               aes(label= ifelse(..count.. > 0, ..count.., ""), 
                   group=typesize,
                   y=1+(..count..))) +
      
      geom_vline(aes(xintercept = 1, color = "Bayesian"), linetype="dashed")+
      geom_vline(aes(xintercept = 0.7, color = "No learning"), linetype="dashed")+

      facet_wrap(~ factor(typesize, levels=get(ntypelist)))+
      #ylim(NA, 15) +
      scale_colour_manual(values = c("blue","red"))+
      scale_x_continuous(breaks=seq(0,1,0.1), limits = c(-0.05, 1.05))+
      ggtitle(
        ifelse(length(gs)>1, 
               paste0("Period 7-12 guess distribution,", ntypelist), 
               paste0("Period 7-12 guess distribution, ",gs,",", ntypelist))
        )+
      
      xlab("Average guess of network") + 
      ylab("Frequency of networks (out of 24)")+
      theme_bw() +
      theme(legend.title = element_blank(), 
            strip.text.x = element_text(size = 13), 
            legend.justification=c(1,1), 
            legend.text = element_text(size=10), 
            legend.position=c(0.2,0.92))

    tname <- ifelse(length(gs)>1, 
                    paste0("../output/guess_distribution-", ntypelist, ".pdf"), 
                    paste0("../output/guess_distribution_",gs,"-", ntypelist, ".pdf"))
    #pdf(tname, width=10, height=7)
    print(final)
    #dev.off()
  }
}

# Regressions (Reg_output) ---------------------------
for (func in c("absolute","indicator_0.1","indicator_0.2","NtNs","NtNs_censored")){
  #sink(paste0("../output/Reg_output_", func ,".txt"))
  #sink()
  
  for (n in c("10","40")){
    for (x in list(c("Goods", "Non-goods"))){
      #for (x in list(c("Goods", "Non-goods"),"Goods", "Non-goods")){
      
      df_reg <- df_pre_reg %>%
        filter(good_signal %in% x, size == n) %>%
        mutate(type = relevel(type, ref = "ER"))
      
      if (func == "NtNs"){
        # Fraction of learning from initial signal (Nt,Ns)
        df_reg <- df_reg %>%
          mutate(cor_con   = (mean-Ns)/(1-Ns),
                 incor_con = (Ns-mean)/(Ns),
                 consensus = pmax(cor_con, incor_con),
                 breakdown = ifelse(consensus <= 0.3, 1, 0),
          )
        lin.1 <- lm.cluster(data=df_reg, formula=consensus ~ type, cluster = "typesizegroup")
        lin.2 <- lm.cluster(data=df_reg, formula=cor_con   ~ type, cluster = "typesizegroup")
        lin.3 <- lm.cluster(data=df_reg, formula=incor_con ~ type, cluster = "typesizegroup")
        logit.4 <- miceadds::glm.cluster(data=df_reg, formula=breakdown ~ type, family = "binomial", cluster="typesizegroup")
        
        #sink(paste0("../output/Reg_output_", func ,".txt"), append = TRUE)
        ifelse(length(x)>1,
               print(paste0("Regression of network size ", n)),
               print(paste0("Regression using partition of ", x, " signals, of network size ", n , collapse = '')))
        print(
          screenreg(list(lin.1, lin.2, lin.3, logit.4), ci.test=0, ci.force.level = 0.95,
                    custom.model.names = c("OLS - Consensus",
                                           "OLS - Correct Consensus",
                                           "OLS - Incorrect Consensus",
                                           "Logit - Breakdown"))
        )
        #sink()
      }
      
      if (func == "NtNs_censored"){
        # Fraction of learning from initial signal censored (Nt,Ns)
        df_reg <- df_reg %>%
          mutate(cor_con   = pmax((mean-Ns),0)/(1-Ns),
                 incor_con = pmax((0.5-mean),0)/(0.5),
                 consensus = pmax(cor_con, incor_con),
                 breakdown = ifelse(consensus <= 0.3, 1, 0),
          )
        tobit.1 <- zelig(data=df_reg, formula=consensus ~ type, model="tobit", robust = TRUE, cluster=typesizegroup, cite=F )
        tobit.2 <- zelig(data=df_reg, formula=cor_con ~ type, model="tobit", robust = TRUE, cluster=typesizegroup, cite=F ) #correct consensus
        tobit.3 <- zelig(data=df_reg, formula=incor_con ~ type, model="tobit", robust = TRUE, cluster=typesizegroup, cite=F )
        logit.4 <- miceadds::glm.cluster(data=df_reg, formula=breakdown ~ type, family = "binomial", cluster="typesizegroup")
        
        #sink(paste0("../output/Reg_output_", func ,".txt"), append = TRUE)
        ifelse(length(x)>1,
               print(paste0("Regression of network size ", n)),
               print(paste0("Regression using partition of ", x, " signals, of network size ", n , collapse = '')))
        print(
          screenreg(list(tobit.1, tobit.2, tobit.3, logit.4), ci.test=0, ci.force.level = 0.95,
                    custom.model.names = c("Tobit - Consensus",
                                           "Tobit - Correct Consensus",
                                           "Tobit - Incorrect Consensus",
                                           "Logit - Breakdown"))
        )
        #sink()
      }
      
      if (func == "absolute"){
        # Absolute function
        df_reg <- df_reg %>%
          mutate(consensus = abs(mean-0.5)*2,
                 cor_con   = ifelse(mean >= 0.5, abs(mean-0.5)*2, 0),
                 incor_con = ifelse(mean <= 0.5, abs(mean-0.5)*2, 0),
                 breakdown = 1- (abs(mean-0.5)*2),
          )
        lin.1 <- lm.cluster(data=df_reg, formula=consensus ~ type, cluster = "typesizegroup") #consensus
        tobit.2 <- zelig(data=df_reg, formula=cor_con ~ type, model="tobit", robust = TRUE, cluster=typesizegroup, cite=F ) #correct consensus
        tobit.3 <- zelig(data=df_reg, formula=incor_con ~ type, model="tobit", robust = TRUE, cluster=typesizegroup, cite=F )
        lin.4 <- lm.cluster(data=df_reg, formula=breakdown ~ type, cluster = "typesizegroup") #consensus
        
        #sink(paste0("../output/Reg_output_", func ,".txt"), append = TRUE)
        ifelse(length(x)>1,
               print(paste0("Regression of network size ", n)),
               print(paste0("Regression using partition of ", x, " signals, of network size ", n , collapse = '')))
        
        print(
          screenreg(list(lin.1, tobit.2, tobit.3, lin.4), ci.test=0, ci.force.level = 0.95,
                    custom.model.names = c("OLS - Consensus", 
                                           "Tobit - Correct Consensus", 
                                           "Tobit - Incorrect Consensus",
                                           "OLS - Breakdown"))
        )
        #sink()
        
        # #Print reg into LaTeX
        # texreg(list(lin.1, tobit.2, tobit.3), ci.test=0, ci.force.level = 0.95,
        #        custom.model.names = c("OLS - Consensus", "Tobit - Correct Consensus", "Tobit - Incorrect Consensus"),
        #        booktabs = TRUE, dcolumn=FALSE, label = "tab:1",
        #        caption= paste0("Regression of network size ",n,", partition of ",x," signals with clustering SE by network, size, group"),
        # )
      }
      
      if (func == "indicator_0.1"){
        # Indicator of width wd
        wd<-0.1
        df_reg <- df_reg %>%
          mutate(consensus = ifelse(mean >=1-wd | mean <= wd, 1,0),
                 cor_con   = ifelse(mean >=1-wd, 1,0),
                 incor_con = ifelse(mean < 0.3, 1,0),
                 breakdown = ifelse(mean < 0.7 & mean >0.3, 1, 0)
          )
        logit.1 <- miceadds::glm.cluster(data=df_reg, formula=consensus ~ type, family = "binomial", cluster="typesizegroup")
        logit.2 <- miceadds::glm.cluster(data=df_reg, formula=cor_con ~ type, family = "binomial", cluster="typesizegroup")
        logit.3 <- miceadds::glm.cluster(data=df_reg, formula=incor_con ~ type, family = "binomial", cluster="typesizegroup")
        logit.4 <- miceadds::glm.cluster(data=df_reg, formula=breakdown ~ type, family = "binomial", cluster="typesizegroup")
        
        #sink(paste0("../output/Reg_output_", func ,".txt"), append = TRUE)
        ifelse(length(x)>1,
               print(paste0("Regression of network size ", n)),
               print(paste0("Regression using partition of ", x, " signals, of network size ", n , collapse = '')))
        print(
          screenreg(list(logit.1, logit.2, logit.3, logit.4), ci.test=0, ci.force.level = 0.95,
                    custom.model.names = c("Logit - Consensus", 
                                           "Logit - Correct Consensus", 
                                           "Logit - Incorrect Consensus", 
                                           "Logit - Breakdown"))
        )
        #sink()
      }
      
      if (func == "indicator_0.2"){
        # Indicator of width wd
        wd<-0.2
        df_reg <- df_reg %>%
          mutate(consensus = ifelse(mean >=1-wd | mean <= wd, 1,0),
                 cor_con   = ifelse(mean >=1-wd, 1,0),
                 incor_con = ifelse(mean < 0.2, 1,0),
                 breakdown = ifelse(mean < 0.8 & mean >0.2, 1, 0)
          )
        logit.1 <- miceadds::glm.cluster(data=df_reg, formula=consensus ~ type, family = "binomial", cluster="typesizegroup")
        logit.2 <- miceadds::glm.cluster(data=df_reg, formula=cor_con ~ type, family = "binomial", cluster="typesizegroup")
        logit.3 <- miceadds::glm.cluster(data=df_reg, formula=incor_con ~ type, family = "binomial", cluster="typesizegroup")
        logit.4 <- miceadds::glm.cluster(data=df_reg, formula=breakdown ~ type, family = "binomial", cluster="typesizegroup")
        
        #sink(paste0("../output/Reg_output_", func ,".txt"), append = TRUE)
        ifelse(length(x)>1,
               print(paste0("Regression of network size ", n)),
               print(paste0("Regression using partition of ", x, " signals, of network size ", n , collapse = '')))
        print(
          screenreg(list(logit.1, logit.2, logit.3, logit.4), ci.test=0, ci.force.level = 0.95,
                    custom.model.names = c("Logit - Consensus", 
                                           "Logit - Correct Consensus", 
                                           "Logit - Incorrect Consensus", 
                                           "Logit - Breakdown"))
        )
        #sink()
      }
    }
  }
}


# Other functions to indicate the 4 dependent variables for robustness ---------------------------
#
# # Squared function
# df_reg <- df_reg %>% #Squared function
#   mutate(consensus = (mean-0.5)^2*4,
#          cor_con   = ifelse(mean >= 0.5, (mean-0.5)^2*4, 0),
#          incor_con = ifelse(mean <= 0.5, (mean-0.5)^2*4, 0),
#   )
#
# # Sigmoid function
# df_reg <- df_reg %>%
#   mutate(consensus = ifelse(mean >= 0.5, (1/(1+exp(-9*(mean-0.5)))-0.5)*2, (1/(1+exp(-9*(0.5-mean)))-0.5)*2),
#          cor_con   = 1/(1+exp(-9*(mean-0.5))),
#          incor_con = 1/(1+exp(-9*(0.5-mean))),
#   )
# lin.1 <- lm.cluster(data=df_reg, formula=consensus ~ type, cluster = "typesizegroup") #consensus
# lin.2 <- lm.cluster(data=df_reg, formula=cor_con ~ type, cluster = "typesizegroup") #consensus
# lin.3 <- lm.cluster(data=df_reg, formula=incor_con ~ type, cluster = "typesizegroup") #consensus
# screenreg(list(lin.1,lin.2,lin.3), ci.test=0, ci.force.level = 0.95,
#           custom.model.names = c("OLS - Consensus", "OLS - Correct Consensus", "OLS - Incorrect Consensus"))


# summary(reg)
# confint(reg) #confidence interval
# wald.test(b = coef(reg), Sigma = vcov(reg), Terms = 3:4) #wald test
# wald.test(b = coef(reg), Sigma = vcov(reg), cbind(0, 1, 0, -1)) #hypothesis testing of beta1=beta2
# exp(cbind(OR = coef(reg), confint(reg))) #odds ratios
# 
# step <- stepAIC(reg_mean, direction="both") # step-wise regression
# step$anova # display results