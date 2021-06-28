# Create permutations of 1,0 for list of agents
binarypermutations <- function(list_of_var, vals = c(0,1)){ 
  tmp <- rep(list(vals), length((list_of_var)))
  tmp <- do.call(expand.grid, tmp)
  colnames(tmp) <- list_of_var
  return(as.matrix(tmp))
}

# Bayesian simulation given network, signals, number of periods
Bayesian <- function(g, s, periods) {
  record <- matrix(nrow=periods, ncol= ncol(g))
  belief <- matrix(data=0.5, nrow=ncol(g), ncol= ncol(g))
  info <- data.frame(period = NA, player = NA, information = NA)
  neigh <- neighborhood(graph_from_adjacency_matrix(as.matrix(g)-diag(ncol(g))),order=1, mode="out") #find out who you observe
  neigh2 <- neighborhood(graph_from_adjacency_matrix(as.matrix(g)-diag(ncol(g))),order=2, mode="out")
  
  for (i in 1:ncol(g)){
    belief[i,i] <- s[i] #first period information is just own signal
    info <- rbind(info, list(period=1, player=i, information=list(i)))
    #Guess majority based on info set. If indifference, guess signal.
    record[1,i] <- ifelse(mean(belief[i,])>0.5, 1, ifelse(mean(belief[i,])<0.5, 0 , s[i]))
  } 
  info <- info[-1,]
  
  
  # looking at their neighbour's information set, and deduce
  for (t in 2:periods){
    for (i in 1:ncol(g)){
      print(paste0("start",i," at period", t))
      I_itx <- filter(info, player==i, period == t-1)$information[[1]] #what is i's info set last period
      
      potential <- data.frame(neighbour = NA, information = NA, newinfo =NA)
      for (j in setdiff(neigh[[i]],i)){ #going over all people i observe, excluding i
        I_jtx <- filter(info, player==j, period == t-1)$information[[1]] #what is j's info set last period
        I_ijtx <- setdiff(I_jtx,I_itx) #how much info of j do i NOT know
        
        ###print(paste0("i=",i,",j=",j,",t=",t,",newinfo=",paste(as.character(I_ijtx),collapse=",")))
        
        #if j knows everything i know and more, maybe i should follow him
        #otherwise, stop caring about j if i know everything j can offer
        if (!is_empty(I_ijtx)) {
          potential <- rbind(potential, list(j, list(I_jtx), list(I_ijtx)))
          }
      }
      potential <- potential[-1,]
      
      if (nrow(potential)==0){ #if i has nothing to learn from anyone
        info <- rbind(info, list(period=t, player=i, information=list(I_itx)))
        belief[i,] <- belief[i,] #i's belief unchanged
      }
      else{
        x <- length(unique(unlist(potential$newinfo))) # total number of possible signals for i to learn
        if (x <=7){
          xx <- unique(unlist(potential$newinfo))
        }
        else {
          xx <- sample(unique(unlist(potential$newinfo)),7)
          x <- 7
        }
        # create a list of possible states given i's info
        tempbelief <- matrix(rep(belief[i,],2^x) ,nrow=2^x, ncol= ncol(g), byrow= TRUE)
        y <- binarypermutations(xx, vals = c(0,1)) #create permutations
        for (k in xx){
          tempbelief[,k] <- y[,paste0(k)] #replace permutations to get list of possible states given i's info
        }
        
        avgbelief <- matrix(nrow=0, ncol= ncol(g))
        for (m in 1:nrow(tempbelief)){ #for every plausible state, check all neighbour j's action consistent with tempbelief given j's knowledge
          a <- c()
          for (j in setdiff(neigh[[i]],i)){
            z <- setdiff(1:ncol(g), filter(info, player==j, period == t-1)$information[[1]]) #j has no perfect info on
            tempbeliefj <- tempbelief[m,]
            tempbeliefj[z]<-0.5 #replace i's tempbelief when j doesn't know information that i does
            # check if behaviour of j under this tempbelief matches guess of j last period
            a <- c(a, record[t-1,j]== ifelse(mean(tempbeliefj)>0.5, 1, ifelse(mean(tempbeliefj)<0.5, 0, record[t-2,j])))
          }
          if (all(a)){ #if all neighbours action match with tempbelief, tempbelief m is accepted as possible
            avgbelief <- rbind(avgbelief, tempbelief[m,])
          }
        }
        
        if (ncol(g)<=10){
          if (is_empty(avgbelief)) { #if no case matches perfectly, try again to find the most likely case
            print(paste0(c("No case match, ",i," will try guessing everyone, at period", t,", knowing info from", I_itx), collapse=" "))
            x <- length(setdiff(1:ncol(g),I_itx)) # total number of unknown signals from neighbours two degree away
            
            # create a list of possible states given i's info
            tempbelief <- matrix(rep(belief[i,],2^x) ,nrow=2^x, ncol= ncol(g), byrow= TRUE)
            y <- binarypermutations(setdiff(1:ncol(g),I_itx), vals = c(0,1)) #create permutations
            for (k in setdiff(1:ncol(g),I_itx)){
              tempbelief[,k] <- y[,paste0(k)] #replace permutations to get list of possible states given i's info
            }
            
            for (m in 1:nrow(tempbelief)){ #for every plausible state, check all neighbour j's action consistent with tempbelief given j's knowledge
              a <- c()
              for (j in potential$neighbour){ #i just look at the people that I can potential learn something from an induce who they could induce from
                tempbeliefj <- tempbelief[m,]
                a <- c(a, record[t-1,j]== ifelse(mean(tempbeliefj)>0.5, 1, 
                                                 ifelse(mean(tempbeliefj)<0.5, 0, record[t-2,j])))
              }
              if (is_empty(a)) {next}
              else if (all(a)){ #if all neighbours action match with tempbelief, tempbelief m is accepted as possible
                avgbelief <- rbind(avgbelief, tempbelief[m,])
              }
            }
          }
        }
        
        if (is_empty(avgbelief)) {
          avgbelief <- rbind(avgbelief, belief[i,])
          print(paste0("No case match even after trying guessing everyone, ",i," will use belief from last period."))
        }
        
        belief[i,] <- colMeans(avgbelief) #average over all possible states
        info <- rbind(info, list( period=t, player=i, information=list(which( colMeans(avgbelief)==1|colMeans(avgbelief)==0 )) )) #record what info i learnt
      }
      record[t,i] <- ifelse(mean(belief[i,])>0.5, 1, ifelse(mean(belief[i,])<0.5, 0 , record[t-1,i])) #i make guess based on belief
    }
    
    # check if t and t-1 info changed for all agents
    b <- 0
    for (i in 1:ncol(g)){
      # if info for i changed, then break this loop and carry on with next period
      if (!setequal(filter(info,period==t, player==i)$information[[1]],filter(info,period==t-1, player==i)$information[[1]])){
        b <- 0
        break
      }
      else { b <- b+1 }
    }
    # if t and t-1 info no change for all agents, then stop process
    if (b == ncol(g)){ break }
  }
  for (tt in (t+1):periods){
    record[tt,] <- record[tt-1,]
  }
  return(list(record, t))
}