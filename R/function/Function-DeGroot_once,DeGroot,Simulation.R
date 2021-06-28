# DeGroot function (graph, signal, discrete or continuous DeGroot) ---------------------------

DeGroot_once <- function(g, s, w=0, tremble=0) {
  guess <- s #set initial guess as first signal
  # Recode network where weighted link of i to j correspond to outdegree of j
  g <- w*t(t(g)*rowSums(g))+(1-w)*g
  avg<- g%*%guess/rowSums(g) #calculate observed average
  
  # DeGroot process
  for (i in 1:length(guess)) {
    if       (is.nan(avg[i]))          {guess[i]=guess[i]}
    else if  (avg[i]>0.5)              {guess[i]=rbinom(1,1,1-as.numeric(tremble))}
    else if  (avg[i]<0.5)              {guess[i]=rbinom(1,1,0+as.numeric(tremble))}
    else if  (avg[i]==0.5)
    {guess[i]=0.5}
    # {
    # if (s[i]==0.5)                   {guess[i]=rbinom(1,1,0.5)}
    # else                             {guess[i]=s[i]}
    # }
  }
  return(guess)
}

# DeGroot simulation function, and repeated simulations ---------------------------
DeGroot <- function(g, s, periods, 
                    w=0, tremble=0, k=1, c=0, 
                    record=matrix(nrow=0, ncol= ncol(g))) {
  guess <- s #set initial guess as first signal
  record<- rbind(record, guess) #record copy of guess at previous periods
  t <- nrow(record)
  
  # Recode network where weighted link of i to j correspond to outdegree of j
  g <- w*t(t(g)*rowSums(g))+(1-w)*g 
  # Recode network where weighted link of i to j correspond to local clustering coefficient of j
  if (c !=0){
    g <- c*t(t(g)*transitivity(graph_from_adjacency_matrix(as.matrix(g)-diag(ncol(g))), type="local"))+(1-c)*g
  }
  
  avg <- g%*%guess/rowSums(g) #calculate observed average
  
  # DeGroot process
  if (k==1){
    for (i in 1:length(guess)) {
      if       (is.nan(avg[i]))          {guess[i]=guess[i]}
      else if  (avg[i]>0.5)              {guess[i]=rbinom(1,1,1-as.numeric(tremble))}
      else if  (avg[i]<0.5)              {guess[i]=rbinom(1,1,0+as.numeric(tremble))}
      else if  (avg[i]==0.5)
      {guess[i]=rbinom(1,1,0.5)}
      # {
      # if (s[i]==0.5)                   {guess[i]=rbinom(1,1,0.5)}
      # else                             {guess[i]=s[i]}
      # }
    }
  }
  
  else if (k>=0 & k<1){
    x <- k/(1-k)
    p <- 1/(1+exp(-x*(avg-0.5)))
    for (i in 1:length(guess)) {
      if       (is.nan(avg[i]))          {guess[i]=guess[i]}
      else     {guess[i] <- rbinom(1,1,p[i])}
    }
  }
  
  # Stop process after x periods
  if (t >= periods-1)
  {record<- rbind(record, guess)
  return(record)}
  
  # repeat function using guess as new signal
  else if (t< periods-1)
  {return(DeGroot(g,guess,periods,w,tremble,k,c,record))}
}

# Simulation function (#repetition, graph, probability) ---------------------------
Simulation <- function(rep, periods, g, p, w=0, tremble=0, k=1, c=0){
  
  output <- data.frame(con_time=NA, mean=NA)
  
  # Generate 1000 sets of signals of size n
  balls <- matrix(data= rbinom(rep*ncol(g),1, p), ncol = ncol(g))
  
  for (r in c(1:rep))
  {
    if (mean(balls[r,])<0.5) { balls[r,] <- 1-balls[r,] }
    sig <- balls[r,]
    history <- DeGroot(g, sig, periods, w, tre, k, c)
    
    # time till convergence
    t<-NA
    for (period in 1:(nrow(history)-1)){
      if (all( abs((history[period,]-history[period+1,])) < 0.0001) ){
        t<-period
        break
      }
    }
    # mean at t
    meanconvg <-mean(history[nrow(history),])
    output <- rbind(output, list(t, meanconvg))
  }
  output <- output[-1,]
  return(output)
}