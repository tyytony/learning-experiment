is_r_cohesive <- function(g, vlist) {
  for (i in c(1:length(vlist))){
    if (degree(induced_subgraph(g, vlist),i)*2 < degree(g,vlist[i]))
      return(FALSE)
  }
  return(TRUE)
}

r_cohesiveness <- function(g,size) {
  n<-vcount(g)
  ds<-0
  dg<-0
  df <- matrix(NA)
  #dg <- nCr(n, 2)+nCr(n, 3)+nCr(n, 4)+nCr(n, 5)+nCr(n, 6)
  for (k in c(2:size)){
    for (vlist in combn(c(1:n),k, simplify = FALSE)){
      if (is.connected(induced_subgraph(g,vlist))){
        #check whether cohesive
        if (is_r_cohesive(g,vlist)) {
          df <- rbind(df,list(vlist))
          ds <- ds+1
          dg <- dg+1
        }
        else {
          dg <- dg+1
        }
      }
    }
  }
  df <- as.data.frame(df)[-1,]
  return(list(ds/dg,df))
}