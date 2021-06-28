Influence_v <- function(g){
  if (is.igraph(g)){
    g<-as.matrix(as_adjacency_matrix(g))
  }
  colnames(g) <- 1:ncol(g)
  g<-g/rowSums(g)
  
  #left unit eigenvector = influence vector
  eigenvalue <- abs(Re(eigen(t(g))$vectors[,1]))
  eigenvalue <- eigenvalue/sum(eigenvalue)
  return(eigenvalue)
}

Influence_ratio <- function(g){
 infl_list<-sort(Influence_v(g))
 return(infl_list[length(infl_list)]/infl_list[1])
}
