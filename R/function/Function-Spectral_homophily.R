Homophily <- function(g){
  if (is.igraph(g)){
    g<-as.matrix(as_adjacency_matrix(g))
  }
  colnames(g) <- 1:ncol(g)
  if (g[1,1]==1){g<-g-diag(x=1,nrow=ncol(g),ncol= ncol(g))}
  m<-ncol(g)/5
  mat <- matrix(data=NA ,nrow=m,ncol=m)
  for (i in 1:m){
    for (j in 1:m){
      gg<-g[(5*(i-1)+1):(5*(i-1)+5),(5*(j-1)+1):(5*(j-1)+5)]
      mat[i,j]<-sum(rowSums(gg))
    }
  }
  
  mat<-mat/rowSums(mat)
  homophily<-sort(eigen(mat)$values, decreasing = TRUE)[2]
  return(homophily)
  #left unit eigenvector = influence vector
  # eigenvalue <- abs(Re(eigen(t(g))$vectors[,1]))
  # eigenvalue <- eigenvalue/sum(eigenvalue)
  # return(eigenvalue)
}
