# Functions for enterotyping

dist.JSD <- function(inMatrix, pseudocount = 0.000001, ...) {
  KLD <- function(x, y) sum(x * log(x/y))
  JSD <- function(x, y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix, 1:2, function(x) ifelse (x == 0, pseudocount, x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i, j] <- JSD(as.vector(inMatrix[ ,i]),
                                 as.vector(inMatrix[ ,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix) -> resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}


pam.clustering <- function(x, k) {
  clust <- as.vector(pam(as.dist(x), k, diss = TRUE)$clustering)
  return(clust)
}


noise.removal <- function(dataframe, percent=0.01, top=NULL){
  dataframe -> Matrix
  bigones <- rowSums(Matrix) * 100 / (sum(rowSums(Matrix))) > percent 
  Matrix_1 <- Matrix[bigones, ]
  print(percent)
  return(Matrix_1)
}