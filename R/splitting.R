#' Function SAR
#'
#' Splitting function from the SAR algorithm proposed by Peña and Tiao (2003)
#' This function calls to function SAR but also remove discriminators and group with sizes smaller than minsize parameter.
#' @param datos Data frame, currently sagra supports only numeric variables.
#' @param minsize The minimum size of a cluster
#' @param r order of the considered discriminators
#' @return A cluster solution
#' @export
functionSAR <- function(datos,minsize=ncol(datos)+log(nrow(datos)-ncol(datos)), r=ncol(datos)+log(nrow(datos)-ncol(datos))-1){
  n <- nrow(datos) #number of observations
  p <- ncol(datos) #number of variables
  original <- datos # we save original data for later
  datos <- as.data.frame(datos) #in previous steps it could be matrix
  Ed <- Effi_d(datos) # We calculate the distance to the convex hull
  d <- Ed$d2 # And we extract the distances
  discriminadores <- Ed$discriminadores # Vector with the convex hull members
  x <- max.col(d, ties.method="last") # what member of convex hull maximizes distance
  discriminadores <- discriminadores[sort(unique(x))] # Discriminators are subset of convex hull
#   if(p>=5){ #if dimension is bigger than 5, we need to use r-discriminators. (See Peña and Tiao, 2003)
#     discriminadores <- as.numeric(names(table(x)[table(x)>r]))  
#   }
  x <- factor(x) #We treat them as factor
  levels(x) <- 1:length(levels(x)) #And we re order the levels
  data <- data.frame(datos,x)
  #Paso1: Ordenar de mayor a menor. 
  tb <- table(x)
  grupo_ordenado <- factor(x,
                           levels = names(tb[order(tb, decreasing = TRUE)])) 
  levels(grupo_ordenado) <- 1:length(levels(grupo_ordenado))
  data2 <- data.frame(datos,grupo_ordenado)
  data2$grupo_ordenado <- factor(grupo_ordenado,ordered=TRUE)
  
  # We clean outliers
  data2 <- data2[order(data2$grupo_ordenado),]
  N <- table(data2$grupo_ordenado)
  N <- as.vector(N)
  N <- rep(N,N)
  data2 <- cbind(data2,N)
  n0 <- minsize
  dataclean <- data2[data2$N>=n0,] #Only those groups bigger than minsize
  outliers <- data2[data2$N<n0,] #And those smaller are grouped as "outliers". They will be recombined later
  dataclean <- dataclean[,-ncol(dataclean)] 
  results <- list(dataclean,outliers,discriminadores)
  names(results) <- c("data","outliers","discriminadores")
  if (dim(results$data)[1]==0){
    results$data <- cbind(original,rep(1,nrow(original)))
    names(results$data) <- c(names(original),"grupo_ordenado")
    results$outliers <- as.data.frame(matrix(ncol=p+1,nrow=0))
    names(results$outliers) <- c(names(original),"grupo_ordenado")
  }
  return(results)
}

# Function SAR2
#
# Splitting function from the SAR algorithm proposed by Peña and Tiao
# This function calls to function SAR but also remove discriminators and group with sizes smaller than minsize parameter.
functionSAR2 <- function(data,minsize=ncol(data)+log(nrow(data)-ncol(data)),r=ncol(data)+log(nrow(data)-ncol(data))){
  
  n <- nrow(data)
  p <- ncol(data)
  if(nrow(data)==minsize){
    data2 <- cbind(data,rep(1,minsize))
    names(data2) <- c(names(data),"grupo_ordenado")
    return(data2)
  }
  resultsSAR <- functionSAR(data,minsize,r)
  newdata <- resultsSAR$data
  discriminadores <- resultsSAR$discriminadores
  newdata <- newdata[rownames(newdata)%in%setdiff(rownames(newdata),discriminadores),]
  newdata[,p+1] <- factor(newdata[,p+1])
  newdata <- elimino(newdata,minsize)
  return(newdata)
}

# Effi_d 
#
# "Efficient" calculation of distances between elements and discriminators
# Is necessary to standardize the data first.
# 
Effi_d=function(data){
  X <- as.matrix(data) #We will need some algebra here
  n <- nrow(data) #number of observations
  p <- ncol(data) #number of variables
  UNO <- rep(1,n) # vector of ones
  P <- diag(n)-(UNO%*%t(UNO))/n 
  Xc <- P%*%X
  # S <- cov(X) # Some times in small groups we can have constant variables, or LD.
  S <- corpcor::make.positive.definite(cov(X))
  E <- eigen(S)
  Sraiz <- E$vectors%*%diag(sqrt(E$values))
  Sraiz2 <- Sraiz%*%t(E$vectors)
  Y <- Xc%*%solve(Sraiz2)
  #})
  Y <- Y/sqrt(n-p) #
  if (p<5){ #Convex hull works for small dimensions
    k <- geometry::convhulln(data, options="Pp") #
    discriminadores <- unique(as.vector(k)) #
    discriminadores <- sort(discriminadores) #
  } else {
    discriminadores <- 1:nrow(data) #So, for big dimensions we just consider all data as possible discriminators
    #Aquí va a estar el cambio. Como soy gilipollas y no tengo tiempo de calcular la distancia de mahalanobis fraccional,
    #lo que voy a hacer es definir el nuevo "convex hull" basado en las distancias fraccionales por pares
    #profit ! Discrimina mucho mejor que considerarlos a todos
    # discriminadores <- convex_hd(Y)
  }
  Dn <- Y%*%t(Y[discriminadores,])
  Dn <- Dn+(1/n)
  Dn <- Dn*Dn
  Dd <- matrix(rep(diag(Y[discriminadores,]%*%t(Y[discriminadores,])),each=n),nrow=n)
  Dd <- (n/(n-1))-Dd
  d2 <- Dn/Dd
  diag(d2[discriminadores,]) <- rep(0,length(discriminadores))
  return(list(d2=d2,discriminadores=discriminadores)) #
}


<<<<<<< HEAD
# convex_hd <- function(Y, p2=0.1){
#   #Given a HD matrix calculate the semi convex hull
#   #The semi convex hull is the set of all elements
#   #Whose fractional distance is bigger to some point in the data
#   #It depends on the pairwise.diffs function below
#   a <- pairwise.diffs(t(Y))
#   a2 <- abs(a)
#   a3 <- (colSums(a2^p2))^(1/p2)
#   m2 <- matrix(0, nrow=nrow(Y), ncol=nrow(Y))
#   m2[lower.tri(m2)] <- a3
#   m3 <- m2+t(m2)
#   # euc.dist <- function(x1, x2, p) (sum(abs(x1 - x2) ^ p))^(1/p)
#   
#   x <- max.col(m3, ties.method="last")
#   x <- sort(unique(x))
#   return(x)
# }
# 
# 
# pairwise.diffs <- function(x)
# {
#   #This function is from Marc Schwartz in the r-help list
#   #https://stat.ethz.ch/pipermail/r-help/2004-August/055324.html
#   stopifnot(is.matrix(x))
#   
#   # create column combination pairs
#   prs <- cbind(rep(1:ncol(x), each = ncol(x)), 1:ncol(x))
#   col.diffs <- prs[prs[, 1] < prs[, 2], , drop = FALSE]
#   
#   # do pairwise differences 
#   result <- x[, col.diffs[, 1]] - x[, col.diffs[, 2], drop = FALSE]
#   
#   # set colnames
#   if(is.null(colnames(x)))
#     colnames(x) <- 1:ncol(x)
#   
#   colnames(result) <- paste(colnames(x)[col.diffs[, 1]], ".vs.",
#                             colnames(x)[col.diffs[, 2]], sep = "")
#   result
# }
=======
convex_hd <- function(Y, p2=0.1){
  #Given a HD matrix calculate the semi convex hull
  #The semi convex hull is the set of all elements
  #Whose fractional distance is bigger to some point in the data
  #It depends on the pairwise.diffs function below
  a <- pairwise.diffs(t(Y))
  a2 <- abs(a)
  a3 <- (colSums(a2^p2))^(1/p2)
  m2 <- matrix(0, nrow=nrow(Y), ncol=nrow(Y))
  m2[lower.tri(m2)] <- a3
  m3 <- m2+t(m2)
  # euc.dist <- function(x1, x2, p) (sum(abs(x1 - x2) ^ p))^(1/p)
  
  x <- max.col(m3, ties.method="last")
  x <- sort(unique(x))
  return(x)
}


pairwise.diffs <- function(x)
{
  #This function is from Marc Schwartz in the r-help list
  #https://stat.ethz.ch/pipermail/r-help/2004-August/055324.html
  stopifnot(is.matrix(x))
  
  # create column combination pairs
  prs <- cbind(rep(1:ncol(x), each = ncol(x)), 1:ncol(x))
  col.diffs <- prs[prs[, 1] < prs[, 2], , drop = FALSE]
  
  # do pairwise differences 
  result <- x[, col.diffs[, 1]] - x[, col.diffs[, 2], drop = FALSE]
  
  # set colnames
  if(is.null(colnames(x)))
    colnames(x) <- 1:ncol(x)
  
  colnames(result) <- paste(colnames(x)[col.diffs[, 1]], ".vs.",
                            colnames(x)[col.diffs[, 2]], sep = "")
  result
}
>>>>>>> 6c6a35454257aa5e31eb0e86b6c1116ba966dcf6


# elimino
#
# Function to remove groups with sizes < minsize
elimino <- function(data,minsize){
  guardo <- data
  p <- dim(data)[2]-1
  a <- as.matrix(table(data[,p+1]))
  b <- which(a<minsize)
  if(length(b)!=0){
    data <- data[data[,p+1]%in%setdiff(data[,p+1],b),]
    data[,p+1] <- factor(data[,p+1])
    if(nrow(data)==0){
      data <- guardo
      data[,p+1] <- factor(rep(1,nrow(data)))
    }
  }
  return(data)
}

