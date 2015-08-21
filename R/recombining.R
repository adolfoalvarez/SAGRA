# recombining
#
# This function uses a Bayes Factor to recombine groups
# 
joining  <-  function(data) {
  #data <- results$data
  names(data)[length(names(data))] <- "grupo"
  n <- dim(data)[1]
  p <- dim(data)[2]
  grupos <- data.frame(matrix(nrow=0,ncol=p))
  gruposcounter <- 0
  z <- NULL
  i <- 0
  repeat {
    if (max(data[,p])==1){
      grupos <- data
      break}
    i <- i+1
    # print(i)
    data <- ordering(data)
    #x11()
    data[,p] <- as.numeric(data[,p])
    resulttest <- testing_function(data)
    if (resulttest >= 0.01){
      #data=siuno(data)
      data[(data[,ncol(data)]==1)|(data[,ncol(data)]==2),p] <- rep(1,length(data[(data[,ncol(data)]==1)|(data[,ncol(data)]==2),p]))
      data[(data[,ncol(data)]!=1),p] <- data[(data[,ncol(data)]!=1),p]-1
      data[(data[,ncol(data)]==1),]
      z <- c(z,as.numeric(rownames(data[(data[,ncol(data)]==1),])))
    } else {
      #data=nouno(data)
      elquesaco <- data[(data[,ncol(data)]==1),]
      elquesaco[,ncol(elquesaco)] <- elquesaco[,ncol(elquesaco)]+gruposcounter
      gruposcounter <- gruposcounter+1
      grupos <- rbind(grupos,elquesaco)
      data <- data[(data[,ncol(data)]!=1),]
      data[(data[,ncol(data)]!=1),p] <- data[(data[,ncol(data)]!=1),p]-1
    }
    if (dim(table(data[,p])) == 1) {
      elquesaco <- data[(data[,ncol(data)]==1),]
      elquesaco[,ncol(elquesaco)] <- elquesaco[,ncol(elquesaco)]+gruposcounter
      gruposcounter <- gruposcounter+1
      grupos <- rbind(grupos,elquesaco)
      break
    }
  }
  #return(grupos)
  return(list(grupos,z))
}

# recombining2
#
# This function calls to the previous recombining function and applies until the result is the same
# 
joining2 <- function(data){
  data2 <- joining(data)
  z <- unique(data2[[2]])
  dataold <- data2
  repeat{
    datanew <- joining(dataold[[1]])
    znew <- c(z,unique(datanew[[2]]))
    z <- unique(znew)
    a <- BIC_cluster(datanew[[1]])
    b <- BIC_cluster(dataold[[1]])
    if (identical(a,b) ){
      break
    }
    dataold<-datanew
  }
  fin <- dataold[[1]]
  return(list(fin,z))
}

# lastcombining
#
# This is the Last combining step in the SAGRA algorithm
# It assigns isolated observations to a given cluster solution based on Mahalanobis' distances
# 
lastcombining <- function(data,resultados){
  data <- as.data.frame(data)
  n <- dim(data)[1]
  p <- dim(data)[2]
  indices <- as.numeric(row.names(resultados))
  originalindex <- as.numeric(row.names(data))
  noindices <- setdiff(originalindex,indices)
  newdata <- as.data.frame(data[row.names(data)%in%noindices,])
  listamedias <- by(resultados[,1:p],resultados[,p+1],colMeans)
  listavarianzas <- by(resultados[,1:p],resultados[,p+1],cov)
  g <- max(as.numeric(resultados[,dim(resultados)[2]]))
  if (g==1){
    resultados <- cbind(data,1)
    names(resultados)[length(names(resultados))] <- "grupo_ord2"
    return(resultados)
  }
  rrr <- matrix(ncol=g,nrow=nrow(newdata))
  for (j in 1:g){    
    rrr[,j] <- mahalanobis(newdata,
                           listamedias[[j]],
                           chol2inv(chol(
                             corpcor::make.positive.definite(
                             listavarianzas[[j]]))),inverted=TRUE)
  }
  rrr <- data.frame(rrr)
  #k <- 0
  themin <- which(rrr==min(rrr), arr.ind=TRUE)
  minrow <- themin[1]
  mincol <- themin[2]
  obs <- newdata[minrow,]
  names(resultados)[length(names(resultados))] <- "grupo_ord2"
  obs$grupo_ord2 <- mincol
  resultados <- rbind(resultados,obs)
  newdata <- newdata[-minrow,]
  rrr <- rrr[-minrow,]
  i <- nrow(resultados)
  while (nrow(newdata)>0){
    # while (i>0){
    # print(max(as.numeric(rownames(resultados))))
    # print(nrow(newdata))
    listamedias[[mincol]] <- colMeans(resultados[resultados[p+1]==mincol,-(p+1)])
    listavarianzas[[mincol]] <- cov(resultados[resultados[p+1]==mincol,-(p+1)])
    rrr[,mincol] <- mahalanobis(newdata,
                                listamedias[[mincol]],
                                chol2inv(chol(
                                  corpcor::make.positive.definite(
                                  listavarianzas[[mincol]]))),inverted=TRUE)
    themin <- which(rrr==min(rrr), arr.ind=TRUE)
    minrow <- themin[1,1]
    mincol <- themin[1,2]
    obs <- newdata[minrow,]
    obs$grupo_ord2 <- mincol
    # resultados <- rbind(resultados,obs) #rbind causes problems with rownames.
    #I know this is not efficient. We need to preallocate, or using dplyr (without losing rownames!)
    resultados[i+1,] <- obs
    newdata <- newdata[-minrow,]
    rrr <- rrr[-minrow,]    
    i <- i+1
  }
  return(resultados)
}

# laststep
#
# This is the last step in the SAGRA algorithm
# After isolated observations are assigned to clusters, is necessary one last step to recombine those groups
# 
laststep <- function(data){
  n <- nrow(data)
  p <- ncol(data)-1
  lastlist <- list()
  data <- ordering(data)
  data[,p+1] <- as.numeric(data[,p+1])
  repeat{
    test <- testfunction(data)  
#     print(test)
    if(test>0.01){
      lastlist[[length(lastlist)+1]] <- data
      break
    } else{
      lastlist[[length(lastlist)+1]] <- data[data[,p+1]==1,]
      data <- data[data[,p+1]!=1,]
      data[,p+1] <- data[,p+1]-1
    }    
  }
  for (i in 1:length(lastlist)){
    lastlist[[i]][,p+1] <- i
  }
  result <- do.call("rbind", lastlist)
  return(result)
}

laststep2 <- function(data){
  n <- nrow(data)
  p <- ncol(data)-1
  results <- laststep(data)
  if(max(data[,p+1])==max(results[,p+1])) return(results)
  repeat{    
    results <- ordering(results)
    results2 <- laststep(results)
    if(max(results2[,p+1])==max(results[,p+1])) break
    results <- results2
  }
  return(results)
}

testing_function  <-  function(data){
  n <- nrow(data)
  p <- ncol(data)
  pairs <- c(1,2)
  testdata <- data[(data[,p]==pairs[1])|(data[,p]==pairs[2]),]
  testresult <- testfunction(testdata)
  rm(n,p,pairs)
  return(testresult)
}

#' testfunction
#'
#'
#' Performs a Bayes Factor test for a given cluster solution
#'
#' @param data Data frame, where the last column is the cluster solution
#' @return Probability of $H_0$
#' @export
testfunction <-  function(data){
  n <- dim(data)[1]
  w1 <- rep(1,n)
  w <- data[,dim(data)[2]]
  w <- as.vector(w)
  w <- as.factor(w)
  data <- data[,-dim(data)[2]]
  r <- pifunction(w)+mfunction(data,w)-mfunction(data,w1)
  return(1/(1+exp(r)))
}

#Funcion pi(w)
pifunction <- function(w){
  #pitatoria de gamma(n_i) sobre gamma (n)
  (sum(lgamma(table(w))))-lgamma(sum(table(w)))
}

mfunction <- function(data,w){
  n <- nrow(data)
  p <- ncol(data)
  medias <- by(data,w,colMeans)
  varianzas <- by(data,w,cov)
  varianzas <- lapply(varianzas, corpcor::make.positive.definite)
  mylist <- split(data, f=w)
  z <- mapply(mahalanobis, x=mylist, center=medias, cov=varianzas)
  if(!is.list(z)) z <- list(z)
  myfun <- function(x){
    ni <- length(x)
    return(log(1+((ni*x)/((ni+1)*(ni-1)))))
  }
  phi3 <- lapply(z, myfun)
  mydet <- sapply(varianzas,det)
  myn <- sapply(mylist,nrow)
  ni <- myn
  phi2 <- lgamma(ni/2)-lgamma((ni-p)/2)-(p/2)*log(ni-1)-(1/2)*log(mydet)
  phi1 <- (p/2)*log(ni/(pi*(ni+1)))
  total <- sum(phi1*ni)+sum(phi2*ni)+sum(sapply(phi3,sum)*(-ni/2))
  return(total)
}
