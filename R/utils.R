#' BIC function
#'
#' Calculates the BIC for a given cluster solution
#' @param data Data frame, last column is the cluster solution
#' @return BIC of the given cluster solution
#' @export
BIC_cluster <- function(data){
  X <- data[,1:ncol(data)-1]
  Etit <- as.numeric(as.character(data[,ncol(data)]))
  Eti <- Etit
  index <- Etit
  n <- dim(X)[1]
  p <- dim(X)[2]
  meti <- max(Eti)
  auxl <- matrix(nrow=meti,ncol=n)
  param <- 0
  
  for (i in 1:meti){
    print(i)
    Rr <- chol(corpcor::make.positive.definite(var(X[index==i,])))
    Xz <- X-matrix(rep(colMeans(X[index==i,]),nrow(X)),ncol=p,byrow=T)
    E <- t(solve(t(Rr))%*%t(Xz))
    ni <- nrow(X[index==i,])
    auxl[i,] <- ((ni/n)*(1/sqrt(2*pi))/prod(diag(Rr)))*exp(-0.5*colSums(t(E*E)))
    param <- param+(p*(p+1)/2)+p
  }
  BIC <- 2*sum(log(colSums(auxl)))-param*log(n)
  return(BIC)
}

# foutliers
# 
# This function detect (and remove) outliers from 
# inside basic groups.
# @param resultados Data frame, last column contains the cluster solution
# @return A data frame without the outliers of each cluster
# @export
foutliers <- function(resultados){
  resultados <- ordering(resultados)
  original <- resultados
  p <- dim(original)[2]-1
  listaz <- split(resultados,resultados[,(p+1)])
  rnames <- vector()
  for (i in 1:length(listaz)){
    busco <- listaz[[i]][,-(p+1)]
    listaz[[i]] <- data.frame()
    listaz[[i]] <- busco[uni.plot(busco)$outliers==F,]
    listaz[[i]] <- cbind(listaz[[i]],rep(i,nrow(listaz[[i]])))
    names(listaz[[i]]) <- names(original)
    rnames2 <- rownames(listaz[[i]])
    rnames <- c(rnames,rnames2)    
  }
  names(listaz) <- NULL
  znew <- do.call("rbind", listaz)
  rownames(znew) <- rnames
  return(znew)
}

#'
#'funciongrafico
#'
#' Permite graficar el resultado de un algoritmo de cluster
#' (It needs to be improved)
#' @param result Data frame, last column contains the cluster solution.
#' @return A plot with the given cluster solution
#' @export
funciongrafico <- function(result, main=NULL) {
  p <- ncol(result)
  #x11()
  par(ask=F)
  plot(result[,1],result[,2], col="white", xlab="", ylab="", main=main)
  text(result[,1], result[,2], labels = as.character(result[,3]))
  #print(grafo)
  #x11()
  #xyplot(result[,2] ~ result[,1]|result[,p])
}

#'
#' uni.plot
#' 
#' This is a modification of the uni.plot function from the mvoutlier package
#' It calculates outliers without plotting
#' @param x Matrix or data frame.
#' @return A list with outliers and distances
#' @export
uni.plot <- function (x, symb = FALSE, quan = 1/2, alpha = 0.025, ...) 
{
  if (!is.matrix(x) && !is.data.frame(x)) 
    stop("x must be matrix or data.frame")
  if (ncol(x) < 2) 
    stop("x must be at least two-dimensional")
  #if (ncol(x) > 10) 
  #  stop("x should not be more than 10-dimensional")
  rob <- robustbase::covMcd(x, alpha = quan)
  # xarw <- mvoutlier::arw(x, rob$center, rob$cov, alpha = alpha) #Sometimes cov is singular
  xarw <- mvoutlier::arw(x, rob$center, corpcor::make.positive.definite(rob$cov), alpha = alpha)
  dist <- mahalanobis(x, center = rob$center, cov = corpcor::make.positive.definite(rob$cov))
  sx <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
  o <- (sqrt(dist) > min(sqrt(xarw$cn), sqrt(qchisq(0.975, 
                                                    dim(x)[2]))))
  l <- list(outliers = o, md = sqrt(dist))
  return(l)
}

#' ordering 
#'
#' Give ordered labels to the groups
#' The bigger group will be labelled as "1", and the rest are labelled by distance to "1" from closer to farest.
#' @param data Data frame, the last column contains the cluster solution.
#' @return A data frame with the cluster labels ordered.
#' @export
ordering  <-  function(data){
  n <- nrow(data)
  p <- ncol(data)
  x <- data[,p]
  tb <- table(x)
  grupo_ordenado <- factor(x,
                           levels = names(tb[order(tb, decreasing = TRUE)])) 
  levels(grupo_ordenado) <- 1:length(levels(grupo_ordenado))
  data2 <- data.frame(data[,-p],grupo_ordenado)
  data2$grupo_ordenado <- factor(grupo_ordenado,ordered=TRUE)
  data2 <- data2[order(data2$grupo_ordenado),]
  data <- data2
  distancias <- HDMD::pairwise.mahalanobis(data[,-p],data[,p])  
  grupo_ordenado <- data[,p]
  tb2 <- table(grupo_ordenado)
  grupo_ord2 <- factor(grupo_ordenado,
                       levels = names(tb2[sort.list(distancias$distance[,1])])) 
  levels(grupo_ord2) <- 1:length(levels(grupo_ord2))
  data3 <- data.frame(data[,-p],grupo_ord2)
  return(data3)
}
