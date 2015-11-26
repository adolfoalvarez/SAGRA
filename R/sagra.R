#' Splitting And Group Recombining Algorithm (SAGRA)
#'
#'
#' This is a cluster algorithm proposed by Alvarez and Peña (2014)
#' Based on Bayes Factors
#'
#' @param data Data frame, currently sagra supports only numeric variables.
#' @param minsize The minimum size of a cluster
#' @return A cluster solution
#' @export
sagra  <- function(data,minsize=ncol(data)+log(nrow(data)-ncol(data))){
  original <- data
  n <- dim(data)[1]
  p <- dim(data)[2]
  r <- ncol(data)+log(nrow(data)-ncol(data))-1
  grupos_candidatos <- list()
  pasos <- 0
  datos <- functionSAR2(data,minsize,r)
  repeat{
    datos <- ordering(datos)
    listasubgrupost <- split(datos,datos[,p+1])
    porgrupos <- list()
    resultados <- data.frame()
    for (i in 1:length(listasubgrupost)){
      porgrupos[[i]] <- functionSAR2(listasubgrupost[[i]][,-(p+1)],minsize,r)
    }
    tgrupos <- vector()
    #system.time({
    for (i in 1:length(porgrupos)){
      tgrupos[i] <- testfunction(porgrupos[[i]])
    }
    porgrupos2 <- list()
    for (i in 1:length(tgrupos)){
      if (tgrupos[i]>=0.01) {#This is the critical ratio of the Bayes Factor
        porgrupos2[[i]] <- listasubgrupost[[i]]
      } else {
        porgrupos2[[i]] <- porgrupos[[i]]
      }
      names(porgrupos2[[i]])[length(names(porgrupos2[[i]]))] <- "group"
    }
    myfun <- function(x){
      y <- length(unique(as.numeric(x[,ncol(x)])))
      return(y)
    }
    if (length(porgrupos2)==1) {
      grupos_candidatos[[length(grupos_candidatos)+1]] <- data
      break}
    a <- sapply(porgrupos2,myfun)==1
    grupos_candidatos <- c(grupos_candidatos,porgrupos2[a])
    porgrupos2 <- porgrupos2[!a]
    if (length(porgrupos2)==0) break
    porgrupos2[[1]][,p+1] <- as.numeric(as.vector(porgrupos2[[1]][,p+1]))
    if (length(porgrupos2)>1) {
      for (i in 2:length(porgrupos2)) {
        porgrupos2[[i]][,p+1] <- (as.numeric(as.vector(porgrupos2[[i]][,p+1]))/min(as.numeric(as.vector(porgrupos2[[i]][,p+1]))))+max(porgrupos2[[i-1]][,p+1])
      }}
    datos <- do.call("rbind", porgrupos2) 
    is.null(datos)
    if(is.null(datos)) {
      break
    }
  }

  #################################################################
  #Recombining Step
  #################################################################
  for (i in 1:length(grupos_candidatos)) {
    grupos_candidatos[[i]][,p+1] <- rep(i,nrow(grupos_candidatos[[i]]))
    names(grupos_candidatos[[i]])[length(names(grupos_candidatos[[i]]))] <- "grupo"
  }
  datos <- do.call("rbind", grupos_candidatos) 
  datos <- foutliers(datos)
  #Aqui esta el problema
  #Si hay grupos menores al tamaño minimo hay que eliminarlos directamente
  datos <- datos[!datos[,p+1]%in%as.numeric(names(table(datos[,p+1]))[table(datos[,p+1])<minsize]),]
  datos <- ordering(datos)
  datos <- joining2(datos)[[1]]
  data <- original
  datos <- lastcombining(data,datos)
  datos <- joining2(datos)[[1]]
  datos <- laststep2(datos) #This is the only difference respect to previous official version
  return(datos)
}
