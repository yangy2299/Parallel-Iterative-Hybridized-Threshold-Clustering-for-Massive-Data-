library(SparkR)
library('mvtnorm')
library('parallel')
library('Itis')
library('distances')
library('scclust')
library('pdist')

n=10^8
d=2; t=2; m=3; k=3
x0<-cbind(rbind(rmvnorm(0.5*n,mean=c(1,2),sigma=matrix(c(1,0,0,0.5),ncol=2,byrow=T)),
         rmvnorm(0.3*n,mean=c(7,8),sigma=matrix(c(2,0,0,1),ncol=2,byrow=T)),
         rmvnorm(0.2*n,mean=c(3,5),sigma=matrix(c(3,0,0,4),ncol=2,byrow=T))),
         c(rep(1,0.5*n),rep(2,0.3*n),rep(3,0.2*n)))
x0<- x0[sample.int(n),]
x<- x0[,1:2]
y<- x0[,3]
rm(x0)

#########################################
dostep1 <- function(chunks) {
  require(distances)
  require(scclust)
  mydat <- as.matrix(x[unlist(chunks),])
  Itis::Itis_fct_nomid(mydat, t, m, d, chunksize) 
}

dostep2 <- function(chunks) {
  mclu<- unlist(mres[unlist(chunks)])
  ll<- sum(nprotos[1:(unlist(chunks))])+1
  ul<- sum(nprotos[1:(unlist(chunks)+1)])
  ft <- Itis::fastJoin2(mclu, ktag[ll:ul])
  return(ft)
}
#########################################
sparkR.session()
test <- function(round){
  #chunksize=n/round
  ichunk1=splitIndices(n, ceiling(round))
  mres <- spark.lapply(ichunk1, dostep1)
  protos<- matrix(unlist(lapply(mres, function(u) t(u[1][[1]]) ) ),ncol=2,byrow=T)
  nprotos <<- c(0, unlist(lapply(mres, function(u) length(unlist(u[1]))/2)))

  kfit <- kmeans(protos, centers = k, nstart = 10)
  ktag <<- kfit$cluster 
  mres <<- lapply(mres, function(u) u[-1] )

  ichunk2=splitIndices(round, ceiling(round))
  mres2<- spark.lapply(ichunk2,dostep2)
  mclusterlabel<- cbind(unlist(lapply(mres2, function(u) t(u))),y)
  match <- rep(0,k)
  for (i in 1:k){
    match[i] <- which.max(as.matrix(table(c(1:k,mclusterlabel[mclusterlabel[,2]==i,1]))))
  }
  print(match)
  mclusterlabel[,1] <- 0 - mclusterlabel[,1] 
  for (i in 1:k){ 
    mclusterlabel[mclusterlabel[,1]==-match[i],1] <- i # replace the predicted label with matching label
  }
  ctab=table(mclusterlabel[,1],mclusterlabel[,2],dnn=c("Predicted","True"))
  acc=sum(diag(ctab))/n
  print(acc)
  return(acc)
}

round=50
chunksize=n/round
tim0=proc.time()
acc1=test(round)
(proc.time()-tim0)
print(round)
