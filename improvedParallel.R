library('mvtnorm')
library('parallel')
library('Itis')
library('distances')
library('scclust')
library('pdist')
# n: sample size
# round: number of batches/subgroups
# chunksize: subsample size in each batch/subgroup
# c: number of cores
n=10^9
round=50
chunksize=n/round
c <- 16
# dostep1: a function to do threshold clustering
# ichunk: a list of sequence of index, 
#         each core get one sequence of index from the list at one time 
dostep1 <- function(ichunk) {  
  require(distances)
  require(scclust)
  mydat <- as.matrix(x[ichunk,])
  Itis::Itis_fct_nomid(mydat, t, m, d, chunksize) 
}

# dostep2: matching the labels of all points to those of prototypes 
dostep2 <- function(ichunk2) {
  mclu<- unlist(mres[ichunk2])
  ll<- sum(nprotos[1:(ichunk2)])+1
  ul<- sum(nprotos[1:(ichunk2+1)])
  ft <- Itis::fastJoin2(mclu, ktag[ll:ul])
  return(ft)
}

# test: main function
test <- function(cls, n, d, t, m, k,chunksize) { 
  ########step1 (parallelized)
  # paralle the threshold clustering step
  ichunk=splitIndices(n, ceiling(n/chunksize)) # create a list of sequence of index
  clusterExport(cls,c("x","dostep1","n","d","t","m","k","chunksize"))  # export variables to worker environment
  time0=proc.time()
  mres <- clusterApplyLB(cls,ichunk,dostep1) # execute dostep1 in each cores, mres is a list of b list(b is number of batches), each list contains of prototypes and labels in this batch
  print(proc.time()-time0)# time of parallelizabel time in step1
  protos<- matrix(unlist(lapply(mres, function(u) t(u[1][[1]]) ) ),ncol=2,byrow=T)# prototype coordinates
  nprotos<<- c(0, unlist(lapply(mres, function(u) length(unlist(u[1]))/2))) # number of prototypes in each batch  
  ########

  ########step2 (sequential)
  # run the k-means with prototypes
  kfit <- kmeans(protos, centers = k, nstart = 10)
  ktag<<- kfit$cluster # get kmeans labels for prototypes
  mres<<- lapply(mres, function(u) u[-1] )#remove prototype coordinates to save memory,only keep labels of all points according to prototypes, because we will pass mres to workers later
  #write.table(mres,file='mres.txt', row.names=F)
  #write.table(ktag,file='ktag.txt', row.names=F)
  #write.table(nprotos,file='nprotos.txt', row.names=F)
  ########

  ########step3 (parallelized)
  # match kmeans labels of prototypes to all points
  clusterExport(cls,c("mres","ktag","nprotos")) # export variables to worker environment
  nchunk=n/chunksize # nchunk is just number of batches
  ichunk2=splitIndices(nchunk, ceiling(nchunk))# create a list of sequence of index, each list has 1 number(index)
  time0=proc.time()
  mres2<- clusterApplyLB(cls,ichunk2,dostep2) # execute dostep2 in each cores
  print(proc.time()-time0)# time of parallelizabel time in step3
  mclusterlabel<- cbind(unlist(lapply(mres2, function(u) t(u) )),y)#combine labels of kmeans and true labels of all data points
  #write.table(mclusterlabel,file='mclusterp1.txt', row.names=F)
  ########   
  rm(list=setdiff(ls(),c("mclusterlabel","k")))
  
  ########step4 (sequential)
  # match the true lables and labels from k-means
  match <- rep(0,k)# k is the number of different true labels, k is 3 in our simulation
  for (i in 1:k){
    match[i] <- which.max(as.matrix(table(c(1:k,mclusterlabel[mclusterlabel[,2]==i,1]))))# for points with true label of i, find the most frequent predicted label
  }
  print(match)
  mclusterlabel[,1] <- 0 - mclusterlabel[,1] # change the predicted label to negative value, to avoid error when replacing the labels in next loop
  for (i in 1:k){ # k is the number of different true labels, k is 3 in our simulation
    mclusterlabel[mclusterlabel[,1]==-match[i],1] <- i # replace the predicted label with matching label
  }
  ########

  ######## rest part (sequential)
  #write.table(mclusterlabel,file='clusterlabelp1.txt', row.names=F)
  ctab=table(mclusterlabel[,1],mclusterlabel[,2],dnn=c("Predicted","True")) # create a confusion matrix
  #write.table(ctab,file='ctabp1.txt')
  acc=sum(diag(ctab))/n # calculate the accuracy
  ########
  return(acc)
}

# simulation
d=2; t=2; m=3; k=3
x0<-cbind(rbind(rmvnorm(0.5*n,mean=c(1,2),sigma=matrix(c(1,0,0,0.5),ncol=2,byrow=T)),
         rmvnorm(0.3*n,mean=c(7,8),sigma=matrix(c(2,0,0,1),ncol=2,byrow=T)),
         rmvnorm(0.2*n,mean=c(3,5),sigma=matrix(c(3,0,0,4),ncol=2,byrow=T))),
         c(rep(1,0.5*n),rep(2,0.3*n),rep(3,0.2*n)))
x0<- x0[sample.int(n),]
x<- x0[,1:2]
y<- x0[,3]
rm(x0)# remove x0 to save memory
cls=makeCluster(c)# make a cluster of cores

# call the main function and record time
tim0=proc.time()
acc=test(cls, n, d, t, m, k, chunksize)
(proc.time()-tim0)# total time(data simulation not included)


