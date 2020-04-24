#################
##  Rdsm Itis  ##
#################
library('mvtnorm')
library('Rdsm')
library('parallel')
library('Rcpp')
library('distances')
library('scclust')
library('pdist')
#############################
##    parallel function    ##
#############################
# The function dostep1 generates the prototypes in shared memory parallelization using package 'Rdsm' and 'parallel'.
# The function dostep2 matches the cluster labels of original points to the labels of prototypes by k-means clustering. 
# The output includes mclusterlabel and proto in main function.
# mclusterlabel is a n by 4 shared matrix, 
#               1st column is an index for tasks, 
#               2nd column is index for myproto in its own chunk of prototypes,
#               3rd column is cluster label assigned by simulation,
#               4th column is the predicted label. 
# proto is a f by 3 matrix from mapreduce, 
#       f is the number of prototypes,
#       the first 2 columns are prototypes,
#       3rd column is an index for tasks.
dostep1 <- function(n,t,d,m,dat,mclusterlabel,cumm) {  
  require(parallel)
  require(Rcpp)
  require(distances)
  require(scclust)
  require(RcppAlgos)
  cppFunction('NumericMatrix fastAgg(NumericMatrix orgMeans, IntegerVector cats) {
            
              //orgMeans is the original dataset
              //cats is the categories, numbered 0 to n-1
            
              //Store dimensons of the for loop
              long catLeng = cats.length();
              int numCols = orgMeans.ncol();
            
              //initialize number of categories, number of observations in each category, aggregated means
              long numCats = max(cats);
              IntegerVector catSize(numCats+1);
              NumericMatrix aggMeans(numCats+1,numCols);
            
              for(int j = 0; j < numCols; j++){
                for(long i = 0; i < catLeng; i++ ){
                  //Different instructions if j = 0 and if j notequal zero
                  if(j == 0){
            
                  //Increase the category total
                  catSize[cats[i]]++;
            
                  //Update the means
                  aggMeans(cats[i],j) = (double)(catSize[cats[i]]-1)/catSize[cats[i]]*aggMeans(cats[i],j) + (double)1/catSize[cats[i]]*orgMeans(i,j);
                  }else{
                  aggMeans(cats[i],j) = (double)aggMeans(cats[i],j) + (double)1/catSize[cats[i]]*orgMeans(i,j);
                  }
                }
              }
              return aggMeans;
              }')
  cppFunction('IntegerVector fastJoin2(IntegerVector mer1, IntegerVector mer2) {

              //mer1 is original data matrix
              //mer2 in final cluster
              //This function want to perform inner join
            
              //Store dimensons of the for loop
             
              long numN = mer1.length();
              IntegerVector mer3(numN);         
            
              for(long j = 0; j < numN; j++){
                mer3(j)=mer2( mer1(j) );
              }
              return mer3;
              }')
  myidxs <- getidxs(n)
  mclusterlabel[myidxs,3] <- dat[myidxs,3]
  mydat <- as.matrix(dat[myidxs,1:2])
  n1 <- length(myidxs)   
  
  clusterlabel <- rep(NA, n1)
  for(i in 1:m){
    my_dist <- distances(mydat)
    my_clustering_new <- sc_clustering(my_dist, t)
    aggdata_old <- fastAgg(as.matrix(mydat), my_clustering_new)
    if(i == 1){
      clusterlabel <- as.integer(my_clustering_new)
    }
    else{
      clusterlabel <- fastJoin2(my_clustering_old, my_clustering_new)
    }
    mydat <- aggdata_old
    my_clustering_old <- clusterlabel
  }
  mclusterlabel[myidxs,2] <- clusterlabel
  myproto <- mydat
  mynpro <- nrow(myproto)
  rdsmlock('lck')
  cumm[1,1]<- cumm[1,1] + 1
  cumm[1,2]<- cumm[1,2] + mynpro
  cum=cumm[1,1]
  rdsmunlock('lck')
  mclusterlabel[myidxs,1] <- cum
  myproto <- cbind(myproto,rep(cum,mynpro))
  return(myproto)
}

dostep2 <- function(klabel,mclusterlabel) {
  require(parallel)
  require(Rcpp)
  cppFunction('IntegerVector fastJoin2(IntegerVector mer1, IntegerVector mer2) {
              //mer1 is original data matrix
              //mer2 in final cluster
              //This function want to perform inner join            
              //Store dimensons of the for loop             
              long numN = mer1.length();
              IntegerVector mer3(numN);                     
              for(long j = 0; j < numN; j++){
                mer3(j)=mer2( mer1(j) );
              }
              return mer3;
              }')
  mywok <- myinfo$id
  mklabel <- klabel[klabel[,1]==mywok,2]
  mcluster <- mclusterlabel[mclusterlabel[,1]==mywok,2]
  ft <- fastJoin2(mcluster, mklabel)
  mclusterlabel[mclusterlabel[,1]==mywok,4] <- ft
  return(0)
}
#############################
##      main function      ##
#############################
# cumm is a 1 by 2 shared matrix, 
#      1st number is accumulated number of tasks, 
#      2nd number is accumulated number of prototypes.
# klabel is a f by 2 shared matrix,
#        f is the number of prototypes,
#        1st column is an index for tasks,
#        2nd column is labels by kmeans clustering.
# ctab is a tabel for cluster labels, the rows are simulated label, the columns are predicted labels.
#      is exported in ctab7.txt.
# acc is the prediction accuracy.
#     is exported in ctab7.txt.
test <- function(datt,cls, n, d, t, m, k) {
  library('parallel')
  library('Rdsm')
  library('Rcpp')
  library('distances')
  library('scclust')
  library('pdist')
  # make shared variables and lock in multiple cores
  mgrinit(cls)  # initial the cores
  makebarr(cls) # shared bar
  mgrmakevar(cls,"dat",n,d+1) # shared variable, must be a matrix
  mgrmakevar(cls,"mclusterlabel",n,4)
  mgrmakevar(cls,"cumm",1,2)
  mgrmakelock(cls,"lck") # shared lock
  dat[,]<- datt # initial value of shared varibales
  cumm[,]<- 0
  clusterExport(cls,c("dostep1","dostep2","n","d","t","m","k"))  # export variables to worker environment
  mproto <- clusterEvalQ(cls,dostep1(n,t,d,m,dat,mclusterlabel,cumm)) # execute dostep1 in each cores
  proto <- matrix(unlist(lapply(mproto, function(u) c(t(u)) ) ), ncol=d+1, byrow=T) # mapreduce
  print(cumm[,])
  kfit <- kmeans(proto[,1:2], centers = k, nstart = 10) # k-means clustering
  mgrmakevar(cls,"klabel",cumm[1,2],2)
  klabel[,1] <- proto[,3]
  klabel[,2] <- kfit$cluster
  clusterEvalQ(cls,dostep2(klabel,mclusterlabel)) # execute dostep2 in each cores   
  #write.table(mclusterlabel[,],file='clusterlabelbef1.txt', row.names=F)
  #ctab0=table(mclusterlabel[,4],mclusterlabel[,3],dnn=c("Predicted","True"))
  #write.table(ctab0,file='ctab0.txt')
  match <- rep(0,k) # match the lables of simulation and k-means
  for (i in 1:k){
    match[i] <- which.max(as.matrix(table(c(1,2,3,mclusterlabel[mclusterlabel[,3]==i,4]))))
  }
  print(match)
  mclusterlabel[,4] <- 0 - mclusterlabel[,4]
  for (i in 1:k){
    mclusterlabel[mclusterlabel[,4]==-match[i],4] <- i
  }
  #write.table(mclusterlabel[,],file='clusterlabelaft1.txt', row.names=F)
  ctab=table(mclusterlabel[,4],mclusterlabel[,3],dnn=c("Predicted","True"))
  #write.table(ctab,file='ctab7.txt')
  acc=sum(diag(ctab))/n
  print(acc)
  return(acc)
}
#############################
##  simulate, test, result ##
#############################
# simulation
d=2; t=2; m=3; k=3
n=10^8
c1 <- 4
c2 <- 8
c3 <- 16
c4 <- 32
x<-cbind(rbind(rmvnorm(0.5*n,mean=c(1,2),sigma=matrix(c(1,0,0,0.5),ncol=2,byrow=T)),
         rmvnorm(0.3*n,mean=c(7,8),sigma=matrix(c(2,0,0,1),ncol=2,byrow=T)),
         rmvnorm(0.2*n,mean=c(3,5),sigma=matrix(c(3,0,0,4),ncol=2,byrow=T))),
         c(rep(1,0.5*n),rep(2,0.3*n),rep(3,0.2*n)))
datt<- x[sample.int(n),]

rm(list=setdiff(ls(),c("datt","n","d","t","m","k","c1","c2","c3","c4","dostep1","dostep2","test")))
cls1=makeCluster(c1)
tim0=proc.time()
acc=test(datt,cls1, n, d, t, m, k)
(tim1=proc.time()-tim0)
stopCluster(cls1)

rm(list=setdiff(ls(),c("datt","n","d","t","m","k","c1","c2","c3","c4","dostep1","dostep2","test")))
cls2=makeCluster(c2)
tim0=proc.time()
acc=test(datt,cls2, n, d, t, m, k)
(tim1=proc.time()-tim0)
stopCluster(cls2)

rm(list=setdiff(ls(),c("datt","n","d","t","m","k","c1","c2","c3","c4","dostep1","dostep2","test")))
cls3=makeCluster(c3)
tim0=proc.time()
acc=test(datt,cls3, n, d, t, m, k)
(tim1=proc.time()-tim0)
stopCluster(cls3)
