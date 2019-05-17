get_data<-function(path="/mnt/powervault/jonhall/Desktop/Popstruct/Save/Pemberton/SGEM/single_start/loglik_",ext=".txt",nfiles=100,method="SGEM"){
# extract logLikelihood scores from individual files nrep number
# by Jon Ahlinder
  

  for(i in 1:nfiles){
   # cat(sprintf('i = %d',i))
    filename<-paste(path,i,ext,sep="")
    data<-read.table(filename)
    l1<-length(data$V1)
   # cat(sprintf("l1 = %d\n",l1))
    it<-seq(1,l1)
 #   print(head(it))
    tmp<-cbind(data,it,rep(method,l1)) 
  #  print(head(tmp))
    if(i>1){ 
      dataf<-rbind(dataf,tmp)
    }
    else{
      dataf<-tmp
    }
  }

  colnames(dataf)<-c("logL","Iteration","Method")
  return(dataf)

}
