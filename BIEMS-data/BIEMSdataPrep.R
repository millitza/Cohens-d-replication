BIEMSdataprep <- function(data){
  N <- table(data[,2])
  data <- cbind(data[,1],c(rep(1,N[1]),rep(0,N[2])),c(rep(0,N[1]),rep(1,N[2])))
}