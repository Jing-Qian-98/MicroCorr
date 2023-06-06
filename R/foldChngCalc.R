foldChngCalc <-
function(data, nSampG1,nSampG2){
g1fold<-apply((data[,1:nSampG1]), 1, weightedMean)
g2fold<-apply((data[,(nSampG1+1):(nSampG1+nSampG2)]),1, weightedMean)
foldChange<-g2fold-g1fold
return(foldChange)
}
