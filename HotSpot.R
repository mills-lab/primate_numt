#Made by: Christian Brion 2019/12/08

#Read the insertion events (Numts) from files
ref<-read.table(file = "reference.txt",header = T,sep = "\t",stringsAsFactors = F)
poly<-read.table(file = "polymorphic.txt",header = T,sep = "\t",stringsAsFactors = F)
chrINFO<-read.table(file = "human_karyotype.txt",header = F,sep = "\t",stringsAsFactors = F)
chrINFO$cumulsum<-cumsum(as.numeric(chrINFO$V3))
chrINFO$cumulsum2<-c(0,cumsum(as.numeric(chrINFO$V3))[-nrow(chrINFO)])
deltabtchr<-50000000
maxchr<-max(chrINFO$cumulsum)+deltabtchr*nrow(chrINFO)
ref$start<-ref$start-1
poly$start<-poly$start-1

#correct the last bin for each chromosome to be the right size (and not 25mb)
ref$start2<-(0:(nrow(ref)-1))*25000000
ref$end2<-(1:(nrow(ref)))*25000000
poly$start2<-(0:(nrow(ref)-1))*25000000
poly$end2<-(1:(nrow(ref)))*25000000
for (n in 1:nrow(ref)) {
  if (ref$end[n]>chrINFO$V3[chrINFO$V4==ref$chr[n]]) {
    ref$end[n]<-chrINFO$V3[chrINFO$V4==ref$chr[n]]
    poly$end[n]<-chrINFO$V3[chrINFO$V4==ref$chr[n]]
  }
  add<-chrINFO$cumulsum2[chrINFO$V4==ref$chr[n]]
  ref$start2[n]<-ref$start[n]+add
  ref$end2[n]<-ref$end[n]+add
  poly$start2[n]<-ref$start[n]+add
  poly$end2[n]<-ref$end[n]+add
}

#general sums or mean per species
npoly<-colSums(poly[,4:9])
nref<-colSums(ref[,4:9])
xpoly<-colMeans(poly[,4:9])
xref<-colMeans(ref[,4:9])

#normalized standart deviation to the expectation based on the permutation (frequency data used in the paper)
poly$HSvalue<-sqrt(rowMeans((poly[,4:9]-xpoly)^2))
ref$HSvalue<-sqrt(rowMeans((ref[,4:9]-xref)^2))

#normalized mean
poly$density<-rowMeans(poly[,4:9]/npoly)*sum(npoly)/6
ref$density<-rowMeans(ref[,4:9]/nref)*sum(nref)/6

#genome size for permutation
minGeno<-min(ref$start2)
maxGeno<-max(ref$end2)

#####permutation 1000 (to run only one time) Uncomment for doing the permutations 
# permref<-ref[,c(1:3,10,11)]
# permpoly<-ref[,c(1:3,10,11)]
# permDref<-ref[,c(1:3,10,11)]
# permDpoly<-ref[,c(1:3,10,11)]
# temppoly<-poly
# tempref<-ref
# for (i in 1:1000) {
#   for (j in 1:6) {
#     insP<-runif(npoly[j],minGeno,maxGeno)
#     insR<-runif(nref[j],minGeno,maxGeno)
#     for (n in 1:nrow(ref)) {
#       temppoly[n,j+3]<-sum(insP>ref$start2[n] & insP<=ref$end2[n])
#       tempref[n,j+3]<-sum(insR>ref$start2[n] & insR<=ref$end2[n])
#     }
#   }
#   tempHSVpoly<-sqrt(rowMeans((temppoly[,4:9]-xpoly)^2))
#   tempHSVref<-sqrt(rowMeans((tempref[,4:9]-xref)^2))
#   permref<-cbind(permref,tempHSVref)
#   permpoly<-cbind(permpoly,tempHSVpoly)
#   tempDpoly<-rowMeans(temppoly[,4:9]/npoly)*sum(npoly)/6
#   tempDref<-rowMeans(tempref[,4:9]/nref)*sum(nref)/6
#   permDref<-cbind(permDref,tempDref)
#   permDpoly<-cbind(permDpoly,tempDpoly)
#   if (i %in% ((1:20)*50)) {
#     print(paste(i,date()))
#   }
# }
# save(permref,permpoly,permDref,permDpoly,file = "permutHS.RData")


load("permutHS.RData")

#calculate the local 5% significant threshold based on the permutation and the Zscore
permref$quant05<-permref$start
permpoly$quant05<-permpoly$start
permDref$quant05<-permDref$start
permDpoly$quant05<-permDpoly$start
permDref$quant05L<-permDref$start
permDpoly$quant05L<-permDpoly$start
zscore_ref<-permDpoly$start
zscore_poly<-permDpoly$start
for (n in 1:nrow(ref)) {
  permref$quant05[n]<-quantile(permref[n,6:1005],0.95)[[1]]
  permpoly$quant05[n]<-quantile(permpoly[n,6:1005],0.95)[[1]]
  permDref$quant05[n]<-quantile(permDref[n,6:1005],0.95)[[1]]
  permDpoly$quant05[n]<-quantile(permDpoly[n,6:1005],0.95)[[1]]
  permDref$quant05L[n]<-quantile(permDref[n,6:1005],0.05)[[1]]
  permDpoly$quant05L[n]<-quantile(permDpoly[n,6:1005],0.05)[[1]]
  zscore_ref[n]<-(ref$density[n]-mean(t(permDref[n,6:1005])))/sd(t(permDref[n,6:1005]))
  zscore_poly[n]<-(poly$density[n]-mean(t(permDpoly[n,6:1005])))/sd(t(permDpoly[n,6:1005]))
}


########ploting

pdf("HS-newmit.pdf",width = 8, height = 5.5)
plot(0,0,typ="n",xlim=c(0,maxchr),ylim=c(0,30),ylab="HS score",xlab="Primate genome",main="in references")
abline(h=mean(permref$quant05),lty=2)
for (i in 1:nrow(chrINFO)) {
  lines(c(0,chrINFO$V3[i])+chrINFO$cumulsum2[i]+deltabtchr*(i-1),c(0,0),col="blue",lwd=2)
  lines(rowMeans(ref[ref$chr==chrINFO$V4[i],2:3])+chrINFO$cumulsum2[i]+deltabtchr*(i-1),ref$HSvalue[ref$chr==chrINFO$V4[i]])
  #lines(rowMeans(ref[ref$chr==chrINFO$V4[i],2:3])+chrINFO$cumulsum2[i]+deltabtchr*(i-1),permref$quant05[ref$chr==chrINFO$V4[i]],col="red")
}

plot(0,0,typ="n",xlim=c(0,maxchr),ylim=c(0,5),ylab="HS score",xlab="Primate genome",main="polymorphic")
abline(h=mean(permpoly$quant05),lty=2)
for (i in 1:nrow(chrINFO)) {
  lines(c(0,chrINFO$V3[i])+chrINFO$cumulsum2[i]+deltabtchr*(i-1),c(0,0),col="blue",lwd=2)
  lines(rowMeans(poly[poly$chr==chrINFO$V4[i],2:3])+chrINFO$cumulsum2[i]+deltabtchr*(i-1),poly$HSvalue[poly$chr==chrINFO$V4[i]])
  #lines(rowMeans(poly[poly$chr==chrINFO$V4[i],2:3])+chrINFO$cumulsum2[i]+deltabtchr*(i-1),permpoly$quant05[poly$chr==chrINFO$V4[i]],col="red")
}

plot(0,0,typ="n",xlim=c(0,maxchr),ylim=c(0,30),ylab="density",xlab="Primate genome",main="in references")
abline(h=mean(permDref$quant05),lty=2)
for (i in 1:nrow(chrINFO)) {
  lines(c(0,chrINFO$V3[i])+chrINFO$cumulsum2[i]+deltabtchr*(i-1),c(0,0),col="blue",lwd=2)
  lines(rowMeans(ref[ref$chr==chrINFO$V4[i],2:3])+chrINFO$cumulsum2[i]+deltabtchr*(i-1),ref$density[ref$chr==chrINFO$V4[i]])
  #lines(rowMeans(ref[ref$chr==chrINFO$V4[i],2:3])+chrINFO$cumulsum2[i]+deltabtchr*(i-1),permDref$quant05[ref$chr==chrINFO$V4[i]],col="red")
}

plot(0,0,typ="n",xlim=c(0,maxchr),ylim=c(0,5),ylab="density",xlab="Primate genome",main="polymorphic")
abline(h=mean(permDpoly$quant05),lty=2)
for (i in 1:nrow(chrINFO)) {
  lines(c(0,chrINFO$V3[i])+chrINFO$cumulsum2[i]+deltabtchr*(i-1),c(0,0),col="blue",lwd=2)
  lines(rowMeans(poly[poly$chr==chrINFO$V4[i],2:3])+chrINFO$cumulsum2[i]+deltabtchr*(i-1),poly$density[poly$chr==chrINFO$V4[i]])
  #lines(rowMeans(poly[poly$chr==chrINFO$V4[i],2:3])+chrINFO$cumulsum2[i]+deltabtchr*(i-1),permDpoly$quant05[poly$chr==chrINFO$V4[i]],col="red")
}

plot(0,0,typ="n",xlim=c(0,maxchr),ylim=c(-15,15),ylab="zscore",xlab="Primate genome",main="in references")
abline(h=c(-2,2),lty=2,col="grey")
abline(h=c(-5,5),lty=2,col="black")
for (i in 1:nrow(chrINFO)) {
  lines(c(0,chrINFO$V3[i])+chrINFO$cumulsum2[i]+deltabtchr*(i-1),c(0,0),col="blue",lwd=2)
  lines(rowMeans(ref[ref$chr==chrINFO$V4[i],2:3])+chrINFO$cumulsum2[i]+deltabtchr*(i-1),zscore_ref[ref$chr==chrINFO$V4[i]])
  #lines(rowMeans(ref[ref$chr==chrINFO$V4[i],2:3])+chrINFO$cumulsum2[i]+deltabtchr*(i-1),permDref$quant05[ref$chr==chrINFO$V4[i]],col="red")
}


plot(0,0,typ="n",xlim=c(0,maxchr),ylim=c(-15,15),ylab="zscore",xlab="Primate genome",main="polymorphic")
abline(h=c(-2,2),lty=2,col="grey")
abline(h=c(-5,5),lty=2,col="black")
for (i in 1:nrow(chrINFO)) {
  lines(c(0,chrINFO$V3[i])+chrINFO$cumulsum2[i]+deltabtchr*(i-1),c(0,0),col="blue",lwd=2)
  lines(rowMeans(ref[ref$chr==chrINFO$V4[i],2:3])+chrINFO$cumulsum2[i]+deltabtchr*(i-1),zscore_poly[ref$chr==chrINFO$V4[i]])
  #lines(rowMeans(ref[ref$chr==chrINFO$V4[i],2:3])+chrINFO$cumulsum2[i]+deltabtchr*(i-1),permDref$quant05[ref$chr==chrINFO$V4[i]],col="red")
}
dev.off()
