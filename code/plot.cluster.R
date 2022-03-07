load("insulin/stem.cluster/STEM.Clusters.July2021.RData")
library(scales)
pval=c(0.0, 0.0, 3e-252, 1e-159, 5e-140, 1e-35, 5e-19, 2e-14, 2e-8, 3e-7, 4e-4, 6e-4)
profileid =c(35,25,32,45,29,19,42,26,10,41,37,40)


# alternatively, plot only
group1 <- c(4,5,10)
group2 <- c(2,3,12) 
group3 <- c(6,9)

pdf("insulin/stem.cluster/STEM.clusters.grouped.Jan2022.pdf")
par(mfrow=c(3,3),mar=c(2, 3, 1, 0.2),oma=c(6,4,4,1)+0.1)
for(j in 1:length(group1)){
  #yrange <- range(clusters[[group1[j]]]) + c(-0.05,0.05)
  yrange <- c(-1.0,1.5)
  plot(x=NULL,y=NULL,xlim=c(-0.1,5.1),ylim=yrange,xlab="",ylab="log2FC")
  
  abline(h=0,col="gray")
  for(i in 1:nrow(clusters[[group1[j]]])){
    lines(0:5,clusters[[group1[j]]][i,],type = "l",col=alpha("blue",alpha = 0.1))
  }
  avg <- apply(clusters[[group1[j]]],2,mean)
  lines(0:5,avg,col="red",lwd=2)
  mtext(text=paste0("p-value=",pval[group1[j]]),side = 3,line = 0,adj = 1,cex=0.8)
  mtext(text=paste0("# ",group1[j]),side=3,line=-1,adj=0,cex=0.8)
  if(j==1){mtext(text=paste0("Group",1),side = 1,line = -1,adj = 0)
      mtext(text="log2 FC",side=2,line=2.5)
    }
}

for(j in 1:length(group2)){
  #yrange <- range(clusters[[group2[j]]]) + c(-0.05,0.05)
  yrange <- c(-1.0,1.5)
  plot(x=NULL,y=NULL,xlim=c(-0.1,5.1),ylim=yrange,xlab="",ylab="log2FC")
  abline(h=0,col="gray")
  for(i in 1:nrow(clusters[[group2[j]]])){
    lines(0:5,clusters[[group2[j]]][i,],type = "l",col=alpha("blue",alpha = 0.1))
  }
  avg <- apply(clusters[[group2[j]]],2,mean)
  lines(0:5,avg,col="red",lwd=2)
  mtext(text=paste0("p-value=",pval[group2[j]]),side = 3,line = 0,adj = 1,cex=0.8)
  mtext(text=paste0("# ",group2[j]),side=3,line=-1,adj=0,cex=0.8)
  if(j==1){mtext(text=paste0("Group",2),side = 1,line = -1,adj = 0)
    mtext(text="log2 FC",side=2,line=2.5)
    }
}


for(j in 1:length(group3)){
  yrange <- c(-1.0, 1.5)
  plot(x=NULL,y=NULL,xlim=c(-0.1,5.1),ylim=yrange,xlab="",ylab="log2FC")
  abline(h=0,col="gray")
  for(i in 1:nrow(clusters[[group3[j]]])){
    lines(0:5,clusters[[group3[j]]][i,],type = "l",col=alpha("blue",alpha = 0.1))
  }
  avg <- apply(clusters[[group3[j]]],2,mean)
  lines(0:5,avg,col="red",lwd=2)
  mtext(text=paste0("p-value=",pval[group3[j]]),side = 3,line = 0,adj = 1,cex=0.8)
  mtext(text=paste0("# ",j=group3[j]),side=3,line=-1,adj=0,cex=0.8)
  if(j==1){mtext(text=paste0("Group",3),side = 1,line = -1,adj = 0)
    mtext(text="log2 FC",side=2,line=2.5)
  }
  if(j==2){ mtext(text="Time Points",side=1,line=3, cex=1.2)
  }
}

dev.off()

