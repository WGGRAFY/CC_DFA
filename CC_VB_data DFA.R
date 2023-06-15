rm(list=ls())

library(dplyr)

# Load VB cohort estimates
load("C:/MY_FILES/ANALYSES/Growth analysis - California current/vbpar_2021.Rdata")

# Loop
for(region in c("ECV","Monterey","Conception")) { # region <- "ECV"

  alldat <- data.frame(array(NA,dim=c(0,10)))
  colnames(alldat) <- c("linf","se.linf","p.linf","k","se.k","p.k","cohort","region","sp","stt")
  
  for(sp in c("Pacific hake","petrale sole","Pacific sanddab","sablefish","lingcod","shortbelly rockfish","darkblotched rockfish")) { # sp <- "Pacific hake"  
  
    print(paste(region, sp, sep=","))
    
    vb <- data.frame(vbpar[,,sp,region])
    vb <- vb[!is.na(vb$p.linf),]
    vb <- vb[vb$p.linf<0.05,]
    vb <- vb[vb$p.k<0.05,]
    
    if(dim(vb)[1]>0){
      vb$cohort <- rownames(vb)
      vb$region <- region
      vb$sp <- sp
      rownames(vb) <- c()
      vb$stt <- (vb[,"linf"]-mean(vb$linf,na.rm=T))/sd(vb$linf,na.rm=T) 
      alldat <- rbind(alldat,vb)
    }
    
    rm(vb)
    
  } # Species loop
  
  save(alldat, file=paste("C:/MY_FILES/ANALYSES/Growth analysis - California current/VB DFA data 2021/alldat_",region,".Rdata",sep=''))
  
  dfa <- data.frame(array(NA,dim=c(length(min(alldat$cohort):max(alldat$cohort)),7)))
  colnames(dfa) <- c("Pacific hake","petrale sole","Pacific sanddab","sablefish","lingcod","shortbelly rockfish","darkblotched rockfish")
  rownames(dfa) <- c(min(alldat$cohort):max(alldat$cohort))
  
  for (i in c("Pacific hake","petrale sole","Pacific sanddab","sablefish","lingcod","shortbelly rockfish","darkblotched rockfish")) {
    for (j in min(alldat$cohort):max(alldat$cohort)) {
      if (length(alldat[alldat$sp==i & alldat$cohort==j,"stt"])>0) { dfa[as.character(j),i] <- alldat[alldat$sp==i & alldat$cohort==j,"stt"] }
    }
  }
  
  save(dfa, file=paste("C:/MY_FILES/ANALYSES/Growth analysis - California current/VB DFA data 2021/dfa_",region,".Rdata",sep=''))
  
} # Region loop

