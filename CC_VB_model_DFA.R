rm(list=ls())

library(MARSS);library(xtable);require(broom);library(dplyr); library(naniar)

# Load covariates

load("C:/MY_FILES/ANALYSES/Growth analysis - California current/AR1_temperature_regions_2023.Rdata")
temp <- data.frame(tempest)
colnames(temp) <- colnames(tempest)
  
# Plot temp 
#windows()
jpeg("C:/MY_FILES/ANALYSES/Growth analysis - California current/DFA/temp.jpeg", width=10, height=3, units="in", res=250)
par(mfrow=c(1,3))
plot(rownames(temp),temp$`Conception_<=183m`,ylim=c(min(na.omit(temp)),max(na.omit(temp))),type="b",ylab="",xlab="",main="Conception")
lines(rownames(temp),temp$`Conception_184-550m`,type="b")
lines(rownames(temp),temp$`Conception_>=550m`,type="b")

plot(rownames(temp),temp$`Monterey_<=183m`,ylim=c(min(na.omit(temp)),max(na.omit(temp))),type="b",ylab="",xlab="",main="Monterey")
lines(rownames(temp),temp$`Monterey_184-550m`,type="b")
lines(rownames(temp),temp$`Monterey_>=550m`,type="b")

plot(rownames(temp),temp$`ECV_<=183m`,ylim=c(min(na.omit(temp)),max(na.omit(temp))),type="b",ylab="",xlab="",main="ECV")
lines(rownames(temp),temp$`ECV_184-550m`,type="b")
lines(rownames(temp),temp$`ECV_>=550m`,type="b")
dev.off()

# Create a fishing step-change covariate
fishing <- data.frame(array(NA,dim=c(length(1977:2018),1)))
rownames(fishing) <- 1977:2018
colnames(fishing) <- "fishing"
fishing[as.character(1977:2000),"fishing"] <- 0.8
fishing[as.character(2001:2018),"fishing"] <- 0.2
# fishing[as.character(1977:2000),"fishing"] <- jitter(fishing[as.character(1977:2000),"fishing"],factor = 2)
# fishing[as.character(2001:2018),"fishing"] <- jitter(fishing[as.character(2001:2018),"fishing"],factor = 2)
#save(fishing, file = "C:/MY_FILES/ANALYSES/Growth analysis - California current/fishing.Rdata")

### For each region

for(region in c("ECV","Monterey","Conception")) { 
#region <- "ECV" 
#region <- "Monterey"
#region <- "Conception"
  
# Load data generated
load(paste("C:/MY_FILES/ANALYSES/Growth analysis - California current/VB DFA data 2021/dfa_",region,".Rdata",sep=''))

# Only keep lingcod in ECV  
if(region!="ECV") {dfa<-subset(dfa, select = -c(lingcod))}
    
# Remove columns with no data (NA only) if needed  
dfa <- dfa[,which(unlist(lapply(dfa, function(x) !all(is.na(x)))))]

# Keep only 1977 onwards (temp data available)
if(min(as.numeric(rownames(dfa)))<1977) {
dfa <- dfa[as.character(1977:rownames(dfa)[length(rownames(dfa))]),] }
    
# Transpose the dataset so time series are in columns
dfa <- t(dfa)

# Plot the data
jpeg(paste("C:/MY_FILES/ANALYSES/Growth analysis - California current/DFA/",region,"_DFA data.jpeg",sep=''), width=15, height=5, units="in", res=250)
spp = rownames(dfa)
par(mfcol=c(2,4), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
for(i in spp){
  plot(colnames(dfa),dfa[i,],xlab="",ylab="Linf", bty="L", pch=16, col="blue", type="b")
  title(i)
}
dev.off()

# Number of time series
N.ts = dim(dfa)[1]
# Length of time series
TT = dim(dfa)[2]


### Test with different numbers of trends

# set new control parameters
cntl.list = list(minit=200, maxit=1000, allow.degen=FALSE)

# set up forms of R matrices
levels.R = c("diagonal and equal","diagonal and unequal","equalvarcov") #,"diagonal and unequal","unconstrained","equalvarcov"

### HOW MANY TRENDS DO YOU WANT TO TEST? ###

n <- 3 # number substracted to the max number of trends  

model.data = data.frame()
# fit lots of models & store results
# NOTE: this will take a long time to run!!!
for(R in levels.R) {
  for(m in 1:(N.ts-n)) {
    dfa.model = list(A="zero", R=R, m=m)
    mm = MARSS(dfa, model=dfa.model, control=cntl.list, form="dfa", z.score=TRUE)
    model.data = rbind(model.data,
                       data.frame(R=R,
                                  m=m,
                                  logLik=mm$logLik,
                                  K=mm$num.params,
                                  AICc=mm$AICc,
                                  stringsAsFactors=FALSE))
    assign(paste("mm", m, R, sep="."), mm)
  }} # end R loop


### Produce model comparison table

# AIC
model.data$AIC = model.data$AICc
# calculate delta-AICc
model.data$delta.AICc = model.data$AICc - min(model.data$AICc)
# calculate Akaike weights
wt = exp(-0.5*model.data$delta.AICc)
model.data$Ak.wt = wt/sum(wt)
# sort results
model.tbl = model.data[order(model.data$AICc),-4]
# drop AICc from table
# calculate cumulative wts
model.tbl$Ak.wt.cum = cumsum(model.tbl$Ak.wt)
model.tbl = model.tbl[,-4]
print(model.tbl)

write.csv(model.tbl, file=paste("C:/MY_FILES/ANALYSES/Growth analysis - California current/DFA/",region,"_model table.csv",sep=''))

### Get the best model

best.model = model.tbl[1,]
fitname = paste("mm",best.model$m,best.model$R,sep=".")
best.fit = get(fitname)


### Varimax (the varimax rotation seeks a rotation matrix H that creates the largest difference between factor loadings)

Z.est = coef(best.fit, type="matrix")$Z

# get the inverse of the rotation matrix
H.inv = 1
if(ncol(Z.est)>1) H.inv = varimax(coef(best.fit, type="matrix")$Z)$rotmat


### Rotations

# rotate factor loadings
Z.rot = Z.est %*% H.inv
# rotate trends
trends.rot = solve(H.inv) %*% best.fit$states


### Get the best model with 1 trend (if applicable)

if (best.model$m>1) {
  best.model1 <- model.tbl[model.tbl$m==1,]; best.model1 <- best.model1[1,]
  fitname1 = paste("mm",best.model1$m,best.model1$R,sep=".")
  best.fit1 = get(fitname1)
  Z.est1 = coef(best.fit1, type="matrix")$Z
  H.inv1 = 1
  Z.rot1 = Z.est1 %*% H.inv1
  trends.rot1 = solve(H.inv1) %*% best.fit1$states
  ts.trends1 = t(trends.rot1)
}


### Plot of the trends

jpeg(paste("C:/MY_FILES/ANALYSES/Growth analysis - California current/DFA/",region,"_DFA trends.jpeg",sep=''), width=8, height=3, units="in", res=500)

#windows()
# get ts of trends
ts.trends = t(trends.rot)
par(mfcol=c(ceiling(dim(ts.trends)[2]),2), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
# loop over each trend
for(i in 1:dim(ts.trends)[2]) {
  # set up plot area
  plot(ts.trends[,i],
       ylim=c(-1.1,1.1)*max(abs(ts.trends)), 
       type="n", lwd=2, bty="L", 
       xlab="", ylab="", xaxt="n", yaxt="n")
  # draw zero-line
  abline(h=0, col="gray")
  # plot trend line
  par(new=TRUE)
  plot(ts.trends[,i],
       ylim=c(-1.1,1.1)*max(abs(ts.trends)), 
       type="l", lwd=2, bty="L", 
       xlab="", ylab="", xaxt="n")
  # add panel labels
  mtext(paste("Trend",i,sep=" "), side=3, line=0.5)
  axis(1,(0:dim(dfa)[2])+1,min(na.omit(as.numeric(colnames(dfa))))+0:dim(dfa)[2])
  # add best model with one trend if applicable
  if (best.model$m>1) {
    lines(ts.trends1[,1], col="gray")
  }
} # end i loop (trends)


### Plot the factor loadings

spp = rownames(dfa)
minZ = 0.00         # minZ = 0.05 # choose min level here 
m=dim(trends.rot)[1]
ylims = c(-1.1*max(abs(Z.rot)), 1.1*max(abs(Z.rot)))
#par(mfrow=c(ceiling(m/2),2), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
for(i in 1:m) {
  plot(c(1:N.ts)[abs(Z.rot[,i])>minZ], as.vector(Z.rot[abs(Z.rot[,i])>minZ,i]),
       type="h", lwd=2, xlab="", ylab="", xaxt="n", ylim=ylims, xlim=c(0,N.ts+1))
  for(j in 1:N.ts) {
    if(Z.rot[j,i] > minZ) {text(j, -0.05, spp[j], srt=90, adj=1, cex=.7)}
    if(Z.rot[j,i] < -minZ) {text(j, 0.05, spp[j], srt=90, adj=0, cex=.7)}
    abline(h=0, lwd=1, col="gray")
  } # end j loop
  mtext(paste("Factor loadings on trend",i,sep=" "),side=3,line=.5)
} # end i loop

dev.off()

# Save trends
save(ts.trends, file=paste("C:/MY_FILES/ANALYSES/Growth analysis - California current/DFA/",region,"_trends.Rdata",sep=''))

# Save factor loadings
save(Z.rot, file=paste("C:/MY_FILES/ANALYSES/Growth analysis - California current/DFA/",region,"_loadings.Rdata",sep=''))


### Examine model fits

# create getDFAfits function
getDFAfits <- function(MLEobj, alpha=0.05, covariates=NULL) { # MLEobj <- best.fit ; alpha=0.05 ; covariates=NULL
  fits <- list()
  Ey <- MARSShatyt(MLEobj) # for var() calcs
  ZZ <- coef(MLEobj, type="matrix")$Z # estimated Z
  nn <- nrow(ZZ) # number of obs ts
  mm <- ncol(ZZ) # number of factors/states
  TT <- ncol(Ey$ytT)  # number of time steps
  ## check for covars
  if(!is.null(covariates)) {
    DD <- coef(MLEobj, type="matrix")$D
    cov_eff <- DD %*% covariates
  } else {
    cov_eff <- matrix(0, nn, TT)
  }
  ## model expectation
  fits$ex <- ZZ %*% MLEobj$states + cov_eff
  ## Var in model fits
  VtT <- MARSSkfss(MLEobj)$VtT
  VV <- NULL
  for(tt in 1:TT) {
    RZVZ <- coef(MLEobj, type="matrix")$R - ZZ%*%VtT[,,tt]%*%t(ZZ)
    SS <- Ey$yxtT[,,tt] - Ey$ytT[,tt,drop=FALSE] %*% t(MLEobj$states[,tt,drop=FALSE])
    VV <- cbind(VV,diag(RZVZ + SS%*%t(ZZ) + ZZ%*%t(SS)))
  }
  SE <- sqrt(VV)
  ## upper & lower (1-alpha)% CI
  fits$up <- qnorm(1-alpha/2)*SE + fits$ex
  fits$lo <- qnorm(alpha/2)*SE + fits$ex
  return(fits)
}

# get the fits
fit.b = getDFAfits(best.fit)

# plot the fits
jpeg(paste("C:/MY_FILES/ANALYSES/Growth analysis - California current/DFA/",region,"_DFA fits.jpeg",sep=''), width=15, height=5, units="in", res=250)
ylbl = rownames(dfa)
w_ts = seq(dim(dfa)[2])
par(mfcol=c(2,4), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
for(i in 1:N.ts) {
  up <- fit.b$up[i,]
  mn <- fit.b$ex[i,]
  lo <- fit.b$lo[i,]
  plot(w_ts, mn, xlab="", ylab=ylbl[i], xaxt="n", type="n", cex.lab=1.2,
       ylim = c(min(lo),max(up)))
  axis(1,(0:dim(dfa)[2])+1,min(na.omit(as.numeric(colnames(dfa))))+0:dim(dfa)[2])
  points(w_ts, dfa[i,], pch=16, col="blue")
  lines(w_ts, up, col="darkgray")
  lines(w_ts, mn, col="black", lwd=2)
  lines(w_ts, lo, col="darkgray")
}
dev.off()


### Add covariates

# Run the best model with the added temp covariate 

if (region == "ECV") {

temp1 = t(matrix(temp[colnames(dfa),"ECV_<=183m"]))
temp2 = t(matrix(temp[colnames(dfa),"ECV_184-550m"]))
temp3 = t(matrix(temp[colnames(dfa),"ECV_>=550m"]))
fishing. = t(matrix(fishing[colnames(dfa),"fishing"]))
model.list=list(m=best.model$m, R=best.model$R)
mm.temp1 = MARSS(dfa, model=model.list, z.score=TRUE, form="dfa", control=cntl.list, covariates=temp1)
mm.temp2 = MARSS(dfa, model=model.list, z.score=TRUE, form="dfa", control=cntl.list, covariates=temp2)
mm.temp3 = MARSS(dfa, model=model.list, z.score=TRUE, form="dfa", control=cntl.list, covariates=temp3)
mm.temp = MARSS(dfa, model=model.list, z.score=TRUE, form="dfa", control=cntl.list, covariates=rbind(temp1,temp2,temp3))
mm.fish = MARSS(dfa, model=model.list, z.score=TRUE, form="dfa", control=cntl.list, covariates=fishing.)
mm.both1 = MARSS(dfa, model=model.list, z.score=TRUE, form="dfa", control=cntl.list, covariates=rbind(temp1,fishing.))
mm.both2 = MARSS(dfa, model=model.list, z.score=TRUE, form="dfa", control=cntl.list, covariates=rbind(temp2,fishing.))
mm.both3 = MARSS(dfa, model=model.list, z.score=TRUE, form="dfa", control=cntl.list, covariates=rbind(temp3,fishing.))
mm.both = MARSS(dfa, model=model.list, z.score=TRUE, form="dfa", control=cntl.list, covariates=rbind(temp1,temp2,temp3,fishing.))

# Compare AICc between models with and without covariate   
cov.tbl <- cbind(model=c("no covars",
                         "temp <=183m",
                         "temp 184-550m",
                         "temp >=550m",
                         "all temp",
                         "fishing",
                         "fishing & temp <=183m",
                         "fishing & temp 184-550m",
                         "fishing & temp >=550m",
                         "all covars"),
                 AICc=round(c(best.fit$AICc,
                              mm.temp1$AICc,
                              mm.temp2$AICc,
                              mm.temp3$AICc,
                              mm.temp$AICc,
                              mm.fish$AICc,
                              mm.both1$AICc,
                              mm.both2$AICc,
                              mm.both3$AICc,
                              mm.both$AICc)))
}

if (region == "Monterey") {

temp1 = t(matrix(temp[colnames(dfa),"Monterey_<=183m"]))
temp2 = t(matrix(temp[colnames(dfa),"Monterey_184-550m"]))
temp3 = t(matrix(temp[colnames(dfa),"Monterey_>=550m"]))
fishing. = t(matrix(fishing[colnames(dfa),"fishing"]))
model.list=list(m=best.model$m, R=best.model$R)
mm.temp1 = MARSS(dfa, model=model.list, z.score=TRUE, form="dfa", control=cntl.list, covariates=temp1)
mm.temp2 = MARSS(dfa, model=model.list, z.score=TRUE, form="dfa", control=cntl.list, covariates=temp2)
mm.temp3 = MARSS(dfa, model=model.list, z.score=TRUE, form="dfa", control=cntl.list, covariates=temp3)
mm.temp = MARSS(dfa, model=model.list, z.score=TRUE, form="dfa", control=cntl.list, covariates=rbind(temp1,temp2,temp3))
mm.fish = MARSS(dfa, model=model.list, z.score=TRUE, form="dfa", control=cntl.list, covariates=fishing.)
mm.both1 = MARSS(dfa, model=model.list, z.score=TRUE, form="dfa", control=cntl.list, covariates=rbind(temp1,fishing.))
mm.both2 = MARSS(dfa, model=model.list, z.score=TRUE, form="dfa", control=cntl.list, covariates=rbind(temp2,fishing.))
mm.both3 = MARSS(dfa, model=model.list, z.score=TRUE, form="dfa", control=cntl.list, covariates=rbind(temp3,fishing.))
mm.both = MARSS(dfa, model=model.list, z.score=TRUE, form="dfa", control=cntl.list, covariates=rbind(temp1,temp2,temp3,fishing.))

# Compare AICc between models with and without covariate   
cov.tbl <- cbind(model=c("no covars",
                         "temp <=183m",
                         "temp 184-550m",
                         "temp >=550m",
                         "all temp",
                         "fishing",
                         "fishing & temp <=183m",
                         "fishing & temp 184-550m",
                         "fishing & temp >=550m",
                         "all covars"),
                 AICc=round(c(best.fit$AICc,
                              mm.temp1$AICc,
                              mm.temp2$AICc,
                              mm.temp3$AICc,
                              mm.temp$AICc,
                              mm.fish$AICc,
                              mm.both1$AICc,
                              mm.both2$AICc,
                              mm.both3$AICc,
                              mm.both$AICc)))

}

if (region == "Conception") {
 
temp1 = t(matrix(temp[colnames(dfa),"Conception_<=183m"]))
temp2 = t(matrix(temp[colnames(dfa),"Conception_184-550m"]))
temp3 = t(matrix(temp[colnames(dfa),"Conception_>=550m"]))
fishing. = t(matrix(fishing[colnames(dfa),"fishing"]))
model.list=list(m=best.model$m, R=best.model$R)
mm.temp1 = MARSS(dfa, model=model.list, z.score=TRUE, form="dfa", control=cntl.list, covariates=temp1)
mm.temp2 = MARSS(dfa, model=model.list, z.score=TRUE, form="dfa", control=cntl.list, covariates=temp2)
mm.temp3 = MARSS(dfa, model=model.list, z.score=TRUE, form="dfa", control=cntl.list, covariates=temp3)
mm.temp = MARSS(dfa, model=model.list, z.score=TRUE, form="dfa", control=cntl.list, covariates=rbind(temp1,temp2,temp3))
mm.fish = MARSS(dfa, model=model.list, z.score=TRUE, form="dfa", control=cntl.list, covariates=fishing.)
mm.both1 = MARSS(dfa, model=model.list, z.score=TRUE, form="dfa", control=cntl.list, covariates=rbind(temp1,fishing.))
mm.both2 = MARSS(dfa, model=model.list, z.score=TRUE, form="dfa", control=cntl.list, covariates=rbind(temp2,fishing.))
mm.both3 = MARSS(dfa, model=model.list, z.score=TRUE, form="dfa", control=cntl.list, covariates=rbind(temp3,fishing.))
mm.both = MARSS(dfa, model=model.list, z.score=TRUE, form="dfa", control=cntl.list, covariates=rbind(temp1,temp2,temp3,fishing.))

# Compare AICc between models with and without covariate   
cov.tbl <- cbind(model=c("no covars",
                         "temp <=183m",
                         "temp 184-550m",
                         "temp >=550m",
                         "all temp",
                         "fishing",
                         "fishing & temp <=183m",
                         "fishing & temp 184-550m",
                         "fishing & temp >=550m",
                         "all covars"),
                 AICc=round(c(best.fit$AICc,
                              mm.temp1$AICc,
                              mm.temp2$AICc,
                              mm.temp3$AICc,
                              mm.temp$AICc,
                              mm.fish$AICc,
                              mm.both1$AICc,
                              mm.both2$AICc,
                              mm.both3$AICc,
                              mm.both$AICc)))
}
  
print(cov.tbl, quote=FALSE)
write.csv(cov.tbl, file=paste("C:/MY_FILES/ANALYSES/Growth analysis - California current/DFA/",region,"_covariates table.csv",sep=''))

} # Region loop
