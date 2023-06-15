rm(list=ls())

library(dplyr);library(mapdata);library(scales);library(marmap);library(minpack.lm)

# Data from Paul Spencer (already filtered)
load("C:/MY_FILES/ANALYSES/Growth analysis - California current/WareHouse.All.Ages.Env_filteredFeb2021.Rdata")
dat<- WareHouse.All.Ages.Env_filtered_Feb2021 # rename for ease
dat$cohort <- dat$year-dat$age_years # add cohort for VB modelling

# Storage of the VB parameter values
vbpar <- array(NA, dim=c(length(min(dat$cohort):max(dat$cohort)),9,7,3),dimnames=list(c(as.character(min(dat$cohort):max(dat$cohort))),
                                                                                      c("linf", "se.linf", "p.linf", "k", "se.k", "p.k", "to", "se.to", "p.to"),
                                                                                      c("Pacific hake","petrale sole","Pacific sanddab","sablefish","lingcod","shortbelly rockfish","darkblotched rockfish"),
                                                                                      c("ECV","Monterey","Conception")))

# Start loop

for(region in c("ECV","Monterey","Conception")) { # region <- "ECV"
  
  for(sp in c("Pacific hake","petrale sole","Pacific sanddab","sablefish","lingcod",
              "shortbelly rockfish","darkblotched rockfish")) { # sp <- "Pacific hake"
    
    for(cohort in c(min(dat$cohort):max(dat$cohort))) { # cohort <- 1977
      
      print(paste(region, sp, cohort, sep=","))
      
      # Select species & cohort
      sdat <- dat[dat$area2==region & dat$common_name==sp & dat$cohort==cohort,]
      
      if(length(sdat[,1])>0) {
        
        # Plot 
        jpeg(paste("C:/MY_FILES/ANALYSES/Growth analysis - California current/VB plots 2021/",region,"/",sp,"/",region , sp, cohort,".jpeg",sep=''), width=5, height=5, units="in", res=250)
        plot(jitter(sdat$age_years, factor = 0.5), jitter(sdat$length_cm, factor = 1), 
             type = "p", pch = ".", cex = 2,
             xlim = c(0,max(na.omit(sdat$age_years))), ylim = c(0, max(na.omit(sdat$length_cm))),
             xlab = "Age (years)", ylab = "Length (cm)",main = paste(region, sp, cohort, sep=","))
        
        # Fit nonlinear least-squares VB
        
        vonbf.minpack.0 <- function(wk.p) {
          wk.linf <- wk.p[1]
          wk.k <- wk.p[2]
          wk.t0 <- wk.p[3]
          len <- wk.linf * (1 - exp(-wk.k * (sdat$age_years - wk.t0))) 
          len.diff <- sdat$length_cm - len # Residuals NOT squared!
          len.diff
        }
        
        linf.0 <- 70
        k.0 <- 0.1
        t0.0 <- -5
        params.0 <- c(linf.0, k.0, t0.0)
        
        mqdt.control.0 <- list(ftol = 0.00001, ptol = 0.00001, gtol = 0, 
                               diag = numeric(), factor = 100, maxfev = 100 * (length(params.0) + 1), 
                               nprint = 1, maxiter = 1000)
        
        try(vb.minpack.0 <- nls.lm(params.0, fn = vonbf.minpack.0, control = mqdt.control.0))
        
        # Skip loop to next cohort if needed
        #if(vb.minpack.0$message!="Conditions for `info = 1' and `info = 2' both hold.") {dev.off()}
        #if(vb.minpack.0$message!="Conditions for `info = 1' and `info = 2' both hold.") {next}
        
        vonbf.plot.0 <- function(age, Linf, K, t0) {
          Linf * (1 - exp(-K * (age - t0))) 
        }
        
        vb.x <- seq(0, max(sdat$age_years), length = 100)
        vb.minpack.y.0 <- vonbf.plot.0(vb.x, vb.minpack.0$par[1], vb.minpack.0$par[2], vb.minpack.0$par[3])
        lines(vb.x, vb.minpack.y.0, lty = 1, lwd = 2, col = 4)
        dev.off()
        
        # Save VB parameters
        try(vbpar[as.character(cohort),"linf",sp,region] <- vb.minpack.0$par[1])
        try(vbpar[as.character(cohort),"p.linf",sp,region] <- summary(vb.minpack.0)$coefficients[1,"Pr(>|t|)"])
        try(vbpar[as.character(cohort),"se.linf",sp,region] <- summary(vb.minpack.0)$coefficient[1,"Std. Error"])
        try(vbpar[as.character(cohort),"k",sp,region] <- vb.minpack.0$par[2])
        try(vbpar[as.character(cohort),"p.k",sp,region] <- summary(vb.minpack.0)$coefficients[2,"Pr(>|t|)"])
        try(vbpar[as.character(cohort),"se.k",sp,region] <- summary(vb.minpack.0)$coefficient[2,"Std. Error"])
        try(vbpar[as.character(cohort),"to",sp,region] <- vb.minpack.0$par[3])
        try(vbpar[as.character(cohort),"p.to",sp,region] <- summary(vb.minpack.0)$coefficients[3,"Pr(>|t|)"])
        try(vbpar[as.character(cohort),"se.to",sp,region] <- summary(vb.minpack.0)$coefficient[3,"Std. Error"])
        
        # Remove model fits before next loop
        rm(vb.minpack.0)
        
      } # Condition  
      
    } # Cohort loop
    
  } # Species loop
  
} # Region loop

save(vbpar, file="C:/MY_FILES/ANALYSES/Growth analysis - California current/vbpar_2021.Rdata")

