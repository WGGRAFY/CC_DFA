rm(list=ls())

library(dplyr)
library(mapdata)


# Data from Christine (UPDATED BY PAUL)
load("C:/MY_FILES/ANALYSES/Growth analysis - California current/WareHouse.All.Ages.Env_filteredFeb2021.Rdata")
ldat <- WareHouse.All.Ages.Env_filtered_Feb2021
ldat$region <- ldat$area2

# Load data to get which age to investigate
A50 <- read.csv("C:/MY_FILES/ANALYSES/Growth analysis - California current/A50.csv",row.names = 1)
maxage <- read.csv("C:/MY_FILES/ANALYSES/Growth analysis - California current/maxage.csv",row.names = 1)


## LOOP TO SAVE MEAN LENGTH-AT-AGE

# Get mean length at age 
mlat <- data.frame(ldat %>% group_by(region,common_name,year,age_years) %>% summarise(ML = mean(length_cm, na.rm=T)))

# Pick 3 age classes: juveniles, A50, and matures
age_juv <- maxage; age_juv[,] <- NA
age_A50 <-  maxage; age_A50[,] <- NA
age_mat <- maxage; age_mat[,] <- NA

for (aa in colnames(age_juv)) { # aa <- "ECV"
  for (sp in rownames(age_juv)) { # sp <- "darkblotched rockfish"
    age_juv[sp,aa] <- round(A50[sp,"A50_average"]/2) # half of A50
    age_A50[sp,aa] <- round(A50[sp,"A50_average"]) # A50
    if (maxage[sp,aa] - A50[sp,"A50_average"]>0) {
      age_mat[sp,aa] <-  round(A50[sp,"A50_average"] + ((maxage[sp,aa] - A50[sp,"A50_average"])/2)) # Half way between A50 and max age observed in region
    }
  }
}
rm(aa,sp)

age_ <- age_juv; save(age_, file = "C:/MY_FILES/ANALYSES/Growth analysis - California current/age_juv.Rdata")
age_ <- age_A50; save(age_, file = "C:/MY_FILES/ANALYSES/Growth analysis - California current/age_A50.Rdata")
age_ <- age_mat; save(age_, file = "C:/MY_FILES/ANALYSES/Growth analysis - California current/age_mat.Rdata")

# Collate and save data for DFA

for (aa in unique(mlat$region)) { # aa <- "ECV"
  for (age in c("age_juv","age_A50","age_mat")) { # age <- "age_juv"
    
    load(paste("C:/MY_FILES/ANALYSES/Growth analysis - California current/",age,".Rdata",sep=''))
    
    ml <- data.frame(array(NA,dim=c(length(min(mlat$year):max(mlat$year)),length(unique(mlat$common_name)))))
    colnames(ml) <- unique(mlat$common_name)
    rownames(ml) <- c(min(mlat$year):max(mlat$year))
    
    for (sp in unique(mlat$common_name)) { # sp <- "darkblotched rockfish"
      
      if(!is.na(age_[sp,aa])) {
      
        for (yy in mlat[mlat$region==aa & mlat$common_name==sp & mlat$age_years==age_[sp,aa],"year"]) { # yy <- 2003
        
          ml[as.character(yy),sp] <- (mlat[mlat$region==aa & mlat$common_name==sp & mlat$age_years==age_[sp,aa] & mlat$year==yy,"ML"] - mean(mlat[mlat$region==aa & mlat$common_name==sp & mlat$age_years==age_[sp,aa],"ML"])) / sd(mlat[mlat$region==aa & mlat$common_name==sp & mlat$age_years==age_[sp,aa],"ML"])

        }
      }
    } 
    save(ml,file=paste("C:/MY_FILES/ANALYSES/Growth analysis - California current/mean-length-at-age/data/ml_",aa,"_",age,".Rdata",sep=''))
  }
}

