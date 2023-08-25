# --------------------
# Clean light data, write to file, extract light summary, light bins 
# D.W. 23/08/23
# --------------------
# Setup --------------------
rm(list = ls())
library(readr)
library(Rcpp)
library(imputeTS)
library(parsedate)

# Functions ----------
# -> Rolling window indices --------
rollingWindowInd <- function(t=c(),
                             window=c(),
                             step=c()
){
  if (length(t)==0) {
    stop("Error: Specify t")
  }
  if (length(window)==0) {
    stop("Error: Specify window")
  }
  if (length(step)==0) {
    stop("Error: Specify step")
  }
  
  tLen <- t[length(t)]- t[1] # Length of time series
  
  if (window >= tLen) {
    stop("Error: Window length must not be longer than time series length")
  }
  if (step > window) {
    stop("Error: Step length must not be longer than window length")
  }
  
  tAbs <- t - t[1] # Time vector starting at zero
  nWin <- ceiling((tLen-window)/step+1) # Number of windows
  
  stInd <- rep(NA, nWin); enInd <- rep(NA, nWin);
  for (jj in 1:nWin) {
    
    tSt <- step*(jj-1)
    tEn <- tSt + window
    stInd[jj] <- which.min(abs(tAbs-tSt))
    
    if (jj == nWin){
      enInd[jj] <- length(t)
    } else {
      enInd[jj] <- which.min(abs(tAbs-tEn))
    }
    
  }
  
  if (stInd[length(stInd)] == enInd[length(enInd)]) {
    stInd <- stInd[-length(stInd)]
    enInd <- enInd[-length(enInd)]
  }
  
  naBl <- (t[enInd] - t[stInd]) < .8*window | (t[enInd] - t[stInd]) > 1.2*window # For cases where windows are too long / short
  stInd[naBl] <- NA; enInd[naBl] <- NA
  
  returnLi <- list(stInd, enInd)
  return(returnLi)
}

# -> Longest interval indices ---------------
# Find indicies of longest continuous string of 1s in a binary vector
longestInterval <- function(binaryIn=c()
){
  if (sum(binaryIn == 1) == length(binaryIn)){
    # If the vector is all 1s
    maxStIn <- 1
    maxEnIn <- length(binaryIn)
    
    returnLi <- list(maxStIn, maxEnIn)
    return(returnLi)
    
  } else {
    oneSt <- which(diff(binaryIn) == 1) + 1
    if (length(oneSt)>0){
      oneDf <- data.frame(ind = oneSt, trans=1)
    } else {
      oneDf <- data.frame()
    }
    
    zSt <- which(diff(binaryIn) == -1) + 1
    if (length(zSt)>0){
      zDf <- data.frame(ind = zSt, trans=0)
    } else {
      zDf <- data.frame()
    }
    
    trV <- rbind(oneDf,zDf)
    trV <- trV[order(trV$ind),]
    if (trV$ind[1]>1){
      top <- data.frame(ind=1,trans=abs(trV$trans[1]-1))
    } else {
      top <- data.frame()
    }
    bottom <- data.frame(ind=(length(binaryIn)+1),trans=NA)
    
    trV2 <- rbind(top,trV,bottom)
    diffIn <- diff(trV2$ind)
    diffIn[trV2$trans[-length(trV2$trans)] == 0] <- 0
    
    maxStIn <- trV2$ind[which.max(diffIn)]
    maxEnIn <- trV2$ind[which(trV2$ind == maxStIn) + 1] - 1
    
    returnLi <- list(maxStIn, maxEnIn)
    return(returnLi)
  }
}


# --> Split linear/quadratic function ------------
quadlin_fun <- function(x){
  y <- rep(NA, length(x))
  y[x <= 2.4 & !is.na(x)] <- (21.55*x - 4.54*x^2)[x <= 2.4 & !is.na(x)]
  y[x > 2.4 & !is.na(x)] <- (10.71*x)[x > 2.4 & !is.na(x)]
  return(y)
}

# --------------------
# Options ---------
folder <- 247 # Specify folder if using local system

readdir <- paste("C:/Users/dan_t/Documents/R/Biobank/Light_Data/",folder,"_output",sep="") # Specify directory if using local system
run_n <- 31
lgtsmryname <- paste0("lgtsummary", run_n, ".csv")
lgtbinname <- paste0("lgtbins", run_n, ".csv")
paramsname <- paste0("params", run_n, ".csv")

SRIname <- "SRI_4.csv"
SWVname <- "SWV_4.csv"
sleepvarname <- "SleepVar_3.csv"

excl_nonWearSWV <- TRUE # T = write NA to light data for all NA epochs in SWV (non-wear defined by GGIR + our own implementation of GGIR non-wear)

smooth_lgt <- TRUE # T = smooth light data based on a rolling-window average
zero_correct <- TRUE # T = correct all convlux values such that smallest value == 0
lux_transform <- TRUE # T = transform AX3 data to approximated 'real' lux, determined by device testing/luxmeter ### SET FALSE IF USING IMPUTATION
transform_type <- "quad_linear" # use "quad_linear" or "interpolation" 

write_lightsummary <- TRUE # T = write summary of light cleaning steps per ppt
write_lightbins <- TRUE # T = to bin light data and write to file

write_params <- TRUE # Write a list of parameters used for this run

SWVtype <- "mclonoff_WASOnap"

# --------------------
# Create light summary file to write to --------
if (write_lightsummary){
  lgtsmryfile <- paste0(readdir,"/",lgtsmryname)
  lgtsmryheader <- c("loc", "id", "slp.excl", "SRI.excl", "cont.excl", "rpt.excl", "val.daysGGIR", "val.daysSWV", "pcLgtInSWV", "dark.thres",
                     "lgtDaysTot", "naPcNW", "naPcNWcov", "impDays", "impWritten", "impNApc")
  write.table(t(lgtsmryheader),lgtsmryfile,sep=",", col.names=FALSE, row.names=FALSE) #--------------------------------------------------------<<
}

# Create light bin file to write to --------
if (write_lightbins){
  lgtbinfile <- paste(readdir,lgtbinname,sep="/")
  lgtbinheader <- c("ID", "binOn", "binOff", "convluxM", "convluxMed", "convluxSD", "Ndata")
  lgtbinNaVec <- rep(NA ,length(lgtbinheader))
  write.table(t(lgtbinheader),lgtbinfile,sep=",", col.names=FALSE, row.names=FALSE) #--------------------------------------------------------<<
}

# --------------------
# Parameters --------
invfr <- 10 # seconds between each light data sample 
smWin <- 600; smStep <- 120 # Window size and step size (seconds) for smoothing
lightBinHrs <- 0.5 # Size of light data bins 

# --------------------
# Create d.f. summarizing details of light cleaning per file ---------
loc <- list.files(readdir, recursive=T, pattern = "light") # Find all files that contain "light" (equivalent to file_list)
id <- substr(loc,(nchar(loc)-25), (nchar(loc)-9)) # (equivalent to file_list2)
lgtsmry <- data.frame(loc = loc, id = id, 
                      slp.excl = F, 
                      SRI.excl = F, 
                      cont.excl = F, 
                      rpt.excl = F, 
                      val.daysGGIR = NA, val.daysSWV = NA, pcLgtInSWV = NA, dark.thres = NA,
                      lgtDaysTot = NA, naPcNW = NA, naPcNWcov = NA,
                      impDays = NA, impWritten = FALSE, impNApc = NA)

# --------------------
# Read in light data from file, fill in gaps with 10sec increment time-stamps and light=NA ----------
lightDatalist <- list() # Initialize empty list
for (i in 1:nrow(lgtsmry)){
  # -------------
  # Working data -----------
  pptID <- lgtsmry$id[i] # Extract ppt ID
  pptdata <- read.csv(paste0(readdir,"/",lgtsmry$loc[i])) # Read light data for this ppt
  # Check for missing light data, fill in blank times and lux values with NA ----------
  pptdata$time <- floor(pptdata$time/10)*10 # Round all time values to the nearest 10 - won't work for invfr != 10
  undif <- unique(diff(pptdata$time))
  
  negs <- undif[undif < 0] # Remove any negative time difference values
  if (length(negs) > 0) {
    for (ng in 1:length(negs)){
      pptdata$time[which(diff(pptdata$time) == negs[ng])+1] <- NA
    }
  }
  pptdata <- pptdata[!is.na(pptdata$time),]
  
  undif <- unique(diff(pptdata$time))
  undif2 <- undif
  
  j=1
  while (length(undif) > 1){ # If there are differences between timestamps  > 10sec
    if(undif2[j] != invfr){
      mi.ind <- which(diff(pptdata$time) %in% undif2[j])
      fdf.li <- list()
      for (k in 1:length(mi.ind)){
        fst <- pptdata$time[mi.ind[k]] + invfr
        fen <- pptdata$time[(mi.ind[k]+1)] # - invfr
        fv <- seq(from=fst,to=fen,by=invfr)
        fv <- fv[-length(fv)]
        fdf <- data.frame(time = fv, lux = NA)
        fdf.li[[k]] <- fdf
      }
      for (k in 1:length(fdf.li)){
        pptdata <- rbind(pptdata[1:mi.ind[k],], fdf.li[[k]], pptdata[(mi.ind[k]+1):nrow(pptdata),])
        mi.ind <- mi.ind + nrow(fdf.li[[k]])
      }
      undif <- unique(diff(pptdata$time))
    }
    j <- j+1
    if (j > 1000){
      break
    }
  }
  # Convert lux, write to list ----------
  convlux <- 10^(pptdata$lux/341) # Extract converted lux value
  abstime <- pptdata$time - pptdata$time[1]
  pptIDdata <- cbind(pptdata,convlux,abstime) # Combine into one d.f.
  lightDatalist[[i]] <- pptIDdata # Assign ppt d.f. to list
  names(lightDatalist)[i] <- pptID # Re-name list segment as ppt ID
  # -------------
}
# --------------------
# Record which files have GGIR sleep data -----------
nightsumlist <- paste(readdir,"/output_",lgtsmry$id,".csv/results/QC/part4_nightsummary_sleep_full.csv",sep="") # List of nightsummary file locations
lgtsmry$slp.excl <- !file.exists(nightsumlist)

# Read in GGIR sleep data ----------
pptSleepLi <- list()
nsl_exists <- nightsumlist[!lgtsmry$slp.excl]
nsl_exists_pt <- lgtsmry$id[!lgtsmry$slp.excl]

for (i in 1:length(nsl_exists)){
  pptSleepLi[[i]] <- read.csv(nsl_exists[i])
  names(pptSleepLi)[i] <- nsl_exists_pt[i]
}

# ----------------------------------------------------------------------------- CHANGE 'SRI_3' DEPENDING UPON THE RELEVANT/RECENT FILE (E.G., SRI3)
# --------------------
# Record which files have SRI scores (>4 days of sleep data) ----------
SRIdf <- read.csv(paste0(readdir,"/",SRIname),header=FALSE) # Read in SRI .csv data
SRIdf_narm <- SRIdf[!is.na(SRIdf[,3]),][-1,]; names(SRIdf_narm) <- SRIdf[1,] # Clean / rearrange it
pptlist_SRI <- paste(SRIdf_narm$participant, SRIdf_narm$tag, sep="_")  # Get list of all participants with an SRI score (> 4 days sleep data)
lgtsmry$SRI.excl <- !(lgtsmry$id %in% pptlist_SRI)

# --------------------
# Extract number of valid days of data (based on GGIR output), list of valid days per ppt, list of >=5 of activity data, midsleep time ----------
incl <- rep(TRUE, length(nightsumlist)) # Create an empty vector for TRUE/FALSE values
inclnights.cons <- vector(mode = "list", length = length(nightsumlist))
inclnights.all <- inclnights.cons
nightcount <- rep(NA,length(nightsumlist))
midsleepdf <- data.frame(ID = lgtsmry$id, avmidsleep = NA)
for (i in 1:length(nightsumlist)){
  if (file.exists(nightsumlist[i])){
    # --------------------
    # Read -------------
    tempdf <- read.csv(nightsumlist[i])
    # Count number of valid nights per participant, calculate midsleep time, extract valid nights to list ----------
    midsleepdf$avmidsleep[i] <- mean((tempdf$wakeup + tempdf$sleeponset)/2) # Average midsleep times (for CBTmin initial condition, used later)
    nightcount[i] <- nrow(tempdf) # Number of nights of activity / sleep data
    inclnights.all[[i]] <- tempdf$night
    
    # Output list of 5 or more consecutive days of list data -----------------
    if (nrow(tempdf) == 5){
      inclnights.cons[[i]] <- tempdf$night
      cons <- sum(diff(tempdf$night)) + 1
      if (cons != nrow(tempdf)) {
        incl[i] <- FALSE
      }
    } else if (nrow(tempdf) == 6){
      inclnights.cons[[i]] <- tempdf$night
      cons <- sum(diff(tempdf$night)) + 1
      if (cons != nrow(tempdf)) {
        dbldiff <- sum(diff(diff(tempdf$night)))
        if (dbldiff == -1) {
          inclnights.cons[[i]] <- c(3,4,5,6,7)
          nightcount[i] <- 5
        } else if (dbldiff == 1) {
          inclnights.cons[[i]] <- c(1,2,3,4,5)
          nightcount[i] <- 5
        } else {
          incl[i] <- FALSE
        }
      }
    }
    # --------------------
  }
}
midsleepdf <- midsleepdf[!is.na(midsleepdf$avmidsleep),]
midsleepdf$avmidsleep[midsleepdf$avmidsleep >= 24] <- midsleepdf$avmidsleep[midsleepdf$avmidsleep >= 24] - 24 # Reduce midsleep values to 0-24 range

lgtsmry$cont.excl <- !incl
lgtsmry$val.daysGGIR <- nightcount

# --------------------
# Record which participants have repeat measures --------
lgtsmry$rpt.excl <- as.numeric(substr(lgtsmry$id,15,15)) != 0

# --------------------
# Read in SWVs ----------
f3col <- as.data.frame(read_csv(paste0(readdir,"/",SWVname))) # Read in first 3 columns only
scol <- strsplit(f3col$SRItype,",")
for (i in 1:length(scol)){
  f3col$SRItype[i] <- scol[[i]][1]
}

reqind <- (1:nrow(f3col))[f3col$SRItype == SWVtype] # Extract indices of all participant SWVs 
reqind <- reqind[!is.na(reqind)]
pptSWVli <- list() # Initiate empty list
for (i in 1:length(reqind)){
  datr <- as.character(read.csv(paste0(readdir,"/",SWVname), nrows=1, skip=reqind[i], header=FALSE)) # Read row
  datrSO <- as.numeric(datr[3:length(datr)]) # Sleep/wake data only
  pptSWVli[[i]] <- data.frame(trans = datrSO[1:(length(datrSO)/2)], 
                              t = datrSO[(length(datrSO)/2+1):length(datrSO)]) # Add transition / timestamp data to list
}
names(pptSWVli) <- substr(f3col$file,1,17)[reqind] # Names for the list

# Account for daylight-saving time in all SWVs ----------
daylon <- c(954205200, 985741200, 1017277200, 1048813200,  1080435600, 1111971600, 1143507600, 1175043600, 1206666000, 
            1238202000, 1269738000, 1301274000, 1332896400, 1364432400, 1395968400, 1427504400, 1459126800, 1490662800, 
            1522198800, 1553734800, 1585357200, 1616893200, 1648429200, 1679965200, 1711587600, 1743123600) # Define daylight savings on-off times in UNIX
dayloff <- c(972957600, 1004493600, 1036029600, 1067565600, 1099188000, 1130724000, 1162260000, 1193796000, 1225418400,
             1256954400, 1288490400, 1320026400, 1351648800, 1383184800, 1414720800, 1446256800, 1477879200, 1509415200,
             1540951200, 1572487200, 1604109600, 1635645600, 1667181600, 1698717600, 1730340000, 1761876000)

for (i in 1:length(pptSWVli)){
  tist <- pptSWVli[[i]]$t[1]
  if (is.na(tist)){
    tist <- pptSWVli[[i]]$t[2]
  }
  tien <- pptSWVli[[i]]$t[nrow(pptSWVli[[i]])]
  if(sum(tist > daylon & tist < dayloff) == 1 & sum(tien > daylon & tien < dayloff) == 1){
    pptSWVli[[i]]$t <- pptSWVli[[i]]$t - 60*60
  } # If both start and end times of recording are within daylight savings time, participants SWV values back 1hr
}

# --------------------
for (i in 1:length(lightDatalist)){
  tryCatch({ # Start of code to catch errors in each iteration, write error, and skip to next iteration
    # ---------------------------
    # Extract ppt data ----------
    pptID <- names(lightDatalist[i])
    pptLD <- lightDatalist[[i]]
    lgtstp <- pptLD$time[2] - pptLD$time[1] # Length of time between light data points (sec)
    pptSWV.in <- which(names(pptSWVli) == pptID) # Index of ppt SWV within ppt SWV list
    
    # ---------------------------
    # Re-zero participant data ----------
    if (zero_correct){
      pptLD$convlux <- pptLD$convlux - min(pptLD$convlux, na.rm = TRUE)
    }
    
    # ---------------------------
    # Extract total length of participant data  -----------
    lgtsmry$lgtDaysTot[lgtsmry$id == pptID] <- round(pptLD$abstime[length(pptLD$abstime)]/60/60/24, 3)
    
    # ---------------------------
    if (length(pptSWV.in) == 1) { # If there is a sleep-wake vector for this participant
      # ----------------------------
      # Re-create binary SWV -------
      pptSWV.sum <- pptSWVli[[pptSWV.in]] # Participant sleep-wake vector summary d.f.
      ppttim <- seq(pptSWV.sum$t[1], pptSWV.sum$t[nrow(pptSWV.sum)], lgtstp) # Times

      pptSWV <- rep(NA, length(ppttim)) # Create empty binary vector
      for (j in 1:(nrow(pptSWV.sum)-1)){
        st <- which(pptSWV.sum$t[j] == ppttim)
        en <- which(pptSWV.sum$t[j+1]-lgtstp == ppttim)
        pptSWV[st:en] <- pptSWV.sum$trans[j]
      }
      
      # ----------------------------
      # Get original length of light and SWV vectors -----------
      pptSWV_l <- (ppttim[length(ppttim)] - ppttim[1])/60/60/24
      pptLD_l <- (pptLD$time[length(pptLD$time)] - pptLD$time[1])/60/60/24
      
      # ----------------------------
      # Match time vectors for light and SWV data -------
      ld1 <- as.double(substr(ppttim[1],nchar(ppttim[1]),nchar(ppttim[1]))) # last digit of SWV time vector, same for light data time vector
      ld2 <- as.double(substr(pptLD$time[1],nchar(pptLD$time[1]),nchar(pptLD$time[1]))) 
      ppttim <- ppttim+ld2-ld1 # Adjust SWV time vector by a few seconds to match light data time vector
      
      pptLD <- pptLD[pptLD$time %in% ppttim,] # Match time vectors for light and SWV 
      if (nrow(pptLD) == 0){ # If no overlapping (for some reason time vectors don't overlap...?)
        print(paste0("Error: Timestamps for light data and sleep-wake vector do not overlap, i = ", i))
        next
      }
      pptLD$abstime <- pptLD$abstime - min(pptLD$abstime) # Re-zero abstime
      pptSWV <- pptSWV[ppttim %in% pptLD$time]
      ppttim <- ppttim[ppttim %in% pptLD$time]
      
      # ----------------------------
      # Get percentage light data in SWV -----------
      pptLD_l2 <- (pptLD$time[length(pptLD$time)] - pptLD$time[1])/60/60/24
      lgtsmry$pcLgtInSWV[lgtsmry$id == pptID] <- round(pptLD_l2/pptLD_l, 3)
      
      # ----------------------------
      # Remove all light epochs matching with NAs in SWVs (non-wear by GGIR and custom method) ----------
      if (excl_nonWearSWV){
        pptLD$convlux[is.na(pptSWV)] <- NA # Call all light periods overlapping with NA SWV periods as NA
      }
      
      # ----------------------------
      # Get NA percentage of light data after non-wear exclusion ----------
      lgtsmry$naPcNW[lgtsmry$id == pptID] <- round(sum(is.na(pptLD$convlux))/length(pptLD$convlux),3)
      
      # ----------------------------
      # Get NA percentage of light data after non-wear exclusion and coverage exclusion ----------
      lgtsmry$naPcNWcov[lgtsmry$id == pptID] <- round(sum(is.na(pptLD$convlux))/length(pptLD$convlux),3)
      
      # ----------------------------
      # Extract number of valid days of data, based on SWV ---------
      if (nrow(pptSWV.sum) >= 2){
        lgtsmry$val.daysSWV[lgtsmry$id == pptID] <- sum(diff(pptSWV.sum$t)[!is.na(pptSWV.sum$trans[-nrow(pptSWV.sum)])])/60/60/24
      }
      
      # ----------------------------
    } else {
      pptSWV <- NA
    }
    # ---------------------------
    # Smooth ppt data ----------
    if (smooth_lgt){
      winOnOff <- rollingWindowInd(t = pptLD$time, window = smWin, step = smStep)
      smPptLD <- data.frame(time = rep(NA,length(winOnOff[[1]])), convlux=NA, abstime=NA)
      for (j in 1:nrow(smPptLD)){
        smPptLD$convlux[j] <- median(pptLD$convlux[winOnOff[[1]][j]:winOnOff[[2]][j]], na.rm=TRUE)
        smPptLD$time[j] <- mean(pptLD$time[winOnOff[[1]][j]:winOnOff[[2]][j]], na.rm=TRUE)
      }
      smPptLD$abstime <- smPptLD$time - smPptLD$time[1]
      pptLD <- smPptLD
      lgtstp <- pptLD$time[2] - pptLD$time[1]
    }
    
    # ---------------------------
    # Convert light data using AX3 to luxmeter conversion function -----------
    if (lux_transform){
      wd <- quadlin_fun(pptLD$convlux)
      pptLD$convlux <- wd
    }

    # ---------------------------
    # Extract light metrics
    # -> Convert UNIX times to clock time and numeric (0-24) ---------
    cltim.st <- substr(as.POSIXct(pptLD$time, origin="1970-01-01",tz="UTC"),12,19) # Convert UNIX to clock time
    cltim.nu <- as.numeric(substr(cltim.st,1,2)) + as.numeric(substr(cltim.st,4,5))/60 + as.numeric(substr(cltim.st,7,8))/3600 # Numeric 24h
    
    # -> Define a 'day' as data start + 24 ---------
    st24 <- seq(from=0, by=24*60*60, length.out=ceiling(max(pptLD$abstime)/(24*60*60))) # Generate a vector of start of each 24h interval
    st.v <- rep(NA,length(st24))
    en.v <- rep(NA,length(st24))
    for (j in 1:length(st24)){ # For each 24hr interval in the dataset (i.e., 0-24,24-48,48-72), find the start index of the interval (hrdim length) with lowest light lvls
      st.v[j] <- which(pptLD$abstime == st24[j])
      en <- which(pptLD$abstime == (st24[j]+24*60*60-10))
      if (length(en) == 0){
        en.v[j] <- length(pptLD$abstime)
      } else {
        en.v[j] <- en
      }
    }
    # ---------
    # Extract light data as bins ----------
    if (write_lightbins){
      cuts <- seq(from=0, to=24, by=lightBinHrs) # Define cuts
      if (cuts[length(cuts)] != 24){
        cuts[length(cuts)+1] <- 24
      }
      
      cuts1 <- cuts[-length(cuts)]
      cuts2 <- cuts[-1]
      
      lgtBinD <- data.frame(ID = rep(pptID, length(cuts1)), 
                            binOn=cuts1, binOff=cuts2, 
                            convluxM=NA, convluxMed=NA, convluxSD=NA,
                            Ndata=NA) # Empty d.f.
      for (j in 1:length(cuts1)){
        binBi <- (cltim.nu >= cuts1[j]) & (cltim.nu < cuts2[j])
        if (j == length(cuts1)){
          binBi <- (cltim.nu >= cuts1[j]) & (cltim.nu <= cuts2[j])
        }
        lgtBinD$convluxM[j] <- round(mean(pptLD$convlux[binBi], na.rm=TRUE),3)
        lgtBinD$convluxMed[j] <- round(median(pptLD$convlux[binBi], na.rm=TRUE),3)
        lgtBinD$convluxSD[j] <- round(sd(pptLD$convlux[binBi], na.rm=TRUE),3)
        lgtBinD$Ndata[j] <- sum(!is.na(pptLD$convlux[binBi]))
      }
    }
    
    # ---------------------------
    # Write light bins to file ----------
    if (write_lightbins){
      write.table(lgtBinD, lgtbinfile, sep = ",", col.names = !file.exists(lgtbinfile), append = TRUE, row.names=FALSE)
    }
    # ---------------------------
  }, error = function(e) {
    print(e)
    print(i)
    
    if (write_lightbins){
      writebinvec <- c(pptID, lgtbinNaVec[-1])
      write.table(t(writebinvec), lgtbinfile, sep = ",", col.names = !file.exists(lgtbinfile), append = TRUE, row.names=FALSE)
    }
    
  }, finally = {
    next # Skip to next loop iteration 
  }) # Other half of the error-catch function
}

# --------------------
# Write light summary data --------
if (write_lightsummary){
  write.table(lgtsmry, lgtsmryfile, sep = ",", col.names = !file.exists(lgtsmryfile), append = TRUE, row.names=FALSE)
}
  
# --------------------
# Write parameters ----------
if (write_params){
  name <- c("run_n",
             "lgtsmryname",
             "lgtbinname",
             "SRIname",
             "SWVname",
             "excl_nonWearSWV",
             "smooth_lgt",
             "zero_correct",
             "lux_transform", 
             "write_lightsummary",
             "write_lightbins",
             "SWVtype",
             "invfr",
             "smWin",
             "lightBinHrs")
  value <- c(run_n,
             lgtsmryname,
             lgtbinname,
             SRIname,
             SWVname,
             excl_nonWearSWV,
             smooth_lgt,
             zero_correct,
             lux_transform,
             write_lightsummary,
             write_lightbins,
             SWVtype,
             invfr,
             smWin,
             lightBinHrs
             )
  params <- data.frame(name=name, value=value)
  paramsfile <- paste0(readdir,"/", paramsname)
  write.csv(params,paramsfile,row.names=FALSE)
}

# --------------------
