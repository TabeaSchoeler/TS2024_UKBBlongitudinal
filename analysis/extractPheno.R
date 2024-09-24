#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)
HOME=args[1]
source(paste0(HOME, "/analysis/input.R"))
source(paste0(HOME, "/analysis/functions.R"))
library("psych") 

start_time <- Sys.time()


con <- file(paste0(HOME, "/output/log/extractUKBB_R.log"))
sink(con, append=TRUE, type="output")

print("Remove old file")
file.remove(paste0(HOME,"/output/rds/UKBB.rds"))

print("Read in individual-level data")
redIn=c(677009:677016, 
        676964,
        677171, 677172,
        677173, 677174,
        672734,
        21069,
        6881,
        676527, 676872)

UKBBList=lapply(redIn, function(x) fread(paste0(UKBB, "/pheno/ukb", x, ".csv"))) 
names(UKBBList)=redIn

print("Select variables")
dfAll=readVarNames()
dfExtract=readVarNames(do="extract")


print("Process data")
UKBBL=lapply(dfExtract$label, function(x) extractPheno(var=x, varFile=dfAll, dfL=UKBBList, return="df"))
UKBBdf=Reduce(function(x,y) merge(x = x, y = y, by = "eid", all=T ),  UKBBL)


print("Add polution data")
pollVar=subset(dfAll, ID %in% c(24016, 24017, 24018, 24003))
pullDat=subset(UKBBList[[selectDat(listIn=UKBBList, ID="24016")]], select=c("eid",paste0(pollVar$ID, "-0.0")))
colnames(pullDat)=c("eid", pollVar$label)
# furst
m=as.data.frame(matrix(nrow=NROW(pullDat), ncol=length(pollVar$label)))
mFirst=m
colnames(mFirst)=paste0("first_", pollVar$label)
mFirst$first_pollution_0=0
mFirst$first_pollution_1=1
mFirst$first_pollution_2=2
mFirst$first_pollution_3=3

# date
mDate=m
colnames(mDate)=paste0("date_", pollVar$label)
mDate$date_pollution_0=as.Date("2005-01-01", origin = "1900-01-01")
mDate$date_pollution_1=as.Date("2006-01-01", origin = "1900-01-01")
mDate$date_pollution_2=as.Date("2007-01-01", origin = "1900-01-01")
mDate$date_pollution_3=as.Date("2010-01-01", origin = "1900-01-01")

# assess
mAssess=m
colnames(mAssess)=paste0(pollVar$label, "_assessment")
mAssess[] <- lapply(mAssess, function(x) "pollution")
# combine
dfPoll=cbind(pullDat, mFirst, mAssess, mDate)
head(dfPoll)


UKBBp=merge(UKBBdf, dfPoll, by="eid", all.x=T)

print("Exclude participants withdrawing from the study")
exclude=read.csv(paste0(HOME, "/data/w16389_2023-04-25.csv"), header=F)
UKBBind=subset(UKBBp, !eid %in% exclude$V1)
print("Excluded individuals (consent withdrawn):")
print(NROW(UKBBp)-NROW(UKBBind))

print("Save selected data")
UKBBour=subset(UKBBind, age_0 >= 40 & age_0 <= 69)
saveRDS(UKBBour, paste0(HOME, "/output/rds/UKBB.rds"))


end_time <- Sys.time()

print("======= RUNNING TIME =======")
print(end_time - start_time)
print("============================")

print("All output saved")
uploadLog(file=paste0("extractUKBB_R.log"))
uploadLog(file=paste0("exractUKBB_C.log"))

print("All data saved on cluster")
