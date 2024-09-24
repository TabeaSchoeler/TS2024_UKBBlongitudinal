#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
HOME=args[1]
source(paste0(HOME, "/analysis/input.R"))
source(paste0(HOME, "/analysis/functions.R"))
UKBBdir=UKBB
library(dplyr)

con <- file(paste0(HOME, "/output/log/processPheno_R.log"))
sink(con, append=TRUE, type="output")

print("Read in raw data (selected variables)")
UKBB=readRDS(paste0(HOME, "/output/rds/UKBB.rds"))

print("Remove old files")
file.remove(paste0(HOME,"/output/rds/UKBBfu.rds")) 
file.remove(paste0(HOME,"/output/rds/dateMiss.rds")) 
file.remove(paste0(HOME,"/output/rds/corL.rds")) 
file.remove(paste0(HOME,"/output/rds/corA.rds")) 

print("Select variables")
varAll=readVarNames()
varAll$additionalFU=ifelse(is.na(varAll$additionalFU)==T, "no", varAll$additionalFU)
varAll$additionalFU=ifelse(varAll$additionalFU=="NA", "no", varAll$additionalFU)
varAll$sign=ifelse(is.na(varAll$sign)==T, "no", varAll$sign)
varAll$sign=ifelse(varAll$sign=="NA", "no", varAll$sign)
varInc=subset(varAll, extract=="yes")
long=readVarNames(do="inLong")
baseline=subset(varAll, baseline=="yes")

print(" ========================================================= ")
print(" ================== Data dictionary ====================== ")
print(" ========================================================= ")
checkCoding(ID=varInc$ID)

print(" ========================================================= ")
print(" ===================== Recode variables ================== ")
print(" ========================================================= ")
recodeDF=subset(varAll, baseline=="yes" | inLong=="yes")
recodeL=lapply(recodeDF$label, function(x) recodeVar(df=UKBB, var=x, varInfo=varAll))
UKBBr=Reduce(function(x,y) merge(x = x, y = y, by = "eid", suffixes=c("", "_rem"), all=T ),  recodeL)

print(" ========================================================= ")
print(" ===================== Missing data ====================== ")
print(" ========================================================= ")
print("Summary of data availability and lastest follow up")
declineCor=subset(varAll, dimension %in% c( "cognition") & inLong=="yes")$label
declinePhys=subset(varAll, dimension %in% c( "physical") & inLong=="yes")$label
declineVar=c(declineCor, declinePhys)
dateMissL=lapply(declineVar, function(x) getDesc(df=UKBBr, var=x, return="desc"))
dateMiss=do.call(rbind, dateMissL)
saveOutput(object=dateMiss, label="dateMiss", upload="yes")

print("Get data for correlation matrix between different time points of the same phenotype")
corL=lapply(declineVar, function(x) getDesc(df=UKBBr, var=x, return="cor"))
saveOutput(object=corL, label="corL", upload="yes")

print("Get correlations across all pheno")
corA=lapply(list(declineCor, declinePhys), function(x) corMat(df=UKBBr, var=x))
names(corA)=c("cognition", "physical")
saveOutput(object=corA, label="corA", upload="yes")


print(" ========================================================= ")
print(" ==================== Get follow up data ================= ")
print(" ========================================================= ")
print("get fu data")
UKBBfuL=lapply(long$label, function(x) deriveFU(var=x, df=UKBBr, varDF=varAll))
UKBBfur=Reduce(function(x,y) merge(x = x, y = y, by = "eid" , suffixes=c("", "_rem"), all=T),  UKBBfuL)

print("additionl baseline data")
blVec=baseline$label
blVecS=blVec[! blVec %in% long$label]

UKBBblL=lapply(blVecS, function(x) getBL(var=x, dfFU=UKBBfur, dfBL=UKBBr))
UKBBbl=Reduce(function(x,y) merge(x = x, y = y, by = "eid", suffixes=c("", "_rem") ),  UKBBblL)

UKBBfubl=merge(UKBBbl, UKBBfur, by = "eid", , suffixes=c("", "_rem"))
UKBBfu=merge(subset(UKBB, select=c(eid)),UKBBfubl, by = "eid", all.x=T, suffixes=c("", "_rem"))


print(" ========================================================= ")
print(" ====================== Save output ====================== ")
print(" ========================================================= ")

print("Save file")
saveRDS(UKBBfu , paste0(HOME,"/output/rds/UKBBfu.rds"))

print("All output saved")
uploadLog(file=paste0("processPheno_R.log"))
uploadLog(file=paste0("processPheno_C.log"))

