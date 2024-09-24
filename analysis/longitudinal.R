#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)
HOME=args[1]
source(paste0(HOME, "/analysis/input.R"))
source(paste0(HOME, "/analysis/functions.R"))
library(lm.beta)
library(readr)
library(relaimpo)
library(formattable)
UKBBdir=UKBB
library( tidyverse)

con <- file(paste0(HOME, "/output/log/longitudinal_R.log"))
sink(con, append=TRUE, type="output")

library(gridExtra)
library(ggplot2)
library(cowplot)
library(grid)
colCognition="aquamarine4"
colPhysical="cornflowerblue"

print("Reading variable names")
varAll=readVarNames()
varAll$sign=ifelse(is.na(varAll$sign)==T, "no", varAll$sign)
varAll$sign=ifelse(varAll$sign=="NA", "no", varAll$sign)
fuvar=subset(varAll, inLong=="yes")$label
physicalVar=subset(varAll, inLong=="yes" & dimension %in% c("physical"))$label
cognitiveVar=subset(varAll, inLong=="yes" & dimension %in% c("cognition"))$label
varSUM=c(cognitiveVar, physicalVar)

print("Read in data")
UKBBFU=readRDS( paste0(HOME,"/output/rds/UKBBfu.rds"))
print(paste0("Data created at ", file.info(paste0(HOME,"/output/rds/UKBBfu.rds"))$mtime))


print("Get n with complete longitudinal data")
missFUL=lapply(varSUM, function(x) getMISSfu(df=UKBBFU, var=x))
missFU=do.call(rbind, missFUL)
saveOutput(object=missFU, label="missFU", upload="yes")

print("Mean years for time 0 phenotype since first assessment")
blDiffL=lapply(varSUM, function(x) mean(UKBBFU[[paste0("bldiff_",x,"_0")]], na.rm=T))
blDiff=data.frame(pheno=varSUM, mean=round(do.call(rbind,blDiffL), 1) )

# select phenotypes with substantial correlations across time
corL=readRDS(paste0(HOME,"/output/rds/corL.rds"))
corML=lapply(corL, function(x) withinCor(x))
corMdf=do.call(rbind, corML)
corMn=merge(corMdf, subset(missFU, time==1), by.x="var", by.y="pheno", all=T)

print("Select phenotyps")
minSample=40000
corMs=subset(corMn, corM>0.4 & n>minSample)
varSUM=corMs$var
print("Selected variables")
print(varSUM)

varPhy=physicalVar[physicalVar %in% varSUM]
varCog=cognitiveVar[cognitiveVar %in% varSUM]

print(" ========================================================= ")
print(" ================ GET SUMMARY SCORE ====================== ")
print(" ========================================================= ")
print("Prepare data for change scores")
changeDat=lapply(varSUM, function(x) changeFunc(var=x, data=UKBBFU, return="data", varlist=varAll))
names(changeDat)=varSUM

print("Get summary of follow-up data")
changeSum=lapply(varSUM, function(x) changeFunc(var=x, data=UKBBFU, return="sum", varlist=varAll))
changeSumDF=do.call(rbind, changeSum)
saveOutput(object=changeSumDF, label="changeSum", upload="yes")

modelType=c("changeRes", "changeDiff", "changeLog") 
modelTypeS=c("changeDiff", "changeLog")
physChange=lapply(modelType, function(x) meanSum(df=changeDat, label="physical", select=varPhy, model=x, ID="eid"))
physChangeDF=Reduce(function(x,y) dplyr::full_join(x = x, y = y, by = "eid"), physChange)
cogsChange=lapply(modelType, function(x) meanSum(df=changeDat, label="cognition", select=varCog, model=x, ID="eid"))
cogsChangeDF=Reduce(function(x,y) merge(x = x, y = y, by = "eid", all=T), cogsChange)

print(" ========================================================= ")
print(" ============= SCATTER PLOTS ==================== ")
print(" ========================================================= ")
library(ggpubr)
box_PL=lapply(varPhy, function(x) boxP(var=x, df=changeDat[[x]], varList=varAll, col=colPhysical, sum=changeSumDF)) 
box_P=ggarrange(plotlist=box_PL, nrow=length(varPhy), ncol=1)
ggsave(file=paste0(HOME, "/output/plot/box_PHY.pdf"), plot = box_P, width = 24, height = 27, units = "cm", bg='transparent', limitsize = FALSE)
drop_upload(paste0(HOME, "/output/plot/box_PHY.pdf"), path = paste0(LOCAL, "/results/figures/"), mode = "overwrite")

box_CL=lapply(varCog, function(x) boxP(var=x, df=changeDat[[x]], varList=varAll, col=colCognition, sum=changeSumDF)) 
box_C=ggarrange(plotlist=box_CL, nrow=length(varCog), ncol=1)
ggsave(file=paste0(HOME, "/output/plot/box_COG.pdf"), plot = box_C, width = 24, height = 40, units = "cm", bg='transparent', limitsize = FALSE)
drop_upload(paste0(HOME, "/output/plot/box_COG.pdf"), path = paste0(LOCAL, "/results/figures/"), mode = "overwrite")


print(" ========================================================= ")
print(" ============= Data for BL analyses ====================== ")
print(" ========================================================= ")

print("Scale physical baseline data")
PhyL=lapply(varPhy, function(x) baselineRes(var=x, df=UKBBFU))
blPhy=Reduce(function(x,y) dplyr::full_join(x = x, y = y, by = "eid"), PhyL)
PhyDF=meanBL(df=blPhy, select=varPhy, label="physical", ID="eid")

print("Scale cognitive data")
CogL=lapply(varCog, function(x) baselineRes(var=x, df=UKBBFU))
blCog=Reduce(function(x,y) dplyr::full_join(x = x, y = y, by = "eid"), CogL)
CogDF=meanBL(df=blCog, select=varCog, label="cognition", ID="eid")

print("Combine")
UKBBb=readRDS(paste0(HOME, "/output/rds/UKBB.rds"))
UKBBb$eid=as.character(UKBBb$eid)

datUse=list(subset(UKBBb, select=c("eid" ,"sex_0")),
                           blPhy, blCog, # scaled and residualized (for age) baseline phenotypes
                           physChangeDF,cogsChangeDF, # scaled and residualized (for age) change scores
                           PhyDF, CogDF ) # baseline summary scores
UKBBSUM=Reduce(function(x,y) merge(x = x, y = y, by = "eid" , suffixes=c("", "_rem"), all.x=T),  datUse)

print(" ========================================================= ")
print(" ====================== Correlations ===================== ")
print(" ========================================================= ")
print("Get correlations across all indexes of change")
corChange=lapply(list(c("cognition", varCog), c("physical", varPhy)), function(x) corMat(df=UKBBSUM, var=x, s=c(modelType))) # "0" (if to include baseline pheno)
names(corChange)=c("cognition", "physical")
saveOutput(object=corChange, label="corChange", upload="yes")


print(" ========================================================= ")
print(" ========Prepare phenotype file for regenie ============== ")
print(" ========================================================= ")
sumVar=c("cognition", "physical")
gwaMod=c(modelType, "0", "FU0", "MISS" )
varSUMinc=c(varCog, varPhy)
cogphyGWA=do.call(c, lapply(gwaMod, function(x) paste0(c(sumVar, varSUMinc), "_", x)))

gwaModL=c(modelType, "0", "FU0", "MISS" )
cogphyGWAL=do.call(c, lapply(gwaModL, function(x) paste0(c(sumVar, varSUMinc), "_", x)))

print('Get baseline data for subset of individuals with complete follow-up data')
compFUL=lapply(c(sumVar, varSUMinc), function(x) baselineFU(df=UKBBSUM, var=x, suffix="FU0", run="FU"))
compFU=do.call(cbind, compFUL)

print('Get baseline data for subset of individuals with missing follow-up data')
compMISSL=lapply(c(sumVar, varSUMinc), function(x) baselineFU(df=UKBBSUM, var=x, suffix="MISS", run="MISS"))
compMISS=do.call(cbind, compMISSL)
changeGWAFU=cbind(UKBBSUM, compFU, compMISS)


print("Save REGENIE files")
prepGWAdat(df=changeGWAFU, varIn=cogphyGWAL, saveFile="gwa")

write.table(data.frame(pheno=cogphyGWA),
            file= paste0(HOME, "/output/gwas/batchName"),
            sep="\t",
            row.names = FALSE,
            col.names=F,
            quote=F)


gwaModS=c(modelTypeS, "0")
cogphyGWAS=do.call(c, lapply(gwaModS, function(x) paste0(c(sumVar, varSUMinc), "_", x)))
write.table(data.frame(pheno=cogphyGWAS),
            file= paste0(HOME, "/output/gwas/batchNameS"),
            sep="\t",
            row.names = FALSE,
            col.names=F,
            quote=F)

print("Phenotypes included in GWA:")
print(cogphyGWA)



print("Read in data (unweighted)")
UKBBFU=readRDS( paste0(HOME,"/output/rds/UKBBfu.rds"))
UKBBFU$eid=as.character(UKBBFU$eid)

print(" ========================================================= ")
print(" ========== Phenotypic predictors of change ============== ")
print(" ========================================================= ")

print("=========== GET REGRESSION RESULTS =============")
ageVar=do.call(c, lapply(c( "changeLog" ), function(x) paste0(c("physical","cognition"), "_", x)))
blVar=subset(varAll, baseline=="yes")$label
UKBBFU$eid=as.character(UKBBFU$eid)
UKBBFUC=merge(UKBBFU, subset(changeGWAFU, select=c("eid", ageVar)))

REG_L=list()
for ( i in 1:length(ageVar) ) {
    outIn=ageVar[i]
    outList=lapply(blVar, function(x) runREG(df=UKBBFUC, predictor=x, outcome=outIn))
    REG_L[[i]]=do.call(rbind, outList)
}
REG=do.call(rbind, REG_L)


print(" ========================================================= ")
print(" ======================== Age effects ==================== ")
print(" ========================================================= ")
print("Add UKBB participation weights")
weightsPa=readRDS(paste0(HOME,"/data/datHSEUKBB.rds"))
weightsP=subset(weightsPa, sampleName=="UKBB", select=c(eid,  propensity.weight.normalized))
weightsP$eid=as.character(weightsP$eid)
UKBBFUp=merge(UKBBFU, weightsP, by="eid", all.x=T)
UKBBFUp$w=UKBBFUp$propensity.weight.normalized

print("Get age coefficients")
outcome=c(cognitiveVar, physicalVar)
ageEffL=lapply(outcome , function(x) ageModel(var=x, datAge=UKBBFUp, return="age"))
ageEff=do.call(rbind, ageEffL)
ageEff$sample="ukbb"

print("Get means per age bin")
ageMeanL=lapply(outcome, function(x) ageModel(var=x, datAge=UKBBFUp, return="data"))
ageMean=do.call(rbind, ageMeanL)

print(" ========================================================= ")
print(" ====================== Save output ====================== ")
print(" ========================================================= ")

print("Upload output")
saveOutput(object=changeGWAFU, label="changeGWA", upload="no")
saveOutput(object=REG, label="lmPred", upload="yes") # Phenotypic predictors of change
saveOutput(object=ageEff, label="ageEff", upload="yes") # age effects
saveOutput(object=ageMean, label="ageMean", upload="yes") # age effects across bins
saveOutput(object=cogphyGWA, label="cogphyGWA", upload="yes")
saveOutput(object=varSUMinc, label="declineVar", upload="yes")

print("All output saved")
uploadLog(file=paste0("longitudinal_R.log"))
uploadLog(file=paste0("longitudinal_C.log"))




