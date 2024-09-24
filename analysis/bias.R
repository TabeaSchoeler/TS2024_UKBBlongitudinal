#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)
HOME=args[1]
task=args[2]

source(paste0(HOME, "/analysis/input.R"))
source(paste0(HOME, "/analysis/functions.R"))

con <- file(paste0(HOME, "/output/log/bias_",task,"_R.log"))
sink(con, append=TRUE, type="output")

print("Read in UKBB data harmonized with HSE")
UKBBHSE=readRDS( paste0(HOME,"/data/UKBBHSE.rds"))
UKBBHSE$ID=paste0(UKBBHSE$eid, "_", UKBBHSE$sampleName)


print("Read in UKBB data with available FU data")
changeUKBB=readRDS( paste0(HOME,"/output/rds/changeGWA.rds"))

print("Get number of individuals with available summary scores")
sumN=data.frame(nCOG=NROW(na.omit(subset(changeUKBB, select=c(cognition_changeDiff)))),
           nPHY=NROW(na.omit(subset(changeUKBB, select=c(physical_changeDiff)))))
saveOutput(object=sumN, label="sumN", upload="yes")

changeUKBB$cog=ifelse(is.na(changeUKBB$cognition_changeDiff)==T, 0, 1)
changeUKBB$phy=ifelse(is.na(changeUKBB$physical_changeDiff)==T, 0, 1)
chUKBB=subset(changeUKBB, select=c(eid, cog, phy))
chUKBB$missFU=chUKBB$cog + chUKBB$phy

print("Read in UKBB longitudinal data")
UKBBfu=readRDS(paste0(HOME, "/output/rds/UKBBfu.rds"))
UKBBfu$eid=as.character(UKBBfu$eid)
UKBB=merge(UKBBfu, chUKBB, by="eid", all.x=T)

print("Dichotomize loss to follow up")
UKBB$sample=ifelse(UKBB$missFU==0, 0, NA )
UKBB$sample=ifelse(UKBB$missFU>0, 1, UKBB$sample )
UKBB$ID=paste0(UKBB$eid, "_UKBB")


descFU=list(info="Data generated in bias.R", 
            date=Sys.Date(),
            ntotUKBB=NROW(UKBB), 
            nFU=NROW(UKBB[which(UKBB$sample == 1),]),
            nCog=NROW(UKBB[which(UKBB$cog == 1),]),
            nPhy=NROW(UKBB[which(UKBB$phy == 1),]))

print("Upload descriptives")
saveOutput(object=descFU, label="descFU", upload="yes")

print("Sample with available FU data, compare with HSE")
UKBBHSEFU=subset(UKBBHSE, !ID %in% subset(UKBB, missFU==0)$ID)
print(table(UKBBHSEFU$sample))

print("Sample with available FU data, compare with droput UKBB sample")
UKBBHSEinc=subset(UKBBHSE, sampleName=="UKBB") # baseline reference sample
UKBBlost=subset(UKBBHSEinc, ID %in% subset(UKBB, missFU==0)$ID)
UKBBlost$sample=0
UKBBlost$sampleName="UKBBlost"
UKBBret=subset(UKBBHSEinc, ID %in% subset(UKBB, missFU>0)$ID) # retained samples
UKBBret$sample=1
UKBBret$sampleName="UKBBretained"
UKBBh=rbind(UKBBlost, UKBBret)
UKBBh$weight_individual=1
print(table(UKBBh$sample))



library("readxl")
varList <- read_excel(paste0(HOME, "/data/variableListWeighting.xlsx"))

print("read in variables for weighting")
vars=readRDS(paste0(HOME, "/data/weightVariables.rds"))

if(task=="getweights"){

    print("================= Start participation (baseline versus HSE) LASSO ======================")
    lassoPB=runLasso(data=UKBBHSE, iteration="participation", varInc=vars, variableList=varList)
    print("================= Finnished participation LASSO ======================")
    
    print("================= Start participation (FU versus HSE ) LASSO ======================")
    lassoPBFU=runLasso(data=UKBBHSEFU, iteration="participationFU", varInc=vars, variableList=varList)
    print("================= Finnished participation LASSO ======================")
    
    print("Combine results")
    lassoDataOut=dplyr::full_join(x = lassoPB[["data"]], y = lassoPBFU[["data"]], by = "eid", suffix=c("", "") )
    lassoSumOut=rbind(lassoPB[["lassoSum"]], lassoPBFU[["lassoSum"]])


    print("Upload data")
    saveOutput(object=lassoSumOut, label="lassoSum", upload="yes")
    saveOutput(object=lassoDataOut, label="lassoWeights", upload="yes")
    print("===== All models uploaded - END ========")
}


if(task=="getpredictors"){
    
    print("Get predictors of participation behaviours")

    print("Read in sampling weights")
    weights=readRDS(paste0(HOME, "/output/rds/lassoWeights.rds"))
    weights$eid=as.character(weights$eid)
    weights$ID=paste0(weights$eid, "_UKBB")

    print("Check with original participation weights")
    weightsPa=readRDS(paste0(HOME,"/data/datHSEUKBB.rds"))
    weightsPa$ID=paste0(weightsPa$eid, "_", weightsPa$sampleName)

    UKBBHSEw=merge(weightsPa, weights, by="ID", all=T, suffixes=c("", "_rem")  )
    wCor=cor.test(UKBBHSEw$propensity.weight.normalized, UKBBHSEw$PW_participation)
    print('Correlation between original PW and newly generated')
    print(wCor)

    print("Get estimates for predictors")
    print("================= Start participation (baseline predictors) ======================")
    print('Unweighted analysis')
    predPB=weightedDist(data=UKBBHSEw, iteration="participation", varInc=vars, variableList=varList, w="weight_individual")
    predPB$iteration="unweighted"

    print('Weighted analysis')
    predPBw=weightedDist(data=UKBBHSEw, iteration="participation", varInc=vars, variableList=varList, w="propensity.weight.normalized")
    predPBw$iteration="weighted"
    print("================= Finnished participation (predictors) ======================")
    
    print("================= Start participation (FU predictors) ======================")
    print('Unweighted analysis')
    UKBBHSEwfu=subset(UKBBHSEw, !ID %in% subset(UKBB, missFU==0)$ID)
    predPBFU=weightedDist(data=UKBBHSEwfu, iteration="participationFU", varInc=vars, variableList=varList, w="weight_individual")
    predPBFU$iteration="unweighted"

    print('Weighted analysis')
    UKBBHSEwfu$PW_participationFU=ifelse(UKBBHSEwfu$sampleName=="HSE", UKBBHSEwfu$weight_individual, UKBBHSEwfu$PW_participationFU)
    predPBFUw=weightedDist(data=UKBBHSEwfu, iteration="participationFU", varInc=vars, variableList=varList, w="PW_participationFU")
    predPBFUw$iteration="weighted"
    print("================= Finnished participation (predictors) ======================")

    print("================= Start FU dropout (predictors) ======================")
    predFU=weightedDist(data=UKBBh, iteration="FU", varInc=vars, variableList=varList, w="weight_individual")
    predFU$iteration="unweighted"
    print("================= Finnished reporting error (predictors) ======================")
    
    predOut=rbind(predPB, predPBw, predPBFU, predPBFUw, predFU)
    saveOutput(object=predOut, label="predPB", upload="yes")


    print("Baseline predictors only included in UKBB")
    varAll=readVarNames()
    blvar=subset(varAll, participation=="yes")$label

    blPredL=lapply(paste0(blvar, "_0"), function(x) glmUni(x, outcome="sample", df=UKBB, weights="no", weightInc="", fam=gaussian))
    blPredc=cleanString(df=do.call(rbind, blPredL), label=blvar)
    blPred=extractRegression(model=blPredc, label=blvar, n= blPredc$n)
    saveOutput(object=blPred, label="baselineP", upload="yes")
}

if(task=="applyweights"){
    print("Reading variable names")
    varSUM=readRDS( paste0(HOME,"/output/rds/declineVar.rds"))

    expV=c("sex_0", "education_length_0", "SBP_0", "pollution_0")
    outV=c(paste0(c("physical","cognition", varSUM), "_changeDiff"),
           paste0(c("physical","cognition", varSUM), "_changeLog"),
           paste0(c("physical","cognition", varSUM), "_changeRes"))

    print("Read in weigts")
    weights=readRDS(paste0(HOME, "/output/rds/lassoWeights.rds"))
    weights$eid=as.character(weights$eid)

    print("Merge UKBB data with weights")
    UKBW=merge(changeUKBB, weights, by="eid", all.x=T)
    UKBW$PW_none=1

    print("Add baseline phenotypes")
    UKBBfu=readRDS(paste0(HOME, "/output/rds/UKBBfu.rds"))
    blPheno=subset(UKBBfu, select=c("eid", expV))
    blPheno$eid=as.character(UKBW$eid)
    UKBWb=merge(UKBW, blPheno, by="eid", all.x=T, suffixes=c("" ,"_rem"))

    print("Get weighted correlations")
    library(lm.beta)
    library(survey)
    UKBWbC <-  UKBWb[complete.cases(UKBWb$PW_participationFU),]
    UKBBcorB=getCorrs(data=UKBWbC, weights="PW_participation", out=outV, exp=expV, model="baselineW")
    UKBBcorFU=getCorrs(data=UKBWbC, weights="PW_participationFU", out=outV, exp=expV, model="fuW")
    UKBBcorNone=getCorrs(data=UKBWbC, weights="PW_none", out=outV, exp=expV, model="none")
    wCor=rbind(UKBBcorB, UKBBcorFU, UKBBcorNone)
    saveOutput(object=wCor, label="wCor", upload="yes")
}

print("All output saved")
uploadLog(file=paste0("bias_",task,"_R.log"))
uploadLog(file=paste0("bias_",task,"_C.log"))
