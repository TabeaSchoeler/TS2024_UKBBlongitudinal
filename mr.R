args <- commandArgs(trailingOnly = TRUE)
HOME=args[1]
outcome=args[2]
source(paste0(HOME, "/analysis/input.R"))
source(paste0(HOME, "/analysis/functions.R"))

con <- file(paste0(HOME, "/output/log/mr_", outcome,"_R.log"))
sink(con, append=TRUE, type="output")

print(paste0(" =========== Start MR for ", outcome, " ==============="))

gwaU=read.table(paste0(HOME, "/output/gwas/batchNameS"), head = F)$V1


if(outcome!="process"){
library(TwoSampleMR)
system(paste0("chmod 755 ",HOME, "/programs/plink"))

print("Select gwa sumstats")
gwaListall=loadGWA()
varGWA=readVarNames(do="gwa")
gwaL=subset(gwaListall, label %in% varGWA$labelGWA)
exposVar=gwaL$label

print(paste0("Read in ",NROW(exposVar), " exposure data"))
library(parallel)
cores=detectCores()
gwa=lapply(exposVar, function(x) readIN(file=paste0(GWA, "/data/GWAssumstatsClean/", x), select="yes", c=cores))
names(gwa)=exposVar

print(paste0("Read in ",NROW(gwaU), " regenie sumstats"))
gwaUL=lapply(gwaU, function(x) readIN(file=paste0(HOME, "/output/gwas/processed/",x), select="snp", c=cores))
names(gwaUL)=gwaU

print("Combine the two sets of GWA")
gwaC <- c(gwa, gwaUL)
print("Read in outcome data")
factor=lapply(outcome, function(x) fread(file= paste0(HOME, "/output/gwas/processed/",x)))
names(factor)=outcome


print(" ==============  START MR (causal effects) ========================= ")
mrOut=funMR(exposVar=names(gwaC), outIn=outcome,  datExp=gwaC, datOut=factor)
mrOut$direction="causal"
print(" ==============  FINISHED MR (causal effects) ====================== ")


mrOut$pval=2*pnorm(-abs(mrOut$b/mrOut$se))
mrOut$uCI=mrOut$b + 1.96 * mrOut$se
mrOut$lCI=mrOut$b - 1.96 * mrOut$se
print(mrOut)

print("All done!")
saveRDS(mrOut , paste0(HOME,"/output/mr/mr_",outcome,".rds"))

print(paste0(" =========== Finished MR for ", outcome, " ==============="))

}

if(outcome=="process"){
    print("Combine MR results")
    mrL=lapply(gwaU, function(x) readRDS( paste0(HOME,"/output/mr/mr_",x,".rds")))
    mrC <- mrL[sapply(mrL, function(df) NCOL(df) > 1)]
    mr=do.call(rbind, mrC)
    print(mr)
    saveOutput(object=mr, label="mr", upload="yes")
}
print("All output saved")
uploadLog(file=paste0("mr_", outcome,"_R.log"))
uploadLog(file=paste0("mr_", outcome,"_C.log"))








