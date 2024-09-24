
args = commandArgs(trailingOnly=TRUE)
HOME=args[1]

source(paste0(HOME, "/analysis/input.R"))
source(paste0(HOME, "/analysis/functions.R"))

con <- file(paste0(HOME, "/output/log/ldak_R.log"))
sink(con, append=TRUE, type="output")

declineVar=readRDS(paste0(HOME, "/output/rds/declineVar.rds"))
declineVar=c(declineVar, "cognition", "physical")
  
f <- list.files(paste0(HOME, "/output/ldak/"), include.dirs = F, full.names = T, recursive = T)
file.remove(f)

print("Save SNP data")
snps=as.data.frame(readRDS(paste0(HOME, "/output/rds/gwaClump.rds")))
snps$ID=paste0(snps$CHR, ":", snps$POS, "_", snps$A2, "_", snps$A1) # A1 flipped in regenie
CHR=levels(as.factor(snps$CHR))
extractSNP=levels(as.factor(snps$ID))
extractDF=data.frame(SNP=extractSNP)
write.table(extractDF, paste0(HOME, "/output/ldak/extractSNP"), col.names = FALSE, row.names = FALSE, quote = FALSE)
    
print("Get genotype data")
gwaCol=read.table(paste0(HOME, "/output/gwas/gwa"), header=TRUE)


print("Start LDAK for difference score")
cDiffL=lapply(paste0(declineVar, "_changeDiff"), function(x) runLDAK(pheno=x, gwaCol=gwaCol, w="PW_participationFU"))
cDiff=do.call(rbind, cDiffL)

print("Start LDAK for log-difference score")
cLogL=lapply(paste0(declineVar, "_changeLog"), function(x) runLDAK(pheno=x, gwaCol=gwaCol, w="PW_participationFU"))
cLog=do.call(rbind, cLogL)

print("Start LDAK for baseline score")
cBLL=lapply(paste0(declineVar, "_0"), function(x) runLDAK(pheno=x, gwaCol=gwaCol, w="PW_participation"))
cBL=do.call(rbind, cBLL)

wGWA=rbind(cDiff, cLog, cBL)
wGWAs=subset(wGWA, is.na(SNP)==FALSE)

print("Save output")
saveOutput(object=wGWAs, label="wGWA", upload="yes")

print("All output saved")
uploadLog(file=paste0("ldak_R.log"))
uploadLog(file=paste0("ldak_C.log"))

