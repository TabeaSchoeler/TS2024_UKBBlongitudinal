#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
HOME=args[1]
pheno=args[2]

source(paste0(HOME, "/analysis/input.R"))
source(paste0(HOME, "/analysis/functions.R"))

library(GenomicSEM)
print("Process GWA")
con <- file(paste0(HOME, "/output/log/processGWA.",pheno,"_R.log"))
sink(con, append=TRUE, type="output")

if(pheno!="process"){
    print('================ Start processing REGENIE output =========================')
    GWAout=readGWA(model=pheno)
    saveGWA(df=GWAout, name=pheno)
    print("All GWA processed and saved on cluster")        
    nSNP=NROW(GWAout)

    print("================ Clump data =========================")
    gwaClump=clumpData(clump_input=subset(GWAout, P<5e-8), name=pheno, p="gwa") # select gwa (not mr) clumping parameters
    gwaClump$nSNP=nSNP

    if(ncol(gwaClump) <= 4){
        cols=c("CHR", "POS", "SNP", "A1", "A2", "BETA", "SE", "P", "INFO", "N", "EAF", "MAF", "pheno", "nSNP")
        gwaClump=as.data.frame(matrix(ncol=length(cols), nrow=0))
        colnames(gwaClump)=cols
    }
    
    saveOutput(object=gwaClump, label=paste0("clump_",pheno), upload="no")
    print("================ Munge data =========================")
    mungeData(name=pheno, path=paste0(HOME, "/output/gwas/processed/"))

}
        

if(pheno=="process"){
       print("================ upload clumped SNPs =========================")
    readIn=read.table(paste0(HOME, "/output/gwas/batchName"), head = F)$V1

    gwaC=lapply(readIn, function(x) readRDS(paste0(HOME,"/output/rds/clump_",x,".rds")))
    gwaClumpDF=do.call(rbind, gwaC)
    saveOutput(object=gwaClumpDF, label="gwaClump", upload="yes")

    extractDF=data.frame(snp=levels(as.factor(gwaClumpDF$SNP)))
    write.table(extractDF, paste0(HOME, "/output/gwas/plink/extractPlink"), col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    print("Extract beta estimates for all gwa sig hits")
    BetaSNPL=lapply(readIn, function(x) getBetaSNP(x ,extract=levels(as.factor(gwaClumpDF$SNP))))
    BetaSNP=do.call(rbind, BetaSNPL)
    saveOutput(object=BetaSNP, label="BetaSNP", upload="yes")

    print("Get SNP heritability estimates")
    gwaListall=loadGWA()
    gwaL=readVarNames(do="gwa")
    addPCA=as.data.frame(matrix(ncol=NCOL(gwaL), nrow=length(readIn)))
    colnames(addPCA)=colnames(gwaL)
    addPCA$label=readIn
    addPCA$path= paste0(HOME, "/output/gwas/processed/",addPCA$label, ".sumstats.gz")
    
    print("Get SNP heritability estimates")
    h2L=lapply(readIn, function(x) ldscRun(df=addPCA, trait1=x, trait2=x, return="h2") )
    h2df=do.call(rbind, h2L)
    h2df$pval=2*pnorm(-abs(h2df$h2/h2df$h2_se))
    h2df$uCI=h2df$h2 + 1.96 * h2df$h2_se
    h2df$lCI=h2df$h2 - 1.96 * h2df$h2_se
    print('Finnished SNP heritability analyses - save output')
    print(h2df)
    saveOutput(object=h2df, label="h2" , upload="yes")


    sumVar=c("cognition", "physical")
    varSUM =readRDS(paste0(HOME, "/output/rds/declineVar.rds"))
    nameLabel=c(sumVar,varSUM)

    print("Get SNP effects for dropout sample versus FU sample")
    missFUeff=compSNPs(t1="MISS", t2="FU0", nameLabel=nameLabel) 
    saveOutput(object=missFUeff, label="missFUeff", upload="yes")
    
    print("Get SNP effects for baseline versus change")
    SNPdiff=compSNPs(t1="0", t2="changeDiff", nameLabel=nameLabel) 
    saveOutput(object=SNPdiff, label="SNPdiff", upload="yes")
    
    print("Get SNP effects for baseline versus change (log scale)")
    SNPlog=compSNPs(t1="0", t2="changeLog", nameLabel=nameLabel) 
    saveOutput(object=SNPlog, label="SNPlog", upload="yes")
    
    print("Get SNP effects for baseline versus change (residual change score)")
    SNPres=compSNPs(t1="0", t2="changeRes", nameLabel=nameLabel) 
    saveOutput(object=SNPres, label="SNPres", upload="yes")

    print("Get SNP effects for baseline (FU sample) versus change")
    SNPFUeff=compSNPs(t1="FU0", t2="changeDiff", nameLabel=nameLabel) 
    saveOutput(object=SNPFUeff, label="SNPFUeff", upload="yes")
    
    print("Get SNP effects for baseline (Dropout sample) versus change")
    SNPmisseff=compSNPs(t1="MISS", t2="changeDiff", nameLabel=nameLabel) 
    saveOutput(object=SNPmisseff, label="SNPmisseff", upload="yes")

 

}


print("All output saved")
uploadLog(file=paste0("processGWA.",pheno,"_R.log"))
uploadLog(file=paste0("processGWA.",pheno,"_C.log"))

