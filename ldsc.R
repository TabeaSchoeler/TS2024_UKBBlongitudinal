args = commandArgs(trailingOnly=TRUE)
HOME=args[1]
pheno=args[2]
source(paste0(HOME, "/analysis/input.R"))
source(paste0(HOME, "/analysis/functions.R"))
library(GenomicSEM)


con <- file(paste0(HOME, "/output/log/ldsc_",pheno,"_R.log"))
sink(con, append=TRUE, type="output")

print(paste0("Run LDSC for ", pheno))

print("Select gwa sumstats")
gwaListall=loadGWA()
varGWA=readVarNames(do="gwa")
gwaL=subset(gwaListall, label %in% varGWA$labelGWA)
gwaL$path=paste0(GWA, "/data/GWAsumstatsLDSC/",gwaL$label, ".sumstats.gz")


cogphy=read.table(paste0(HOME, "/output/gwas/batchName"), head = F)$V1
cogphyA=as.data.frame(matrix(ncol=NCOL(gwaL), nrow=length(cogphy)))
colnames(cogphyA)=colnames(gwaL)
cogphyA$label=cogphy
cogphyA$path= paste0(HOME, "/output/gwas/processed/",cogphy, ".sumstats.gz")


h2=readRDS(paste0(HOME,"/output/rds/h2.rds"))
h2rem=subset(h2, is.na(h2)==TRUE)
h2=subset(h2, is.na(h2)==FALSE)
use=pheno

gwaIn=rbind( cogphyA)
gwaS=subset(gwaIn, !label %in% h2rem$pheno)

if(pheno!="process"){
file.remove(paste0(HOME,"/output/rds/rg_", pheno, ".rds"))

if(NROW(subset(h2, pheno %in% use))==0){
    print("No heritability to estimate genetic correlatiobs")
}

if(NROW(subset(h2, pheno %in% use))==1){
    print("Get genetic correlations")
    rgL=list()
    rgLP=lapply(gwaS$label, function(x) ldscRun(df=gwaS, trait1=pheno, trait2=x, return="rg") )
    rgdf=do.call(rbind, rgLP)
    rgdf$pval=2*pnorm(-abs(rgdf$rg/rgdf$rg_se))
    rgdf$uCI=rgdf$rg + 1.96 * rgdf$rg_se
    rgdf$lCI=rgdf$rg - 1.96 * rgdf$rg_se
    rgdf=subset(rgdf, factor!=pheno)
    print(paste0('Finnished genetic correlation analyses for ', pheno))
    saveRDS(rgdf , paste0(HOME,"/output/ldsc/rds/rg_",pheno,".rds"))

}
}

if(pheno=="process"){
 rgL=lapply(h2$pheno, function(x) readRDS(  paste0(HOME,"/output/ldsc/rds/rg_",x,".rds")))
 rg=do.call(rbind, rgL)
 saveOutput(object=rg, label="rg", upload="yes")

}

print("All output saved")
uploadLog(file=paste0("ldsc_",pheno,"_R.log"))
uploadLog(file=paste0("ldsc_",pheno,"_C.log"))