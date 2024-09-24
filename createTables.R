
################################################################
# ====================== Create output tables ==================
################################################################

# +++++++++++++++++++++++++ RUN FUNCTIOND AND DEFINE STYLES ++++++++++++++++++++++++++++++++++++++++++++++++++
ColStart=2
RowHeader=2
RowSubheaderStart=3
RowSubheaderEnds=3
RowTable=4

library(openxlsx)

# Create info text
createInfo=function(dataInfoPath){
  datOut=read.csv(dataInfoPath,header=T)
  datOut$X=NULL
  datOut_merged=paste0(datOut[,1],": " ,datOut[,2])
  return(datOut_merged)
}

# Headerstyle
hs1 <- createStyle(halign = "CENTER", textDecoration = "Bold",
                   border = "Bottom", fontColour = "black", fgFill = "white")

h_info <- createStyle(halign = "left", textDecoration = "Bold",
                      border = "Bottom", fontColour = "black", fgFill = "white")

addTable=function(sheet, table){
  writeDataTable(wb, sheet, table, headerStyle=hs1, tableStyle = "TableStyleLight1",
                 startRow = RowTable, startCol = ColStart)
  setColWidths(wb, sheet, cols=2:10, widths = 15)
}

# HEADER
headerFunc=function(TITLE, sheet){
  writeData(wb, sheet = sheet, TITLE,
            colNames = FALSE, rowNames = FALSE,
            startCol = ColStart, startRow = RowHeader)
}
# INFO ROW
InfoFunc=function(TITLE, sheet){
  writeData(wb, sheet = sheet, TITLE,
            colNames = FALSE, rowNames = FALSE,
            startCol = ColStart, startRow = RowSubheaderStart)
}

# Create new workbook
setwd(paste0(HOME,"/results/tables/"))
wb <- openxlsx::createWorkbook()




# ============== Table 1 =================
dateMissS=subset(dateMiss, select=c(label_clean, timeDate, n, meanDate, minDate, maxDate, assB))
dateMissS$meanDate=as.character(dateMissS$meanDate)
dateMissS$minDate=as.character(dateMissS$minDate)
dateMissS$maxDate=as.character(dateMissS$maxDate)
colnames(dateMissS)=c("Aging phenotype", "Time point", "n (complete cases)", "Mean date of assessment", "Earliest date of assessment", "Latest date of assessment", "Type of assessment (face-to-face versus online)")

missDat="sTable 1"
addWorksheet(wb, missDat)
# Add parameters
title_name="sTable 1. Missing data and dates of follow-up for cognitive and physical variables"
sheet=missDat
table=dateMissS
Info_text=""
# Run functions
addTable(sheet, table)
headerFunc(title_name, sheet)
InfoFunc(Info_text, sheet)

# ============== Table 2 =================
# Sumstats included in Mendelian Randomization
setwd(HOME)
setwd("..")
GWA=getwd()
mrList=subset(variables, gwa=="yes")
sumstats <- read.csv(paste0(GWA, "/GWAsumstats//data/sumstatsInfo.csv"), na = c("N/A", "n/a", NA, "NA"))
datMR=merge(mrList, sumstats, by.x="labelGWA", by.y="label", all.x=T)
datMR$label_clean=ifelse(is.na(datMR$label_clean.x)==T, datMR$label_clean.y, datMR$label_clean.x)
datMRs=subset(datMR, select=c(label_clean, Publication, link, N))
datMRs$Publication=ifelse(is.na(datMRs$Publication)==T, "Generated in the UKBB using REGENIE", datMRs$Publication)
colnames(datMRs)=c("Trait", "Publication", "Link to summary statistic file", "N")

mrStats="sTable 2"
addWorksheet(wb, mrStats)
# Add parameters
title_name="sTable 2. Summary statistic files included in Mendelian Randomization analyses"
sheet=mrStats
table=datMRs
Info_text=""
# Run functions
addTable(sheet, table)
headerFunc(title_name, sheet)
InfoFunc(Info_text, sheet)

# ============== Table 3 =================
changeSum=readRDS(paste0(HOME,"/output/rds/changeSum.rds"))
changeSumC=recodeChange(df=changeSum, var="var")
changeSum$var=changeSumC$label_clean
changeSum$meanFU=round(changeSum$meanFU, 2)
changeSum$dimension=changeSumC$dimension
changeSumO=changeSum[order(changeSum$dimension, decreasing = FALSE), ] 
colnames(changeSumO)=c("Aging phenotype", "Follow-up time (minimum, in years)", "Follow-up time (maximum, in years)", "Follow-up time (mean, in years)","Number of removed outliers","n (complete cases with follow-up data)", "Dimension")

FUtime="sTable 3"
addWorksheet(wb, FUtime)
# Add parameters
title_name="sTable 3. Follow-up duration and sample size per longitudinal aging phenotype"
sheet=FUtime
table=changeSumO
Info_text=""
# Run functions
addTable(sheet, table)
headerFunc(title_name, sheet)
InfoFunc(Info_text, sheet)

# ============== Table 4 =================
# age effects
head(ageEffS)
ageEffS$N=round(ageEffS$n, 0)
ageEffS$eff=paste0(round(ageEffS$betaSTD, 3), " (", round(ageEffS$lCISTD, 3), "; ",round(ageEffS$uCISTD, 3), ")")
#ageEffS$Pf=formatC(ageEffS$P, 4)
ageEffsupp=subset(ageEffS, select=c(label_clean, groupC, N, eff))
colnames(ageEffsupp)=c("Aging phenotype", "UK Biobank sample", "(Effective) sample size", "beta (95% CI)")

ageEfft="sTable 4"
addWorksheet(wb, ageEfft)
# Add parameters
title_name="sTable 4. Age effects on aging phenotypes across levels of sample representativeness"
sheet=ageEfft
table=ageEffsupp
Info_text=""
# Run functions
addTable(sheet, table)
headerFunc(title_name, sheet)
InfoFunc(Info_text, sheet)


# ============== Table 5 =================
# weighted gwa
wGWAcomb=rbind(wGWAlog, wGWAdiff, wGWASbl)
wGWAcombM=merge(wGWA, subset(wGWAcomb, select=c(SNP, gene)), by="SNP", all.x=T)
wGWAIDs=subset(wGWAcombM, select=c(pheno_clean, model, SNP, gene,A1, A2, CHR, POS, BETA, SE, P, N, NEFF, type))
colnames(wGWAIDs)=c("Aging phenotype", "Model (unweighted or IP-weighted)", "SNP", "Nearest gene","A1", "A2","CHR", "POS", "BETA", "SE", "P", "N", "NEFF", "Outcome type")

gwaWEI="sTable 5"
addWorksheet(wb, gwaWEI)
# Add parameters
title_name="sTable 5. Standard and Inverse Probability Weighted genetic variant effects"
sheet=gwaWEI
table=wGWAIDs
Info_text="The colume 'outcome type' specifies model that was used to derive the outcome phenotype, including 1) log = log-difference model (relative change), diff = difference score (absolute change) 3) res = residual change score model (conditional change) and 4) baseline = the cross-sectional phenotype assessed at baseline. "
# Run functions
addTable(sheet, table)
headerFunc(title_name, sheet)
InfoFunc(Info_text, sheet)

# ============== Table 6 =================
# REGENIE results
clumpSupp=subset(clump, !typeR %in% c("other", "baselineFU", "missed"))
clumpSupps=subset(clumpSupp, select=c(label_clean, SNP, A1, A2, CHR, POS, BETA, SE, P, N, INFO, EAF, typeR))
colnames(clumpSupps)=c("Aging phenotype", "SNP", "A1", "A2","CHR", "POS", "BETA", "SE", "P", "N", "INFO", "EAF","Outcome type")

gwaCHA="sTable 6"
addWorksheet(wb, gwaCHA)
# Add parameters
title_name="sTable 6. Cross-sectional and longitudinal genetic variant effects"
sheet=gwaCHA
table=clumpSupps
Info_text="The colume 'outcome type' specifies model that was used to derive the outcome phenotype, including 1) log = log-difference model (relative change), diff = difference score (absolute change) 3) res = residual change score model (conditional change) and 4) baseline = the cross-sectional phenotype assessed at baseline. "
# Run functions
addTable(sheet, table)
headerFunc(title_name, sheet)
InfoFunc(Info_text, sheet)

# ============== Table 7 =================
ukbbLMsupp=ukbbLM
ukbbLMsupp$b_CI=paste0(round(ukbbLM$beta, 3), " (", round(ukbbLM$lCI, 3), "; ", round(ukbbLM$uCI, 3), ")")
ukbbLMsupp$p=formatC(ukbbLMsupp$pval, 2)
ukbbLMsupp$label=paste0(ukbbLMsupp$labelOut, " - log-difference score")
ukbbLMS=subset(ukbbLMsupp, select=c(labelExp, label, b_CI, p, n))
colnames(ukbbLMS)=c("Exposure", "Outcome", "beta (95% CI)", "p-value", "Sample size")
phenoAsso="sTable 7"
addWorksheet(wb, phenoAsso)
# Add parameters
title_name="sTable 7. Phenotypoc associations with change (log-difference)"
sheet=phenoAsso
table=ukbbLMS
Info_text=""
# Run functions
addTable(sheet, table)
headerFunc(title_name, sheet)
InfoFunc(Info_text, sheet)

# ============== Table 8 =================
mrSupp$b_CI=paste0(round(mrSupp$b, 3), " (", round(mrSupp$lCI, 3), "; ", round(mrSupp$uCI, 3), ")")
mrSupp$p=formatC(mrSupp$pval, 2)
mrSuppS=subset(mrSupp, select=c(phenoX, phenoY, phenoYtype, b_CI, p, nsnp, Nexp, Nout))
colnames(mrSuppS)=c("Exposure", "Outcome", "Outcome type", "beta (95% CI)", "p-value", "n genetic instruments", "Sample size (exposure)", "Sample size (outcome)")

mrRES="sTable 8"
addWorksheet(wb, mrRES)
# Add parameters
title_name="sTable 8. Mendelian Randomization results"
sheet=mrRES
table=mrSuppS
Info_text="The colume 'outcome type' specifies model that was used to derive the outcome phenotype, including 1) log = log-difference model (relative change), diff = difference score (absolute change) 3) res = residual change score model (conditional change) and 4) baseline = the cross-sectional phenotype assessed at baseline. "
# Run functions
addTable(sheet, table)
headerFunc(title_name, sheet)
InfoFunc(Info_text, sheet)

################################################################
# ====================== Export table ==================
################################################################

# Create new styles
s <- createStyle(fgFill = "#FFFFFF")
h_info <- createStyle(halign = "left",
                      border = "BOTTOM", fontColour = "black", fgFill = "white", fontSize=16, textDecoration = "Bold", numFmt="TEXT", borderColour = "black")
info_info <- createStyle(halign = "left",
                         border = NULL, fontColour = "black", fgFill = "white", fontSize=14, textDecoration = NULL, numFmt="TEXT", wrapText=TRUE)
# Run loop
for(curr_sheet in names(wb)){
  addStyle(wb,sheet = curr_sheet, s, cols=1:40, rows=1:2000, gridExpand = TRUE)
  setColWidths(wb, sheet = curr_sheet, cols=1:40, widths = 20)
  addStyle(wb,sheet = curr_sheet, h_info, cols=ColStart:20, rows=RowHeader, gridExpand = TRUE)
  addStyle(wb,sheet = curr_sheet, info_info, cols=ColStart:5, rows=RowSubheaderStart, gridExpand = TRUE)
  mergeCells(wb,sheet = curr_sheet, cols = 2:100, rows = RowSubheaderStart:RowSubheaderEnds)
}

library( openxlsx)
openxlsx::saveWorkbook(wb, paste0(HOME,"/results/tables/declineUKBB.xlsx"), overwrite = TRUE)
# Open File
openXL(wb)
