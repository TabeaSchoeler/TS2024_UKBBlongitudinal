
#######################################################
# =================== FUNCTIONS =======================
#######################################################
rm(list = ls())
library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..")
HOME=getwd()
source(paste0(HOME, "/analysis/funcLongitudinal.R"))

# =======================================================================
# ======================= RUN SIMULATIONS  ============================
# =======================================================================
source(paste0(HOME,"/analysis/sim.R"))

# =======================================================================
# ======================= Data availability  ============================
# =======================================================================
# Data availability (number of waves)
dateMiss=readRDS(paste0(HOME,"/output/rds/dateMiss.rds"))
dateMiss=subset(dateMiss, n>0)
dateMiss$assB=ifelse(dateMiss$assessment=="date_assessment", "Centre", "Online")
missLabel=recodeChange(df=dateMiss, var="var")
dateMiss$label_clean=missLabel$label_clean
dateMiss$dimension=missLabel$dimension
dateMiss$TA=paste0(dateMiss$assB, " (" , dateMiss$timeDate+1, " waves)")
dateMiss$TA=ifelse(dateMiss$timeDate==0, paste0(dateMiss$assB, " (1 wave)"), dateMiss$TA)

# Data availability (longitudinal)
missFU=readRDS(paste0(HOME,"/output/rds/missFU.rds"))
missFUL=recodeChange(df=missFU, var="pheno")
missFU$label_clean=missFUL$label_clean
missFU$dimension=missFUL$dimension


dateMiss$m=paste0(dateMiss$var, "_", dateMiss$timeDate)
missFU$m=paste0(missFU$pheno, "_", missFU$time)
datUKFU=merge(dateMiss, missFU, by="m", all=T, suffixes = c("", "_FU"))

# plot availability of follow up assessments
dataPlotC=availPlot(df=datUKFU, dim="cognition", c1="darkgreen", c2="magenta4", left=0, right=0,  caption="Cognitive measures")
dataPlotP=availPlot(df=datUKFU, dim="physical", c1="darkblue", c2="darkred", left=0, right=20, caption="Physical measures")
dataPlot=ggarrange(dataPlotC, dataPlotP, labels=c("", ""), nrow=2)
saveFigure(fileName=paste0(HOME,"/results/figures/dataPlot"), plotName=dataPlot, w=60, h=30)

# =======================================================================
# ========== Within-trait correlations across time ======================
# =======================================================================
corL=readRDS(paste0(HOME,"/output/rds/corL.rds"))
corTime=do.call(rbind, corL)
corLabel=recodeChange(df=corTime, var="var")
corTime$dimension=corLabel$dimension
corTime$label=corLabel$label_clean
corTime$Var1=paste0( "W ", corLabel$t1)
corTime$Var2=paste0( "W ", corLabel$t2)
cogCorP=CorPlot(df=corTime, dim="cognition", col="aquamarine4")
saveFigure(fileName=paste0(HOME,"/results/figures/cogCorP"), plotName=cogCorP, w=20, h=40)
phyCorP=CorPlot(df=corTime, dim="physical", col="cornflowerblue")
saveFigure(fileName=paste0(HOME,"/results/figures/phyCorP"), plotName=phyCorP, w=15, h=21)


# =======================================================================
# ==================== Means across age bins ============================
# =======================================================================
ageMean=readRDS(paste0(HOME,"/output/rds/ageMean.rds"))
ageL=data.frame(pheno=levels(as.factor(ageMean$pheno)))
ageC=recodeChange(df=ageL, var="pheno")

plotAgeL=lapply(ageC$pheno, function(x) plotAge(df=ageMean, sel=x))
names(plotAgeL)=ageC$label_clean
agePlotS=ggarrange(plotlist=plotAgeL, common.legend = T, ncol=3, nrow=6)
saveFigure(fileName=paste0(HOME,"/results/figures/plotAge"), plotName=agePlotS, w=40, h=70 )


# =======================================================================
# ================= Age effects in subgroups ============================
# =======================================================================
# get data from UKBB
ageEff=readRDS(paste0(HOME,"/output/rds/ageEff.rds"))
ageEff=subset(ageEff, var %in% fuvar)

reorderAge=c("UKBB (baseline sample, 1-wave)", 
             "UKBB (weighted sample, 1-wave)",
             "UKBB (follow-up sample, 2-waves)",
             "UKBB (follow-up sample, 3-waves)")

ageEff$groupC=revalue(  ageEff$group, c(
  "betweenAll" = reorderAge[1],
  "weighted" = reorderAge[2],
  "betweenFU" = reorderAge[3],
  "between2FU"= reorderAge[4]
))

ageEff$groupC <- factor(ageEff$groupC, levels=reorderAge)


levels(as.factor(ageEff$groupC))
colsAge= c("wheat3", "coral","#AE4371", "pink", "brown" )
ageEffP=recodeChange(df=ageEff, var="var")
ageEffS=ageEffP

# ====== Assess shrinkage in beta coefficents
ageShrinkAll=shrinkDat(df=ageEffS, reference="betweenAll", comp="betweenFU")
ageShrinkW=shrinkDat(df=ageEffS, reference="weighted", comp="betweenFU")
ageShrinkW$nLoss=ageShrinkAll$nLoss
ageShrink=rbind(ageShrinkW, ageShrinkAll)
selBL=unique(subset(ageShrinkAll, n_all>450000)$var)
ageShrink=subset(ageShrink, var %in% selBL)

reorderShrink=c("Reference: weighted sample", 
                "Reference: baseline sample")
ageShrink$sampleC=revalue(  ageShrink$groupS_all, c(
  "weighted_ukbb" = reorderShrink[1],
  "betweenAll_ukbb" = reorderShrink[2]))

ageShrink$sampleC<- factor(ageShrink$sample, levels=reorderShrink)
colsShrink= c( "gray62", "gray23" )

e1=expression(paste( beta[0], ": UKBB (baseline sample); ", beta[1], ": UKBB (2-wave follow-up sample)" ))
e2=expression(paste( beta[0], ": UKBB (weighted baseline sample); ", beta[1], ": UKBB (2-wave follow-up sample)" ))
labelShrink=c(e1, e2)

shrinkPlot <- ggplot(ageShrink, aes(x=nLoss, y=bLoss, colour=sample)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  theme_minimal() +
  xlim(50,100) +
  guides( shape = "none")  +
  scale_colour_manual("",values = colsShrink, label=labelShrink) +
  theme(legend.position="top",
        legend.direction = "vertical",
        panel.border = element_rect(colour = "black", fill=NA, size=0.1),
        plot.margin = margin(t=0, r=0, b=0, l=1, "cm")) +
  labs(x =  "Loss of sample at follow up (%)", 
       caption="",
       y = expression(paste( italic(beta[0]) - italic(beta[1]) )))
print(shrinkPlot)

# ====== Plot beta coefficents of age effects
ageEffS=subset(ageEffP, var %in% fuvar)
ageEffS=subset(ageEffS, var %in% selBL)
ageEffplot=plotAgeComp(df=ageEffS, bottom=0, top=1, col=colsAge)

# ===== Combine all figures
ageExamples=ggarrange(plotAgeL[["Height"]], plotAgeL[["Fluid intelligence score"]], nrow=1, labels=c("A", ""), common.legend = T)
ageSumPlot=ggarrange(ageEffplot , shrinkPlot, nrow=1, labels=c("B", "C"), widths=c(1, 1))

ageEffComb=ggarrange(ageExamples, ageSumPlot, 
                     ncol=1, 
                     hjust=0, 
                     heights=c(1, 1))

saveFigure(fileName=paste0(HOME,"/results/figures/agePlot"), plotName=ageEffComb, w=28, h=28 )




# =======================================================================
# ====== Between baseline/change score correlations  ====================
# =======================================================================
colCognition="aquamarine4"
colPhysical="cornflowerblue"


# ======= Change scores ==========
corChangeL=readRDS(paste0(HOME,"/output/rds/corChange.rds"))
corChangeC=prepMatrix(df=corChangeL[["cognition"]], col=colCognition)
corChangeP=prepMatrix(df=corChangeL[["physical"]], col=colPhysical, left=2)
# combine
corChange=ggarrange(corChangeC, corChangeP, nrow=2, heights=c(1, 0.7))
saveFigure(fileName=paste0(HOME,"/results/figures/corChange"), plotName=corChange, w=30, h=45)

physicalVarsA=recodeChange(df=data.frame(var=rownames(corChangeL[["physical"]])), var="var")
physicalVarsR=levels(as.factor(physicalVarsA$label))
physicalVars=physicalVarsR[ !physicalVarsR == 'physical']

cognitiveVarsA=recodeChange(df=data.frame(var=rownames(corChangeL[["cognition"]])), var="var")
cognitiveVarsR=levels(as.factor(cognitiveVarsA$label))
cognitiveVars=cognitiveVarsR[ !cognitiveVarsR == 'cognition']

# ======= Baseline ==========
corA=readRDS(paste0(HOME,"/output/rds/corA.rds"))
corBLC=prepMatrix(df=corA[["cognition"]], col=colCognition, wave="bl")
corBLP=prepMatrix(df=corA[["physical"]], col=colPhysical, wave="bl", left=1)
# combine
corBL=ggarrange(corBLC, corBLP, nrow=2, heights=c(1, 0.7))
saveFigure(fileName=paste0(HOME,"/results/figures/corBL"), plotName=corBL, w=20, h=25)


# =======================================================================
# ======= Compare predictors participation and FU =======================
# =======================================================================
predPB=readRDS( paste0(HOME,"/output/rds/predPB.rds"))
predPB=subset(predPB, iteration=="unweighted")
predPB=subset(predPB, label!="FU")
predPB$P=predPB[['Pr(>|t|)']]
predPB$labelW=paste0(predPB$label,"_" ,predPB$iteration)

predPB$labelRecoded=recodeHarmVar1(predPB$labelShort)
predPB=subset(predPB, variable!="(Intercept)")
predPB$order=seq(1, NROW(predPB), 1)
colSELps=c('violetred4', 'seagreen', "lightpink", "darkseagreen")
predPB$labelClean=recodeParticipation(predPB$labelW)
predPB$pval_fdr <- do.call(rbind, lapply(predPB$P, function(x) p.adjust(x, "fdr",  NROW(levels(as.factor(predPB$labelShort))) )))
predPB$pval_fdr_sig = ifelse(predPB$pval_fdr <0.05 , "sig", "ns")
levels(predPB$labelClean)

predPBplot <- ggplot(predPB, aes(x = fct_reorder(labelRecoded, -est), 
                                 y = est, 
                                 ymin = lCI, 
                                 ymax = uCI,
                                 fill=labelClean)) +
  geom_bar(stat = "identity", position = "dodge" ) + 
  geom_errorbar(width=.2, position=position_dodge(.9))  +
  scale_fill_manual("",values = colSELps) + 
  theme_classic() + 
  guides( color = "none") +
  theme(legend.position="top",
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        plot.title = element_text(size = 15, face = "bold", hjust = 0),
        plot.caption = element_text(size=11, hjust = 0),
        plot.margin = margin(t=0, r=0, b=0, l=1.5, "cm"),
        axis.text.x=element_text(size=12, angle=45, hjust=1),
        axis.title.y = element_text( size=20), 
        legend.text=element_text(size=15),
        strip.background = element_blank(),
        panel.grid.minor = element_blank())  + 
        guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  scale_x_discrete(position = "bottom") +
  labs(x =  "", 
       title="A. Predictors of UK Biobank participation",
       y = expression(paste( italic(beta[STD]) ))) + themeTransparent 
print(predPBplot)


# =======================================================================
# ======= Baseline UKBB predictors FU participation =====================
# =======================================================================
baselineP=readRDS( paste0(HOME,"/output/rds/baselineP.rds"))
baselineP$label=baselineP$labelShort %>% str_replace("_0", "")
baselineP=recodeChange(df=baselineP, var="label")

baselinePplot <- ggplot(baselineP, aes(x = fct_reorder(label_clean, -est), 
                                       y = est, 
                                       ymin = lCI, 
                                       ymax = uCI)) +
  geom_bar(stat = "identity", position = "dodge"  , fill = "deepskyblue4") + #
  theme_classic() + 
  guides( color = "none") +
  geom_errorbar(width=.2, position=position_dodge(.9))  +
  theme(legend.position="top",
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 15, face = "bold", hjust = 0),
        plot.caption = element_text(vjust = 1, size=14),
        plot.margin = margin(t=0, r=0, b=0, l=1.5, "cm"),
        axis.title.y = element_text( size=20), 
        axis.text.x=element_text(size=12, angle=45, hjust=1))  + 
  scale_x_discrete(position = "bottom") +
  labs(x =  "", 
       title="B. Predictors of participation status (follow-up versus dropout)",
       y = expression(paste( italic(beta[STD]) ))) + themeTransparent 
print(baselinePplot)


# =======================================================================
# ========================== Sampling weights ===========================
# =======================================================================
lassoW=readRDS( paste0(HOME,"/output/rds/lassoWeights.rds"))
lassoRes=readRDS( paste0(HOME,"/output/rds/lassoSum.rds"))
participationFU=subset(lassoRes, model=="participationFU")
participation=subset(lassoRes, model=="participation")

fuText1=TeX(paste0("$N_{EFF}/N=$", comma(participationFU$nEFF), "/", comma(participationFU$n)))
fuText2=TeX(paste0("$\\sigma^2=$ ", round(var(lassoW$PW_participationFU, na.rm=T), 2), " (range = ", round(min(lassoW$PW_participationFU, na.rm=T), 0), " - ", round(max(lassoW$PW_participationFU, na.rm=T), 0) ,")"))

blText1=TeX(paste0("$N_{EFF}/N=$", comma(participation$nEFF), "/", comma(participation$n)))
blText2=TeX(paste0("$\\sigma^2=$ ", round(var(lassoW$PW_participation, na.rm=T), 2), " (range = ", round(min(lassoW$PW_participation, na.rm=T), 0), " - ", round(max(lassoW$PW_participation, na.rm=T), 0) ,")"))


blWL=TeX(paste0("Baseline participation weights"))
fuWL=TeX(paste0("Follow-up participation weights"))

lassoL=data.frame(weights=c(lassoW$PW_participation, lassoW$PW_participationFU), 
                  label=c(rep("bl", NROW(lassoW)), rep("fu", NROW(lassoW)) ))

wDens=ggplot(lassoL,aes(x=weights, colour=label, fill=label)) + 
  geom_density( alpha = 0.5) +
  xlim(0,5) +
  themeTransparent +
  theme_minimal() + 
  scale_colour_manual("",values = c("violetred4", "seagreen" ), guide = "none") + 
  scale_fill_manual("", values = c("violetred4", "seagreen" ), labels=c(blWL, fuWL)) + 
  theme(legend.position="right",
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.title.y = element_text( size=20), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t=0, r=0, b=0, l=0, "cm"), 
        legend.text=element_text(size=20),
        plot.title = element_text(size = 15, face = "bold", hjust = 0)) + #        plot.title = element_text(size = 13,  hjust = 0)
  geom_vline(xintercept=1) +
  labs(x =  "",
       y = "Density",
       title =  "C. Distribution of baseline and follow-up participation weights")  +
  annotate("text",size=5, x=1.2,y=0.6,hjust=0,vjust=0,label=blText1, colour="violetred4") +
  annotate("text",size=5, x=1.25,y=0.5,hjust=0,vjust=0,label=blText2, colour="violetred4") +
  annotate("text",size=5, x=1.5,y=0.3,hjust=0,vjust=0,label=fuText1, colour="seagreen") +
  annotate("text",size=5, x=1.55,y=0.2,hjust=0,vjust=0,label=fuText2, colour="seagreen") 
wDens


# ====== Combine the two plots
fuBiasComb=ggarrange(predPBplot, baselinePplot, wDens, 
                     nrow=3, heights=c(0.9,1,0.6),
                     labels=c("", "", ""))
saveFigure(fileName=paste0(HOME,"/results/figures/fuBiasPlot"), plotName=fuBiasComb, w=37, h=50 )


# =======================================================================
# ============= Dropout versus FU sample (GWA results) ==================
# =======================================================================
missFU=readRDS( paste0(HOME,"/output/rds/missFUeff.rds"))
decVar=readRDS( paste0(HOME,"/output/rds/declineVar.rds"))
colsFU=c(rep(colPhysical, length(physicalVars)), rep(colCognition, length(cognitiveVars)))
colsDF=data.frame(var=c(physicalVars,cognitiveVars), col=colsFU)
missFUp=lapply(c(physicalVars,cognitiveVars), function(x) SNPPlot(df=missFU[[x]], g=x, title="", test="difference", nmin=0, addT="", colTitle=subset(colsDF, var==x)$col))

ggexport(plotlist = missFUp, 
         filename =paste0(HOME,"/results/figures/missFUs.pdf"), 
         nrow = length(missFUp)/2, 
         ncol = 2,
         width = 14, height = 16)


# =======================================================================
# ========= SNP effects: Cross-sectional versus change ==================
# =======================================================================
# clumped SNPs
clump=readRDS(paste0(HOME,"/output/rds/gwaClump.rds"))
clump$typeR=recodeChange(df=clump, var="pheno", return="type")
clump$phenoR=recodeChange(df=clump, var="pheno", return="pheno")
clumpC=recodeChange(df=clump, var="phenoR")
clump$label_clean=clumpC$label_clean
clump$dimension=clumpC$dimension
clump$uCI=clump$BETA + 1.96 * clump$SE
clump$lCI=clump$BETA - 1.96 * clump$SE

# Phewas plot
clumpChange=levels(as.factor(subset(clump, typeR %in% c("log", "diff") & P < 5e-8)$SNP))
snpPWL=lapply(clumpChange, function(x) pheWAplot(snp=x, allele=clump))
names(snpPWL)=clumpChange
snpPW=ggarrange(plotlist=snpPWL, ncol = 2, nrow=4)
legendPW=pheWAplot(snp="rs429358", allele=clump, return = "legend")
legendPWo=as_ggplot(legendPW)
snpPWl=ggarrange(legendPWo, snpPW, nrow=2, heights=c(0.2,1))
saveFigure(fileName=paste0(HOME,"/results/figures/snpPW"), plotName=snpPWl, w=30, h=50)


# ================= Compare beta estimnates ====================================
cols=recodeDelta(var=NULL, return="colour")

# Log difference score
SNPlog=readRDS( paste0(HOME,"/output/rds/SNPlog.rds"))
SNPlogCog=SNPPlot(df=SNPlog[["cognition"]], g="cognition", addT="LOG",  title="no", nmin=0, colTitle=cols[1])
SNPlogPhy=SNPPlot(df=SNPlog[["physical"]], g="physical",  addT="LOG",  title="no", nmin=0, colTitle=cols[1])

# Difference score
SNPdiff=readRDS( paste0(HOME,"/output/rds/SNPdiff.rds"))
SNPdiffCog=SNPPlot(df=SNPdiff[["cognition"]], g="cognition",  title="no",addT="DIFF", nmin=0, colTitle=cols[2])
SNPdiffPhy=SNPPlot(df=SNPdiff[["physical"]], g="physical", title="no",addT="DIFF",  nmin=0, colTitle=cols[2])

# Combine all plots
cogDEC=ggarrange(SNPlogCog, SNPdiffCog, ncol = 1)
phyDEC=ggarrange(SNPlogPhy, SNPdiffPhy, ncol = 1)
SNPeffPlot=ggarrange(cogDEC, phyDEC, ncol = 2,  hjust=-0.7, font.label = list(size = 18, face="plain"),
                     labels=c("Cognitive decline", "Physical decline"))

saveFigure(fileName=paste0(HOME,"/results/figures/SNPeffPlot"), plotName=SNPeffPlot, h=20, w=25)

# ================ SUPPLEMENT PLOTS: LONGITUDINAL VERSUS CROSS-SECTIONAL ====================
SNPres=readRDS( paste0(HOME,"/output/rds/SNPres.rds"))
snpPlotCombined(trait=names(SNPlog))


# =======================================================================
# =========================== Weighted GWA ==============================
# =======================================================================
wGWA=readRDS(paste0(HOME,"/output/rds/wGWA.rds"))
wGWA$model=wGWA$label
wGWA$type=recodeChange(df=wGWA, var="pheno", return="type")
wGWA$label=recodeChange(df=wGWA, var="pheno", return="pheno")
phenoW=recodeChange(df=wGWA, var="label")
wGWA$pheno_clean=phenoW$label_clean
wGWA$dimension=phenoW$dimension

# log difference
wGWAlog=weightedGWA(df=do.call(rbind, SNPlog), dfW=wGWA, model="log")
# difference sc0re
wGWAdiff=weightedGWA(df=do.call(rbind, SNPdiff), dfW=wGWA, model="diff")

# plots
wGWAlogP=forestWeighted(df=wGWAlog)
wGWAdiffP=forestWeighted(df=wGWAdiff) 
wGWAchange=ggarrange(wGWAlogP + guides(shape = "none"), wGWAdiffP + guides(shape = "none"), nrow=1)

# baseline weighted effects
wGWASbl=weightedGWA(df=subset((do.call(rbind, SNPdiff)), label %in% c("cognition_0", "physical_0")), dfW=wGWA, model="baseline")
wGWASblp=forestWeighted(df=wGWASbl,  m="baseline")
wGWAplot=ggarrange(wGWAchange, wGWASblp, nrow=2, heights=c(0.8, 1))

# ==== Over-versus under-estimation 
sIPW_change=sumIPW(df=rbind(subset(wGWAlog, model=="weighted"), subset(wGWAdiff, model=="weighted")), m="change")
sIPW_BL=sumIPW(df=rbind(subset(wGWASbl, model=="weighted"), subset(wGWAdiff, model=="weighted")), m="baseline")
sIPW=rbind(sIPW_change, sIPW_BL)
sIPW$ID=seq(1,NROW(sIPW), 1)
colVar=c(recodeDelta("baseline", return="colmatch"), "darkgreen")
labelsIPW=c(TeX(paste0("$P_0$")), TeX(paste0("$\\Delta_{DIFF/LOG}$")))
sIPW$text=paste0("k=", sIPW$k)
#saveRDS(sIPW, paste0(HOME,"/output/rds/sIPW.rds"))


#
sIPWp <- ggplot(data=sIPW, aes(x=ID, y=meanChange, ymin=lCI, ymax=uCI, colour=model, shape = data, label=text)) +
  geom_text(hjust=1.3, show.legend = FALSE) + 
  geom_pointrange(lwd = 0.4, position = position_dodge(width = 0.5), size=1) + 
  geom_hline(yintercept=0, lty=2) + 
  scale_shape_manual("",  values = c(15,16,17), labels=c("Both", "Overestimation only", "Underestimation only")) +
  theme_minimal() +
  facet_grid(~model , scales="free", space = "free") +
  theme(legend.position="top", 
        legend.text = element_text(size = 15), 
        legend.box = "vertical",
        legend.key.size = unit(0.5, 'cm'),
        legend.direction = "horizontal",
        strip.text.x = element_blank(),
        plot.margin = margin(t=0, r=2, b=0, l=2, "cm"),
        axis.text.x=element_text(size=8, angle=45, hjust=1),
        plot.title = element_text(size = 13,  hjust = 0,  face = "bold"),
        plot.caption = element_text(vjust = 1))  + 
  scale_x_discrete(position = "bottom") +
  scale_colour_manual(name="",values = colVar, labels= labelsIPW) +
  labs(x =  "", 
       caption="",
       y = TeX(paste0("mean( $\\frac{ | \\beta| - | \\beta_w|} {|\\beta|}$ )"))   ) 
sIPWp
# Combine both figure
wGWAplotC=ggarrange(wGWAplot, sIPWp, nrow=2, heights=c(1, 0.5), labels=c("A", "B"))
saveFigure(fileName=paste0(HOME,"/results/figures/wGWAplot"), plotName=wGWAplotC, w=22, h=32 )


# =======================================================================
# =================== Mendelian Randomization ===========================
# =======================================================================
mr=readRDS(paste0(HOME,"/output/rds/mr.rds"))
mr$labelE=recodeChange(df=mr, var="exposure", return="pheno")
mrexp=recodeChange(df=mr, var="labelE")
mr$phenoX=mrexp$label_clean
mr$labelRaw=mrexp$labelRaw
mr$phenoXdim=mrexp$dimension
mr=subset(mr, is.na(phenoXdim)==F & !exposure %in% paste0(fuvar, "_RG"))
mr$phenoXtype=recodeChange(df=mr, var="exposure", return="type")
mr=subset(mr, phenoXtype %in% c("baseline", "other") )
mr$label=recodeChange(df=mr, var="outcome", return="pheno")
mrout=recodeChange(df=mr, var="label")
mr$phenoY=mrout$label_clean
mr$outcomeRaw=mr$outcome
mr$phenoYdim=mrout$dimension
mr$phenoYtype=recodeChange(df=mr, var="outcome", return="type")


mrSupp=mr
mrSum=forestRG(df=mr, beta="b", exposure="phenoX", outcome="phenoY", vars=c("cognition", "physical"), save="full")
saveFigure(fileName=paste0(HOME,"/results/figures/mrSum"), plotName=mrSum, h=45, w=28)

mrSumPhy=forestRG(df=mr, beta="b", exposure="phenoX", outcome="phenoY", vars=physicalVar[physicalVar %in% declineVar], save="suppP")
saveFigure(fileName=paste0(HOME,"/results/figures/mrSuppPhy"), plotName=mrSumPhy, h=45, w=30)

mrSumCog=forestRG(df=mr, beta="b", exposure="phenoX", outcome="phenoY", vars=cognitionVar[cognitionVar %in% declineVar], save="suppC")
saveFigure(fileName=paste0(HOME,"/results/figures/mrSuppCog"), plotName=mrSumCog, h=45, w=40)

# =======================================================================
# ================== Phenotypic correltions ============================
# =======================================================================
# read in results
ukbbLM=readRDS( paste0(HOME,"/output/rds/lmPred.rds"))
# Recode outcome
ukbbLM$typeOut=recodeChange(df=ukbbLM, var="outcome", return="type")
ukbbLM$label=recodeChange(df=ukbbLM, var="outcome", return="pheno")
ukbbLMout=recodeChange(df=ukbbLM, var="label")
ukbbLM$labelOut=ukbbLMout$label_clean

# Recode predictor
ukbbLMexp=recodeChange(df=ukbbLM, var="predictor")
ukbbLM$labelExp=ukbbLMexp$label_clean
ukbbLM$dimExp=ukbbLMexp$dimension

# harmonize column names
colsNew=c("exposure", "outcome", "beta", "lCI", "uCI","pval", "dimExp", "typeOut", "label" )
lmSel=subset(ukbbLM, select=c(labelExp, labelOut, beta, lCI, uCI, pval, dimExp, typeOut, label), typeOut=="log")
colnames(lmSel)=colsNew
lmSel$model="lm"
lmSel$m=paste0(lmSel$exposure, "_", lmSel$label)
mrSel=subset(mr, select=c(phenoX, phenoY, b, lCI, uCI, pval, phenoXdim, phenoYtype, label), label %in% c("cognition", "physical") & phenoYtype=="log")
colnames(mrSel)=colsNew
mrSel$model="mr"
mrSel$m=paste0(mrSel$exposure, "_", mrSel$label)
lmmrWide=merge(lmSel, subset(mrSel, select=c(m, beta, pval, lCI, uCI)), by="m", suffixes = c("_lm", "_mr") )
lmmrWide=na.omit(lmmrWide)

# Scatter plot
scatterRE_COG=scatterRE(df=lmmrWide, select="cognition")
scatterRE_PHY=scatterRE(df=lmmrWide, select="physical")  
rePlotSUM=ggarrange(scatterRE_PHY, scatterRE_COG, nrow=2, common.legend = TRUE)
saveFigure(fileName=paste0(HOME,"/results/figures/rePlot"), plotName=rePlotSUM, h=28, w=17)


# =======================================================================
# ============================ Create table =============================
# =======================================================================
source(paste0(HOME, "/analysis/createTables.R"))


