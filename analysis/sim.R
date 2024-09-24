
#######################################################
# =================== FUNCTIONS =======================
#######################################################
rm(list = ls())
library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..")
HOME=getwd()
source(paste0(HOME, "/analysis/funcLongitudinal.R"))
source(paste0(HOME, "/analysis/simFunc.R"))

library(lm.beta)
library(patchwork)


# =============== RELIABILITY OF CHANGE SCORES ========================
p1L=lapply(seq(0, 1, 0.05), function(x) simPower(rel=x, cor=0))
p3L=lapply(seq(0, 1, 0.05), function(x) simPower(rel=x, cor=0.3))
p6L=lapply(seq(0, 1, 0.05), function(x) simPower(rel=x, cor=0.6))
p9L=lapply(seq(0, 1, 0.05), function(x) simPower(rel=x, cor=0.9))

pDF=rbind(do.call(rbind, p1L), do.call(rbind, p3L), do.call(rbind, p6L), do.call(rbind, p9L))
pDF$cor=as.factor(pDF$cor)

changeRplot= ggplot(pDF, aes(x=rel, y=rel_change, colour=cor)) +
  geom_abline(intercept = 0, slope = 1, colour="lightgrey") +
  geom_point(size=3) + 
  geom_line() +
  theme_minimal() +
  xlim(0,1) +
  ylim(0,1) +
  scale_colour_manual(expression(paste(r(P["t0"],P["t1"]))) ,values = c("steelblue", "darkblue","red4", "orange3")) +
  theme(legend.position="top") +
  labs( y=expression(paste("Difference score (D) reliability (", rho, ")")), 
        x=expression(paste("Phenotype (P) reliability (", rho, ")"))) 

saveFigure(fileName=paste0(HOME,"/results/figures/changeRplot"), plotName=changeRplot, h=10, w=15)

# =============== BIAS IN CHANGE SCORES ========================
textS=1
nrep=10000; gammaS=0.5; betaS=1; deltaS=1
#gammaS=0.5; betaS=0.5; deltaS=1
s1=simSum(rep=nrep, alpha=0.00, gamma=gammaS, beta=betaS, delta=deltaS)
s2=simSum(rep=nrep, alpha=0.05, gamma=gammaS, beta=betaS, delta=deltaS) 
s3=simSum(rep=nrep, alpha=0.1, gamma=gammaS, beta=betaS, delta=deltaS) 

mDF=rbind(s1, s2, s3)
max=max(c(mDF$snpChange_o, mDF$alpha), na.rm=T)
min=min(c(mDF$snpChange_o, mDF$alpha), na.rm=T)

relP1=plotSim(df=s1, max=max, min=min, title=TeX(paste0("$\\alpha = 0$")))
relP2=plotSim(df=s2, max=max, min=min, title=TeX(paste0("$\\alpha = 0.05$")))
relP3=plotSim(df=s3, max=max, min=min, title=TeX(paste0("$\\alpha = 0.1$")))

relSImPlot=ggarrange( relP1, relP2,relP3, 
                      nrow=1, 
                      common.legend = TRUE, 
                      legend="top", 
                      font.label = list(size = 12),
                      labels=c("A1", "A2", "A3"))
relSImPlotA=annotate_figure(relSImPlot, bottom = textGrob(expression(paste("Measurement reliability (", rho, ")")), gp = gpar(cex = textS)))


# Get estimates across different levels of causal effects
relSImC=lapply(seq(0,0.2,0.01), function(x) simSum(alpha=x, secSim=c(0.3, 0.6, 1), rep=nrep, beta=betaS, delta=deltaS, gamma=gammaS))
relSC=do.call(rbind, relSImC)

snpCplot1=plotRel(df=subset(relSC, rel==1)) 
snpCplot2=plotRel(df=subset(relSC, rel==0.6))
snpCplot3=plotRel(df=subset(relSC, rel==0.3))
snpCplot=ggarrange(snpCplot1,snpCplot2,snpCplot3, nrow=1, font.label = list(size = 12), labels=c("B1", "B2", "B3"))

snpCplotA=annotate_figure(snpCplot, bottom = textGrob(expression(paste("Simulated baseline genetic effect (", alpha, ") on ", P[t])), gp = gpar(cex = textS)))

# Combine plots
emptyP=cowplot::ggdraw() +
  cowplot::draw_label('')
simPlot=ggarrange(relSImPlotA, emptyP, snpCplotA, ncol=1, common.legend = TRUE, legend="top", heights=c(1,0.05,0.8))
simPlotA=annotate_figure(simPlot, left = textGrob(expression(paste("Estimated baseline genetic effects on ", Delta , " (", widehat(pi), ")")), rot = 90, vjust = 1, gp = gpar(cex = 1.5)))


# =============== risk of bias when using wrong model========================
# false positives
nullM=lapply(0.1, function(x) compModel(gamma=gammaS, beta=betaS, alpha=x, delta=deltaS, rel=0.8, N = 10000, rep=nrep)) #seq(0,0.1,0.05)
mN=do.call(rbind, lapply(nullM, function(x) x$sum)) 
mN$type="(false positive rate)"
# power of detecting GxE
powerM=lapply(0.1, function(x) compModel(alpha=x, beta=betaS, gamma=gammaS, delta=deltaS, rel=0.8, N = 10000, rep=nrep, return="interaction")) #seq(0,0.1,0.05)
mP=do.call(rbind, lapply(powerM, function(x) x$sum)) 
mP$type="(power)"

mS=rbind(mN, mP)
compModelN <-plotE(df=mS) # plot false positives + power
#compModelP <-plotE(df=mP, label="(power)")
compBN=plotB(df=do.call(rbind, lapply(nullM, function(x) x$df))) +  ggtitle(TeX(paste0("$\\alpha = 0.1$")))

compModelN <- compModelN + guides(fill = guide_legend(nrow = 2))
# plot false positives beta estimates
compMSlop=ggarrange(compModelN, compBN, nrow=1, labels=c("C", "D"), common.legend = T, legend="bottom", widths=c(0.9,1.2))



# Combine all plots
simPlotComb=ggarrange(simPlotA, emptyP, compMSlop, ncol=1, heights=c(1, 0.05,0.7))


# =======================================================================
# ============= Weighted correlates of change  ==========================
# =======================================================================
wCor=readRDS( paste0(HOME,"/output/rds/wCor.rds"))
wCor$outcomeD=recodeChange(df=wCor, var="outcome", return="pheno")
wCor$type=recodeChange(df=wCor, var="outcome", return="type")
wCor$typeC=droplevels(wCor$type)
wCor$typeC=recodeDelta(wCor$typeC, return="factor")
wCorL=recodeChange(df=wCor, var="outcomeD")
wCor$outcomeC=wCorL$label_clean
wCor$exposureD=recodeChange(df=wCor, var="exposure", return="pheno")
wCorP=recodeChange(df=wCor, var="exposureD")
wCor$exposureC=wCorP$label_clean
wCor$model <- factor(wCor$model, levels = c("none", "baselineW", "fuW"))
wCor$dimension=wCorL$dimension
wCor=subset(wCor, model!="baselineW")


# Plot for cognitive and physical decline
wCorP=subset(wCor, outcomeD %in% c("physical", "cognition"))
wCorP$outcomeC=recodeTriang(wCorP$outcomeD)
reorderVar=recodeDelta(return="label")

wCorPlot <- ggplot(wCorP, aes(x = reorder(exposureC, -as.numeric(beta)), 
                              y = beta, 
                              ymin = lCI, 
                              ymax = uCI,
                              col=typeC,
                              shape=model)) +
  scale_colour_manual("",values = recodeDelta(return="colour"), labels=reorderVar ) +
  scale_shape_manual("" , values=c(21,22,23), labels=c("No weighting", "Follow-up participation weighted") ) +
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2)+ 
  geom_linerange(aes( ymin = lCI,
                      ymax = uCI),  lwd = .2, position = position_dodge(width = 0.5)) + 
  geom_pointrange(aes(ymin = lCI,
                      ymax = uCI), lwd = 0.4, position = position_dodge(width = 0.5)) + 
  theme_minimal() +
  theme(legend.position="top", 
        #legend.direction = "vertical",
        axis.text.x=element_text(size=10, angle=45, hjust=1, face="bold"),
        plot.title = element_text(size = 13, face = "bold"), # hjust = -1, 
        plot.caption = element_text(vjust = 1),
        axis.title.y = element_text(angle = 0, vjust = 0.5,  face = "bold"),
        plot.margin = margin(t=1, r=4, b=0, l=4, "cm"))  + 
  scale_x_discrete(position = "bottom") +
  facet_wrap(vars(outcomeC), scales = "free",  nrow = 1) +
  labs(x =  "", 
       caption="",
       #subtitle="E",
       y = expression(paste( italic(r) ))) + guides(colour = "none") +
  ggtitle("E")

wCorPlot

simPlotEmp=ggarrange(simPlotComb, emptyP, wCorPlot, ncol=1, heights=c(1,0.03,0.4), 
                     font.label = list(size = 18, color = "black", face="italic"),
                     labels=c("Simulation results", "", "Real data application"))

saveFigure(fileName=paste0(HOME,"/results/figures/simPlot"), plotName=simPlotEmp, h=45, w=30)








