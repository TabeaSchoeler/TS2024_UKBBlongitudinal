

simMesErr=function(s, N = 10000, rel, alpha, beta, gamma, delta, return="baseline"){
    m=0.1
    sd=1

  set.seed(s)
  G0=rnorm(N, mean=0, sd=1)
  E0=rnorm(N, mean=0, sd=1)
  GE=rnorm(N, mean=0, sd=1)
  E1= delta*E0 + rnorm(N, mean=m, sd=sd)

  # baseline 
  intercept=170
  t0True=intercept + alpha*G0 + beta*E0 + gamma*GE*E0 
  t0Error=rnorm(N, sd = sqrt((1 - rel) * var(t0True))) 
  t0Obs <-  t0True + t0Error
  
  # follow up
  t1True= intercept + alpha*G0 + beta*E1 + gamma*GE*E1 
  t1Error=rnorm(N, sd = sqrt((1 - rel) * var(t1True)))
  t1Obs <- t1True + t1Error
  
  t1TrueE = exp(t1True) 
  t0TrueE = exp(t0True) 
  t1ObsE = exp(t1Obs) 
  t0ObsE = exp(t0Obs)
  
  if(return=="interaction"){
    G=GE
  } else{
    G=G0
  }

  
  # exponential: true model
  trueExpM=summary(lm(log(t1ObsE)-log(t0ObsE) ~  G ))
  trueExp_DIFF_B=trueExpM$coefficients[2,1]
  trueExp_DIFF_SE=trueExpM$coefficients[2,2]
  trueExp_DIFF_P=trueExpM$coefficients[2,4]
  
  # exponential: false model
  falseExpM=summary(lm(log(t1Obs)-log(t0Obs) ~  G ))
  falseExp_DIFF_B=falseExpM$coefficients[2,1]
  falseExp_DIFF_SE=falseExpM$coefficients[2,2]
  falseExp_DIFF_P=falseExpM$coefficients[2,4]
  
  # linear: true model
  trueLinearM=summary(lm(t1Obs-t0Obs ~  G ))
  trueLin_DIFF_B=trueLinearM$coefficients[2,1]
  trueLin_DIFF_SE=trueLinearM$coefficients[2,2]
  trueLin_DIFF_P=trueLinearM$coefficients[2,4]
  
  # linear: false model
  falseLinearM=summary(lm(t1ObsE-t0ObsE ~  G ))
  falseLin_DIFF_B=falseLinearM$coefficients[2,1]
  falseLin_DIFF_SE=falseLinearM$coefficients[2,2]
  falseLin_DIFF_P=falseLinearM$coefficients[2,4]
  
  # residual change: true model
  trueLinearMr=summary(lm(t1Obs-t0Obs ~  G + t0Obs ))
  trueLin_RES_B=trueLinearMr$coefficients[2,1]
  trueLin_RES_SE=trueLinearMr$coefficients[2,2]
  trueLin_RES_P=trueLinearMr$coefficients[2,4]

  # residual change: false model
  falseLinearMr=summary(lm(t1ObsE-t0ObsE ~  G + t0ObsE))
  falseLin_RES_B=falseLinearMr$coefficients[2,1]
  falseLin_RES_SE=falseLinearMr$coefficients[2,2]
  falseLin_RES_P=falseLinearMr$coefficients[2,4]
  
  dfOut=data.frame(alpha, beta, gamma, delta, intercept, N, rel, 
                   trueExp_DIFF_B, trueExp_DIFF_P, trueExp_DIFF_SE,
                   falseExp_DIFF_B, falseExp_DIFF_P, falseExp_DIFF_SE,
                   trueLin_DIFF_B, trueLin_DIFF_P, trueLin_DIFF_SE,
                   trueLin_RES_B, trueLin_RES_P, trueLin_RES_SE,
                   falseLin_DIFF_B, falseLin_DIFF_P, falseLin_DIFF_SE,
                   falseLin_RES_B, falseLin_RES_P, falseLin_RES_SE)

  return(dfOut)
}

rep=1000
alpha=0.1
delta=1
beta=0.5
gamma=0.5
return="interaction"

compModel=function(alpha, beta, gamma, delta, rep, rel, N=10000, return="baseline"){
  print(paste0("G on T1: ", delta))

  simOutL=lapply(1:rep, function(x) simMesErr(s=x, rel=rel, alpha=alpha, delta=delta, beta=beta, gamma=gamma, N=N, return=return))
  simOut=do.call(rbind, simOutL)

  pt=0.05
  signB=mean(sign(simOut$delta))
  

    tLsig=subset(simOut, trueLin_DIFF_P<pt)
    fLsig=subset(simOut, falseLin_DIFF_P<pt)
    tLsigR=subset(simOut, trueLin_RES_P<pt)
    fLsigR=subset(simOut, falseLin_RES_P<pt)
    tEsig=subset(simOut, trueExp_DIFF_P<pt)
    fEsig=subset(simOut, falseExp_DIFF_P<pt)

  pT1_trueLinearDIFF=NROW(tLsig)/NROW(simOut)
  pT1_falseLinearDIFF=NROW(fLsig)/NROW(simOut)
  pT1_trueLinearRES=NROW(tLsigR)/NROW(simOut)
  pT1_falseLinearRES=NROW(fLsigR)/NROW(simOut)
  pT1_trueLOG=NROW(tEsig)/NROW(simOut)
  pT1_falseLOG=NROW(fEsig)/NROW(simOut)

  
  dfOut=data.frame(pT1_trueLinearDIFF,
                   pT1_falseLinearDIFF, 
                   pT1_trueLinearRES,
                   pT1_falseLinearRES, 
                   pT1_falseLOG, 
                   pT1_trueLOG, 
                   e0=mean(simOut$beta), 
                   e1=mean(simOut$gamma), 
                   g0=mean(simOut$alpha),
                   g1=mean(simOut$delta),
                   Ip0=mean(simOut$intercept),
                   #theta=mean(simOut$theta),
                   N=mean(simOut$N),
                   rep=rep,
                   rel=rel)
  
  
  listOut=list(dfOut, simOut)
  names(listOut)=c("sum", "df")
  return(listOut)
  
}




plotModel=function(df, out){
  
  dfP=data.frame(model=c(rep("true", NROW(df)), rep("false", NROW(df))),
                 out=c(df[[paste0("true", out)]], df[[paste0("false", out)]]),
                 type=out,
                 trueSim=  rep(df$gonT1, 2),
                 m=as.character(1:NROW(df)))
  
  
  dfP$model <- factor( dfP$model, levels = c("true", "false"))
  
  labels=c("Correctly specified model", "Incorrectly specified model")
  if(out=="Linear"){
    title="Linear relationship"
  }
  if(out=="Exp"){
    title="Exponential relationship"
  }

  dfP$trueSim=as.numeric(as.character(dfP$trueSim))
  dfP <- dfP[order(dfP$m),]
  dfP$ID=seq(1,NROW(dfP), 1)
  dfP$trueSimF=as.factor(dfP$trueSim)
  
  plot=ggplot(data=dfP, aes(x=trueSimF, y=out, colour=model)) +
    scale_colour_manual("",values = c(  "orange3", "red4" ), labels) +
    #stat_summary(data=dfM, geom = "hpline", width = 0.5, size = 0.1) +
    #geom_vpline(aes(x = trueSim), size = 1.5, height = 0.7, color = "#D55E00") +
    geom_hpline(aes(y = trueSim), colour="grey", stat = "summary", width = 0.6) +
    geom_point(size=3, position=position_dodge(0.5))+
    ylim(-0.05, 0.2) +
    labs( x="Simulated causal effect on change",
          y="Estimated effect on change") +
    ggtitle(title) +
    theme_classic() +
    theme(legend.position="none",
          plot.margin = margin(t=0, r=0, b=0, l=0, "cm")) 
  plot
  
  return(plot)
  
}


simSum=function(alpha, secSim=seq(0.1, 1, 0.1), rep, N=10000, delta, gamma,beta ){
  relL=list()
  
  for ( i in 1:length(secSim) ) {
    print(paste0("Simulate reliability of ", i/10))
    rel=secSim[i]
  
    simOutL=lapply(1:rep, function(x) simMesErr(s=x, rel=rel, alpha=alpha, gamma=gamma, delta=delta, beta=beta, N=N))
    simOut=combRes(listIn=simOutL)

    relL[[i]]=simOut
  }
  
  relSIm=do.call(rbind, relL)
  return(relSIm)
}



simPower=function(rel, cor ){
  corOB = cor * rel # observed correlation
  rel_change = (rel  - corOB) / (1-corOB)
  dfOut=data.frame(rel, cor, corOB, rel_change)
  return(dfOut)
}



combRes=function(listIn){
  simOut=do.call(rbind, listIn)

  trueLin_DIFF_B=mean(simOut$trueLin_DIFF_B, na.rm=T)
  trueLin_RES_B=mean(simOut$trueLin_RES_B, na.rm=T)
  trueExp_DIFF_B=mean(simOut$trueExp_DIFF_B, na.rm=T)

  dfDIFF_abs=data.frame(var="DIFF", snpChange_o=trueLin_DIFF_B) 
  dfRES_abs=data.frame(var="RES", snpChange_o=trueLin_RES_B) 
  dfDIFF_rel=data.frame(var="LOG", snpChange_o=trueExp_DIFF_B) 
  
  dfOut=rbind(dfDIFF_abs, dfRES_abs, dfDIFF_rel)
  
  dfOut$alpha=simOut$alpha[1]
  dfOut$beta=simOut$beta[1]
  dfOut$gamma=simOut$gamma[1]
  dfOut$delta=simOut$delta[1]
  dfOut$rel=simOut$rel[1]
  dfOut$intercept=simOut$intercept[1]
  dfOut$expI=simOut$expI[1]
  dfOut$N=simOut$N[1]
  
  
  return(as.data.frame(dfOut))
}



          
plotRel=function(df, max=0.12, min=-0.05, maxX=0.5){
  maxX=max(df$alpha)
  
  df$var_c=recodeDelta(var=tolower(df$var), return="factor")
  cols=recodeDelta(var=df$var_c, return="colour")
  df$size=recodeDelta(var=  df$var_c, return="size")
  
  snpCplot= ggplot(df, aes(y=snpChange_o, x=alpha, color=var_c)) +
    geom_point(size= df$size)+
    theme_classic() +
    scale_colour_manual("",values = cols) + 
    ylim(min,max) + 
    xlim(0,maxX) +
    theme(legend.position="none",
          plot.margin = margin(t=0, r=0, b=0, l=0, "cm")) +
    labs( y="", 
          x="") +
    ggtitle(TeX(paste0("$\\rho=$", df$rel[1]))) 
  print(snpCplot)
  return(snpCplot)
}


plotSim=function(df, max, min, title){
  alpha=df$alpha
  df$var_c=recodeDelta(var=tolower(df$var), return="factor")
  
  cols=recodeDelta(var=df$var_c, return="colour")
  reorderVar=recodeDelta(var=  df$var_c, return="label")
  df$size=recodeDelta(var=  df$var_c, return="size")
  
  plot=ggplot(data=df, aes(x=rel, y=snpChange_o, colour=var_c)) +
    scale_colour_manual("",values = cols, labels=reorderVar) +
    geom_hline(yintercept = alpha, colour = "black", linetype="dotted") +
    geom_point(size=df$size)+
    scale_x_reverse() +
    annotate("text", x=0.30, y=alpha+0.005, 
             label=  expression(paste("Simulated baseline genetic effect (", alpha, ") on ", P[t])), 
             colour = gray(1/2), size=2) +
    ylim(min, max) +
    labs(y="", 
          x="") +
    theme_classic() +
    theme(legend.position="top",
          legend.text = element_text(size = 20), 
          plot.margin = margin(t=0.5, r=0, b=0, l=0, "cm")) +
    guides(color = guide_legend(override.aes = list(size=10))) +
    ggtitle(title) 
  print(plot)
  return(plot)
  
}

# Define a function to pivot the data
pivot_data <- function(df, cols) {
  df %>%
    pivot_longer(
      cols = cols,
      names_to = "out",
      values_to = "value",
      values_transform = ~ as.numeric(gsub(",", "", .x))
    )
}



plotE=function(df){
  library(latex2exp)
  library(cowplot)
  colsL <- c("pT1_trueLinearDIFF", "pT1_falseLinearDIFF", "pT1_trueLinearRES", "pT1_falseLinearRES", "pT1_trueLOG", "pT1_falseLOG")

  # Apply the function to both sets of columns
  dfL <- as.data.frame(pivot_data(df, colsL))
  dfL$modelFunc =  factor(   dfL$out, levels = colsL)
  
  #labels=c(expression(paste("Correctly specified linear model (", P[1]^"*" - P[0]^"*" == pi * G[0] + epsilon, ")")),
  #         expression(paste("Incorrectly specified linear model (", e^P[1]^"*" - e^P[0]^"*" == pi * G[0] + epsilon, ")")),
  #         expression(paste("Correctly specified linear model (", P[1]^"*" - P[0]^"*" == pi * G[0] + P[0]^"*"+ epsilon, ")")),
  #         expression(paste("Incorrectly specified linear model (", e^P[1]^"*" - e^P[0]^"*" == pi * G[0] + e^P[0]^"*" + epsilon, ")")),
  #         expression(paste("Correctly specified exponential model (", Log(e^P[1]^"*") - Log(e^P[0]^"*") == pi * G[0] + epsilon, ")")),
  #         expression(paste("Incorrectly specified exponential model (", Log(P[1]^"*") - Log(P[0]^"*") == pi * G[0] + epsilon, ")")))
  
  labels=c("Correctly specified linear model", 
           "Incorrectly specified linear model",
           "Correctly specified linear model",
           "Incorrectly specified linear model",
           "Correctly specified exponential model",
           "Incorrectly specified exponential model")
  
  #dfL$trueSimF=label
  dfL$value= dfL$value*100

  cols=c( "blue3", "darkblue" , "cornsilk4", "gray43", "orange" ,"orange3")
  
  p<-ggplot(data=dfL, 
            aes(x=type, y=value, fill=modelFunc, label = value, group=modelFunc)) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_fill_manual("",values = cols, labels=labels) +
    labs( x= "",
          y=expression(paste( "Percentage (%) of ", italic(p)<0.05 ))) +
    theme_classic() +
    ylim(0,100) +
    theme(legend.position="bottom",
          plot.margin = margin(t=0.5, r=0, b=0, l=0.5, "cm")) +
    geom_text(vjust = 1, position = position_dodge(.9), colour="white") +
    guides(fill = guide_legend(nrow = 6)) +
    geom_abline(intercept = 5, slope = 0, colour="grey", linetype="dotted") 
  p 
  
  return(p)
  
}




plotB=function(df){

  dfS=df
  b=c("trueLin_DIFF", "falseLin_DIFF","trueLin_RES", "falseLin_RES","trueExp_DIFF" ,"falseExp_DIFF")
  dfO <- lapply(b, function(x) orderDat(df=dfS, var=x))
  dfP=do.call(rbind, dfO)
  dfP$modelR=ifelse(grepl("true",   dfP$model)==T, "true",   "false")
  dfP$modelR=  factor(   dfP$modelR, levels = c("true", "false"))
  head(dfP)
  
  dfP$pMC=str_remove(dfP$model, "true")
  dfP$pMC=str_remove(dfP$pMC, "false")

  cols=c(  "blue3" , "darkblue", "cornsilk4", "gray43", "orange" ,"orange3")
  
  dfP$pM =  factor(   dfP$model, levels = paste0(b, "_B"))
  
  dfP$beta=dfP$beta/dfP$SE
  dfPs=subset(dfP, p=="sig")
  
  p=ggplot(dfP, aes(x=pM, y=beta, colour=pM)) +
    scale_colour_manual("", values = cols) +
    geom_violin() +
    geom_boxplot(width=0.01, outlier.shape = NA) +
    geom_jitter(data=dfPs, shape=16, position=position_jitter(0.2), size=0.2, alpha=0.2) +
    stat_summary(fun.data=mean_sdl,
                   geom="pointrange") +
    labs( x="",
          y=expression(paste("Estimated baseline genetic effects on ", Delta, " (",widehat(pi)/SE(widehat(pi)), ")" ))) +
    #facet_wrap(model ~ modelR, scales = "free",  nrow = 2) +
    facet_wrap(modelR ~ pMC, scales = "free",  nrow = 2) +
    ylim(-5, 5) +
    theme(legend.position="none",
          plot.margin = margin(t=0.5, r=0.5, b=0, l=0.5, "cm"),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          panel.grid.major.y = element_line( size=.1, color="darkgrey" ) ,
          panel.background = element_rect(fill = "white",
                                          colour = "white",
                                          size = 0.5, linetype = "solid"),
          panel.grid.major.x = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
    p
  return( p)
  
}


orderDat=function(var, df){
  print(var)
  dfOut=data.frame(beta=df[[paste0(var, "_B")]])
  dfOut$model=paste0(var, "_B")
  dfOut$p=ifelse(df[[paste0(var, "_P")]]<0.05, "sig", "ns")
  dfOut$SE=df[[paste0(var, "_SE")]]
  return(dfOut)
}


plotM=function(df, lim=1, seqS, title){
  
  df$diffL=round(ifelse((df$x %in% seqS)==T, df$diff, NA ),2)
  df$p=ifelse((df$x %in% seqS)==T, df$x, NA )
  intercept=signif(subset(df, x==0 & g==1)$y, digits = 1)

  plot=ggplot(df, aes(x = x, y = y,  color =as.factor( g), label=diffL)) + #
    geom_abline(intercept = 0, slope = 1, color="grey", linetype="solid", size=.5) +
    geom_line(size=2) +
    labs(x =expression(paste(P[0] )), y = expression(paste(P[1] )), title = title) +
    theme_minimal() +
    scale_color_manual("Genotype", values = c( "forestgreen", "lightseagreen"))  + # Set color for growth and decay
    scale_x_continuous(breaks = round(seq(min(df$x), max(df$x), by = 0.2),1)) +
    scale_y_continuous(breaks = round(seq(min(df$y), max(df$y), by = 0.2),1)) +
    theme(legend.position = "top",
          plot.title = element_text( size=16, hjust=0.5),
          axis.title.y = element_text(angle = 360, vjust = 0.5, size=16),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.title.x = element_text(size=16)) 
 
  print(plot)
  return(plot)
  
}


