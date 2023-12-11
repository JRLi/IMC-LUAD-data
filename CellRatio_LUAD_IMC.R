#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Down-sampling
rm(list=ls())
library(survival)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(ggtext) 
library(glmnet)
load('IMC_LUAD.RData')

pnn=100 # down-sampling times
q.thr = 0.05 # FDR threshold
d.size = 0.8 # 80%
Fig1 = 'Frequency_downsampling_samples.pdf'

# Clinical
time = info$Survival.time
event = info$Progression
info = cbind(time, event, info)
comxx = intersect(row.names(info), row.names(data))
data = data[comxx,]
info = info[comxx,]
raw.data = data
# Fraction
imm = data[,grepl('FRAC.ALL_\\$_',colnames(data))]
colnames(imm) = gsub('FRAC.ALL_\\$_','',colnames(imm))
# Fraction IMM normalized
imm2 = data[,grepl('FRAC.IMM_\\$_',colnames(data))]
colnames(imm2) = gsub('FRAC.IMM_\\$_','',colnames(imm2))
# Density
den = data[,grepl('INTEN_\\$_',colnames(data))]
colnames(den) = gsub('INTEN_\\$_','',colnames(den))
# Ratio
rat = data[,grepl('RATIO_\\$_',colnames(data))]
colnames(rat) = gsub('RATIO_\\$_','',colnames(rat))
write.csv(imm, 'Cell_Fraction.csv',quote = F)
write.csv(imm2, 'Cell_Fraction_rescaled_after_removed_cancer_endothelial.csv',quote = F)
write.csv(den, 'Cell_Density.csv',quote = F)
write.csv(rat, 'CellCell_Ratio.csv',quote = F)
raw.info = info # clinical
raw.rat = rat # ratio
raw.imm = imm # proportion
raw.den = den # density ####

rat.List = imm.List = den.List = list(NULL) ####

for(p in 1:pnn)
{
  cat("\r", p)
  rat = raw.rat # ratio-pair: 240
  se = sample(1:nrow(rat), nrow(rat)*d.size)
  rat = rat[se,]
  info = raw.info[se,]
  imm = raw.imm[se,]
  den = raw.den[se,]  #####
  PV1 = PV2 =  rep(0, ncol(rat))
  HR1 = HR2 =  rep(0, ncol(rat))
  PV3 = PV4 =  rep(0, ncol(imm))
  HR3 = HR4 =  rep(0, ncol(imm))
  PV5 = PV6 =  rep(0, ncol(den))  ###
  HR5 = HR6 =  rep(0, ncol(den))  ###
  
  for(k in 1:ncol(rat))
  {
    mytf = as.numeric(rat[,k])
    xx = cbind(mytf, info)
    xx = xx[xx[, "time"]>0,]
    mycox = coxph(Surv(time, event)~mytf, xx) 
    mycox = summary(mycox)
    PV1[k] = mycox$coefficients[5]
    tmp = mycox$conf.int
    HR1[k] = tmp[1]
    mycox = coxph(Surv(time, event)~mytf+Age.75+Sex.F+Smoking.Status.non1+Stage.Late1, xx) 
    mycox = summary(mycox)
    PV2[k] = mycox$coefficients["mytf",5]
    tmp = mycox$conf.int
    HR2[k] = tmp["mytf", 1]
  }
  tmp = data.frame(HR1, PV1, HR2, PV2)
  row.names(tmp) = colnames(rat)
  rat.List[[p]] = tmp
  
  for(k in 1:ncol(imm))
  {
    mytf = as.numeric(imm[,k])
    xx = cbind(mytf, info)
    xx = xx[xx[, "time"]>0,]
    mycox = coxph(Surv(time, event)~mytf, xx) 
    mycox = summary(mycox)
    PV3[k] = mycox$coefficients[5]
    tmp = mycox$conf.int
    HR3[k] = tmp[1]
    mycox = coxph(Surv(time, event)~mytf+Age.75+Sex.F+Smoking.Status.non1+Stage.Late1, xx) 
    mycox = summary(mycox)
    PV4[k] = mycox$coefficients["mytf",5]
    tmp = mycox$conf.int
    HR4[k] = tmp["mytf", 1]
  }
  tmp = data.frame(HR3, PV3, HR4, PV4)
  row.names(tmp) = colnames(imm)
  imm.List[[p]] = tmp
  
  for(k in 1:ncol(den)) ###
  {
    mytf = as.numeric(den[,k])
    xx = cbind(mytf, info)
    xx = xx[xx[, "time"]>0,]
    mycox = coxph(Surv(time, event)~mytf, xx) 
    mycox = summary(mycox)
    PV5[k] = mycox$coefficients[5]
    tmp = mycox$conf.int
    HR5[k] = tmp[1]
    mycox = coxph(Surv(time, event)~mytf+Age.75+Sex.F+Smoking.Status.non1+Stage.Late1, xx) 
    mycox = summary(mycox)
    PV6[k] = mycox$coefficients["mytf",5]
    tmp = mycox$conf.int
    HR6[k] = tmp["mytf", 1]
  } ###
  tmp = data.frame(HR5, PV5, HR6, PV6) ##
  row.names(tmp) = colnames(den) ##
  den.List[[p]] = tmp ##
}

#---------------------------- # Density ###
nsig5 = nsig6 = rep(0, pnn)
for(k in 1:pnn)
{
  xx = den.List[[k]]
  qval = p.adjust(xx$PV5, method = "BH")
  nsig5[k] = sum(qval<=q.thr)
  qval = p.adjust(xx$PV6, method = "BH")
  nsig6[k] = sum(qval<=q.thr)
}
sum(nsig5>0) # Univariate
sum(nsig6>0) # Multivariate

mean(nsig5)
mean(nsig6) # 0.11

#---------------------------- # proportion
nsig3 = nsig4 = rep(0, pnn)
for(k in 1:pnn)
{
  xx = imm.List[[k]]
  qval = p.adjust(xx$PV3, method = "BH")
  nsig3[k] = sum(qval<=q.thr)
  qval = p.adjust(xx$PV4, method = "BH")
  nsig4[k] = sum(qval<=q.thr)
}
sum(nsig3>0)
sum(nsig4>0)	

mean(nsig3)
mean(nsig4)

#---------------------------- # ratio
nsig1 = nsig2 = rep(0, pnn)
for(k in 1:pnn)
{
  xx = rat.List[[k]]
  qval = p.adjust(xx$PV1, method = "BH")
  nsig1[k] = sum(qval<=q.thr)
  qval = p.adjust(xx$PV2, method = "BH")
  nsig2[k] = sum(qval<=q.thr)
}
sum(nsig1>0)
sum(nsig2>0)

mean(nsig1)
mean(nsig2)


# Plot
df1 = as.data.frame(table(nsig2))
df2 = as.data.frame(table(nsig4))
df3 = as.data.frame(table(nsig6))
df1$Type = 'Cell Ratio'
df2$Type = 'Cell Fraction'
df3$Type = 'Cell Density'
colnames(df1)[1] = colnames(df2)[1] =colnames(df3)[1] ='SigCell'

df = rbind(df1,df2,df3)
df$SigCell = factor(df$SigCell, levels = sort(as.numeric(levels(df$SigCell))))
cols2 =c('#f94040', '#437fc7', '#00c000')

pdf(Fig1, width=5, height=3.5)
ggplot(df, aes(x=SigCell,y=Freq, fill = Type)) + 
  geom_bar(stat='identity', color='black',position = 'identity', alpha = 0.6, width=0.8)+
  scale_y_continuous(expand = expansion(mult = c(0.0, 0.0)))+
  scale_fill_manual(values = cols2)+
  theme_classic()+
  labs(y = "Frequency",
       title = "Down-sampling analysis of samples")+
  theme(axis.text=element_text(size=10, colour = 'black'),
        axis.title.x = element_blank(),
        legend.position = c(0.7,0.8),
        legend.box.spacing = unit(0, "pt"),
        legend.text = element_text(size=13.5),
        legend.title = element_text(size=14, face="bold"),
        plot.title = element_text(size=16, face="bold", hjust = 0.5), 
        axis.title=element_text(size=16,face="bold"))
dev.off()



#-------------------------------------------------------------------------------
# Survival analysis
load('IMC_LUAD.RData')
raw.data = data
#-----------------------------------
se = grep("FRAC.ALL", colnames(raw.data))
data = raw.data[,se]
colnames(data) = gsub("FRAC.ALL_\\$_", "", colnames(data))
data = data*100
coxph.pval1 = coxph.pval2 = coxph.pval3 = rep(0, ncol(data))
hr1 = hr2 = hr3 = rep(0, ncol(data))
for(k in 1:ncol(data))
{
  cat("\r", k)
  mytf = as.numeric(data[,k])
  mytag = ifelse(mytf>=median(mytf, na.rm=T), 1, 0)
  xx = cbind(mytf, mytag, info)
  xx = xx[xx[, "time"]>0,]
  mycox = coxph(Surv(time, event)~mytag, xx) 
  mycox = summary(mycox)
  coxph.pval1[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hr1[k] = tmp[1]
  mycox = coxph(Surv(time, event)~mytf, xx) 
  mycox = summary(mycox)
  coxph.pval2[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hr2[k] = tmp[1]
  mycox = coxph(Surv(time, event)~mytf+Age.75+Sex.F+Smoking.Status.non1+Stage.Late1, xx) 
  mycox = summary(mycox)
  coxph.pval3[k] = mycox$coefficients["mytf",5]
  tmp = mycox$conf.int
  hr3[k] = tmp["mytf", 1]
}
name = colnames(data)
qval2 = p.adjust(coxph.pval2, method="BH")
qval3 = p.adjust(coxph.pval3, method="BH")
res = data.frame(name, uHR=hr2, uPV = coxph.pval2, uQV= qval2, mHR=hr3, mPV=coxph.pval3, mQV= qval3)
res = res[order(res$uPV), ]
res.all = res

#-----------------------------------
se = grep("FRAC.IMM", colnames(raw.data))
data = raw.data[,se]
colnames(data) = gsub("FRAC.IMM_\\$_", "", colnames(data))
data = data*100
coxph.pval1 = coxph.pval2 = coxph.pval3 = rep(0, ncol(data))
hr1 = hr2 = hr3 = rep(0, ncol(data))
for(k in 1:ncol(data))
{
  cat("\r", k)
  mytf = as.numeric(data[,k])
  mytag = ifelse(mytf>=median(mytf, na.rm=T), 1, 0)
  xx = cbind(mytf, mytag, info)
  xx = xx[xx[, "time"]>0,]
  mycox = coxph(Surv(time, event)~mytag, xx) 
  mycox = summary(mycox)
  coxph.pval1[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hr1[k] = tmp[1]
  mycox = coxph(Surv(time, event)~mytf, xx) 
  mycox = summary(mycox)
  coxph.pval2[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hr2[k] = tmp[1]
  mycox = coxph(Surv(time, event)~mytf+Age.75+Sex.F+Smoking.Status.non1+Stage.Late1, xx) 
  mycox = summary(mycox)
  coxph.pval3[k] = mycox$coefficients["mytf",5]
  tmp = mycox$conf.int
  hr3[k] = tmp["mytf", 1]
}
name = colnames(data)
qval2 = p.adjust(coxph.pval2, method="BH")
qval3 = p.adjust(coxph.pval3, method="BH")
res = data.frame(name, uHR=hr2, uPV = coxph.pval2, uQV= qval2, mHR=hr3, mPV=coxph.pval3, mQV= qval3)
res = res[order(res$uPV), ]
res.imm = res

#-----------------------------------
se = grep("INTEN", colnames(raw.data))
data = raw.data[,se]
colnames(data) = gsub("INTEN_\\$_", "", colnames(data))

coxph.pval1 = coxph.pval2 = coxph.pval3 = rep(0, ncol(data))
hr1 = hr2 = hr3 = rep(0, ncol(data))
for(k in 1:ncol(data))
{
  cat("\r", k)
  mytf = as.numeric(data[,k])
  mytag = ifelse(mytf>=median(mytf, na.rm=T), 1, 0)
  xx = cbind(mytf, mytag, info)
  xx = xx[xx[, "time"]>0,]
  mycox = coxph(Surv(time, event)~mytag, xx) 
  mycox = summary(mycox)
  coxph.pval1[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hr1[k] = tmp[1]
  mycox = coxph(Surv(time, event)~mytf, xx) 
  mycox = summary(mycox)
  coxph.pval2[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hr2[k] = tmp[1]
  mycox = coxph(Surv(time, event)~mytf+Age.75+Sex.F+Smoking.Status.non1+Stage.Late1, xx) 
  mycox = summary(mycox)
  coxph.pval3[k] = mycox$coefficients["mytf",5]
  tmp = mycox$conf.int
  hr3[k] = tmp["mytf", 1]
}
name = colnames(data)
qval2 = p.adjust(coxph.pval2, method="BH")
qval3 = p.adjust(coxph.pval3, method="BH")
res = data.frame(name, uHR=hr2, uPV = coxph.pval2, uQV= qval2, mHR=hr3, mPV=coxph.pval3, mQV= qval3)
res = res[order(res$uPV), ]
res.int = res

#-----------------------------------
se = grep("RATIO", colnames(raw.data))
data = raw.data[,se]
colnames(data) = gsub("RATIO_\\$_", "", colnames(data))
myavg = apply(data,2, mean, na.rm=T)
mystd = apply(data,2, sd, na.rm=T)
for(k in 1:ncol(data))
{
  data[,k] = (data[,k]-myavg[k])/mystd[k]
}
coxph.pval1 = coxph.pval2 = coxph.pval3 = rep(0, ncol(data))
hr1 = hr2 = hr3 = rep(0, ncol(data))
for(k in 1:ncol(data))
{
  cat("\r", k)
  mytf = as.numeric(data[,k])
  mytag = ifelse(mytf>=median(mytf, na.rm=T), 1, 0)
  xx = cbind(mytf, mytag, info)
  xx = xx[xx[, "time"]>0,]
  mycox = coxph(Surv(time, event)~mytag, xx) 
  mycox = summary(mycox)
  coxph.pval1[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hr1[k] = tmp[1]
  mycox = coxph(Surv(time, event)~mytf, xx) 
  mycox = summary(mycox)
  coxph.pval2[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hr2[k] = tmp[1]
  mycox = coxph(Surv(time, event)~mytf+Age.75+Sex.F+Smoking.Status.non1+Stage.Late1, xx) 
  mycox = summary(mycox)
  coxph.pval3[k] = mycox$coefficients["mytf",5]
  tmp = mycox$conf.int
  hr3[k] = tmp["mytf", 1]
}
name = colnames(data)
qval2 = p.adjust(coxph.pval2, method="BH")
qval3 = p.adjust(coxph.pval3, method="BH")
res = data.frame(name, uHR=hr2, uPV = coxph.pval2, uQV= qval2, mHR=hr3, mPV=coxph.pval3, mQV= qval3)
res = res[order(res$uPV), ]
res.rat = res

# Survival analysis result of fraction of all cell types
res.all
# Survival analysis result of fraction of all 14 re-scaled immune cell types 
res.imm
# Survival analysis result of density of all cell types
res.int
# Survival analysis result of cell-cell ratio of all cell-cell pairs
res.rat

