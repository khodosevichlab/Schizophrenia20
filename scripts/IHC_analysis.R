### Mixed effects models for analysing IHC data
### LMM necessary because there were multiple measurements from the same sample - nonindependent data



library(multcomp)
library(lme4)
library(nlme)
library(ggthemes)
library(ggplot2)

## data table for CR, PV
sumd <- read.table("alldatag.csv", header=TRUE, sep=",", dec=".")
str(sumd)
sumd$diag   <- as.factor(sumd$diag)
sumd$gender <- as.factor(sumd$gender)
#sumd$gendernum <- as.factor(sumd$gendernum)
#sumd$gendernum <- as.numeric(sumd$gendernum)
sumd$Layer  <- as.factor(sumd$Layer)
str(sumd)

#subsets 
  L1<-sumd[sumd$Layer=="1",]
  L2<-sumd[sumd$Layer=="2",]
  L3<-sumd[sumd$Layer=="3",]
  L4<-sumd[sumd$Layer=="4",]
  L5<-sumd[sumd$Layer=="5",]
  L6<-sumd[sumd$Layer=="6",]
    c<-sumd[sumd$diag=="ctr",]
    s <- sumd[sumd$diag=="sch",]
    c1<-c[c$Layer=="1",]
    c2<-c[c$Layer=="2",]
    c3<-c[c$Layer=="3",]
    c4<-c[c$Layer=="4",]
    c5<-c[c$Layer=="5",]
    c6<-c[c$Layer=="6",]
    s1<-s[s$Layer=="1",]
    s2<-s[s$Layer=="2",]
    s3<-s[s$Layer=="3",]
    s4<-s[s$Layer=="4",]
    s5<-s[s$Layer=="5",]
    s6<-s[s$Layer=="6",]

# data table for CB, NPY, SMI, NISSL densities
sum6 <- read.table("6sch6ctrlayerwiseg.csv", header=TRUE, sep=",", dec=".")
str(sum6)
sum6$diag   <- as.factor(sum6$diag)
sum6$gender <- as.factor(sum6$gender)
sum6$Layer  <- as.factor(sum6$Layer)
str(sum6)
  
  L1sum <-sum6[sum6$Layer=="1",]
  L2sum <-sum6[sum6$Layer=="2",]
  L3sum <-sum6[sum6$Layer=="3",]
  L4sum <-sum6[sum6$Layer=="4",]
  L5sum <-sum6[sum6$Layer=="5",]
  L6sum <-sum6[sum6$Layer=="6",]
    csum <-sum6[sum6$diag=="ctr",]
    ssum <- sum6[sum6$diag=="sch",]
      c1sum <-csum[csum$Layer=="1",]
      c2sum <-csum[csum$Layer=="2",]
      c3sum <-csum[csum$Layer=="3",]
      c4sum <-csum[csum$Layer=="4",]
      c5sum <-csum[csum$Layer=="5",]
      c6sum <-csum[csum$Layer=="6",]
      s1sum <-ssum[ssum$Layer=="1",]
      s2sum <-ssum[ssum$Layer=="2",]
      s3sum <-ssum[ssum$Layer=="3",]
      s4sum <-ssum[ssum$Layer=="4",]
      s5sum <-ssum[ssum$Layer=="5",]
      s6sum <-ssum[ssum$Layer=="6",]

## PARVALBUMIN

parv = lme(density ~ diag*Layer, 
           random = ~1|ID, 
           data=sumd, 
           na.action = na.exclude)
summary(parv)
qqnorm(resid(parv))
qqline(resid(parv))
reff = ranef(parv)
qqnorm(reff[["(Intercept)"]])
qqline(reff[["(Intercept)"]])
hist(reff[["(Intercept)"]]) 

#L2
KontP =       rbind("Layer=2, diag=ctr" = c(1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=2, diag=sch" = c(1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0))
KontP

Kont_kulP2 =  rbind("Layer=2, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0))

kont_LP2 = glht(parv, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L3
KontP =       rbind("Layer=3, diag=ctr" = c(1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=3, diag=sch" = c(1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0))
KontP

Kont_kulP2 =  rbind("Layer=3, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0))

kont_LP2 = glht(parv, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L4
KontP =       rbind("Layer=4, diag=ctr" = c(1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=4, diag=sch" = c(1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0))
KontP

Kont_kulP2 =  rbind("Layer=4, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0))

kont_LP2 = glht(parv, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L5
KontP =       rbind("Layer=5, diag=ctr" = c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0), 
                    "Layer=5, diag=sch" = c(1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0))
KontP

Kont_kulP2 =  rbind("Layer=5, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0))

kont_LP2 = glht(parv, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L6
KontP =       rbind("Layer=6, diag=ctr" = c(1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0), 
                    "Layer=6, diag=sch" = c(1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1))
KontP

Kont_kulP2 =  rbind("Layer=6, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1))

kont_LP2 = glht(parv, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

# no need for additional correction for multiple comp this way

# L2 correlations

ggplot(L2, aes(x=pmi, y=density, color=diag)) +
  ylim(0, 12000)+
  xlim(0, 100)+
  scale_color_manual(values=c("#00CED1", "#FA8072"))+
  geom_smooth(method="lm", fill =NA)+
  geom_point(size=3)+
  theme_bw()+
  ggtitle(" PMI ~ Parvalbumin density in  Layer 2") +
  xlab("PMI") + ylab("Density (cell/cm^2)")+
  theme(text = element_text(size=20))

shapiro.test(c2$density)
shapiro.test(s2$density)
shapiro.test(L2$density)
cor.test(c2$density, c2$pmi, method=c("pearson"))
cor.test(s2$density, s2$pmi, method=c("pearson"))
cor.test(L2$density, L2$pmi, method=c("spearman"))


ggplot(L2, aes(x=age, y=density, color=diag)) +
  ylim(0, 12000)+
  xlim(45, 100)+
  scale_color_manual(values=c("#00CED1", "#FA8072"))+
  geom_smooth(method="lm", fill =NA)+
  geom_point(size=3)+
  theme_bw()+
  ggtitle(" Age ~ Parvalbumin density in  Layer 2") +
  xlab("Age") + ylab("Density (cell/cm^2)")+
  theme(text = element_text(size=20))

cor.test(c2$density, c2$age, method=c("pearson"))
cor.test(s2$density, s2$age, method=c("pearson"))
cor.test(L2$density, L2$age, method=c("spearman"))

ggplot(L2, aes(x=gendernum, y=density, color=diag)) +
  ylim(0, 11000)+
  xlim(-0.5, 1.5)+
  geom_smooth(method="lm", fill =NA)+
  scale_color_manual(values=c("#00CED1", "#FA8072"))+
  geom_point(size=3)+
  theme_bw()+
  ggtitle(" Gender ~ PV density in Layer 2") +
  xlab("Gender") + ylab("PV cell density (cell/cm^2)")+
  theme(text = element_text(size=20))

cor.test(c2$density, c2$gendernum, method=c("pearson"))
cor.test(s2$density, s2$gendernum, method=c("pearson"))
cor.test(L2$density, L2$gendernum, method=c("spearman"))


## **CALRETININ**

calr = lme(calretinin ~ diag*Layer,
           random = ~1|ID,
           data=sumd,
           na.action = na.exclude)
summary(calr)
qqnorm(resid(calr))
qqline(resid(calr))
reff = ranef(calr)
qqnorm(reff[["(Intercept)"]])
qqline(reff[["(Intercept)"]])
hist(reff[["(Intercept)"]]) 

#L2
KontP =       rbind("Layer=2, diag=ctr" = c(1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=2, diag=sch" = c(1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0))
KontP

Kont_kulP2 =  rbind("Layer=2, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0))

kont_LP2 = glht(calr, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L3
KontP =       rbind("Layer=2, diag=ctr" = c(1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=2, diag=sch" = c(1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0))
KontP

Kont_kulP2 =  rbind("Layer=2, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0))

kont_LP2 = glht(calr, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L4
KontP =       rbind("Layer=2, diag=ctr" = c(1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=2, diag=sch" = c(1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0))
KontP

Kont_kulP2 =  rbind("Layer=2, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0))

kont_LP2 = glht(calr, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L5
KontP =       rbind("Layer=2, diag=ctr" = c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0), 
                    "Layer=2, diag=sch" = c(1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0))
KontP

Kont_kulP2 =  rbind("Layer=2, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0))

kont_LP2 = glht(calr, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L6
KontP =       rbind("Layer=6, diag=ctr" = c(1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0), 
                    "Layer=6, diag=sch" = c(1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1))
KontP

Kont_kulP2 =  rbind("Layer=6, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1))

kont_LP2 = glht(calr, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)


# CR L2 correlations
ggplot(L2, aes(x=pmi, y=calretinin, color=diag)) +
  ylim(0, 30000)+
  xlim(0, 100)+
  scale_color_manual(values=c("#00CED1", "#FA8072"))+
  geom_smooth(method="lm", fill =NA)+
  geom_point(size=3)+
  theme_bw()+
  ggtitle(" PMI ~ Calretinin density in  Layer 2") +
  xlab("PMI") + ylab("Density (cell/cm^2)")+
  theme(text = element_text(size=20))

shapiro.test(c2$calretinin) # normal
shapiro.test(s2$calretinin) # normal
shapiro.test(L2$calretinin) # norm
cor.test(c2$calretinin, c2$pmi, method=c("pearson"))
cor.test(s2$calretinin, s2$pmi, method=c("pearson"))
cor.test(L2$calretinin, L2$pmi, method=c("pearson"))

ggplot(L2, aes(x=age, y=calretinin, color=diag)) +
  ylim(0, 30000)+
  xlim(45, 100)+
  scale_color_manual(values=c("#00CED1", "#FA8072"))+
  geom_smooth(method="lm", fill =NA)+
  geom_point(size=3)+
  theme_bw()+
  ggtitle(" Age ~ Calretinin density in  Layer 2") +
  xlab("Age") + ylab("Density (cell/cm^2)")+
  theme(text = element_text(size=20))

cor.test(c2$calretinin, c2$age, method=c("pearson"))
cor.test(s2$calretinin, s2$age, method=c("pearson"))
cor.test(L2$calretinin, L2$age, method=c("pearson"))

ggplot(L2, aes(x=gendernum, y=calretinin, color=diag)) +
  ylim(0, 30000)+
  xlim(-0.5, 1.5)+
  geom_smooth(method="lm", fill =NA)+
  scale_color_manual(values=c("#00CED1", "#FA8072"))+
  geom_point(size=3)+
  theme_bw()+
  ggtitle(" Gender ~ Calretinin density in Layer 2") +
  xlab("Gender") + ylab("PV cell density (cell/cm^2)")+
  theme(text = element_text(size=20))

cor.test(c2$calretinin, c2$gendernum, method=c("pearson"))
cor.test(s2$calretinin, s2$gendernum, method=c("pearson"))
cor.test(L2$calretinin, L2$gendernum, method=c("pearson"))


## **NPY**

npy = lme(npy ~ diag*Layer,
          random = ~1|ID,
          data=sum6)
summary(npy)
qqnorm(resid(npy))
qqline(resid(npy))
reff = ranef(npy)
qqnorm(reff[["(Intercept)"]])
qqline(reff[["(Intercept)"]])
hist(reff[["(Intercept)"]])

#L2
KontP =       rbind("Layer=2, diag=ctr" = c(1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=2, diag=sch" = c(1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0))
KontP

Kont_kulP2 =  rbind("Layer=2, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0))

kont_LP2 = glht(npy, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L3
KontP =       rbind("Layer=3, diag=ctr" = c(1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=3, diag=sch" = c(1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0))
KontP

Kont_kulP2 =  rbind("Layer=3, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0))

kont_LP2 = glht(npy, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L4
KontP =       rbind("Layer=4, diag=ctr" = c(1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=4, diag=sch" = c(1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0))
KontP

Kont_kulP2 =  rbind("Layer=4, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0))

kont_LP2 = glht(npy, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L5
KontP =       rbind("Layer=5, diag=ctr" = c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0), 
                    "Layer=5, diag=sch" = c(1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0))
KontP

Kont_kulP2 =  rbind("Layer=5, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0))

kont_LP2 = glht(npy, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L6
KontP =       rbind("Layer=6, diag=ctr" = c(1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0), 
                    "Layer=6, diag=sch" = c(1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1))
KontP

Kont_kulP2 =  rbind("Layer=6, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1))

kont_LP2 = glht(npy, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#NPY L4 correlations
ggplot(L4sum, aes(x=pmi, y=npy, color=diag)) +
  ylim(0, 750)+
  xlim(0, 30)+
  scale_color_manual(values=c("#00CED1", "#FA8072"))+
  geom_smooth(method="lm", fill =NA)+
  geom_point(size=3)+
  theme_bw()+
  ggtitle(" PMI ~ NPY density in  Layer 4") +
  xlab("PMI") + ylab("Density (cell/cm^2)")+
  theme(text = element_text(size=20))

shapiro.test(c4sum$npy) 
shapiro.test(s4sum$npy) 
shapiro.test(L4sum$npy)
cor.test(c4sum$npy, c4sum$pmi, method=c("pearson"))
cor.test(s4sum$npy, s4sum$pmi, method=c("spearman"))
cor.test(L4sum$npy, L4sum$pmi, method=c("spearman"))

ggplot(L4sum, aes(x=age, y=npy, color=diag)) +
  ylim(0, 750)+
  xlim(45, 100)+
  scale_color_manual(values=c("#00CED1", "#FA8072"))+
  geom_smooth(method="lm", fill =NA)+
  geom_point(size=3)+
  theme_bw()+
  ggtitle(" Age ~ NPY density in  Layer 4") +
  xlab("Age") + ylab("Density (cell/cm^2)")+
  theme(text = element_text(size=20))

cor.test(c4sum$npy, c4sum$age, method=c("pearson"))
cor.test(s4sum$npy, s4sum$age, method=c("spearman"))
cor.test(L4sum$npy, L4sum$age, method=c("spearman"))

ggplot(L4sum, aes(x=gendernum, y=npy, color=diag)) +
  ylim(0, 1000)+
  xlim(-0.5, 1.5)+
  geom_smooth(method="lm", fill =NA)+
  scale_color_manual(values=c("#00CED1", "#FA8072"))+
  geom_point(size=3)+
  theme_bw()+
  ggtitle(" Gender ~ NPY density in Layer 2") +
  xlab("Gender") + ylab("Density (cell/cm^2)")+
  theme(text = element_text(size=20))

cor.test(c4sum$npy, c4sum$gendernum, method=c("pearson"))
cor.test(s4sum$npy, s4sum$gendernum, method=c("spearman"))
cor.test(L4sum$npy, L4sum$gendernum, method=c("spearman"))

## CALBINDIN
cb = lme(calbindin ~ diag*Layer, 
         random = ~1|ID, 
         data=sum6)
summary(cb)
qqnorm(resid(cb))
qqline(resid(cb))
reff = ranef(cb)
qqnorm(reff[["(Intercept)"]])
qqline(reff[["(Intercept)"]])
hist(reff[["(Intercept)"]])

#L2
KontP =       rbind("Layer=2, diag=ctr" = c(1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=2, diag=sch" = c(1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0))
KontP

Kont_kulP2 =  rbind("Layer=2, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0))

kont_LP2 = glht(cb, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L3
KontP =       rbind("Layer=3, diag=ctr" = c(1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=3, diag=sch" = c(1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0))
KontP

Kont_kulP2 =  rbind("Layer=3, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0))

kont_LP2 = glht(cb, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L4
KontP =       rbind("Layer=4, diag=ctr" = c(1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=4, diag=sch" = c(1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0))
KontP

Kont_kulP2 =  rbind("Layer=4, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0))

kont_LP2 = glht(cb, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L5
KontP =       rbind("Layer=5, diag=ctr" = c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0), 
                    "Layer=5, diag=sch" = c(1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0))
KontP

Kont_kulP2 =  rbind("Layer=5, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0))

kont_LP2 = glht(cb, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L6
KontP =       rbind("Layer=6, diag=ctr" = c(1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0), 
                    "Layer=6, diag=sch" = c(1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1))
KontP

Kont_kulP2 =  rbind("Layer=6, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1))

kont_LP2 = glht(cb, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

# CB L2 correlations
ggplot(L2sum, aes(x=pmi, y=calbindin, color=diag)) +
  ylim(0, 7500)+
  xlim(0, 30)+
  scale_color_manual(values=c("#00CED1", "#FA8072"))+
  geom_smooth(method="lm", fill =NA)+
  geom_point(size=3)+
  theme_bw()+
  ggtitle(" PMI ~ Calbindin density in  Layer 4") +
  xlab("PMI") + ylab("Density (cell/cm^2)")+
  theme(text = element_text(size=20))

shapiro.test(c2sum$calbindin) #  normal
shapiro.test(s2sum$calbindin) #  normal
shapiro.test(L2sum$calbindin)
cor.test(c2sum$calbindin, c2sum$pmi, method=c("pearson"))
cor.test(s2sum$calbindin, s2sum$pmi, method=c("pearson"))
cor.test(L2sum$calbindin, L2sum$pmi, method=c("pearson"))

ggplot(L2sum, aes(x=age, y=calbindin, color=diag)) +
  ylim(0, 7500)+
  xlim(45, 100)+
  scale_color_manual(values=c("#00CED1", "#FA8072"))+
  geom_smooth(method="lm", fill =NA)+
  geom_point(size=3)+
  theme_bw()+
  ggtitle(" Age ~ Calbindin density in  Layer 2") +
  xlab("Age") + ylab("Density (cell/cm^2)")+
  theme(text = element_text(size=20))

cor.test(c2sum$calbindin, c2sum$age, method=c("pearson"))
cor.test(s2sum$calbindin, s2sum$age, method=c("pearson"))
cor.test(L2sum$calbindin, L2sum$age, method=c("pearson"))


ggplot(L2sum, aes(x=gendernum, y=calbindin, color=diag)) +
  ylim(0, 7500)+
  xlim(-0.5, 1.5)+
  geom_smooth(method="lm", fill =NA)+
  scale_color_manual(values=c("#00CED1", "#FA8072"))+
  geom_point(size=3)+
  theme_bw()+
  ggtitle(" Gender ~ Calbindin density in Layer 2") +
  xlab("Gender") + ylab("Density (cell/cm^2)")+
  theme(text = element_text(size=20))

cor.test(c2sum$calbindin, c2sum$gendernum, method=c("pearson"))
cor.test(s2sum$calbindin, s2sum$gendernum, method=c("pearson"))
cor.test(L2sum$calbindin, L2sum$gendernum, method=c("pearson"))

# *SMI 311*

smi = lme(SMI ~ diag*Layer, random = ~1|ID, data=sum6)
summary(smi)
qqnorm(resid(smi))
qqline(resid(smi))
reff = ranef(smi)
qqnorm(reff[["(Intercept)"]])
qqline(reff[["(Intercept)"]])
hist(reff[["(Intercept)"]]) 

# no significant differences

# *NISSL*

nissl = lme(nissl ~ diag*Layer, random = ~1|ID, data=sumd)
summary(nissl)
qqnorm(resid(nissl))
qqline(resid(nissl))
reff = ranef(nissl)
qqnorm(reff[["(Intercept)"]])
qqline(reff[["(Intercept)"]])
hist(reff[["(Intercept)"]])

#L2
KontP =       rbind("Layer=2, diag=ctr" = c(1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=2, diag=sch" = c(1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0))
KontP

Kont_kulP2 =  rbind("Layer=2, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0))

kont_LP2 = glht(nissl, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L3
KontP =       rbind("Layer=2, diag=ctr" = c(1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=2, diag=sch" = c(1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0))
KontP

Kont_kulP2 =  rbind("Layer=2, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0))

kont_LP2 = glht(nissl, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L4
KontP =       rbind("Layer=2, diag=ctr" = c(1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=2, diag=sch" = c(1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0))
KontP

Kont_kulP2 =  rbind("Layer=2, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0))

kont_LP2 = glht(nissl, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L5
KontP =       rbind("Layer=2, diag=ctr" = c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0), 
                    "Layer=2, diag=sch" = c(1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0))
KontP

Kont_kulP2 =  rbind("Layer=2, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0))

kont_LP2 = glht(nissl, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L6
KontP =       rbind("Layer=6, diag=ctr" = c(1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0), 
                    "Layer=6, diag=sch" = c(1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1))
KontP

Kont_kulP2 =  rbind("Layer=6, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1))

kont_LP2 = glht(nissl, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

# NISSL L2 correlations

shapiro.test(c2sum$nissl) 
shapiro.test(s2sum$nissl) 
shapiro.test(L2sum$nissl)
cor.test(c2sum$nissl, c2sum$pmi, method=c("pearson"))
cor.test(s2sum$nissl, s2sum$pmi, method=c("pearson"))
cor.test(L2sum$nissl, L2sum$pmi, method=c("pearson"))


cor.test(c2sum$nissl, c2sum$age, method=c("pearson"))
cor.test(s2sum$nissl, s2sum$age, method=c("pearson"))
cor.test(L2sum$nissl, L2sum$age, method=c("pearson"))


cor.test(c2sum$nissl, c2sum$gendernum, method=c("pearson"))
cor.test(s2sum$nissl, s2sum$gendernum, method=c("pearson"))
cor.test(L2sum$nissl, L2sum$gendernum, method=c("pearson"))

## NISSL L4 correlations

shapiro.test(c4sum$nissl) #  normal
shapiro.test(s4sum$nissl) #  normal
shapiro.test(L4sum$nissl)
cor.test(c4sum$nissl, c4sum$pmi, method=c("pearson"))
cor.test(s4sum$nissl, s4sum$pmi, method=c("pearson"))
cor.test(L4sum$nissl, L4sum$pmi, method=c("pearson"))


cor.test(c4sum$nissl, c4sum$age, method=c("pearson"))
cor.test(s4sum$nissl, s4sum$age, method=c("pearson"))
cor.test(L4sum$nissl, L4sum$age, method=c("pearson"))


cor.test(c4sum$nissl, c4sum$gendernum, method=c("pearson"))
cor.test(s4sum$nissl, s4sum$gendernum, method=c("pearson"))
cor.test(L4sum$nissl, L4sum$gendernum, method=c("pearson"))


## CRCRHVIP RNAscope data


neu <- read.table("recountedall_withfirst12102020g.csv", 
                  header=TRUE, 
                  sep=",", 
                  dec=".")
neu$Layer   <- as.character(neu$Layer)
neu$Layer   <- as.factor(neu$Layer)
neu$diag    <- as.factor(neu$diag)
neu$gender  <- as.factor(neu$gender)
str(neu)
#neu$diag <- relevel(neu$diag, "ctr")

  L1crh<-neu[neu$Layer=="1",]
  L2crh<-neu[neu$Layer=="2",]
  L3crh<-neu[neu$Layer=="3",]
  L4crh<-neu[neu$Layer=="4",]
  L5crh<-neu[neu$Layer=="5",]
  L6crh<-neu[neu$Layer=="6",]
  ccrh<-neu[neu$diag=="ctr",]
  scrh <- neu[neu$diag=="sch",]
    c1crh<-ccrh[ccrh$Layer=="1",]
    c2crh<-ccrh[ccrh$Layer=="2",]
    c3crh<-ccrh[ccrh$Layer=="3",]
    c4crh<-ccrh[ccrh$Layer=="4",]
    c5crh<-ccrh[ccrh$Layer=="5",]
    c6crh<-ccrh[ccrh$Layer=="6",]
    s1crh<-scrh[scrh$Layer=="1",]
    s2crh<-scrh[scrh$Layer=="2",]
    s3crh<-scrh[scrh$Layer=="3",]
    s4crh<-scrh[scrh$Layer=="4",]
    s5crh<-scrh[scrh$Layer=="5",]
    s6crh<-scrh[scrh$Layer=="6",]

mod = lme(crcrhvip_first ~ diag*Layer, 
              random = ~1|ID, 
              data=neu, 
              na.action = na.exclude)
summary(mod)
qqnorm(resid(mod))
qqline(resid(mod))
reff = ranef(mod)
qqnorm(reff[["(Intercept)"]])
qqline(reff[["(Intercept)"]])
hist(reff[["(Intercept)"]]) 

#L2

Kont =        rbind("Layer=2, diag=ctr" = c(1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=2, diag=sch" = c(1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0))
Kont

Kont_kul2 =   rbind("Layer=2, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0))

kont_L2 = glht(mod, linfct=Kont_kul2)
summary(kont_L2)

#L3

Kont =        rbind("Layer=3, diag=ctr" = c(1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=3, diag=sch" = c(1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0))
Kont

Kont_kul2 =   rbind("Layer=3, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0))

kont_L2 = glht(mod, linfct=Kont_kul2)
summary(kont_L2)

#L4
KontP =       rbind("Layer=4, diag=ctr" = c(1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=4, diag=sch" = c(1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0))
KontP

Kont_kulP2 =  rbind("Layer=4, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0))

kont_LP2 = glht(mod, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L5
KontP =       rbind("Layer=5, diag=ctr" = c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0), 
                    "Layer=5, diag=sch" = c(1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0))
KontP

Kont_kulP2 =  rbind("Layer=5, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0))

kont_LP2 = glht(mod, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L6
KontP =       rbind("Layer=6, diag=ctr" = c(1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0), 
                    "Layer=6, diag=sch" = c(1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1))
KontP

Kont_kulP2 =  rbind("Layer=6, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1))

kont_LP2 = glht(mod, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

## CRCRHVIP Rnasdcope data: PMI, age, gender korrel

# L2 
ggplot(L2crh, aes(x=pmi, y=crcrhvip_first, color=diag)) +
  ylim(0, 20000)+
  xlim(0, 30)+
  scale_color_manual(values=c("#00CED1", "#FA8072"))+
  geom_smooth(method="lm", fill =NA)+
  geom_point(size=3)+
  theme_bw()+
  ggtitle(" PMI ~ CRCRHVIP density in  Layer 2") +
  xlab("PMI") + ylab("Density (cell/cm^2)")+
  theme(text = element_text(size=20))

shapiro.test(c2crh$crcrhvip_first) # not normal
shapiro.test(s2crh$crcrhvip_first) # not normal
shapiro.test(L2crh$crcrhvip_first) # not norm
cor.test(c2crh$crcrhvip_first, c2crh$pmi, method=c("spearman"))
cor.test(s2crh$crcrhvip_first, s2crh$pmi, method=c("spearman"))
cor.test(L2crh$crcrhvip_first, L2crh$pmi, method=c("spearman"))

ggplot(L2crh, aes(x=age, y=crcrhvip_first, color=diag)) +
  ylim(0, 20000)+
  xlim(45, 100)+
  scale_color_manual(values=c("#00CED1", "#FA8072"))+
  geom_smooth(method="lm", fill =NA)+
  geom_point(size=3)+
  theme_bw()+
  ggtitle(" Age ~ CRCRHVIP density in  Layer 2") +
  xlab("Age") + ylab("Density (cell/cm^2)")+
  theme(text = element_text(size=20))

cor.test(c2crh$crcrhvip_first, c2crh$age, method=c("spearman"))
cor.test(s2crh$crcrhvip_first, s2crh$age, method=c("spearman"))
cor.test(L2crh$crcrhvip_first, L2crh$age, method=c("spearman"))

ggplot(L2crh, aes(x=gendernum, y=crcrhvip_first, color=diag)) +
  ylim(0, 20000)+
  xlim(-0.5, 1.5)+
  geom_smooth(method="lm", fill =NA)+
  scale_color_manual(values=c("#00CED1", "#FA8072"))+
  geom_point(size=3)+
  theme_bw()+
  ggtitle(" Gender ~ CRCRHVIP density in Layer 2") +
  xlab("Gender") + ylab("Density (cell/cm^2)")+
  theme(text = element_text(size=20))

cor.test(c2crh$crcrhvip_first, c2crh$gendernum, method=c("spearman"))
cor.test(s2crh$crcrhvip_first, s2crh$gendernum, method=c("spearman"))
cor.test(L2crh$crcrhvip_first, L2crh$gendernum, method=c("spearman"))


# L3

ggplot(L3crh, aes(x=pmi, y=crcrhvip_first, color=diag)) +
  ylim(0, 20000)+
  xlim(0, 30)+
  scale_color_manual(values=c("#00CED1", "#FA8072"))+
  geom_smooth(method="lm", fill =NA)+
  geom_point(size=3)+
  theme_bw()+
  ggtitle(" PMI ~ CRCRHVIP density in  Layer 3") +
  xlab("PMI") + ylab("Density (cell/cm^2)")+
  theme(text = element_text(size=20))

shapiro.test(c3crh$crcrhvip_first) 
shapiro.test(s3crh$crcrhvip_first) 
shapiro.test(L3crh$crcrhvip_first)
cor.test(c3crh$crcrhvip_first, c3crh$pmi, method=c("pearson"))
cor.test(s3crh$crcrhvip_first, s3crh$pmi, method=c("pearson"))
cor.test(L3crh$crcrhvip_first, L3crh$pmi, method=c("pearson"))

ggplot(L3crh, aes(x=age, y=crcrhvip_first, color=diag)) +
  ylim(0, 20000)+
  xlim(45, 100)+
  scale_color_manual(values=c("#00CED1", "#FA8072"))+
  geom_smooth(method="lm", fill =NA)+
  geom_point(size=3)+
  theme_bw()+
  ggtitle(" Age ~ CRCRHVIP density in  Layer 3") +
  xlab("Age") + ylab("Density (cell/cm^2)")+
  theme(text = element_text(size=20))

cor.test(c3crh$crcrhvip_first, c3crh$age, method=c("pearson"))
cor.test(s3crh$crcrhvip_first, s3crh$age, method=c("pearson"))
cor.test(L3crh$crcrhvip_first, L3crh$age, method=c("pearson"))

ggplot(L3crh, aes(x=gendernum, y=crcrhvip_first, color=diag)) +
  ylim(0, 20000)+
  xlim(-0.5, 1.5)+
  geom_smooth(method="lm", fill =NA)+
  scale_color_manual(values=c("#00CED1", "#FA8072"))+
  geom_point(size=3)+
  theme_bw()+
  ggtitle(" Gender ~ CRCRHVIP density in Layer 3") +
  xlab("Gender") + ylab("Density (cell/cm^2)")+
  theme(text = element_text(size=20))

cor.test(c3crh$crcrhvip_first, c3crh$gendernum, method=c("pearson"))
cor.test(s3crh$crcrhvip_first, s3crh$gendernum, method=c("pearson"))
cor.test(L3crh$crcrhvip_first, L3crh$gendernum, method=c("pearson"))


