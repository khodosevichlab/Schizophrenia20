Notebook4: IHC, smFISH, RNAscope
================

``` r
library(multcomp)
library(lme4)
library(nlme)
library(ggthemes)
library(ggplot2)
```

``` r
#load tables with density data
sumd <- read.table("alldatag.csv", header=TRUE, sep=",", dec=".")
sum6 <- read.table("6sch6ctrlayerwiseg.csv", header=TRUE, sep=",", dec=".")
```

### 1\. IHC densities

#### PARVALBUMIN

``` r
#lme
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

#glht
#L2
KontP = rbind("Layer=2, diag=ctr" = c(1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=2, diag=sch" = c(1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0))
Kont_kulP2 =  rbind("Layer=2, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0))
kont_LP2 = glht(parv, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L3
KontP = rbind("Layer=3, diag=ctr" = c(1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=3, diag=sch" = c(1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0))
Kont_kulP2 =  rbind("Layer=3, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0))
kont_LP2 = glht(parv, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L4
KontP = rbind("Layer=4, diag=ctr" = c(1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=4, diag=sch" = c(1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0))
Kont_kulP2 =  rbind("Layer=4, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0))
kont_LP2 = glht(parv, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L5
KontP = rbind("Layer=5, diag=ctr" = c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0), 
                    "Layer=5, diag=sch" = c(1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0))
Kont_kulP2 =  rbind("Layer=5, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0))
kont_LP2 = glht(parv, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L6
KontP = rbind("Layer=6, diag=ctr" = c(1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0), 
                    "Layer=6, diag=sch" = c(1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1))
Kont_kulP2 =  rbind("Layer=6, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1))
kont_LP2 = glht(parv, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)
```

#### PV L2 correlation test

``` r
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
```

#### CALRETININ

``` r
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
KontP = rbind("Layer=2, diag=ctr" = c(1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=2, diag=sch" = c(1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0))
Kont_kulP2 =  rbind("Layer=2, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0))
kont_LP2 = glht(calr, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L3
KontP = rbind("Layer=2, diag=ctr" = c(1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=2, diag=sch" = c(1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0))
Kont_kulP2 =  rbind("Layer=2, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0))
kont_LP2 = glht(calr, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L4
KontP = rbind("Layer=2, diag=ctr" = c(1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=2, diag=sch" = c(1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0))
Kont_kulP2 =  rbind("Layer=2, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0))
kont_LP2 = glht(calr, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L5
KontP = rbind("Layer=2, diag=ctr" = c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0), 
                    "Layer=2, diag=sch" = c(1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0))
Kont_kulP2 =  rbind("Layer=2, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0))
kont_LP2 = glht(calr, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L6
KontP = rbind("Layer=6, diag=ctr" = c(1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0), 
                    "Layer=6, diag=sch" = c(1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1))
Kont_kulP2 =  rbind("Layer=6, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1))
kont_LP2 = glht(calr, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)
```

#### CR L2 correlation test

``` r
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
```

#### NPY

``` r
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
KontP = rbind("Layer=2, diag=ctr" = c(1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=2, diag=sch" = c(1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0))
Kont_kulP2 =  rbind("Layer=2, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0))
kont_LP2 = glht(npy, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L3
KontP = rbind("Layer=3, diag=ctr" = c(1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=3, diag=sch" = c(1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0))
Kont_kulP2 =  rbind("Layer=3, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0))
kont_LP2 = glht(npy, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L4
KontP = rbind("Layer=4, diag=ctr" = c(1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=4, diag=sch" = c(1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0))
Kont_kulP2 =  rbind("Layer=4, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0))
kont_LP2 = glht(npy, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L5
KontP = rbind("Layer=5, diag=ctr" = c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0), 
                    "Layer=5, diag=sch" = c(1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0))
Kont_kulP2 =  rbind("Layer=5, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0))
kont_LP2 = glht(npy, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L6
KontP = rbind("Layer=6, diag=ctr" = c(1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0), 
                    "Layer=6, diag=sch" = c(1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1))
Kont_kulP2 =  rbind("Layer=6, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1))
kont_LP2 = glht(npy, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)
```

### NPY L4 correlation test

``` r
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
```

#### CALBINDIN

``` r
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
KontP = rbind("Layer=2, diag=ctr" = c(1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=2, diag=sch" = c(1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0))
Kont_kulP2 =  rbind("Layer=2, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0))
kont_LP2 = glht(cb, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L3
KontP = rbind("Layer=3, diag=ctr" = c(1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=3, diag=sch" = c(1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0))
Kont_kulP2 =  rbind("Layer=3, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0))
kont_LP2 = glht(cb, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L4
KontP = rbind("Layer=4, diag=ctr" = c(1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=4, diag=sch" = c(1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0))
Kont_kulP2 =  rbind("Layer=4, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0))
kont_LP2 = glht(cb, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L5
KontP = rbind("Layer=5, diag=ctr" = c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0), 
                    "Layer=5, diag=sch" = c(1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0))
Kont_kulP2 =  rbind("Layer=5, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0))
kont_LP2 = glht(cb, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L6
KontP =       rbind("Layer=6, diag=ctr" = c(1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0), 
                    "Layer=6, diag=sch" = c(1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1))
Kont_kulP2 =  rbind("Layer=6, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1))
kont_LP2 = glht(cb, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)
```

#### CB L2 correlation test

``` r
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
```

#### SMI 311

``` r
smi = lme(SMI ~ diag*Layer, random = ~1|ID, data=sum6)
summary(smi)
qqnorm(resid(smi))
qqline(resid(smi))
reff = ranef(smi)
qqnorm(reff[["(Intercept)"]])
qqline(reff[["(Intercept)"]])
hist(reff[["(Intercept)"]]) 

# no significant differences
```

#### NISSL

``` r
nissl = lme(nissl ~ diag*Layer, random = ~1|ID, data=sumd)
summary(nissl)
qqnorm(resid(nissl))
qqline(resid(nissl))
reff = ranef(nissl)
qqnorm(reff[["(Intercept)"]])
qqline(reff[["(Intercept)"]])
hist(reff[["(Intercept)"]])

#L2
KontP = rbind("Layer=2, diag=ctr" = c(1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=2, diag=sch" = c(1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0))
Kont_kulP2 =  rbind("Layer=2, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0))
kont_LP2 = glht(nissl, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L3
KontP = rbind("Layer=2, diag=ctr" = c(1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=2, diag=sch" = c(1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0))
Kont_kulP2 =  rbind("Layer=2, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0))
kont_LP2 = glht(nissl, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L4
KontP = rbind("Layer=2, diag=ctr" = c(1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=2, diag=sch" = c(1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0))
Kont_kulP2 =  rbind("Layer=2, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0))
kont_LP2 = glht(nissl, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L5
KontP = rbind("Layer=2, diag=ctr" = c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0), 
                    "Layer=2, diag=sch" = c(1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0))
Kont_kulP2 =  rbind("Layer=2, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0))
kont_LP2 = glht(nissl, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L6
KontP = rbind("Layer=6, diag=ctr" = c(1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0), 
                    "Layer=6, diag=sch" = c(1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1))
Kont_kulP2 =  rbind("Layer=6, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1))
kont_LP2 = glht(nissl, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)
```

#### NISSL L2 correlation

``` r
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
```

### 2\. RNA scope data

#### import data CRCRHVIP

``` r
neu <- read.table("recountedall_withfirst12102020g.csv", 
                  header=TRUE, 
                  sep=",", 
                  dec=".")
neu$Layer   <- as.character(neu$Layer)
neu$Layer   <- as.factor(neu$Layer)
neu$diag    <- as.factor(neu$diag)
neu$gender  <- as.factor(neu$gender)

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
```

``` r
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

Kont = rbind("Layer=2, diag=ctr" = c(1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=2, diag=sch" = c(1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0))
Kont_kul2 =   rbind("Layer=2, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0))
kont_L2 = glht(mod, linfct=Kont_kul2)
summary(kont_L2)

#L3

Kont = rbind("Layer=3, diag=ctr" = c(1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=3, diag=sch" = c(1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0))
Kont_kul2 =   rbind("Layer=3, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0))
kont_L2 = glht(mod, linfct=Kont_kul2)
summary(kont_L2)

#L4
KontP = rbind("Layer=4, diag=ctr" = c(1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0), 
                    "Layer=4, diag=sch" = c(1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0))
Kont_kulP2 =  rbind("Layer=4, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0))
kont_LP2 = glht(mod, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L5
KontP = rbind("Layer=5, diag=ctr" = c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0), 
                    "Layer=5, diag=sch" = c(1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0))
Kont_kulP2 =  rbind("Layer=5, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0))
kont_LP2 = glht(mod, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)

#L6
KontP = rbind("Layer=6, diag=ctr" = c(1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0), 
                    "Layer=6, diag=sch" = c(1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1))
Kont_kulP2 =  rbind("Layer=6, diag=2-1" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1))
kont_LP2 = glht(mod, linfct=Kont_kulP2)
summary(kont_LP2)
confint(kont_LP2)
```

#### CRCRHVIP Rnascope data: PMI, age, gender correlation

``` r
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
```

### 3\. IHC density plots

#### Import data in separate tables

``` r
#CR
cr_dens <- t(read.csv("./cr_densities.csv", header = F, row.names = 1)) %>% as.data.frame %>% 
  `rownames<-`(., NULL) %>% `colnames<-`(., gsub("Layer ", "L", colnames(.)))
cr_dens <-  melt(cr_dens, 
                 id.vars = c("sample", "condition"), 
                 measure.vars =  c("L6", "L5", "L4", "L3", "L2", "L1"), 
                 variable.name = "Layer", 
                 value.name = "Density")
cr_dens$Layer %<>% factor(., levels = sort(unique(.) %>% as.vector))
cr_dens$Density %<>% as.numeric

#PV
pv_dens <- t(read.csv("./pv_densities.csv", header = F, row.names = 1)) %>% as.data.frame %>% 
  `rownames<-`(., NULL) %>% `colnames<-`(., gsub("Layer ", "L", colnames(.)))
pv_dens <- melt(pv_dens, 
                id.vars = c("sample", "condition"), 
                measure.vars = c("L6", "L5", "L4", "L3", "L2", "L1"), 
                variable.name = "Layer", 
                value.name = "Density")
pv_dens$Layer %<>% factor(., levels = sort(unique(.) %>% as.vector))
pv_dens$Density %<>% as.numeric

#CB
cb_dens <- t(read.csv("./cb_densities.csv", header = F, row.names = 1)) %>% as.data.frame %>% 
  `rownames<-`(., NULL) %>% `colnames<-`(., gsub("Layer ", "L", colnames(.)))
cb_dens <- melt(cb_dens, 
                id.vars = c("sample", "condition"), 
                measure.vars =  c("L6", "L5", "L4", "L3", "L2", "L1"), 
                variable.name = "Layer", 
                value.name = "Density")
cb_dens$Layer %<>% factor(., levels = sort(unique(.) %>% as.vector))
cb_dens$Density %<>% as.numeric

#NISSL
nissl_dens <- t(read.csv("./nissl_densities.csv", header = F, row.names = 1)) %>% as.data.frame %>% 
  `rownames<-`(., NULL) %>% `colnames<-`(., gsub("Layer ", "L", colnames(.)))
nissl_dens <- melt(nissl_dens, 
                   id.vars = c("sample", "condition"), 
                   measure.vars =  c("L6", "L5", "L4", "L3", "L2", "L1"), 
                   variable.name = "Layer", 
                   value.name = "Density")
nissl_dens$Layer %<>% factor(., levels = sort(unique(.) %>% as.vector))
nissl_dens$Density %<>% as.numeric

#NPY
npy_dens <- t(read.csv("./npy_densities.csv", header = F, row.names = 1)) %>% as.data.frame %>% 
  `rownames<-`(., NULL) %>% `colnames<-`(., gsub("Layer ", "L", colnames(.)))
npy_dens <- melt(npy_dens, 
                 id.vars = c("sample", "condition"), 
                 measure.vars = c("L6", "L5", "L4", "L3", "L2", "L1"), 
                 variable.name = "Layer", 
                 value.name = "Density")
npy_dens$Layer %<>% factor(., leve

#SMI
smi_dens <- t(read.csv("./smi_densities.csv", header = F, row.names = 1)) %>% as.data.frame %>% 
  `rownames<-`(., NULL) %>% `colnames<-`(., gsub("Layer ", "L", colnames(.)))
smi_dens <- melt(smi_dens, 
                 id.vars = c("sample", "condition"), 
                 measure.vars = c("L6", "L5", "L4", "L3", "L2", "L1"),
                 variable.name = "Layer", 
                 value.name = "Density")
smi_dens$Layer %<>% factor(., levels = sort(unique(.) %>% as.vector))
smi_dens$Density %<>% as.numeric
```

calretinin density plot

``` r
cr.plot <-
ggplot(data = cr_dens, aes(x = Layer, y = Density, dodge = condition, fill = condition)) +
 geom_point(aes(y = Density, color = condition), position = 
               position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2), 
               size = 1.5, alpha = 0.2, show.legend = F) +
   stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "pointrange",
                 size = 0.4,
                 show.legend = T, position = position_dodge(width = 0.7), 
                 mapping = aes_(group = ~condition, color = ~condition), fill = "grey",
                 shape = 23) +
  geom_signif(annotations = "p = 0.0028",
             xmin = 1.8, 
             xmax = 2.2, 
             y_position = max(cr_dens$Density)*1.07,
             textsize = 5 ) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size = 13),
          axis.text.y=element_text(size = 13),
          axis.title.x = element_text(size = 13),
          text = element_text(family = "Liberation Sans", size = 13),
          legend.title = element_blank(),
          legend.text = element_text(size = 13),
          legend.position = "top",
          plot.subtitle = element_text(size = 13, hjust = 0.5, face = "bold"),
          panel.grid.minor = element_blank(),
          legend.margin =  margin(t = -1, b = -2.5, r = 0, l = 0, unit = "mm")) +
          labs(x = "n=5254", y = bquote("Density," ~ cells/cm^2), subtitle = "Calretinin+ neurons") +
    scale_color_manual(values = palette_45_2[c(8, 15)]) +
    scale_fill_manual(values = palette_45_2[c(8, 15)]) +
  ylim(0, max(cr_dens$Density)*1.15)
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/calretinin.jpg" width="60%" style="display: block; margin: auto;" />

parvalbumin density plot

``` r
pv.plot <-
ggplot(data = pv_dens, aes(x = Layer, y = Density, dodge = condition, fill = condition)) +
  geom_point(aes(y = Density, color = condition), position = 
               position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2), 
               size = 1.5, alpha = 0.2, show.legend = F) +
   stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "pointrange",
                 size = 0.4,
                 show.legend = T, position = position_dodge(width = 0.7), 
                 mapping = aes_(group = ~condition, color = ~condition), fill = "grey",
                 shape = 23) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size = 13),
          axis.text.y=element_text(size = 13),
          axis.title.x = element_text(size = 13),
          text = element_text(family = "Liberation Sans", size = 13),
          legend.title = element_blank(),
          legend.text = element_text(size = 13),
          legend.position = "top",
          plot.subtitle = element_text(size = 13, hjust = 0.5, face = "bold"),
          panel.grid.minor = element_blank(),
          legend.margin =  margin(t = -1, b = -2.5, r = 0, l = 0, unit = "mm")) +
          labs(x = "n=6507", y = bquote("Density," ~ cells/cm^2), subtitle = "Parvalbumin+ neurons") +
    scale_color_manual(values = palette_45_2[c(8, 15)]) +
    scale_fill_manual(values = palette_45_2[c(8, 15)]) +
  ylim(-60, max(pv_dens$Density)*1.15)
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/parvalbumin.jpg" width="60%" style="display: block; margin: auto;" />

calbindin density plot

``` r
cb.plot <-
ggplot(data = cb_dens, aes(x = Layer, y = Density, dodge = condition, fill = condition)) +
  geom_point(aes(y = Density, color = condition), position = 
               position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2), 
               size = 1.5, alpha = 0.2, show.legend = F) +
   stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "pointrange",
                 size = 0.4,
                 show.legend = T, position = position_dodge(width = 0.7), 
                 mapping = aes_(group = ~condition, color = ~condition), fill = "grey",
                 shape = 23) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size = 13),
          axis.text.y=element_text(size = 13),
          axis.title.x = element_text(size = 13),
          text = element_text(family = "Liberation Sans", size = 13),
          legend.title = element_blank(),
          legend.text = element_text(size = 13),
          legend.position = "top",
          plot.subtitle = element_text(size = 13, hjust = 0.5, face = "bold"),
          panel.grid.minor = element_blank(),
          legend.margin =  margin(t = -1, b = -2.5, r = 0, l = 0, unit = "mm")) +
          labs(x = "n=613", y = bquote("Density," ~ cells/cm^2), subtitle = "Calbindin+ neurons") +
    scale_color_manual(values = palette_45_2[c(8, 15)]) +
    scale_fill_manual(values = palette_45_2[c(8, 15)])
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/calbindin.jpg" width="60%" style="display: block; margin: auto;" />

npy density plot

``` r
npy.plot <-
ggplot(data = npy_dens, aes(x = Layer, y = Density, dodge = condition, fill = condition)) +
  geom_point(aes(y = Density, color = condition), position = 
               position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2), 
               size = 1.5, alpha = 0.2, show.legend = F) +
   stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "pointrange",
                 size = 0.4,
                 show.legend = T, position = position_dodge(width = 0.7), 
                 mapping = aes_(group = ~condition, color = ~condition), fill = "grey",
                 shape = 23) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size = 13),
          axis.text.y=element_text(size = 13),
          axis.title.x = element_text(size = 13),
          text = element_text(family = "Liberation Sans", size = 13),
          legend.title = element_blank(),
          legend.text = element_text(size = 13),
          legend.position = "top",
          plot.subtitle = element_text(size = 13, hjust = 0.5, face = "bold"),
          panel.grid.minor = element_blank(),
          legend.margin =  margin(t = -1, b = -2.5, r = 0, l = 0, unit = "mm")) +
          labs(x = "n=208", y = bquote("Density," ~ cells/cm^2), subtitle = "NPY+ neurons") +
    scale_color_manual(values = palette_45_2[c(8, 15)]) +
    scale_fill_manual(values = palette_45_2[c(8, 15)])
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/npy.jpg" width="60%" style="display: block; margin: auto;" />

nissl density plot

``` r
nissl.plot <-
ggplot(data = nissl_dens, aes(x = Layer, y = Density, dodge = condition, fill = condition)) +
  geom_point(aes(y = Density, color = condition), position = 
               position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2), 
               size = 1.5, alpha = 0.2, show.legend = F) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "pointrange",
                 size = 0.4,
                 show.legend = T, position = position_dodge(width = 0.7), 
                 mapping = aes_(group = ~condition, color = ~condition), fill = "grey",
                 shape = 23) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size = 13),
          axis.text.y=element_text(size = 13),
          axis.title.x = element_text(size = 13),
          text = element_text(family = "Liberation Sans", size = 13),
          legend.title = element_blank(),
          legend.text = element_text(size = 13),
          legend.position = "top",
          plot.subtitle = element_text(size = 13, hjust = 0.5, face = "bold"),
          panel.grid.minor = element_blank(),
          legend.margin =  margin(t = -1, b = -2.5, r = 0, l = 0, unit = "mm")) +
          labs(x = "n=45446", y = bquote("Density," ~ cells/cm^2), subtitle = "Nissl stained neurons") +
    scale_color_manual(values = palette_45_2[c(8, 15)]) +
    scale_fill_manual(values = palette_45_2[c(8, 15)]) +
  
ylim(0, max(nissl_dens$Density)*1.15)
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/nissl.jpg" width="60%" style="display: block; margin: auto;" />

smi density plot

``` r
smi.plot <-
ggplot(data = smi_dens, aes(x = Layer, y = Density, dodge = condition, fill = condition)) +
  geom_point(aes(y = Density, color = condition), position = 
               position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2), 
               size = 1.5, alpha = 0.2, show.legend = F) +
   stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "pointrange",
                 size = 0.4,
                 show.legend = T, position = position_dodge(width = 0.7), 
                 mapping = aes_(group = ~condition, color = ~condition), fill = "grey",
                 shape = 23) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size = 13),
          axis.text.y=element_text(size = 13),
          axis.title.x = element_text(size = 13),
          text = element_text(family = "Liberation Sans", size = 13),
          legend.title = element_blank(),
          legend.text = element_text(size = 13),
          legend.position = "top",
          plot.subtitle = element_text(size = 13, hjust = 0.5, face = "bold"),
          panel.grid.minor = element_blank(),
          legend.margin =  margin(t = -1, b = -2.5, r = 0, l = 0, unit = "mm")) +
          labs(x = "n=41310", y = bquote("Density," ~ cells/cm^2), subtitle = "NEFL/NEFM+ neurons") +
    scale_color_manual(values = palette_45_2[c(8, 15)]) +
    scale_fill_manual(values = palette_45_2[c(8, 15)])
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/smi.jpg" width="60%" style="display: block; margin: auto;" />

### 4\. CR CHR VIP PLOTS

#### import data

``` r
cr.p_crh.p_vip.p <- read.csv("./cr.p_crh.p_vip.p.csv")
cr.p_crh.p_vip.n <- read.csv("./cr.p_crh.p_vip.n.csv")
cr.p_crh.n_vip.n <- read.csv("./cr.p_crh.n_vip.n.csv")
```

``` r
#CR CRH VIP
cr.p_crh.p_vip.p <-  melt(cr.p_crh.p_vip.p, 
                          id.vars = c("sample", "condition"), 
                          measure.vars = c("L6", "L5", "L4", "L3", "L2", "L1"), 
                          variable.name = "Layer", 
                          value.name = "Density")
cr.p_crh.p_vip.p$Layer %<>% factor(., levels = sort(unique(.) %>% as.vector))
cr.p_crh.p_vip.p$Density %<>% as.numeric

#CR CRH
cr.p_crh.p_vip.n <-  melt(cr.p_crh.p_vip.n, 
                          id.vars = c("sample", "condition"),
                          measure.vars = c("L6", "L5", "L4", "L3", "L2", "L1"),
                          variable.name = "Layer", 
                          value.name = "Density")
cr.p_crh.p_vip.n$Layer %<>% factor(., levels = sort(unique(.) %>% as.vector))
cr.p_crh.p_vip.n$Density %<>% as.numeric

#CR
cr.p_crh.n_vip.n <- melt(cr.p_crh.n_vip.n, 
                         id.vars = c("sample", "condition"),
                         measure.vars = c("L6", "L5", "L4", "L3", "L2", "L1"),
                         variable.name = "Layer",
                         value.name = "Density")
cr.p_crh.n_vip.n$Layer %<>% factor(., levels = sort(unique(.) %>% as.vector))
cr.p_crh.n_vip.n$Density %<>% as.numeric
```

#### CR+ CRH+ VIP+ neurons

``` r
cr_crh_vip.plot <-
ggplot(data = cr.p_crh.p_vip.p, aes(x = Layer, y = Density, dodge = condition, fill = condition)) +
  geom_point(aes(y = Density, color = condition), position = 
               position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2), 
               size = 1.5, alpha = 0.2, show.legend = F) +
   stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "pointrange",
                 size = 0.4,
                 show.legend = T, position = position_dodge(width = 0.7), 
                 mapping = aes_(group = ~condition, color = ~condition), fill = "grey",
                 shape = 23) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size = 13),
          axis.text.y=element_text(size = 13),
          axis.title.x = element_text(size = 13),
          text = element_text(family = "Liberation Sans", size = 13),
          legend.title = element_blank(),
          legend.text = element_text(size = 13),
          legend.position = "top",
          plot.subtitle = element_text(size = 13, hjust = 0.5, face = "bold"),
          panel.grid.minor = element_blank(),
          legend.margin =  margin(t = -1, b = -2.5, r = 0, l = 0, unit = "mm")) +
          labs(x = NULL, y = bquote("Density," ~ cells/cm^2), subtitle = "CR+ CRH+ VIP+ neurons") +
    scale_color_manual(values = palette_45_2[c(8, 15)]) +
    scale_fill_manual(values = palette_45_2[c(8, 15)]) +
ylim(0, max(cr.p_crh.p_vip.p$Density)*1.15)
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/crchrvip.jpg" width="60%" style="display: block; margin: auto;" />

##### CR+ CRH+ VIP- neurons

``` r
cr_crh.plot <-
ggplot(data = cr.p_crh.p_vip.n, aes(x = Layer, y = Density, dodge = condition, fill = condition)) +
  geom_point(aes(y = Density, color = condition), position = 
               position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2), 
               size = 1.5, alpha = 0.2, show.legend = F) +
   stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "pointrange",
                 size = 0.4,
                 show.legend = T, position = position_dodge(width = 0.7), 
                 mapping = aes_(group = ~condition, color = ~condition), fill = "grey",
                 shape = 23 +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size = 13),
          axis.text.y=element_text(size = 13),
          axis.title.x = element_text(size = 13),
          text = element_text(family = "Liberation Sans", size = 13),
          legend.title = element_blank(),
          legend.text = element_text(size = 13),
          legend.position = "top",
          plot.subtitle = element_text(size = 13, hjust = 0.5, face = "bold"),
          panel.grid.minor = element_blank(),
          legend.margin =  margin(t = -1, b = -2.5, r = 0, l = 0, unit = "mm")) +
          labs(x = NULL, y = bquote("Density," ~ cells/cm^2), subtitle = "CR+ CRH+ VIP- neurons") +
    scale_color_manual(values = palette_45_2[c(8, 15)]) +
    scale_fill_manual(values = palette_45_2[c(8, 15)]) +
      ylim(-190, max(cr.p_crh.p_vip.n$Density)*1.15)
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/crchr.jpg" width="60%" style="display: block; margin: auto;" />

#### CR+ CRH- VIP- neurons

``` r
cr2.plot <-
ggplot(data = cr.p_crh.n_vip.n, aes(x = Layer, y = Density, dodge = condition, fill = condition)) +
  geom_point(aes(y = Density, color = condition), position = 
               position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2), 
               size = 1.5, alpha = 0.2, show.legend = F) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "pointrange",
                 size = 0.4,
                 show.legend = T, position = position_dodge(width = 0.7), 
                 mapping = aes_(group = ~condition, color = ~condition), fill = "grey",
                 shape = 23) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size = 13),
          axis.text.y=element_text(size = 13),
          axis.title.x = element_text(size = 13),
          text = element_text(family = "Liberation Sans", size = 13),
          legend.title = element_blank(),
          legend.text = element_text(size = 13),
          legend.position = "top",
          plot.subtitle = element_text(size = 13, hjust = 0.5, face = "bold"),
          panel.grid.minor = element_blank(),
          legend.margin =  margin(t = -1, b = -2.5, r = 0, l = 0, unit = "mm")) +
          labs(x = NULL, y = bquote("Density," ~ cells/cm^2), subtitle = "CR+ CRH- VIP- neurons") +
    scale_color_manual(values = palette_45_2[c(8, 15)]) +
    scale_fill_manual(values = palette_45_2[c(8, 15)]) +
  ylim(-1000, max(cr.p_crh.n_vip.n$Density)*1.15)
```

<img src="C:/Users/Katarina/Desktop/scznotebooks/cr.jpg" width="60%" style="display: block; margin: auto;" />
