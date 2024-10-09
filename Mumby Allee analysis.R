library(tidyverse) ## for data wrangling
library(cowplot) ## for multi-panel plots using ggplot
library(lme4) ## for (generalised) mixed-effects models
library(lmerTest) ## for getting p-values from lmer models
library(emmeans) ## for posthoc multiple comparisons
library(readxl) ## to import xl
library(MASS)
library(forecast)
library(r2glmm)
library(broom)
library(glmmTMB)


dat <- read.csv("spawn.csv") # for most analyses
dat$floatgroup=as.factor(dat$floatgroup)
dat$colony=as.factor(dat$colony)
dat$night=as.factor(dat$night)

## using beta-binomial gives much better model fit as suggested by the analysis of model distribution earlier
m3<-glmmTMB(cbind(fert,100) ~ nn+night+(1|colony),betabinomial(link="logit"),data = dat)
m3sum=summary(m3)
capture.output(m3sum, file = "betabinom_m3_sumMAY.txt")
summary(m3)
qqnorm(resid(m3))
qqline(resid(m3))
plot(emmeans(m3, specs= ~nn, cov.keep=2,at=list(nn=c(0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20)),type = "response"), xlab = "Fertilisation success", ylab="Nearest neightbour (m)", horizontal = F)

m3rg<-ref_grid(m3,at=list(nn=c(0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20),night=c("1","2")))
m3rgoutput=plot(m3rg, type="response")

m3mean=m3rgoutput$data$the.emmean
m3LCL=m3rgoutput$data$asymp.LCL
m3UCL=m3rgoutput$data$asymp.UCL
capture.output(m3mean,file="m3meanoutput.txt")
capture.output(m3LCL,file="m3lcl.txt")
capture.output(m3UCL,file="m3ucl.txt")

## Run the original m3 but with interaction between night and nn so that slope can vary between nights
m4<-glmmTMB(cbind(fert,100) ~ nn*night+(1|colony),betabinomial(link="logit"),data = dat)
summary(m4)

## run the different density neighbourhoods
dat3 <- read.csv("fert_neigh_density2_5_10.csv") # for most analyses
dat3$colony=as.factor(dat$colony)
dat3$night=as.factor(dat$night)

m13<-glmmTMB(cbind(fert,100) ~ d2m+night+(1|colony),betabinomial(link="logit"),data = dat3)
summary(m13)
m14<-glmmTMB(cbind(fert,100) ~ d5m+night+(1|colony),betabinomial(link="logit"),data = dat3)
summary(m14)
m15<-glmmTMB(cbind(fert,100) ~ d10m+night+(1|colony),betabinomial(link="logit"),data = dat3)
summary(m15)

## Use Weighted colony area within 10 m
m6<-glmmTMB(cbind(fert,100) ~ nn_weighted+night+(1|colony),betabinomial(link="logit"),data = dat)
summary(m6)
qqnorm(resid(m6))
qqline(resid(m6))


# responses on night 1 by floatgroup. No effect over and above the nn
datn1<-dat%>%
  filter(night=='1')
datn1$floatgroup[8]=2
datn1$floatgroup[15]=2
datn1$floatgroup[12]=1

m5<-glmmTMB(cbind(fert,100) ~ nn+floatgroup,betabinomial(link="logit"),data = datn1)
summary(m5)


## M11 add the time of spawning release as an additional factor, as minutes after first spawning per evening

dat <- read.csv("spawn3.csv") # only when adding some alternative metrics
dat$floatgroup=as.factor(dat$floatgroup)
dat$colony=as.factor(dat$colony)
dat$night=as.factor(dat$night)


m11<-glmmTMB(cbind(fert,100) ~ nn+night/poly(time,2)+(1|colony),betabinomial(link="logit"),data = dat)
summary(m11)
qqnorm(resid(m11))
qqline(resid(m11))

dat2<-dat[-5,]## remove 5th row: testing for outlier but no impact
m11b<-glmmTMB(cbind(fert,100) ~ nn+night/poly(time,2)+(1|colony),betabinomial(link="logit"),data = dat2)
summary(m11b)


#### Quantify intensity of full spawning on each night

spawn <- read.csv("spawnintensity.csv") # for most analyses
spawn$Night=as.factor(spawn$Night)
sp1<-glmmTMB(Fullspawn ~ Night,binomial(link="logit"),data = spawn)
summary(sp1)
sp1sum=summary(sp1)
capture.output(sp1sum, file = "binomial_sp1_sum.txt")
qqnorm(resid(m16))
qqline(resid(m16))
plot(emmeans(m16, specs= ~night, cov.keep=2,at=list(night=c(1,2)),type = "response"), xlab = "Probability of damage", ylab="Night", horizontal = F)

### Examine whether correlation between synchrony (greater when time difference closer to zero) and proximity between colonies but allowing for 
# differences between nights and whether there are two different people doing the measurement (2 vs 1).

syndist <- read.csv("syndist.csv") # for most analyses
synm1<-glmmTMB(cbind(syn,100) ~ distance * diff_person + (1|night), family = binomial(link = "logit"), data = syndist)
summary(synm1)
synm1sum=summary(synm1)
capture.output(synm1sum, file = "synch_by_dist.txt")
qqnorm(resid(synm1))
qqline(resid(synm1))
plot(emmeans(synm1, specs= ~distance, cov.keep=2,at=list(distance=c(0,5,10,15,20),diff_person=c("1","2")),type = "response"), xlab = "Difference in spawning time (min)", ylab="Distance between colonies (m)", horizontal = F)

synm1rg<-ref_grid(synm1,at=list(distance=c(0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20),diff_person=c("1","2")))
synm1rgoutput=plot(synm1rg, type="response")
synm1outputs=predict(synm1rg, by="diff_person", type="response",interval=c("confidence") )
capture.output(synm1outputs,file="synm1output.txt")

## analyse damaged deformed and fragmented eggs/embryos
def <- read.csv("def.csv") # for most analyses
m16<-glmmTMB(damage ~ night+(1|colony),binomial(link="logit"),data = def)
m16<-glm(damage ~ night,family=quasibinomial(link="logit"),data = def)# use this one

summary(m16)
m16sum=summary(m16)
capture.output(m16sum, file = "damagedembryosum.txt")
qqnorm(resid(m16))
qqline(resid(m16))
m16rg<-ref_grid(m16,at=list(night=c(1,2)))
m16rgoutput=plot(m16rg, type="response")
m16output=predict(m16rg, by="night", type="response",interval=c("confidence") )
capture.output(m16output,file="damageoutput.txt")
                
