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
#install.packages('car')        # To check multicollinearity 
#install.packages("corrplot")   # plot correlation plot
library('car')
library('corrplot')
library('performance')
library('stats')
# if getting conflicts with mass then use
#oo <- options(repos = "https://cran.r-project.org/")
#install.packages("Matrix")
#install.packages("lme4")
#options(oo)

dat <- read.csv("spawn.csv") # for most analyses
##dat <- read.csv("spawn3.csv") # only when adding some alternative metrics
dat$floatgroup=as.factor(dat$floatgroup)
dat$colony=as.factor(dat$colony)
dat$night=as.factor(dat$night)
m1<-glmer(cbind(fert,100) ~ nn + (1|colony), family = binomial(link = "logit"), data = dat)
summary(m1)
qqnorm(resid(m1))
qqline(resid(m1))

m2<-glmer(cbind(fert,100) ~ nn + night+(1|colony), family = binomial(link = "logit"), data = dat)
summary(m2)
qqnorm(resid(m2))
qqline(resid(m2))
plot(emmeans(m2, specs= ~ nn, cov.keep="nn",type = "response"), xlab = "Fertilisation success", ylab="Nearest neightbour (m)", horizontal = F)
check_model(m2)
check_residuals(m2)
check_overdispersion(m2)
check_singularity(m2)
check_distribution(m2)
plot(check_distribution(m2))

# suggests using a negative binomial (zero inflated)
m2<-glmer(cbind(fert,100) ~ nn + night+(1|colony), family = binomial(link = "logit"), data = dat)

r2m2=r2beta(m2,partial=TRUE,method="nsj",data=dat)
capture.output(r2m2, file = "R2_m2.txt")

## using beta-binomial gives much better model fit as suggested by the analysis of model distribution earlier
m3<-glmmTMB(cbind(fert,100) ~ nn+night+(1|colony),betabinomial(link="logit"),data = dat)
m3sum=summary(m3)
capture.output(m3sum, file = "betabinom_m3_sumMAY.txt")

qqnorm(resid(m3))
qqline(resid(m3))
plot(emmeans(m3, specs= ~nn, cov.keep=2,at=list(nn=c(0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20)),type = "response"), xlab = "Fertilisation success", ylab="Nearest neightbour (m)", horizontal = F)

#plot(emmeans(m3, specs= ~nn, cov.keep=2,at=list(nn=c(0,2,4,6,8,10,12,14,16,18,20),night=c("1","2")),type = "response"), xlab = "Fertilisation success", ylab="Nearest neightbour (m)", horizontal = F)
#m3rg<-ref_grid(m3,at=list(nn=c(0,2,4,6,8,10,12,14,16,18,20),night=c("1","2")))
plot(emmeans(m3, specs= ~nn, cov.keep=2,at=list(nn=c(0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20),night=c("1","2")),type = "response"), xlab = "Fertilisation success", ylab="Nearest neightbour (m)", horizontal = F)
m3rg<-ref_grid(m3,at=list(nn=c(0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20),night=c("1","2")))
m3rgoutput=plot(m3rg, type="response")
capture.output(m3rgoutput,file="m3outputMAY.txt")



m3outputs=predict(m3rg, by="night", type="response",interval=c("confidence") )
r2m3=r2beta(m3,partial=TRUE,method="nsj",data=dat)
capture.output(r2m3, file = "R2_m3.txt")

m3mean=m3rgoutput$data$the.emmean
m3LCL=m3rgoutput$data$asymp.LCL
m3UCL=m3rgoutput$data$asymp.UCL
capture.output(m3mean,file="m3meanoutput.txt")
capture.output(m3LCL,file="m3lcl.txt")
capture.output(m3UCL,file="m3ucl.txt")

detach("package:emmeans", unload = TRUE)# this uses predict function in glmmTMB
nd_nn <- data.frame(nn=c(0,1,2,3,4,5,6,7,8,9,10,15,20),night="1",
                     colony=NA)

predict(m3,type = c("response"),newdata=nd_nn)
predict(m2,type = c("response"),newdata=nd_nn)

##outputing values for use in fertilisation potential calculation
datcolnnlist <- read.csv("colonynnlist.csv") # predict fert success using model for night 1 for each colony
m3rg_sample<-ref_grid(m3,at=list(nn=c(0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20),night=c("1","2")))
nnlist=t(datcolnnlist$nn) # t transposes matrix
m3rg_sample<-ref_grid(m3,at=list(nn=c(0,nnlist),night=c("1","2")))
m3outputsamples=predict(m3rg_sample, by="night", type="response",interval=c("confidence") )
capture.output(m3outputsamples,file="m3outputsamples.txt")

# attempt at negative binomial
m3<-glmer.nb(fert ~ nn+night+(1|colony),data = dat)## didn't converge properly
summary(m3)

# responses on night 1 by floatgroup. No effect over and above the nn
datn1<-dat%>%
  filter(night=='1')
datn1$floatgroup[8]=2
datn1$floatgroup[15]=2
datn1$floatgroup[12]=1

m4<-glmmTMB(cbind(fert,100) ~ nn+floatgroup,betabinomial(link="logit"),data = datn1)
summary(m4)

# run with more complex neighbourhood metrics. First sum of all neighbours (area) within 10 m diameter
dat10unw <- read.csv("spawn_10m_unweighted.csv")
dat10unw$floatgroup=as.factor(dat10unw$floatgroup)
dat10unw$colony=as.factor(dat10unw$colony)
dat10unw$night=as.factor(dat10unw$night)
dat10unw
m5<-glmmTMB(cbind(fert,100) ~ nn_weighted+night+(1|colony),betabinomial(link="logit"),data = dat10unw)
summary(m5)
qqnorm(resid(m5))
qqline(resid(m5))
plot(emmeans(m5, specs= ~nn_weighted, cov.keep=2,at=list(nn_weighted=c(0,2,4,6,8,10,12,14,16,18,20)),type = "response"), xlab = "Fertilisation success", ylab="Nearest neightbour unweighted (m)", horizontal = F)

# run with more complex neighbourhood metrics. So first didn't work. Now with weightings diluted by distance

## using beta-binomial gives much better model fit as suggested by the analysis of model distribution earlier
m6<-glmmTMB(cbind(fert,100) ~ nn_weighted+night+(1|colony),betabinomial(link="logit"),data = dat)
summary(m6)
qqnorm(resid(m6))
qqline(resid(m6))
plot(emmeans(m6, specs= ~nn_weighted, cov.keep=2,at=list(nn_weighted=c(0,2,4,6,8,10,12,14,16,18,20)),type = "response"), xlab = "Fertilisation success", ylab="Nearest neightbour (m)", horizontal = F)

m6sum=summary(m6)
capture.output(m6sum, file = "betabinom_m6.txt")
## AICs much worse when using larger neighbourhoods regardless of weighting.


## Now try using the observed NN from the 26 we tracked, though this could be unreliable as there may have been other closer colonies that spawned, that we tagged, but couldn't watch

m7<-glmmTMB(cbind(fert,100) ~ nn_spawn_obs+night+(1|colony),betabinomial(link="logit"),data = dat)
m7sum=summary(m7)
capture.output(m7sum, file = "betabinom_m7sum.txt")

summary(m7)
qqnorm(resid(m7))
qqline(resid(m7))
plot(emmeans(m7, specs= ~nn_spawn_obs, cov.keep=2,at=list(nn_spawn_obs=c(0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20)),type = "response"), xlab = "Fertilisation success", ylab="Nearest neightbour observed (m)", horizontal = F)
# AIC for this model is 245.6 whereas original AIC for nn model3 is 241.2, therefore better


### now try 2 m local density (count of colonies within 2 m radius)These were gravid but didn't necessarily all go both nights
m8<-glmer(cbind(fert,100) ~ density_2m + night+(1|colony), family = binomial(link = "logit"), data = dat)
summary(m8)
qqnorm(resid(m8))
qqline(resid(m8))
plot(emmeans(m8, specs= ~ density_2m, cov.keep="density_2m",type = "response"), xlab = "Fertilisation success", ylab="Nearest neightbour (m)", horizontal = F)
check_model(m8)
check_residuals(m8)
check_overdispersion(m8)
check_singularity(m8)
check_distribution(m8)
plot(check_distribution(m8))

m9<-glmer.nb(fert ~ density_2m+night+(1|colony),data = dat)## didn't converge properly
summary(m9)

m10<-glmmTMB(cbind(fert,100) ~ density_2m+night+(1|colony),betabinomial(link="logit"),data = dat)
m10sum=summary(m10)
capture.output(m10sum, file = "betabinom_m10sum.txt")

summary(m10)
qqnorm(resid(m10))
qqline(resid(m10))
plot(emmeans(m10, specs= ~density_2m, cov.keep=2,at=list(density_2m=c(0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20)),type = "response"), xlab = "Fertilisation success", ylab="Count of colonies within 2m", horizontal = F)
## count of colonies within 2 m is positive and significant (p=0.04) but the CL are great and the fertlisation success at a count of zero are misleading (0.12 or so) because so many colonies have nothing within 2 m though 
## some are truly isolated and others not really. 

## M11 add the time of spawning release as an additional factor, as minutes after first spawning per evening

m11<-glmmTMB(cbind(fert,100) ~ nn+night/poly(time,2)+(1|colony),betabinomial(link="logit"),data = dat)
m11sum=summary(m11)
capture.output(m11sum, file = "betabinom_m11_sumJul.txt")
summary(m11)
qqnorm(resid(m11))
qqline(resid(m11))
plot(emmeans(m11, specs= ~nn, cov.keep=2,at=list(nn=c(0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20)),type = "response"), xlab = "Fertilisation success", ylab="Nearest neightbour (m)", horizontal = F)
# note that the aic of this model is higher (243.3) than m3 (241.2) so that added factor does not add overall benefits
# moreover, the only significant term associated with the time of spawning is the linear component on night 1 but this has limited interpretive value
# as while early releasers may have faired poorly those later in the period did well (cf results in Levitan 04. Tested whether sensitive to one data point
# which is the first early release but didn't change outcome.

dat2<-dat[-5,]## remove 5th row
m11b<-glmmTMB(cbind(fert,100) ~ nn+night/poly(time,2)+(1|colony),betabinomial(link="logit"),data = dat2)
summary(m11b)

m12<-glmmTMB(cbind(fert,100) ~ nn+night/poly(time,2),betabinomial(link="logit"),data = dat)
summary(m12) ## just rerunning model 11 without random effects which are very small. NO sign diff to model structure



## Run the original m3 but with interaction between night and nn so that slope can vary between nights

m13<-glmmTMB(cbind(fert,100) ~ nn*night+(1|colony),betabinomial(link="logit"),data = dat)
m13sum=summary(m13)
capture.output(m13sum, file = "betabinom_m13_sumMAY.txt")
summary(m13)
qqnorm(resid(m13))
qqline(resid(m13))
plot(emmeans(m13, specs= ~nn, cov.keep=2,at=list(nn=c(0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20)),type = "response"), xlab = "Fertilisation success", ylab="Nearest neightbour (m)", horizontal = F)

#plot(emmeans(m3, specs= ~nn, cov.keep=2,at=list(nn=c(0,2,4,6,8,10,12,14,16,18,20),night=c("1","2")),type = "response"), xlab = "Fertilisation success", ylab="Nearest neightbour (m)", horizontal = F)
#m3rg<-ref_grid(m3,at=list(nn=c(0,2,4,6,8,10,12,14,16,18,20),night=c("1","2")))
plot(emmeans(m13, specs= ~nn, cov.keep=2,at=list(nn=c(0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20),night=c("1","2")),type = "response"), xlab = "Fertilisation success", ylab="Nearest neightbour (m)", horizontal = F)
m13rg<-ref_grid(m13,at=list(nn=c(0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20),night=c("1","2")))
m13rgoutput=plot(m13rg, type="response")
capture.output(m13rgoutput,file="m13outputMAY.txt")

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


## run the different density neighbourhoods though 2 already done above
dat3 <- read.csv("fert_neigh_density2_5_10.csv") # for most analyses
##dat <- read.csv("spawn3.csv") # only when adding some alternative metrics
dat3$colony=as.factor(dat$colony)
dat3$night=as.factor(dat$night)
m14<-glmmTMB(cbind(fert,100) ~ d5m+night+(1|colony),betabinomial(link="logit"),data = dat3)
m14sum=summary(m14)
capture.output(m14sum, file = "betabinom_m14sum.txt")

m15<-glmmTMB(cbind(fert,100) ~ d10m+night+(1|colony),betabinomial(link="logit"),data = dat3)
m15sum=summary(m15)
capture.output(m15sum, file = "betabinom_m15sum.txt")


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
                
## testing effect of depth class on spawning time and fertilisation
spawn4 <- read.csv("spawn4.csv") # for most analyses
spawn4$night=as.factor(spawn4$night)
sp4<-lm(time~night+depth,data=spawn4)
summary(sp4)
sp5<-glmmTMB(cbind(fert,100) ~ night+depth+(1|colony),betabinomial(link="logit"),data = spawn4)
summary(sp5)
sp4sum=summary(sp4)
capture.output(sp4sum, file = "depth_timingsum.txt")
sp5sum=summary(sp5)
capture.output(sp5sum, file = "depth_fertsuccess.txt")


qqnorm(resid(m16))
qqline(resid(m16))
plot(emmeans(m16, specs= ~night, cov.keep=2,at=list(night=c(1,2)),type = "response"), xlab = "Probability of damage", ylab="Night", horizontal = F)