# load packages
library(brms)
library(rstan)

# rstan options
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
        
## read and prepare Data

setwd("~/GitHub/PapioHormoneComparison")
papio<- read.csv("papio.csv")

# center and standardize
papio$DateOfCollection <- as.Date(papio$DateOfCollection, format = "%Y-%m-%d")
papio$log_CSF_OXT      <- log(as.numeric(papio$CSF_pg.mL))
papio$log_Plasma_OXT   <- log(as.numeric(papio$Plasma_pg.mL))
papio$log_Urine_OXT    <- log(as.numeric(papio$uOT.Creatinine))
papio$log_CSF_AVP      <- log(as.numeric(papio$CSF_AVP_pg.mL))
papio$log_Weight       <- log(as.numeric(papio$Weight)) 
papio$log_weight.z     <- as.numeric(scale(papio$log_Weight))
papio$age.z            <- as.numeric(scale(papio$AgeAtCollection))
papio$date.z           <- as.numeric(scale(papio$DateOfCollection))
papio$corral_size.z    <- as.numeric(scale(papio$CorralSize))
papio$harem_size.z     <- as.numeric(scale(papio$HaremSize_AtDOC))
papio$nOffspring.z     <- as.numeric(scale(papio$nOffspring))
papio$log_CSF_OXT.c    <- papio$log_CSF_OXT-mean(papio$log_CSF_OXT, na.rm=TRUE)
papio$log_plasma_OXT.c <- papio$log_Plasma_OXT-mean(papio$log_Plasma_OXT, na.rm=TRUE)
papio$log_urine_OXT.c  <- papio$log_Urine_OXT-mean(papio$log_Urine_OXT, na.rm=TRUE)
papio$log_CSF_AVP.c    <- papio$log_CSF_AVP-mean(papio$log_CSF_AVP, na.rm=TRUE)  

# males get 0 gestational day and not lactating
papio$Lactating[papio$Sex=="M"]  <- "Not_Lactating"
papio$GD_DOC[papio$Sex=="M"]  <- 0

# 3 females with infant present = No and pregnant = No get 0 gestational day and not lactating
papio$Lactating[is.na(papio$Lactating)]  <- "Not_Lactating"
papio$GD_DOC[is.na(papio$GD_DOC)]  <- 0

papio$Lactating<- relevel(as.factor(papio$Lactating), ref = "Not_Lactating")


## set up multivariate model

CSF_OXT <- bf(log_CSF_OXT.c | mi() ~ log_weight.z + s(date.z) + corral_size.z + Blood_in_CSF + Lactating + GD_DOC + (1|q|Corral) + Subspecies * Sex)

Plasma_OXT <- bf(log_plasma_OXT.c | mi() ~ log_weight.z + s(date.z) + corral_size.z + Plasma_Lysis + Lactating + GD_DOC + (1|q|Corral) + Subspecies * Sex)

Urinary_OXT <- bf(log_urine_OXT.c | mi() ~ log_weight.z + s(date.z) + corral_size.z + Lactating + GD_DOC + (1|q|Corral) + Subspecies * Sex)

CSF_AVP <- bf(log_CSF_AVP.c | mi() ~ log_weight.z + s(date.z) + corral_size.z + Blood_in_CSF + Lactating + GD_DOC + (1|q|Corral) + Subspecies * Sex)

## fit model
m1 <- brm(CSF_OXT + Plasma_OXT + Urinary_OXT + CSF_AVP + 
                             set_rescor(TRUE),
                           
                           data = papio,
                           
                           prior = c(
                             prior(normal(0,2), class = "Intercept"),
                             prior(normal(0,1), class = "b"),
                             prior(exponential(2), class = "sd", resp = "logCSFOXTc"),
                             prior(exponential(2), class = "sd", resp = "logplasmaOXTc"),
                             prior(exponential(2), class = "sd", resp = "logurineOXTc"),
                             prior(exponential(2), class = "sd", resp = "logCSFAVPc"),
                             prior(lkj(2), class="cor")),
                           
                           chains=3, cores=3, warmup=2000, iter=4000, control = list(max_treedepth = 15, adapt_delta=0.9999))
# 3 divergent transitions - this is not a problem, just indicates slightly inefficient sampling

# check convergence
plot(m1) # good mixing
summary(m1) # all Rhat = 1, all ESS > 1000

# save model object (and read as needed)
saveRDS(m1, file = "m1.rds")  
m1<- readRDS("m1.rds")  
  
# examine associations
conditional_effects(m1) 


## run univariate models with censoring
CSF_OXT_cens <- bf(log_CSF_OXT.c | cens(CSF_OXT_censored) ~ log_weight.z + s(date.z) + corral_size.z + Blood_in_CSF + Lactating + GD_DOC  + age.z + (1|q|Corral) + Subspecies * Sex)

Plasma_OXT_cens <- bf(log_plasma_OXT.c | cens(Plasma_OXT_censored) ~ log_weight.z + s(date.z) + corral_size.z + Plasma_Lysis + Lactating + GD_DOC + age.z + (1|q|Corral) + Subspecies * Sex)

Urinary_OXT_cens <- bf(log_urine_OXT.c | cens(Urine_OXT_censored) ~ log_weight.z + s(date.z) + corral_size.z + Lactating + GD_DOC + age.z + (1|q|Corral) + Subspecies * Sex)

CSF_AVP_cens <- bf(log_CSF_AVP.c | cens(CSF_AVP_censored) ~ log_weight.z + s(date.z) + corral_size.z + Blood_in_CSF + Lactating + GD_DOC + age.z + (1|q|Corral) + Subspecies * Sex)

## set initial values, because censored models don't fit well otherwise (see https://bookdown.org/content/3686/tools-in-the-trunk.html)
inits <- list(Intercept = 0, b = rep(0,9))
inits_list <- list(inits, inits, inits)

## fit models
m2a <- brm(CSF_OXT_cens, data = papio,
            prior = c(
            prior(normal(0,2), class = "Intercept"),
            prior(normal(0,1), class = "b"),
            prior(exponential(2), class = "sd")),
          chains=3, cores=3, warmup=2000, iter=4000, init = inits_list,
          control = list(max_treedepth = 19, adapt_delta=0.9999))
# no divergent transitions
# check convergence
plot(m2a) # good mixing
summary(m2a) # all Rhat = 1, all ESS > 1000
# save model objects (and read as needed)
saveRDS(m2a, file = "m2a.rds")  
m2a<- readRDS("m2a.rds")  

m2b <- brm(Plasma_OXT_cens, data = papio,
           prior = c(
             prior(normal(0,2), class = "Intercept"),
             prior(normal(0,1), class = "b"),
             prior(exponential(2), class = "sd")),
           chains=3, cores=3, warmup=2000, iter=4000, init = inits_list, 
           control = list(max_treedepth = 19, adapt_delta=0.9999))
# no divergent transitions
# check convergence
plot(m2b) # good mixing
summary(m2b) # all Rhat = 1, all ESS > 1000
# save model objects (and read as needed)
saveRDS(m2b, file = "m2b.rds")  
m2b<- readRDS("m2b.rds")  

inits <- list(Intercept = 0, b = rep(0,8))
inits_list <- list(inits, inits, inits)

m2c <- brm(Urinary_OXT_cens, data = papio,
           prior = c(
             prior(normal(0,2), class = "Intercept"),
             prior(normal(0,1), class = "b"),
             prior(exponential(2), class = "sd")),
           chains=3, cores=3, warmup=2000, iter=4000, init = inits_list, 
           control = list(max_treedepth = 19, adapt_delta=0.9999))
# no divergent transitions
# check convergence
plot(m2c) # good mixing
summary(m2c) # all Rhat = 1, all ESS > 1000
# save model objects (and read as needed)
saveRDS(m2c, file = "m2c.rds")  
m2c<- readRDS("m2c.rds")  

inits <- list(Intercept = 0, b = rep(0,9))
inits_list <- list(inits, inits, inits)

m2d <- brm(CSF_AVP_cens, data = papio,
           prior = c(
             prior(normal(0,2), class = "Intercept"),
             prior(normal(0,1), class = "b"),
             prior(exponential(2), class = "sd")),
           chains=3, cores=3, warmup=2000, iter=4000, init = inits_list, 
           control = list(max_treedepth = 19, adapt_delta=0.9999))
# no divergent transitions
# check convergence
plot(m2d) # good mixing
summary(m2d) # all Rhat = 1, all ESS > 1000
# save model objects (and read as needed)
saveRDS(m2d, file = "m2d.rds")  
m2d<- readRDS("m2d.rds")  

