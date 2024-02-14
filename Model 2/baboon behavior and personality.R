####################################################################
####################################################################

#Personality analysis of Dan's baboon data

####################################################################
################# Workspace preparation ############################
####################################################################

#load packages
library(brms)
library(cmdstanr)
library(tidyr)
library(ggplot2)

#set path for cmdstan
set_cmdstan_path("...")

#stan settings
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#load data
setwd("...")
data = read.csv("baboon behavior and hormones.csv")

####################### Data Prep ##################################
####################################################################

#log transform, standardize to 2 SD for comparison with binary effects
data$CSF_pg.mL = scale(log(data$CSF_pg.mL),
                        scale= 2 * sd(log(data$CSF_pg.mL),na.rm=TRUE))
data$CSF_AVP_pg.mL = scale(log(data$CSF_AVP_pg.mL),
                            scale= 2 * sd(log(data$CSF_AVP_pg.mL),na.rm=TRUE))
data$Plasma_pg.mL = scale(log(data$Plasma_pg.mL),
                           scale= 2 * sd(log(data$Plasma_pg.mL), na.rm=TRUE))
data$uOT.Creatinine = scale(log(data$uOT.Creatinine),
                             scale= 2 * sd(log(data$uOT.Creatinine), na.rm=TRUE))
data$std.Age = scale(data$Age, scale= 2 * sd(data$Age))

#force negative duration values to 0 (single value from subtracted score)
data$Prox_duration[data$Prox_duration<0] = 0

#check proportion of zeroes to assess whether
#zero-inflation or hurdle models are necessary
sum(data$Prox_event > 0)/length(data$Prox_event)
sum(data$Close_event > 0)/length(data$Close_event)
sum(data$Groom_event > 0)/length(data$Groom_event)
sum(data$Dominance_tot > 0)/length(data$Dominance_tot)
sum(data$Submission_tot > 0)/length(data$Submission_tot)
sum(data$Self_directed > 0)/length(data$Self_directed)

#add small integer to avoid unnecessary hurdle Gamma for Prox_duration
data$Prox_duration[data$Prox_duration==0]=0+0.001

####################################################################
#descriptive raw individual data plot
{
temp = subset(data, select=c(Subject, Prox_event:Self_directed))
aggr = aggregate(. ~ Subject, temp, mean)
aggr.sd = aggregate(. ~ Subject, temp, sd)
aggr.l=gather(aggr[,-1], value="mean")
aggr.l_sd=gather(aggr.sd[,-1], value="sd")
aggr.l$sd = aggr.l_sd$sd

aggr.l$Subj = rep(c(1:70), 9)
#trait labels
trait = c("Prox_event", "Prox_duration", "Close_event", "Close_duration",
           "Groom_event", "Groom_duration", "Dominance", "Submission", "Self_directed")

aggr.l$trait = factor(rep(trait, each=70),levels=trait)

#plot
raw.plot=
  ggplot(aggr.l, aes(x = Subj, y = mean, color=trait)) +
  geom_point(size=2) + geom_errorbar(aes(ymin = mean - sd, ymax= mean + sd),
                                      width=0, size=1)+
  coord_cartesian(ylim=c(0,25))+
  facet_wrap(. ~ trait)+  
  xlab("\n Subject ID (mean + SD)")+
  theme(plot.title = element_text(size=12, hjust=0.5),
        axis.ticks.y= element_blank(),
        axis.title.y= element_blank(),
        axis.ticks.x= element_blank(),
        axis.title.x= element_text(size=12),
        axis.text.x= element_text(size=10),
        axis.text.y= element_text(size=9),
        plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        strip.text=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(1.5, "lines"),
        strip.text.y.right = element_text(angle=0))+
  guides(color=FALSE)

#save
png(file="Raw data subject means.png", res=600, width=13, height=8.5, units="in")  
plot(raw.plot)
dev.off()
}

####################################################################
################## Multivariate Model  #############################
####################################################################

#Basic model without subspecies-specific correlations

#Model formulas and prior
{
proxe.f = 
  bf(Prox_event ~ std.Age + Sex + Subspecies + TimeofDay + 
       (1|cor|Subject) + (1|Corral) + (1|cor2|Observation)) + poisson()

proxd.f = 
  bf(Prox_duration ~ std.Age + Sex + Subspecies + TimeofDay + 
       (1|cor|Subject) + (1|Corral) + (1|cor2|Observation)) + Gamma(link="log")

close.f = 
  bf(Close_event ~ std.Age + Sex + Subspecies + TimeofDay + 
       (1|cor|Subject) + (1|Corral) + (1|cor2|Observation)) + poisson()

closd.f = 
  bf(Close_duration ~ std.Age + Sex + Subspecies + TimeofDay + 
       (1|cor|Subject) + (1|Corral) + (1|cor2|Observation)) + hurdle_gamma(link="log")

groome.f = 
  bf(Groom_event ~ std.Age + Sex + Subspecies + TimeofDay + 
       (1|cor|Subject) + (1|Corral) + (1|cor2|Observation)) + poisson()

groomd.f = 
  bf(Groom_duration ~ std.Age + Sex + Subspecies + TimeofDay + 
       (1|cor|Subject) + (1|Corral) + (1|cor2|Observation)) + hurdle_gamma(link="log")

dome.f = 
  bf(Dominance_tot ~ std.Age + Sex + Subspecies + TimeofDay + 
       (1|cor|Subject) + (1|Corral) + (1|cor2|Observation)) + poisson()

sube.f = 
  bf(Submission_tot ~ std.Age + Sex + Subspecies + TimeofDay + 
       (1|cor|Subject) + (1|Corral) + (1|cor2|Observation)) + poisson()

selfe.f = 
  bf(Self_directed ~ std.Age + Sex + Subspecies + TimeofDay + 
       (1|cor|Subject) + (1|Corral) + (1|cor2|Observation)) + poisson()

#multivariate prior
prior = c(prior(normal(0,1), class=Intercept, resp="Proxevent"),
          prior(normal(0,1), class=b, resp="Proxevent"),
          prior(cauchy(0,1), class=sd, resp="Proxevent"),
          
          prior(normal(0,1), class=Intercept, resp="Proxduration"),
          prior(normal(0,1), class=b, resp="Proxduration"),
          prior(cauchy(0,1), class=sd, resp="Proxduration"),
          prior(constant(1), class=shape, resp="Proxduration"),
          
          prior(normal(0,1), class=Intercept, resp="Closeevent"),
          prior(normal(0,1), class=b, resp="Closeevent"),
          prior(cauchy(0,1), class=sd, resp="Closeevent"),
          
          prior(normal(0,1), class=Intercept, resp="Closeduration"),
          prior(normal(0,1), class=b, resp="Closeduration"),
          prior(cauchy(0,1), class=sd, resp="Closeduration"),
          prior(constant(1), class=shape, resp="Closeduration"),
          
          prior(normal(0,1), class=Intercept, resp="Groomevent"),
          prior(normal(0,1), class=b, resp="Groomevent"),
          prior(cauchy(0,1), class=sd, resp="Groomevent"),
          
          prior(normal(0,1), class=Intercept, resp="Groomduration"),
          prior(normal(0,1), class=b, resp="Groomduration"),
          prior(cauchy(0,1), class=sd, resp="Groomduration"),
          prior(constant(1), class=shape, resp="Groomduration"),
          
          prior(normal(0,1), class=Intercept, resp="Dominancetot"),
          prior(normal(0,1), class=b, resp="Dominancetot"),
          prior(cauchy(0,1), class=sd, resp="Dominancetot"),
          
          prior(normal(0,1), class=Intercept, resp="Submissiontot"),
          prior(normal(0,1), class=b, resp="Submissiontot"),
          prior(cauchy(0,1), class=sd, resp="Submissiontot"),
            
          prior(normal(0,1), class=Intercept, resp="Selfdirected"),
          prior(normal(0,1), class=b, resp="Selfdirected"),
          prior(cauchy(0,1), class=sd, resp="Selfdirected"),
        
          prior(lkj(3), class=cor))
}

#estimate multivariate model
mv.m = brm(proxe.f + proxd.f + close.f + closd.f + groome.f + groomd.f +
            dome.f + sube.f + selfe.f + set_rescor(FALSE), data=data,
            prior=prior, warmup=1000, iter=3000, chains=8, seed=99, init=0.001,
            control=list(adapt_delta=0.99, max_treedepth=20))

saveRDS(mv.m, "mv_m.RDS")
mv.m = readRDS("mv_m.RDS")
summary(mv.m, prob=0.9)

####################################################################

#Model with subspecies correlations

#Model formulas and prior
{
  proxe.f = 
    bf(Prox_event ~ std.Age + Sex + Subspecies + TimeofDay + 
         (1|cor|gr(Subject, by=Subspecies)) + (1|Corral) + (1|cor2|gr(Observation, by=Subspecies)) ) + poisson()
  
  proxd.f = 
    bf(Prox_duration ~ std.Age + Sex + Subspecies + TimeofDay + 
         (1|cor|gr(Subject, by=Subspecies)) + (1|Corral) + (1|cor2|gr(Observation, by=Subspecies)) ) + Gamma(link="log")
  
  close.f = 
    bf(Close_event ~ std.Age + Sex + Subspecies + TimeofDay + 
         (1|cor|gr(Subject,by=Subspecies)) + (1|Corral) + (1|cor2|gr(Observation, by=Subspecies)) ) + poisson()
  
  closd.f = 
    bf(Close_duration ~ std.Age + Sex + Subspecies + TimeofDay + 
         (1|cor|gr(Subject,by=Subspecies)) + (1|Corral) + (1|cor2|gr(Observation, by=Subspecies)) ) + hurdle_gamma(link="log")
  
  groome.f = 
    bf(Groom_event ~ std.Age + Sex + Subspecies + TimeofDay + 
         (1|cor|gr(Subject,by=Subspecies)) + (1|Corral) + (1|cor2|gr(Observation, by=Subspecies)) ) + poisson()
  
  groomd.f = 
    bf(Groom_duration ~ std.Age + Sex + Subspecies + TimeofDay + 
         (1|cor|gr(Subject,by=Subspecies)) + (1|Corral) + (1|cor2|gr(Observation, by=Subspecies)) ) + hurdle_gamma(link="log")
  
  dome.f = 
    bf(Dominance_tot ~ std.Age + Sex + Subspecies + TimeofDay + 
         (1|cor|gr(Subject,by=Subspecies)) + (1|Corral) + (1|cor2|gr(Observation, by=Subspecies)) ) + poisson()
  
  sube.f = 
    bf(Submission_tot ~ std.Age + Sex + Subspecies + TimeofDay + 
         (1|cor|gr(Subject,by=Subspecies)) + (1|Corral) + (1|cor2|gr(Observation, by=Subspecies)) ) + poisson()
  
  selfe.f = 
    bf(Self_directed ~ std.Age + Sex + Subspecies + TimeofDay + 
         (1|cor|gr(Subject,by=Subspecies)) + (1|Corral) + (1|cor2|gr(Observation, by=Subspecies)) ) + poisson()
  
  #multivariate prior
  prior = c(prior(normal(0,1), class=Intercept, resp="Proxevent"),
            prior(normal(0,1), class=b, resp="Proxevent"),
            prior(cauchy(0,1), class=sd, resp="Proxevent"),
            
            prior(normal(0,1), class=Intercept, resp="Proxduration"),
            prior(normal(0,1), class=b, resp="Proxduration"),
            prior(cauchy(0,1), class=sd, resp="Proxduration"),
            prior(constant(1), class=shape, resp="Closeduration"),
            
            prior(normal(0,1), class=Intercept, resp="Closeevent"),
            prior(normal(0,1), class=b, resp="Closeevent"),
            prior(cauchy(0,1), class=sd, resp="Closeevent"),
            
            prior(normal(0,1), class=Intercept, resp="Closeduration"),
            prior(normal(0,1), class=b, resp="Closeduration"),
            prior(cauchy(0,1), class=sd, resp="Closeduration"),
            prior(constant(1), class=shape, resp="Closeduration"),
            
            prior(normal(0,1), class=Intercept, resp="Groomevent"),
            prior(normal(0,1), class=b, resp="Groomevent"),
            prior(cauchy(0,1), class=sd, resp="Groomevent"),
            
            prior(normal(0,1), class=Intercept, resp="Groomduration"),
            prior(normal(0,1), class=b, resp="Groomduration"),
            prior(cauchy(0,1), class=sd, resp="Groomduration"),
            prior(constant(1), class=shape, resp="Groomduration"),
            
            prior(normal(0,1), class=Intercept, resp="Dominancetot"),
            prior(normal(0,1), class=b, resp="Dominancetot"),
            prior(cauchy(0,1), class=sd, resp="Dominancetot"),
            
            prior(normal(0,1), class=Intercept, resp="Submissiontot"),
            prior(normal(0,1), class=b, resp="Submissiontot"),
            prior(cauchy(0,1), class=sd, resp="Submissiontot"),
            
            prior(normal(0,1), class=Intercept, resp="Selfdirected"),
            prior(normal(0,1), class=b, resp="Selfdirected"),
            prior(cauchy(0,1), class=sd, resp="Selfdirected"),
            
            prior(lkj(3), class=cor))
  
}  

#estimate multivariate model
mv.m_sub = brm(proxe.f + proxd.f + close.f + closd.f + groome.f + groomd.f +
                  dome.f + sube.f + selfe.f + set_rescor(FALSE), data=data,
                  prior=prior, warmup=1000, iter=3000, chains=8, seed=99, init=0.001,
                  control=list(adapt_delta=0.99, max_treedepth=20)) 

saveRDS(mv.m_sub, "mv_m_sub.RDS")
mv.m_sub = readRDS("mv_m_sub.RDS")
summary(mv.m_sub, prob=0.9)

####################################################################
####################################################################
#Model comparison
compare_ic(WAIC(mv.m), WAIC(mv.m_sub))

####################################################################
#Subspecies correlations + hormones measures

#Model formulas and prior
{
  
  #missing data imputation models
  
  CSFpg.mi = bf(CSF_pg.mL|mi() ~ mi(CSF_AVP_pg.mL) + mi(Plasma_pg.mL) + mi(uOT.Creatinine) +
                   std.Age + Sex + Subspecies + TimeofDay +  (1|Corral) ) + gaussian()
  
  CSF_AVP_pg.mi = bf(CSF_AVP_pg.mL|mi() ~ mi(CSF_pg.mL) + mi(Plasma_pg.mL) + mi(uOT.Creatinine) +
                   std.Age + Sex + Subspecies + TimeofDay +  (1|Corral) ) + gaussian()
  
  Plasma_pg.mi = bf(Plasma_pg.mL|mi() ~ mi(CSF_AVP_pg.mL) + mi(CSF_pg.mL) + mi(uOT.Creatinine) +
                   std.Age + Sex + Subspecies + TimeofDay +  (1|Corral) ) + gaussian()
  
  uOT.Creatinine.mi = bf(uOT.Creatinine|mi() ~ mi(CSF_AVP_pg.mL) + mi(CSF_pg.mL) + mi(Plasma_pg.mL) +
                   std.Age + Sex + Subspecies + TimeofDay +  (1|Corral) ) + gaussian()
  
  #behavioral response models
  proxe.f = 
    bf(Prox_event ~ mi(CSF_pg.mL) + mi(CSF_AVP_pg.mL) + mi(Plasma_pg.mL) + mi(uOT.Creatinine) +
         std.Age + Sex + Subspecies + TimeofDay + (1|cor|gr(Subject, by=Subspecies)) + (1|Corral) + (1|cor2|gr(Observation, by=Subspecies)) ) + poisson(link="log")
  
  proxd.f = 
    bf(Prox_duration ~ mi(CSF_pg.mL) + mi(CSF_AVP_pg.mL) + mi(Plasma_pg.mL) + mi(uOT.Creatinine) +
         std.Age + Sex + Subspecies + TimeofDay + (1|cor|gr(Subject, by=Subspecies)) + (1|Corral) + (1|cor2|gr(Observation, by=Subspecies)) ) + Gamma(link="log")
  
  close.f = 
    bf(Close_event ~ mi(CSF_pg.mL) + mi(CSF_AVP_pg.mL) + mi(Plasma_pg.mL) + mi(uOT.Creatinine) +
         std.Age + Sex + Subspecies + TimeofDay + (1|cor|gr(Subject,by=Subspecies)) + (1|Corral) + (1|cor2|gr(Observation, by=Subspecies)) ) + poisson(link="log")
  
  closd.f = 
    bf(Close_duration ~ mi(CSF_pg.mL) + mi(CSF_AVP_pg.mL) + mi(Plasma_pg.mL) + mi(uOT.Creatinine) +
         std.Age + Sex + Subspecies + TimeofDay + (1|cor|gr(Subject,by=Subspecies)) + (1|Corral) + (1|cor2|gr(Observation, by=Subspecies)) ) + hurdle_gamma(link="log")
  
  groome.f = 
    bf(Groom_event ~ mi(CSF_pg.mL) + mi(CSF_AVP_pg.mL) + mi(Plasma_pg.mL) + mi(uOT.Creatinine) +
         std.Age + Sex + Subspecies + TimeofDay + (1|cor|gr(Subject,by=Subspecies)) + (1|Corral) + (1|cor2|gr(Observation, by=Subspecies)) ) + poisson(link="log")
  
  groomd.f = 
    bf(Groom_duration ~ mi(CSF_pg.mL) + mi(CSF_AVP_pg.mL) + mi(Plasma_pg.mL) + mi(uOT.Creatinine) +
         std.Age + Sex + Subspecies + TimeofDay + (1|cor|gr(Subject,by=Subspecies)) + (1|Corral) + (1|cor2|gr(Observation, by=Subspecies)) ) + hurdle_gamma(link="log")
  
  dome.f = 
    bf(Dominance_tot ~ mi(CSF_pg.mL) + mi(CSF_AVP_pg.mL) + mi(Plasma_pg.mL) + mi(uOT.Creatinine) +
         std.Age + Sex + Subspecies + TimeofDay + (1|cor|gr(Subject,by=Subspecies)) + (1|Corral) + (1|cor2|gr(Observation, by=Subspecies)) ) + poisson(link="log")
  
  sube.f = 
    bf(Submission_tot ~ mi(CSF_pg.mL) + mi(CSF_AVP_pg.mL) + mi(Plasma_pg.mL) + mi(uOT.Creatinine) +
         std.Age + Sex + Subspecies + TimeofDay + (1|cor|gr(Subject,by=Subspecies)) + (1|Corral) + (1|cor2|gr(Observation, by=Subspecies)) ) + poisson(link="log")
  
  selfe.f = 
    bf(Self_directed ~ mi(CSF_pg.mL) + mi(CSF_AVP_pg.mL) + mi(Plasma_pg.mL) + mi(uOT.Creatinine) +
         std.Age + Sex + Subspecies + TimeofDay + (1|cor|gr(Subject,by=Subspecies)) + (1|Corral) + (1|cor2|gr(Observation, by=Subspecies)) ) + poisson(link="log")
  
  #multivariate prior
  prior = c(prior(normal(0,1), class=Intercept, resp="CSFpgmL"),
            prior(normal(0,1), class=b, resp="CSFpgmL"),
            prior(cauchy(0,1), class=sd, resp="CSFpgmL"),
            prior(cauchy(0,1), class=sigma, resp="CSFpgmL"),
            
            prior(normal(0,1), class=Intercept, resp="CSFAVPpgmL"),
            prior(normal(0,1), class=b, resp="CSFAVPpgmL"),
            prior(cauchy(0,1), class=sd, resp="CSFAVPpgmL"),
            prior(cauchy(0,1), class=sigma, resp="CSFAVPpgmL"),
            
            prior(normal(0,1), class=Intercept, resp="PlasmapgmL"),
            prior(normal(0,1), class=b, resp="PlasmapgmL"),
            prior(cauchy(0,1), class=sd, resp="PlasmapgmL"),
            prior(cauchy(0,1), class=sigma, resp="PlasmapgmL"),
            
            prior(normal(0,1), class=Intercept, resp="uOTCreatinine"),
            prior(normal(0,1), class=b, resp="uOTCreatinine"),
            prior(cauchy(0,1), class=sd, resp="uOTCreatinine"),
            prior(cauchy(0,1), class=sigma, resp="uOTCreatinine"),
    
            prior(normal(0,1), class=Intercept, resp="Proxevent"),
            prior(normal(0,1), class=b, resp="Proxevent"),
            prior(cauchy(0,1), class=sd, resp="Proxevent"),
            
            prior(normal(0,1), class=Intercept, resp="Proxduration"),
            prior(normal(0,1), class=b, resp="Proxduration"),
            prior(cauchy(0,1), class=sd, resp="Proxduration"),
            prior(constant(1), class=shape, resp="Proxduration"),
            
            prior(normal(0,1), class=Intercept, resp="Closeevent"),
            prior(normal(0,1), class=b, resp="Closeevent"),
            prior(cauchy(0,1), class=sd, resp="Closeevent"),
            
            prior(normal(0,1), class=Intercept, resp="Closeduration"),
            prior(normal(0,1), class=b, resp="Closeduration"),
            prior(cauchy(0,1), class=sd, resp="Closeduration"),
            prior(constant(1), class=shape, resp="Closeduration"),
            
            prior(normal(0,1), class=Intercept, resp="Groomevent"),
            prior(normal(0,1), class=b, resp="Groomevent"),
            prior(cauchy(0,1), class=sd, resp="Groomevent"),
            
            prior(normal(0,1), class=Intercept, resp="Groomduration"),
            prior(normal(0,1), class=b, resp="Groomduration"),
            prior(cauchy(0,1), class=sd, resp="Groomduration"),
            prior(constant(1), class=shape, resp="Groomduration"),
            
            prior(normal(0,1), class=Intercept, resp="Dominancetot"),
            prior(normal(0,1), class=b, resp="Dominancetot"),
            prior(cauchy(0,1), class=sd, resp="Dominancetot"),
            
            prior(normal(0,1), class=Intercept, resp="Submissiontot"),
            prior(normal(0,1), class=b, resp="Submissiontot"),
            prior(cauchy(0,1), class=sd, resp="Submissiontot"),
            
            prior(normal(0,1), class=Intercept, resp="Selfdirected"),
            prior(normal(0,1), class=b, resp="Selfdirected"),
            prior(cauchy(0,1), class=sd, resp="Selfdirected"),
            
            prior(lkj(3), class=cor))
  
}  

#estimate multivariate model
mv.m_subhormones = brm( CSFpg.mi  + CSF_AVP_pg.mi + Plasma_pg.mi + uOT.Creatinine.mi +
                           proxe.f + proxd.f + close.f + closd.f + groome.f + groomd.f +
                           dome.f + sube.f + selfe.f + set_rescor(FALSE), data=data,
                  prior=prior, warmup=1000, iter=3000, chains=8, seed=99, init=0.001,
                  control=list(adapt_delta=0.99, max_treedepth=20)) 

saveRDS(mv.m_subhormones, "mv_m_subhormones.RDS")
mv.m_subhormones = readRDS("mv_m_subhormones.RDS")
summary(mv.m_subhormones, prob=0.9)

####################################################################
####################################################################
#Repeatabilities (pred w/ - pred w/o Subject random effects)
p1 = fitted(mv.m, resp = "Proxevent", summary = F, scale = "linear",
            re_formula = ~ (1|Subject) +  (1|Corral) + (1|Observation))
p1n = fitted(mv.m, resp = "Proxevent", summary = F, scale = "linear",
            re_formula = ~ (1|Corral) + (1|Observation))
r = as.vector(NA)
for(i in 1:nrow(p1)){r[[i]] = (var(p1[i,]) - var(p1n[i,])) / var(p1[i,])}
median(r); quantile(r, c(0.05,0.95))


p1 = fitted(mv.m, resp = "Proxduration", summary = F, scale = "linear",
            re_formula = ~ (1|Subject) +  (1|Corral) + (1|Observation))
p1n = fitted(mv.m, resp = "Proxduration", summary = F, scale = "linear",
            re_formula = ~ (1|Corral) + (1|Observation))
r = as.vector(NA)
for(i in 1:nrow(p1)){r[[i]] = (var(p1[i,]) - var(p1n[i,])) / var(p1[i,])}
median(r); quantile(r, c(0.05,0.95))


p1 = fitted(mv.m, resp = "Closeevent", summary = F, scale = "linear",
            re_formula = ~ (1|Subject) +  (1|Corral) + (1|Observation))
p1n = fitted(mv.m, resp = "Closeevent", summary = F, scale = "linear",
            re_formula = ~ (1|Corral) + (1|Observation))
r = as.vector(NA)
for(i in 1:nrow(p1)){r[[i]] = (var(p1[i,]) - var(p1n[i,])) / var(p1[i,])}
median(r); quantile(r, c(0.05,0.95))

p1 = fitted(mv.m, resp = "Closeduration", summary = F, scale = "linear",
            re_formula = ~ (1|Subject) +  (1|Corral) + (1|Observation))
p1n = fitted(mv.m, resp = "Closeduration", summary = F, scale = "linear",
            re_formula = ~ (1|Corral) + (1|Observation))
r = as.vector(NA)
for(i in 1:nrow(p1)){r[[i]] = (var(p1[i,]) - var(p1n[i,])) / var(p1[i,])}
median(r); quantile(r, c(0.05,0.95))

p1 = fitted(mv.m, resp = "Groomevent", summary = F, scale = "linear",
            re_formula = ~ (1|Subject) +  (1|Corral) + (1|Observation))
p1n = fitted(mv.m, resp = "Groomevent", summary = F, scale = "linear",
            re_formula = ~ (1|Corral) + (1|Observation))
r = as.vector(NA)
for(i in 1:nrow(p1)){r[[i]] = (var(p1[i,]) - var(p1n[i,])) / var(p1[i,])}
median(r); quantile(r, c(0.05,0.95))

p1 = fitted(mv.m, resp = "Dominancetot", summary = F, scale = "linear",
            re_formula = ~ (1|Subject) +  (1|Corral) + (1|Observation))
p1n = fitted(mv.m, resp = "Dominancetot", summary = F, scale = "linear",
            re_formula = ~ (1|Corral) + (1|Observation))
r = as.vector(NA)
for(i in 1:nrow(p1)){r[[i]] = (var(p1[i,]) - var(p1n[i,])) / var(p1[i,])}
median(r); quantile(r, c(0.05,0.95))

p1 = fitted(mv.m, resp = "Submissiontot", summary = F, scale = "linear",
            re_formula = ~ (1|Subject) +  (1|Corral) + (1|Observation))
p1n = fitted(mv.m, resp = "Submissiontot", summary = F, scale = "linear",
            re_formula = ~ (1|Corral) + (1|Observation))
r = as.vector(NA)
for(i in 1:nrow(p1)){r[[i]] = (var(p1[i,]) - var(p1n[i,])) / var(p1[i,])}
median(r); quantile(r, c(0.05,0.95))

p1 = fitted(mv.m, resp = "Selfdirected", summary = F, scale = "linear",
            re_formula = ~ (1|Subject) +  (1|Corral) + (1|Observation))
p1n = fitted(mv.m, resp = "Selfdirected", summary = F, scale = "linear",
            re_formula = ~ (1|Corral) + (1|Observation))
r = as.vector(NA)
for(i in 1:nrow(p1)){r[[i]] = (var(p1[i,]) - var(p1n[i,])) / var(p1[i,])}
median(r); quantile(r, c(0.05,0.95))


####################################################################
####################################################################
#labels
trait = c("Prox_event", "Prox_duration", "Close_event", "Close_duration",
           "Groom_event", "Groom_duration", "Dominance", "Submission", "Self_directed")

#Posterior correlation matrices (aggregated)
cor.post=VarCorr(mv.m)
cor.subj = cor.subj$cor
#among-individual
cor.tot = round(matrix(cor.subj[,1,], nrow=9, ncol=9, dimnames=list(trait,trait)),2)
cor.tot

#among-individual
cor.subj = cor.post$Subject


#Posterior correlation matrices x subspecies
cor.post=VarCorr(mv.m_sub)

#among-individual
cor.subj = cor.post$Subject
cor.subj.anub = cor.subj$cor[c(1:9),,c(1:9)]
cor.subj.hamad = cor.subj$cor[c(10:18),,c(10:18)]
cov.subj.anub = cor.subj$cov[c(1:9),,c(1:9)]
cov.subj.hamad = cor.subj$cov[c(10:18),,c(10:18)]

#within-individual
cor.obs = cor.post$Observation
cor.obs.anub = cor.obs$cor[c(1:9),,c(1:9)]
cor.obs.hamad = cor.obs$cor[c(10:18),,c(10:18)]
cov.obs.anub = cor.obs$cov[c(1:9),,c(1:9)]
cov.obs.hamad = cor.obs$cov[c(10:18),,c(10:18)]

#correlation network plots
{
library(qgraph)
library(psych)
  
#prep correlation matrices

#among-individual
cor.anub = round(matrix(cor.subj.anub[,1,], nrow=9, ncol=9, dimnames=list(trait,trait)),2)
diag(cor.anub) = 0 
cor.hamad = round(matrix(cor.subj.hamad[,1,], nrow=9, ncol=9, dimnames=list(trait,trait)),2)
diag(cor.hamad) = 0 
cor.diff = fisherz(cor.anub) - fisherz(cor.hamad)

#within-individual
cor.obs.anub = round(matrix(cor.obs.anub[,1,], nrow=9, ncol=9, dimnames=list(trait,trait)),2)
diag(cor.obs.anub) = 0 
cor.obs.hamad = round(matrix(cor.obs.hamad[,1,], nrow=9, ncol=9, dimnames=list(trait,trait)),2)
diag(cor.obs.hamad) = 0 
cor.obs.diff = fisherz(cor.obs.anub) - fisherz(cor.obs.hamad)

#network labels
abbr=c("PRX_E","PRX_D","CNT_E","CNT_D","GRM_E","GRM_D", "DOM", "SBM", "SLF")

colnames(cor.anub)=abbr;rownames(cor.anub)=abbr
colnames(cor.hamad)=abbr;rownames(cor.hamad)=abbr
colnames(cor.obs.anub)=abbr;rownames(cor.obs.anub)=abbr
colnames(cor.obs.hamad)=abbr;rownames(cor.obs.hamad)=abbr

#plot
png(file = "Baboon behavior correlations.png", res=600, width=11, height=8, units="in")
{
  layout(matrix(c(1,4,2,5,3,6),nrow=2,ncol=3))
    
  #between-subject
  btw.graph1=qgraph(cor.anub, maximum=1, layout="circle", esize=10,vsize=10,colFactor=0.4,
                     posCol="#35abd6", fade=TRUE, minimum=0.1, label.prop=0.5, 
                     edge.labels=TRUE, node.resolution=300, mar=c(2,2,4,2))
  title("Among-Subject Correlation (Anubis)",adj=1)
  
  btw.graph2=qgraph(cor.hamad, maximum=1, layout="circle", esize=10,vsize=10,colFactor=0.4,
                     posCol="#35abd6", fade=TRUE, minimum=0.1, label.prop=0.5, 
                     edge.labels=TRUE,node.resolution=600,mar=c(2,2,4,2))
  title("Among-Subject Correlation (Hamadryas)",adj=1)
  
  btw.graph3=qgraph(cor.diff, maximum=1, layout="circle", esize=10,vsize=10,colFactor=0.4,
                     posCol="#35abd6", fade=TRUE, minimum=0.1, label.prop=0.5, 
                     edge.labels=TRUE,node.resolution=600,mar=c(2,2,4,2))
  title("Among-Subject Cohen's q (Difference in Fisher's z)\n Anubis - Hamadryas",adj=1)
  
  #within-subject
  btw.graph4=qgraph(cor.obs.anub, maximum=1, layout="circle", esize=10,vsize=10,colFactor=0.4,
                     posCol="#35abd6", fade=TRUE, minimum=0.1, label.prop=0.5, 
                     edge.labels=TRUE, node.resolution=300, mar=c(2,2,4,2))
  title("Within-Subject Correlation (Anubis)",adj=1)
  
  btw.graph5=qgraph(cor.obs.hamad, maximum=1, layout="circle", esize=10,vsize=10,colFactor=0.4,
                     posCol="#35abd6", fade=TRUE, minimum=0.1, label.prop=0.5, 
                     edge.labels=TRUE,node.resolution=600,mar=c(2,2,4,2))
  title("Within-Subject Correlation (Hamadryas)",adj=1)
  
  btw.graph6=qgraph(cor.obs.diff, maximum=1, layout="circle", esize=10,vsize=10,colFactor=0.4,
                     posCol="#35abd6", fade=TRUE, minimum=0.1, label.prop=0.5, 
                     edge.labels=TRUE,node.resolution=600,mar=c(2,2,4,2))
  title("Within-Subject Cohen's q (Difference in Fisher's z)\n Anubis - Hamadryas",adj=1)
  }
  dev.off()
}

####################################################################
#random effect density plots

#organize SDs for ggplot
{

library(tidyr)
post= posterior_samples(mv.m_sub)

#Between-subject
df.subj = post[, grepl( "sd_Subject", colnames(post))]
df.subj.anub = df.subj[, grepl( "SubspeciesAnubis", colnames(df.subj))]
df.subj.hamad = df.subj[, grepl( "SubspeciesHamadryas", colnames(df.subj))]
colnames(df.subj.anub) = trait
colnames(df.subj.hamad) = trait

#Within-subject
df.obs = post[, grepl( "sd_Observation", colnames(post))]
df.obs.anub = df.obs[, grepl( "SubspeciesAnubis", colnames(df.obs))]
df.obs.hamad = df.obs[, grepl( "SubspeciesHamadryas", colnames(df.obs))]
colnames(df.obs.anub) = trait
colnames(df.obs.hamad) = trait

#Corral
df.corral = post[, grepl( "sd_Corral", colnames(post))]
colnames(df.corral) = trait

#wide to long
df.subj.anub  = gather(df.subj.anub, value="sd")
df.subj.anub$subspecies = "Anubis"
df.subj.anub$source = "Among-individual"

df.subj.hamad = gather(df.subj.hamad, value="sd")
df.subj.hamad$subspecies = "Hamadryas"
df.subj.hamad$source = "Among-individual"

df.obs.anub  = gather(df.obs.anub, value="sd")
df.obs.anub$subspecies = "Anubis"
df.obs.anub$source = "Within-individual"

df.obs.hamad = gather(df.obs.hamad, value="sd")
df.obs.hamad$subspecies = "Hamadryas"
df.obs.hamad$source = "Within-individual"

df.corral = gather(df.corral, value="sd")
df.corral$subspecies = "NA"
df.corral$source = "NA"

#gather together
longdf = rbind(df.subj.anub, df.subj.hamad,
                df.obs.anub, df.obs.hamad,
                df.corral)
}

#repeatability estimates


#plot
{
#change factor levels
longdf$subspecies = as.factor(longdf$subspecies)
levels(longdf$subspecies) = c("Anubis", "Hamadryas", "Corral")
longdf$trait = factor(longdf$key, levels=c("Prox_event", "Prox_duration",
                                            "Close_event", "Close_duration",
                                            "Groom_event", "Groom_duration",
                                            "Dominance", "Submission", "Self_directed"))
#plot
dens.plot=
  ggplot(longdf, aes(x = sd, fill = interaction(source),
                    group=interaction(subspecies, source), alpha=0.3)) +
  geom_density(aes(y=..scaled..))+ facet_grid(trait ~ subspecies) + 
  scale_fill_manual(values=c("#e8052b","#7305e8","#057ae8"),
                    labels=c("Among-individual", "Among-corral", "Within-individual"))+
  xlab("\nStandard Deviation of Random Intercepts")+
  ylab("Posterior density\n")+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(xlim=c(0,2.5))+
  theme(plot.title =element_text(size=12, hjust=0.5),
        legend.text =element_text(size=9),
        legend.title=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_text(size=12),
        axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=12),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        strip.text=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(1.5, "lines"),
        strip.text.y.right = element_text(angle=0))+
  guides(alpha=FALSE)

#save
png(file="Random effect sd.png", res=600, width=11, height=12, units="in")  
plot(dens.plot)
dev.off()

}


####################################################################
#fixed effect plots
fe = data.frame(fixef(mv.m_subhormones, robust=TRUE, probs=c(0.05, 0.95)))

#organize effects
{
#labels
trait = c("Prox_event", "Prox_duration", "Close_event", "Close_duration",
             "Groom_event", "Groom_duration", "Dominance", "Submission", "Self_directed")  
  
  
df.sex = fe[grepl( "SexMale", rownames(fe)), ]
df.sex = df.sex[-c(1:4),]
df.sex$trait = trait
df.sex$parameter = "Sex_Male"

df.subspecies = fe[grepl( "Subspecies", rownames(fe)), ]
df.subspecies = df.subspecies[-c(1:4),]
df.subspecies$trait = trait
df.subspecies$parameter = "Subspecies_Hamadryas"

df.time = fe[grepl( "TimeofDay", rownames(fe)), ]
df.time = df.time[-c(1:4),]
df.time$trait = trait
df.time$parameter = "Time of day_Morning"

df.CSF_AVP = fe[grepl( "miCSF_AVP", rownames(fe)), ]
df.CSF_AVP = df.CSF_AVP[-c(1:3),]
df.CSF_AVP$trait = trait
df.CSF_AVP$parameter = "CSF_AVP"

df.CSF_pg = fe[grepl( "miCSF_pg", rownames(fe)), ]
df.CSF_pg = df.CSF_pg[-c(1:3),]
df.CSF_pg$trait = trait
df.CSF_pg$parameter = "CSF_pg"

df.Plasma_pg = fe[grepl( "miPlasma_pg", rownames(fe)), ]
df.Plasma_pg = df.Plasma_pg[-c(1:3),]
df.Plasma_pg$trait = trait
df.Plasma_pg$parameter = "Plasma_pg"

df.uOT.Creatinine = fe[grepl( "uOT.Creatinine", rownames(fe)), ]
df.uOT.Creatinine = df.uOT.Creatinine[-c(1:3),]
df.uOT.Creatinine$trait = trait
df.uOT.Creatinine$parameter = "uOT.Creatinine"


fe.df = rbind(df.sex, df.subspecies, df.time,
               df.CSF_AVP, df.CSF_pg, df.Plasma_pg, df.uOT.Creatinine)  
}

#plot
{

#relevel  
fe.df$trait=factor(fe.df$trait, levels=trait)  
fe.df$parameter =factor(fe.df$parameter,
                         levels=c("CSF_AVP", "CSF_pg", "Plasma_pg", "uOT.Creatinine",
                                  "Time of day_Morning", "Sex_Male", "Subspecies_Hamadryas"))

#plot
fe.plot=
  ggplot(fe.df, aes(x = Estimate, y = parameter, color=trait)) +
  geom_point(size=4) + geom_errorbarh(aes(xmin = Q5, xmax=Q95),
                                          height=0, size=1.5)+
  geom_vline(xintercept=0, linetype="dashed")+
  facet_wrap(. ~ trait)+  
  xlab("\n Log Scale Beta Coefficients (median & 90% CI)")+
  theme(plot.title = element_text(size=12, hjust=0.5),
          axis.ticks.y= element_blank(),
          axis.title.y= element_blank(),
          axis.ticks.x= element_blank(),
          axis.title.x= element_text(size=12),
          axis.text.x= element_text(size=10),
          axis.text.y= element_text(size=9),
          plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          strip.text=element_text(size=12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.spacing = unit(1.5, "lines"),
          strip.text.y.right = element_text(angle=0))+
    guides(alpha=FALSE, color=FALSE)
  
  #save plot
  png(file="Fixed effects_with hormones.png", res=600, width=12, height=11, units="in")  
  plot(fe.plot)
  dev.off()
  
  #data frame
  write.csv(fe, "fixed effects_with hormones.csv")
  
}

#organize effects (just hormones)
{
  #labels
  trait = c("Prox_event", "Prox_duration", "Close_event", "Close_duration",
             "Groom_event", "Groom_duration", "Dominance", "Submission", "Self_directed")  
  
  df.CSF_AVP = fe[grepl( "miCSF_AVP", rownames(fe)), ]
  df.CSF_AVP = df.CSF_AVP[-c(1:3),]
  df.CSF_AVP$trait = trait
  df.CSF_AVP$parameter = "CSF_AVP"
  
  df.CSF_pg = fe[grepl( "miCSF_pg", rownames(fe)), ]
  df.CSF_pg = df.CSF_pg[-c(1:3),]
  df.CSF_pg$trait = trait
  df.CSF_pg$parameter = "CSF_pg"
  
  df.Plasma_pg = fe[grepl( "miPlasma_pg", rownames(fe)), ]
  df.Plasma_pg = df.Plasma_pg[-c(1:3),]
  df.Plasma_pg$trait = trait
  df.Plasma_pg$parameter = "Plasma_pg"
  
  df.uOT.Creatinine = fe[grepl( "uOT.Creatinine", rownames(fe)), ]
  df.uOT.Creatinine = df.uOT.Creatinine[-c(1:3),]
  df.uOT.Creatinine$trait = trait
  df.uOT.Creatinine$parameter = "uOT.Creatinine"
  
  
  fe.df = rbind(df.CSF_AVP, df.CSF_pg, df.Plasma_pg, df.uOT.Creatinine)  
}

#plot
{
  #relevel  
  fe.df$trait=factor(fe.df$trait, levels=trait)  
  fe.df$parameter =factor(fe.df$parameter,
                           levels=c("CSF_AVP", "CSF_pg", "Plasma_pg", "uOT.Creatinine"))
  
  #plot
  fe.plot=
    ggplot(fe.df, aes(x = Estimate, y = parameter, color=trait)) +
    geom_point(size=4) + geom_errorbarh(aes(xmin = Q5, xmax=Q95),
                                        height=0, size=1.5)+
    geom_vline(xintercept=0, linetype="dashed")+
    facet_wrap(. ~ trait)+  
    xlab("\n Log Scale Beta Coefficients (median & 90% CI)")+
    theme(plot.title = element_text(size=12, hjust=0.5),
          axis.ticks.y= element_blank(),
          axis.title.y= element_blank(),
          axis.ticks.x= element_blank(),
          axis.title.x= element_text(size=12),
          axis.text.x= element_text(size=10),
          axis.text.y= element_text(size=9),
          plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          strip.text=element_text(size=12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.spacing = unit(1.5, "lines"),
          strip.text.y.right = element_text(angle=0))+
    guides(alpha=FALSE, color=FALSE)
  
  #save plot
  png(file="Fixed effects_only hormones.png", res=600, width=12, height=11, units="in")  
  plot(fe.plot)
  dev.off()
  
  #data frame
  write.csv(fe, "fixed effects_with hormones.csv")
  
}

#####################################################################

#If you have any questions or concerns please contact
#Jordan Scott Martin (jordan.martin@uzh.ch)

