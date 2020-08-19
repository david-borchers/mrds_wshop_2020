library(Distance)
library(hmltm)

data("beaked.ship")

# need to rename come columns for use with ds
bkdsdat = beaked.ship
names(bkdsdat)[1:4] = c("Region.Label", "Area", "Sample.Label", "Effort")
names(bkdsdat)[which(names(bkdsdat)=="x")] = "distance"

head(bkdsdat)

conversion.factor <- convert_units("meter", "kilometer", "Square kilometer")

bk.hn <- ds(data=bkdsdat, key="hn", adjustment=NULL,convert.units=conversion.factor)
bk.hr <- ds(data=bkdsdat, key="hr", adjustment=NULL,convert.units=conversion.factor)
AIC(bk.hn,bk.hr)

summary(bk.hr)

cutpoints = seq(0,max(na.omit(bkdsdat$distance)),length=21)
plot(bk.hr, breaks=cutpoints, main="Half normal model, beaked shipboardl survey")
gof_ds(bk.hr)

cds.N = c(est=bk.hr$dht$individuals$N$Estimate[3],
             lcl=bk.hr$dht$individuals$N$lcl[3],
             ucl=bk.hr$dht$individuals$N$ucl[3])
cds.N

# Correction factors:
xcuts = seq(0,max(na.omit(beaked.ship$x)),length=21)
ycuts = seq(0,max(na.omit(beaked.ship$y)),length=21)
layout(matrix(1:4,nrow=2))
hist(beaked.ship$x,breaks=xcuts,xlab="Perpendicular Distance",main="")
plot(beaked.ship$x, beaked.ship$y, pch="+",xlab="Perpendicular Distance",ylab="Forward Distance")
plot(0,0,xaxt="n",yaxt="n",bty="n",xlab="",ylab="",type="n")
hist(beaked.ship$y,breaks=xcuts,xlab="Forward Distance",main="")

ymax = 4500

Eu=1580 # mean time UNavailable -per dive cycle
Ea=122 # mean time available per dive cycle
spd = 3.5 # ship speed in m/s
# Calculate time in view, in minutes
tiv = ymax/spd / 60 # minutes in view
tiv
Ea/60 # meand nuber minutes up (available)
Eu/60 # mean number minutes down (unavailable)

# Make availability parameters object
availpars = make.hmm.pars.from.Et(Eu=Eu, Ea=Ea, seEu=134.9, seEa=9.62)

simple.cf = simple.a(beaked.hmm.pars)
mclaren.cf = mclaren.a(beaked.hmm.pars,spd=spd, ymax=ymax)
laake.cf = laake.a(beaked.hmm.pars,spd=spd, ymax=ymax)

simple.cf = simple.a(availpars)
laake.cf = laake.a(availpars,spd=spd, ymax=ymax)
mclaren.cf = mclaren.a(availpars,spd=spd, ymax=ymax)
# put them in a list for convenience
cf = data.frame(simple=simple.cf$mean, 
                mclaren=mclaren.cf$mean, 
                laake=laake.cf$mean) 
cf

# "Corrected" total abundance estimates
bk.hr$dht$individuals$N$Estimate[3]/cf[,"simple"]
bk.hr$dht$individuals$N$Estimate[3]/cf[,"mclaren"]
bk.hr$dht$individuals$N$Estimate[3]/cf[,"laake"]

# Point estimates of total abundance and 95% CI with each correction method
simple.N = cds.N/cf[,"simple"]
mclaren.N = cds.N/cf[,"simple"]
laake.N = cds.N/cf[,"simple"]
Nests = t(data.frame(simple=round(simple.N),mclaren=round(mclaren.N),laake=round(laake.N)))
Nests



# Using Hidden Markov Model
# =========================

#data(beaked.hmm.pars) 

# set survey and fitting parameters
spd=3.5
Wl=0
W=4000
ymax=4500
survey.pars=make.survey.pars(spd=3.5,Wl=0,W=4000,ymax=4500) # specify survey parameters
control.fit=list(hessian=FALSE,nx=64) # fitting parameters
control.opt=list(trace=5,maxit=1000) # optimisation parameters
W.est=4000# set perpendicular truncation distance for Horvitz-Thompson-like estimator

# specify model and starting parameter values and fit model
hfun="h.EP2x.0";models=list(y=NULL,x=NULL) # specify detection hazard model (no covariates)
pars=c(0.45, 0.31, 7.7, 62.8) # detection hazard parameter starting values
# Point estimation:
bkEP2x.null=est.hmltm(beaked.ship,pars,hfun,models,survey.pars,availpars,
                      control.fit,control.opt,W.est=W.est)

# display point estimates of density and abundance:
bkEP2x.null$point$ests

# Look at goodness of fit
bkEP2x.null.fx=fxfit.plot(bkEP2x.null,breaks=xcuts,type="prob",allx=FALSE)
bkEP2x.null.fy=fyfit.plot(bkEP2x.null,breaks=ycuts,allx=FALSE,nys=250)
bkEP2x.null.gof.x=hmmlt.gof.x(bkEP2x.null);bkEP2x.null.gof.x$p.ks
bkEP2x.null.gof.y=hmmlt.gof.y(bkEP2x.null,breaks=ycuts);bkEP2x.null.gof.y$p.ks

# Bootstrap (treating availability parameter as fixed)
set.seed(1)
system.time(bests<-bs.hmltm(bkEP2x.null,B=50,hmm.pars.bs=availpars,bs.trace=0,report.by=10,fixed.avail=TRUE,bsfile=NULL))
bsum = bootsum(bests)
estable = strat.estable(bkEP2x.null$point$ests,bsum$cv)
hmltmNests = estable[3,c("N","N.lo","N.hi")]
names(hmltmNests) = colnames(Nests)
Nests = rbind(Nests,hmltm=hmltmNests)
Nests







data("aerial.survey")

# need to rename for use with ds
bwdsdat = aerial.survey
names(bwdsdat)[1:4] = c("Region.Label", "Area", "Sample.Label", "Effort")
names(bwdsdat)[which(names(bwdsdat)=="x")] = "distance"

head(bwdsdat)
head(wren_lt)

conversion.factor <- convert_units("meter", "kilometer", "Square kilometer")

bw.hn <- ds(data=bwdsdat, key="hn", adjustment=NULL,convert.units=conversion.factor)

summary(bw.hn)

cutpoints = seq(0,max(na.omit(bwdsdat$distance)),length=21)
plot(bw.hn, breaks=cutpoints, main="Half normal model, bowhead aerial survey")


# Correction factors:
data(bowhead.hmm.pars) # get bowhead availability HMM parameters (8 sets of parameters)
spd=46.3
ymax=2200

