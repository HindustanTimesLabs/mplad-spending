#################
#               #
#  HT INSIGHT   #
#  21 FEB 2018  # 
#               #
#################

rm(list = ls())

library(lme4)
library(MASS)

setwd("/Users/nsircar/Dropbox/Papers/HTAnalysis/")

dat <- read.csv("/Users/nsircar/Dropbox/Papers/HTAnalysis/MPLADS_ADR_PRS_15thLS.csv", stringsAsFactors = F)

## Measures of Unspent Funds

unspent <- 1 - dat$Expenditure/dat$AvailablewInterest

unspent[unspent < 0] <- 0

unspent0 <- as.numeric(unspent == 0)  ## All Funds Spent
unspent5 <- as.numeric(unspent < .05 ) ## Spent More than 95%
unspent10 <- as.numeric(unspent < .1)  ## Spent More than 90%
unspent20 <- as.numeric(unspent < .2) ## Spent More than 80%

invlogit <- function(x) 1/(1 + exp(-x))  # Inverse of Logit Function


sc <- as.numeric(dat$Category == "SC")  ## SC Constituency
st <- as.numeric(dat$Category == "ST") ## ST Constituency
pg <- as.numeric(dat$Educational.qualifications %in% c("Post Graduate", "Doctorate")) # Postgraduate Degree
crim <- as.numeric(dat$MajorCriminal == 1) ## Indicator of Major Criminal Cases (ADR)
state <- as.numeric(as.factor(dat$State)) ## State
party <- as.numeric(as.factor(paste(dat$Political.party, dat$State)))
logmoveable <- log10(dat$Moveable) ## Log (base 10) of moveable assets
keepval <- which(logmoveable > 0 & !is.na(logmoveable)) ## Subset of data to analyze
margin <- dat$MOV  ## Margin of Victory


fit1 <- glmer(unspent0[keepval]~logmoveable[keepval] + crim[keepval] + pg[keepval] + (1|party[keepval]), family = binomial(link = "logit")) ## Estimating Probability of Spending ALL Funds
fit2 <- glmer(unspent10[keepval]~logmoveable[keepval] + crim[keepval] + pg[keepval] + (1|party[keepval]), family = binomial(link = "logit")) ## Estimating Probability of Spending >90% Funds
fit3 <- glmer(unspent5[keepval]~logmoveable[keepval] + crim[keepval] + pg[keepval] + (1|party[keepval]), family = binomial(link = "logit")) ## Estimating Probability of Spending >95% Funds

coef.sims <- mvrnorm(1000, summary(fit3)$coef[,1], summary(fit3)$vcov)  ## Simulation from Last Model

## Estimated Probability of Effective Spending -- 10 Lakh and 50 Lakh Moveable Wealth

xval <- logmoveable[keepval][log10(1000000) <= logmoveable[keepval] & log10(5000000) >= logmoveable[keepval]]
xpred <- cbind(rep(1, length(xval)), xval, rep(mean(crim[keepval]), length(xval)), rep(mean(pg[keepval]), length(xval)))

mean(colMeans(invlogit((xpred %*% t(coef.sims)))))

## Estimated Probability of Effective Spending --  Greater than 1 Crore Moveable Wealth

xval <- logmoveable[keepval][log10(10000000) <= logmoveable[keepval] ]
xpred <- cbind(rep(1, length(xval)), xval, rep(mean(crim[keepval]), length(xval)), rep(mean(pg[keepval]), length(xval)))

mean(colMeans(invlogit((xpred %*% t(coef.sims)))))




crimval <- cbind(rep(1,2), rep(mean(logmoveable[keepval]), 2), c(0,1), rep(mean(pg[keepval]), 2)) ## X values for Crim Case Simulation -- Mean moveable assets, mean PG completion
pgval <- cbind(rep(1,2), rep(mean(logmoveable[keepval]), 2), rep(mean(crim[keepval]), 2),  c(0,1)) ## X values for PG degree Simulation -- Mean moveable assets, mean criminality



## GENERATE FIGURE 1


vals <- c(apply(invlogit((pgval %*% t(coef.sims))), 1, mean), apply(invlogit((crimval %*% t(coef.sims))), 1, mean)) ## Predicted Values from Simulations
vals.ci <- cbind(apply(invlogit((pgval %*% t(coef.sims))), 1, quantile, c(.05, .95)), apply(invlogit((crimval %*% t(coef.sims))), 1, quantile, c(.05, .95))) ## 90% intervals from simulations

pdf("/Users/nsircar/Dropbox/Papers/HTanalysis/unspent5_fig1.pdf", width=8, height=6)

par(mar = rep(0,4), mai = rep(0,4))
plot(0,0, xlim = c(-.85,7), ylim = c(-.08, .7), ann = F, axes = F, type = "n")
x.coords <- barplot(vals, col = c(rgb(1,125/255,125/255), rgb(1,125/255,125/255), rgb(102/255,178/255,1), rgb(102/255,178/255,1)), space = c(.5,.5,1,.5), axes = F, add = T )
axis(2, at=seq(0,.7,.1), labels = seq(0,.7,.1), pos = 0)

cats <- c("Graduate\nor Less", "Post\nGraduate", "Not Major\nCriminal", "Major\nCriminal")

cols <- c(rgb(204/255,0,0), rgb(204/255,0,0), rgb(0,0,204/255), rgb(0,0,204/255))
for (i in 1:4){
  
  text(x.coords[i],-.05, cats[i], cex = 1.2)
  arrows(x.coords[i], vals.ci[1,i], x.coords[i], vals.ci[2,i], col = cols[i], code = 3, angle = 90, length = .07, lwd = 2)
  text(x.coords[i]+.3, vals[i], round(vals[i], 2), pos = 3, cex = 1.5)
  
  
}

text(-.83, .35, "Estimated Probability of 95%+ MPLADS Utilisation", srt=90, cex = 1.2)

dev.off()





cutpoints <- quantile(logmoveable[keepval], seq(0,1,.05)) ## Generate Cutpoints for "Bins" every 5 Percentile Points


vec <- seq(5,8, .01) ## Values over which Simulated (Log) Moveable Assets Displayed

assetval <- cbind(rep(1, length(vec)), vec, rep(mean(crim[keepval]), length(vec)), rep(mean(pg[keepval]), length(vec))) ## Mean levels of PG degree completion and criminality

X <- cbind(rep(1, length(unspent5)), logmoveable, crim, pg)[keepval, ] ## Predictor Matrix


## Generating Binned Values. Bins refer to interval on the x axis between 2 cutpoints. 
## An observation is said to be in a bin if its x value falls between cutpoints of that bin.
## Binned X value = midpoint between two cutpoints
## Binned Y value = Model prediction + mean error of prediction of all observations in a bin
## The reason to do binning is that the true Y take 0 and 1 values, this procedure allows us to assess the predicted probability by taking the expectation in a bin


xbin5 <- ybin5 <- error <- rep(NA, length(cutpoints) - 1)  ## xbin5 = Binned X value; ybin5 = Binned Y value
for (i in 2:length(cutpoints)){
  xbin5[i-1] <- mean(cutpoints[(i-1):i])  
  
  
  xval <- logmoveable[keepval][(logmoveable[keepval] > cutpoints[i-1]) & (logmoveable[keepval] <= cutpoints[i]) ]
  yval <- unspent5[keepval][(logmoveable[keepval] > cutpoints[i-1]) & (logmoveable[keepval] <= cutpoints[i]) ]
  estval <- invlogit(X[(logmoveable[keepval] > cutpoints[i-1]) & (logmoveable[keepval] <= cutpoints[i]), ] %*% fixef(fit3))
  
  error[i-1] <- mean(yval - estval)
  
  fitval <- mean(invlogit(assetval[(assetval[,2] > cutpoints[i-1]) & (assetval[,2] <= cutpoints[i]), ] %*% fixef(fit3)))
  
  ybin5[i-1] <- fitval + error[i-1]  
}  


vals <- apply(invlogit((assetval %*% t(coef.sims))), 1, mean) ## Predicted Value from Simulations
vals.ci <- apply(invlogit((assetval %*% t(coef.sims))), 1, quantile, c(.05, .95)) ## 90% Intervals from Simulations

colsline <- rgb(153/255,255/255,153/255)
colsfill <- rgb(0,102/255,0)





pdf("/Users/nsircar/Dropbox/Papers/HTanalysis/unspent5_fig2.pdf", width=8, height=6)

par(mar = rep(0,4), mai = rep(0,4))
plot(0,0, xlim = c(min(vec) - .5, max(vec)+.1), ylim = c(-.2, 1), ann =F, axes = F, type = "n")


polygon(c(vec, rev(vec)), c(vals.ci[1,], rev(vals.ci[2,])), col = colsline,border = NA   )  ## Colored Interval Band
points(vec, vals, lwd = 2.5, type = "l", col = colsfill)  ## Predicted Value Curve
points(xbin5[2:19], ybin5[2:19], pch=19, col = colsfill, cex = 1.5) ## Plotting Points for Binned X and Y values

axis(1, at = seq(5, 8, 1), labels = c("1 Lakh", "10 Lakh", "1 Crore", "10 Crore"),  pos=0)  ## Axis on a Log (base 10) scale
axis(2, at = seq(0,1,.2), pos = 5)

text(6.5, -.18, "Reported Candidate Asset Wealth (Moveable)", cex = 1.2)
text(4.55, .5, "Estimated Probability of 95%+ MPLADS Utilisation", cex = 1.2, srt = 90)

dev.off()


## Unspent Funds by Party

partywise <- aggregate(unspent5[keepval], by = list(dat$Political.party[keepval]), length)
o <- order(partywise[,2])

## Consider parties with 15 seats and above

partyval <- partywise[o,][27:33,1]  

keepval2 <- intersect(keepval, which(dat$Political.party %in% partyval))

partywise2 <- aggregate(unspent5[keepval2], by = list(dat$Political.party[keepval2]), mean)

partywise_mean <- aggregate(unspent5[keepval], by = list(dat$Political.party[keepval]), mean)[o,]

partywise2 ## Percentage of MPs spending effectively by party (> 15 seats)