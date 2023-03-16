# --------------------
# Title: Decay.models.plus.loop.R
# Author: TW modified by DBJ
# Date: 2021-Jan-7
#
# Purpose: Script for fitting two-pool models to exponential decay data
#           + loops

# --------------------
# Load package
library(minpack.lm)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(knitr)
# --------------------

setwd("C:/Users/danab/Box/WhitmanLab/Projects/WoodBuffalo/FireSim2019/data")


### Running model on our own data ##

### Two-pool model ###
# Mt = M1*e^(-k1*t)+M2*e^(-k2*t)  The basic two-pool model
# Mt = 1 = M1+M2    Setting Mt as total initial C to 1, M1+M2 make up total
# Mt = M1*e^(-k1*t)+(1-M1)*e^(-k2*t)    Rearranging equation
# Mt = M1*e^(-k1*t)-M1*e^(-k2*t)+e^(-k2*t)  Rearranging equation

# Describe function form
TwoPoolDecay = function(params,t) params$M1*(exp(-params$k1*t))-params$M1*exp(-params$k2*t)+exp(-params$k2*t)

## Define residual function (our actual value minus the exponential function)
residualsFun.twopool = function(p, Mt, t) Mt - TwoPoolDecay(p, t)

# (k1=0.005,k2=0.00001,M1=0.02) was pretty good except for a few dry burns (3, 8, 12, 13, 17) - too linear
# Trying different start parameters for those ones
# Looks like it's decaying faster than we'd think, so starting with a larger M1 might help.
# StartPar.TwoPool = list(k1=0.06,k2=0.00005,M1=0.01)

# OBJECTIVE: set up input and output files for loop
df <- read.csv("fractional-c-remaining.csv", row.names = 1)
df <- subset(df, select = c(core.id, core.id.incub, site, burn.trtmt, DayFactor, fractional.C.remaining))
df <- subset(df, DayFactor < 35)
head(df)


# Create a vector of Core IDs
# Optionally filter to pull out cores of interest
names.df <- df %>%
  subset(DayFactor == 3) %>%
  #filter(core.id %in% c("19UW-WB-03-02","19UW-WB-08-10","19UW-WB-12-02","19UW-WB-13-03","19UW-WB-17-01"))%>%
  subset(select = c(core.id.incub))

# Convert df to vector
names.vec <- names.df[,"core.id.incub"]
# Check if this worked. 
class(names.vec)

# --------------------
# Create an output df to populate with coefficients later
coefs.df <- data.frame("core.id.incub" = "sample", "k1" = 0, "k2" = 0, "M1" = 0, "M2" = 0)

output.fit.df <- data.frame("core.id" = "sample", 
                            "core.id.incub" = "sample", 
                            'r.squared' = 0,
                            'pvalue' = 0,
                            "site" = 0,
                            "burn.trtmt" = "trtmt",
                            "t" = 0, 
                            "Mt" = 0, "Mt.fit" = 0)



# List start parameters - can set these however
# StartPar.TwoPool = list(k1=0.005,k2=0.00001,M1=0.02)
# Dry: k1=0.005,k2=0.00001,M1=0.002
# Wet: k1=0.005,k2=0.00001,M1=0.02
# Control: k1=0.005,k2=0.00001,M1=0.02

# Run the decay.models.R script in a for loop. 
for (i in seq_along(names.vec)) {
  
  df.x <- df %>%
    rename(t = DayFactor, Mt = fractional.C.remaining) %>%
    subset(core.id.incub == names.vec[[i]])
  
  if (df.x$burn.trtmt[1] == "dry") {
    StartPar.TwoPool = list(k1=0.005,k2=0.00001,M1=0.002)
  } else if (df.x$burn.trtmt[1] == "wet") {
    StartPar.TwoPool = list(k1=0.005,k2=0.00001,M1=0.02)
  } else {
    StartPar.TwoPool = list(k1=0.005,k2=0.00001,M1=0.01)
  }
  
  TwoPool.nls.fit <- nls.lm(par=StartPar.TwoPool, fn = residualsFun.twopool, 
                            Mt = df.x$Mt, 
                            t = df.x$t, 
                            lower = c(0,0,0), upper = c(Inf, Inf, 1),
                            control = nls.lm.control(nprint=1,maxiter=500,ftol=sqrt(.Machine$double.eps)/10))
  
  df.x$Mt.fit <- TwoPoolDecay(as.list(coef(TwoPool.nls.fit)),df.x$t)
  
  df.lm <- lm(Mt ~ Mt.fit, data = df.x)
  df.x$r.squared <- summary(df.lm)$r.squared
  df.x$pvalue <- summary(df.lm)$coef[2,4]
  

  ## summary information on parameter estimates
  coefs <- summary(TwoPool.nls.fit) 
  k1 <- coefs$coefficients[[1]]
  k2 <- coefs$coefficients[[2]]
  M1 <- coefs$coefficients[[3]]
  M2 <- df.x$Mt[1] - M1
  pvalue <- df.x$pvalue
  
  coefs.df <- coefs.df %>%
    add_row("core.id.incub" = names.vec[[i]], "k1" = k1, "k2" = k2, "M1" = M1, "M2" = M2)
  coefs.df
  
  output.fit.df <- rbind(output.fit.df, df.x)
}



# --------------------
### Second iteration:

# Set cutoff for R squared
rsquared = 0.98

# If R-squared is above cut-off, it's good. Otherwise, redo parameters. 
output.fit.df <- merge(output.fit.df, coefs.df)
df.rsquared <- subset(output.fit.df, r.squared >rsquared)
df.redo <- subset(output.fit.df, r.squared < rsquared)

names.vec <- df.redo %>%
  subset(t == 3) %>%
  subset(select = c(core.id.incub))


df <- read.csv("fractional-c-remaining.csv", row.names = 1)
df <- merge(names.vec, df)

df <- df %>%
  separate(core.id.incub, c('year','project','Site','core','incub'), sep = '-', remove = FALSE) %>%
  subset(select = c(core.id.incub,
                    core.id,
                    site,
                    burn.trtmt,
                    DayFactor,
                    fractional.C.remaining))

names.vec = as.vector(names.vec$core.id.incub)

# Create an output df to populate with coefficients later
coefs.df <- data.frame("core.id.incub" = "sample", "k1" = 0, "k2" = 0, "M1" = 0, "M2" = 0)

output.fit.df <- data.frame("core.id" = "sample", 
                            "core.id.incub" = "sample", 
                            'r.squared' = 0,
                            'pvalue' = 0,
                            "site" = 0,
                            "burn.trtmt" = "trtmt",
                            "t" = 0, 
                            "Mt" = 0, "Mt.fit" = 0)


# List start parameters - can set these however
StartPar.TwoPool = list(k1=0.005,k2=0.0001,M1=0.01)
head(df)


# Run the decay.models.R script in a for loop. 
for (i in seq_along(names.vec)) {
  
  df.x <- df %>%
    rename(t = DayFactor, Mt = fractional.C.remaining) %>%
    subset(core.id.incub == names.vec[[i]])
  
  TwoPool.nls.fit <- nls.lm(par=StartPar.TwoPool, fn = residualsFun.twopool, 
                            Mt = df.x$Mt, 
                            t = df.x$t, 
                            lower = c(0,0,0), upper = c(Inf, Inf, 1),
                            control = nls.lm.control(nprint=1,maxiter=500,ftol=sqrt(.Machine$double.eps)/10))
  
  df.x$Mt.fit <- TwoPoolDecay(as.list(coef(TwoPool.nls.fit)),df.x$t)
  
  df.lm <- lm(Mt ~ Mt.fit, data = df.x)
  df.x$r.squared <- summary(df.lm)$r.squared
  df.x$pvalue <- summary(df.lm)$coef[2,4]
  
  
  ## summary information on parameter estimates
  coefs <- summary(TwoPool.nls.fit) 
  k1 <- coefs$coefficients[[1]]
  k2 <- coefs$coefficients[[2]]
  M1 <- coefs$coefficients[[3]]
  M2 <- df.x$Mt[1] - M1
  
  coefs.df <- coefs.df %>%
    add_row("core.id.incub" = names.vec[[i]], "k1" = k1, "k2" = k2, "M1" = M1, "M2" = M2)
  coefs.df
  
  output.fit.df <- rbind(output.fit.df, df.x)
}




# --------------------
### Third iteration

# If R-squared is above cut-off, it's good. Otherwise, redo parameters.
output.fit.df <- merge(output.fit.df, coefs.df)
df.rsquared <- df.rsquared %>%
  rbind(output.fit.df) %>%
  subset(r.squared >rsquared)
df.redo <- subset(output.fit.df, r.squared < rsquared)

names.vec <- df.redo %>%
  subset(t == 3) %>%
  subset(select = c(core.id.incub))

df <- read.csv("fractional-c-remaining.csv", row.names = 1)
df <- merge(names.vec, df)

df <- df %>%
  separate(core.id.incub, c('year','project','Site','core','incub'), sep = '-', remove = FALSE) %>%
  subset(select = c(core.id.incub,
                    core.id,
                    site,
                    burn.trtmt,
                    DayFactor,
                    fractional.C.remaining))

names.vec = as.vector(names.vec$core.id.incub)

# Create an output df to populate with coefficients later
coefs.df <- data.frame("core.id.incub" = "sample", "k1" = 0, "k2" = 0, "M1" = 0, "M2" = 0)

output.fit.df <- data.frame("core.id" = "sample", 
                            "core.id.incub" = "sample", 
                            'r.squared' = 0,
                            'pvalue' = 0,
                            "site" = 0,
                            "burn.trtmt" = "trtmt",
                            "t" = 0, 
                            "Mt" = 0, "Mt.fit" = 0)


# List start parameters - can set these however
StartPar.TwoPool = list(k1=0.05,k2=0.0001,M1=0.001)

# Run the decay.models.R script in a for loop. 
for (i in seq_along(names.vec)) {
  
  df.x <- df %>%
    rename(t = DayFactor, Mt = fractional.C.remaining) %>%
    subset(core.id.incub == names.vec[[i]])
  
  TwoPool.nls.fit <- nls.lm(par=StartPar.TwoPool, fn = residualsFun.twopool, 
                            Mt = df.x$Mt, 
                            t = df.x$t, 
                            lower = c(0,0,0), upper = c(Inf, Inf, 1),
                            control = nls.lm.control(nprint=1,maxiter=500,ftol=sqrt(.Machine$double.eps)/10))
  
  df.x$Mt.fit <- TwoPoolDecay(as.list(coef(TwoPool.nls.fit)),df.x$t)
  
  df.lm <- lm(Mt ~ Mt.fit, data = df.x)
  df.x$r.squared <- summary(df.lm)$r.squared
  df.x$pvalue <- summary(df.lm)$coef[2,4]
  
  ## summary information on parameter estimates
  coefs <- summary(TwoPool.nls.fit) 
  k1 <- coefs$coefficients[[1]]
  k2 <- coefs$coefficients[[2]]
  M1 <- coefs$coefficients[[3]]
  M2 <- df.x$Mt[1] - M1
  pvalue <- df.x$pvalue
  
  coefs.df <- coefs.df %>%
    add_row("core.id.incub" = names.vec[[i]], "k1" = k1, "k2" = k2, "M1" = M1, "M2" = M2)
  coefs.df
  
  output.fit.df <- rbind(output.fit.df, df.x)
}


# If R-squared is above cut-off, it's good. Otherwise, redo parameters. 
output.fit.df <- merge(output.fit.df, coefs.df)

df.rsquared <- df.rsquared %>%
  rbind(output.fit.df) %>%
  subset(r.squared >rsquared)
df.redo <- subset(output.fit.df, r.squared < rsquared)

names.vec <- df.redo %>%
  subset(t == 3) %>%
  subset(select = c(core.id.incub))

# Complete at R^2 = 0.98



head(df.rsquared)
max(df.rsquared$k1)


# --------------------
# Clean up output files
df.rsquared <- df.rsquared %>%
  separate(core.id.incub, c("Year", "Project", "Site", "Core", "incub.trtmt"), sep = "-", remove = FALSE) 

#write.csv(df.rsquared,'../data/2-pool-decay-model-output.csv')
