# --------------------
# Title: FireSim2019 Experiment 2 corncob analysis, temperature cutoff
# Author: Dana Johnson
# Date: 2021-Oct-13
#
# set-up -----

library(corncob)
library(phyloseq)
library(magrittr)
library(dplyr)

#setwd("C:/Users/danab/Box/WhitmanLab/Projects/WoodBuffalo/FireSim2019/")

### CORNCOB TAKES RAW SEQ DATA AS INPUT ###

### Step 1. Import relevant phyloseq data ----
# Raw seqs
ps.raw.full <- readRDS('data/sequence-data/LibCombined/phyloseq-objects/ps.raw.full')
ps.raw.full

ps.norm.full <- readRDS('data/sequence-data/LibCombined/phyloseq-objects/ps.norm.full')
ps.norm.full

# Remove duplicate samples:
ps.raw.full <- prune_samples(sample_data(ps.raw.full)$Full.id != '19UW-WB-06-08-A-SI-duplicate' &
                               sample_data(ps.raw.full)$Full.id != '19UW-WB-07-02-O-SI-duplicate' &
                               sample_data(ps.raw.full)$Full.id != '19UW-WB-11-03-O-SI-duplicate'&
                               sample_data(ps.raw.full)$Full.id != '19UW-WB-19-07-A-SI-duplicate' &
                               sample_data(ps.raw.full)$Full.id != '19UW-WB-08-10-O-SI', ps.raw.full)


ps.norm.full <- prune_samples(sample_data(ps.norm.full)$Full.id != '19UW-WB-06-08-A-SI-duplicate' &
                                sample_data(ps.norm.full)$Full.id != '19UW-WB-07-02-O-SI-duplicate' &
                                sample_data(ps.norm.full)$Full.id != '19UW-WB-11-03-O-SI-duplicate'&
                                sample_data(ps.norm.full)$Full.id != '19UW-WB-19-07-A-SI-duplicate' &
                                sample_data(ps.norm.full)$Full.id != '19UW-WB-08-10-O-SI', ps.norm.full)

# Remove taxa with 0 abundance:
ps.raw.prune <-  prune_samples(sample_data(ps.raw.full)$incub.trtmt %in% c('pb','SI') &
                                 sample_data(ps.raw.full)$DNA.type == 'gDNA', ps.raw.full) 
ps.raw.prune <- prune_taxa(taxa_sums(ps.raw.prune)>0, ps.raw.prune) 


ps.norm.prune <-  prune_samples(sample_data(ps.norm.full)$incub.trtmt %in% c('pb','SI') &
                                  sample_data(ps.norm.full)$DNA.type == 'gDNA', ps.norm.full) 
ps.norm.prune <- prune_taxa(taxa_sums(ps.norm.prune)>0, ps.norm.prune)


# Create phylum dataset:
#ps.raw.phylum <- tax_glom(ps.raw.prune, taxrank = 'Phylum')
#ps.raw.prune <- ps.raw.phylum


### Step 2: Create phyloseq objects with Temp cutoff for three comparison -----

## PB DRY vs. 5-WEEK DRY: ##
# Subset phyloseq object
ps.DvD <- ps.raw.prune %>%
  phyloseq::subset_samples(burn.trtmt %in% c('dry')) 

# Create burn temp cutoff:
df.dry <- data.frame(sample_data(ps.DvD)) %>%
  group_by(site) %>%
  mutate(Dry.burn.max = max(Thermo.mid.max)) 

# Recombine sample data with phyloseq object
ps.df <- sample_data(df.dry)

# Fix phyloseq names
sample_names(ps.df) = df.dry$Full.id
sample_data(ps.DvD) <- ps.df



## PB WET vs. 5-WEEK WET ##
# Subset phyloseq object
ps.WvW <- ps.raw.prune %>%
  phyloseq::subset_samples(burn.trtmt %in% c('wet')) 

# Create burn temp cutoff:
df.wet <- data.frame(sample_data(ps.WvW)) %>%
  group_by(site) %>%
  mutate(Wet.burn.max = max(Thermo.mid.max)) 

# Recombine sample data with phyloseq object
ps.df <- sample_data(df.wet)

# Fix phyloseq names
sample_names(ps.df) = df.wet$Full.id
sample_data(ps.WvW) <- ps.df



### Step 3:  Set up and run corncob analyzes -----
##    Set up corncob O horizon, dry burn, Max temp > 50 ##
ps.corncob.DvD.O.gt50 <-  prune_samples(sample_data(ps.DvD)$horizon == 'O' & 
                                          sample_data(ps.DvD)$Dry.burn.max>50, ps.DvD)

ps.corncob.DvD.O.gt50 <- prune_taxa(taxa_sums(ps.corncob.DvD.O.gt50)>0, ps.corncob.DvD.O.gt50)
ps.corncob.DvD.O.gt50

# Run corncob
dt.DvD.O.gt50 <- differentialTest(formula = ~incub.trtmt+pH,
                                         phi.formula = ~ incub.trtmt+pH,
                                         formula_null = ~ pH,
                                         phi.formula_null = ~incub.trtmt+ pH,
                                         test = 'Wald',boot=FALSE,
                                         data=ps.corncob.DvD.O.gt50,
                                         fdr_cutoff = 0.05)



##    Set up corncob A horizon, dry burn, Max temp > 50
ps.corncob.DvD.A.gt50 <-  prune_samples(sample_data(ps.DvD)$horizon == 'A' & 
                                          sample_data(ps.DvD)$Dry.burn.max>50, ps.DvD)


ps.corncob.DvD.A.gt50 <- prune_taxa(taxa_sums(ps.corncob.DvD.A.gt50)>0, ps.corncob.DvD.A.gt50)
ps.corncob.DvD.A.gt50

sample_names(ps.corncob.DvD.A.gt50)
# Results: 

# Run corncob
dt.DvD.A.gt50 <- differentialTest(formula = ~incub.trtmt+pH,
                                  phi.formula = ~ incub.trtmt+pH,
                                  formula_null = ~ pH,
                                  phi.formula_null = ~incub.trtmt+ pH,
                                  test = 'Wald',boot=FALSE,
                                  data=ps.corncob.DvD.A.gt50,
                                  fdr_cutoff = 0.05)



# ##    Set up corncob O horizon, dry burn, Max temp < 50
# ps.corncob.DvD.O.lt50 <-  prune_samples(sample_data(ps.DvD)$horizon == 'O' &
#                                           sample_data(ps.DvD)$Dry.burn.max<50, ps.DvD)
# 
# 
# ps.corncob.DvD.O.lt50 <- prune_taxa(taxa_sums(ps.corncob.DvD.O.lt50)>0, ps.corncob.DvD.O.lt50)
# ps.corncob.DvD.O.lt50
# 
# sample_names(ps.corncob.DvD.O.lt50)
# # Results: 3
# 
# # Run corncob
# dt.DvD.O.lt50 <- differentialTest(formula = ~incub.trtmt+pH,
#                                   phi.formula = ~ incub.trtmt+pH,
#                                   formula_null = ~ pH,
#                                   phi.formula_null = ~incub.trtmt+ pH,
#                                   test = 'Wald',boot=FALSE,
#                                   data=ps.corncob.DvD.O.lt50,
#                                   fdr_cutoff = 0.05)
# 
# 
# 
# ##    Set up corncob A horizon, dry burn, Max temp < 50
# ps.corncob.DvD.A.lt50 <-  prune_samples(sample_data(ps.DvD)$horizon == 'A' &
#                                           sample_data(ps.DvD)$Dry.burn.max<50, ps.DvD)
# 
# 
# ps.corncob.DvD.A.lt50 <- prune_taxa(taxa_sums(ps.corncob.DvD.A.lt50)>0, ps.corncob.DvD.A.lt50)
# ps.corncob.DvD.A.lt50
# 
# sample_names(ps.corncob.DvD.A.lt50)
# # Results:
# 
# # Run corncob
# dt.DvD.A.lt50 <- differentialTest(formula = ~incub.trtmt+pH,
#                                   phi.formula = ~ incub.trtmt+pH,
#                                   formula_null = ~ pH,
#                                   phi.formula_null = ~incub.trtmt+ pH,
#                                   test = 'Wald',boot=FALSE,
#                                   data=ps.corncob.DvD.A.lt50,
#                                   fdr_cutoff = 0.05)
# 



##    Set up corncob O horizon, wet burn, Max temp > 50
ps.corncob.WvW.O.gt50 <-  prune_samples(sample_data(ps.WvW)$horizon == 'O' & 
                                          sample_data(ps.WvW)$Wet.burn.max>50, ps.WvW)

ps.corncob.WvW.O.gt50 <- prune_taxa(taxa_sums(ps.corncob.WvW.O.gt50)>0, ps.corncob.WvW.O.gt50)
ps.corncob.WvW.O.gt50

# Run corncob
dt.WvW.O.gt50 <- differentialTest(formula = ~incub.trtmt+pH,
                                  phi.formula = ~ incub.trtmt+pH,
                                  formula_null = ~ pH,
                                  phi.formula_null = ~incub.trtmt+ pH,
                                  test = 'Wald',boot=FALSE,
                                  data=ps.corncob.WvW.O.gt50,
                                  fdr_cutoff = 0.05)



##    Set up corncob A horizon, wet burn, Max temp > 50
ps.corncob.WvW.A.gt50 <-  prune_samples(sample_data(ps.WvW)$horizon == 'A' & 
                                          sample_data(ps.WvW)$Wet.burn.max>50, ps.WvW)

ps.corncob.WvW.A.gt50 <- prune_taxa(taxa_sums(ps.corncob.WvW.A.gt50)>0, ps.corncob.WvW.A.gt50)
ps.corncob.WvW.A.gt50

sample_names(ps.corncob.WvW.A.gt50)
# Results: 

# Run corncob
dt.WvW.A.gt50 <- differentialTest(formula = ~incub.trtmt+pH,
                                  phi.formula = ~ incub.trtmt+pH,
                                  formula_null = ~ pH,
                                  phi.formula_null = ~incub.trtmt+ pH,
                                  test = 'Wald',boot=FALSE,
                                  data=ps.corncob.WvW.A.gt50,
                                  fdr_cutoff = 0.05)




# ##    Set up corncob O horizon, wet burn, Max temp < 50
# ps.corncob.WvW.O.lt50 <-  prune_samples(sample_data(ps.WvW)$horizon == 'O' &
#                                           sample_data(ps.WvW)$Wet.burn.max<50, ps.WvW)
# 
# ps.corncob.WvW.O.lt50 <- prune_taxa(taxa_sums(ps.corncob.WvW.O.lt50)>0, ps.corncob.WvW.O.lt50)
# ps.corncob.WvW.O.lt50
# 
# # Run corncob
# dt.WvW.O.lt50 <- differentialTest(formula = ~incub.trtmt+pH,
#                                   phi.formula = ~ incub.trtmt+pH,
#                                   formula_null = ~ pH,
#                                   phi.formula_null = ~incub.trtmt+ pH,
#                                   test = 'Wald',boot=FALSE,
#                                   data=ps.corncob.WvW.O.lt50,
#                                   fdr_cutoff = 0.05)
# 
# 
#
# ##    Set up corncob A horizon, wet burn, Max temp < 50
# ps.corncob.WvW.A.lt50 <-  prune_samples(sample_data(ps.WvW)$horizon == 'A' &
#                                           sample_data(ps.WvW)$Wet.burn.max<50, ps.WvW)
# 
# 
# ps.corncob.WvW.A.lt50 <- prune_taxa(taxa_sums(ps.corncob.WvW.A.lt50)>0, ps.corncob.WvW.A.lt50)
# ps.corncob.WvW.A.lt50
# 
# # Run corncob
# dt.WvW.A.lt50 <- differentialTest(formula = ~incub.trtmt+pH,
#                                   phi.formula = ~ incub.trtmt+pH,
#                                   formula_null = ~ pH,
#                                   phi.formula_null = ~incub.trtmt+ pH,
#                                   test = 'Wald',boot=FALSE,
#                                   data=ps.corncob.WvW.A.lt50,
#                                   fdr_cutoff = 0.05)





### Step 4. Clean up corncob output -----
# Create empty dataframes to fill
m = data.frame('intercept'=0, 'mu.incub.trtmt'=0, 'mu.pH'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
               'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- 0
dT <- dt.DvD.O.gt50
for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'pb.dry.vs.SI.dry'
  m[i,7] = 'fast growth'
  m[i,8] = 'Mid thermo > 50C'
  m[i,9] = 'null'
  m[i,10] = 'O'
  m[i,11] = 'gDNA'
  m[i,12] = 'pH'
}

m.concat = m



# REPEAT FOR A HORIZON:
m = data.frame('intercept'=0, 'mu.incub.trtmt'=0, 'mu.pH'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
               'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- 0
dT <- dt.DvD.A.gt50
for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'pb.dry.vs.SI.dry'
  m[i,7] = 'fast growth'
  m[i,8] = 'Mid thermo > 50C'
  m[i,9] = 'null'
  m[i,10] = 'A'
  m[i,11] = 'gDNA'
  m[i,12] = 'pH'
}

m.concat = rbind(m.concat, m)



#REPEAT FOR O HORIZON, temps < 50C:
m = data.frame('intercept'=0, 'mu.incub.trtmt'=0, 'mu.pH'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
               'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- 0
dT <- dt.DvD.O.lt50

for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'pb.dry.vs.SI.dry'
  m[i,7] = 'fast growth'
  m[i,8] = 'Mid thermo < 50C'
  m[i,9] = 'null'
  m[i,10] = 'O'
  m[i,11] = 'gDNA'
  m[i,12] = 'pH'
}

m.concat = rbind(m.concat, m)



# REPEAT FOR A HORIZON, temps < 50C:
m = data.frame('intercept'=0, 'mu.incub.trtmt'=0, 'mu.pH'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
               'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- 0
dT <- dt.DvD.A.lt50
for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'pb.dry.vs.SI.dry'
  m[i,7] = 'fast growth'
  m[i,8] = 'Mid thermo < 50C'
  m[i,9] = 'null'
  m[i,10] = 'A'
  m[i,11] = 'gDNA'
  m[i,12] = 'pH'
}

m.concat = rbind(m.concat, m)



# Create empty dataframes to fill
m = data.frame('intercept'=0, 'mu.incub.trtmt'=0, 'mu.pH'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
               'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- 0
dT <- dt.WvW.O.gt50
for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'pb.wet.vs.SI.wet'
  m[i,7] = 'fast growth'
  m[i,8] = 'Mid thermo > 50C'
  m[i,9] = 'null'
  m[i,10] = 'O'
  m[i,11] = 'gDNA'
  m[i,12] = 'pH'
}

m.concat = rbind(m.concat, m)



# REPEAT FOR A HORIZON:
# m = data.frame('intercept'=0, 'mu.incub.trtmt'=0, 'mu.pH'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
               'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- 0
dT <- dt.WvW.A.gt50
for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'pb.wet.vs.SI.wet'
  m[i,7] = 'fast growth'
  m[i,8] = 'Mid thermo > 50C'
  m[i,9] = 'null'
  m[i,10] = 'A'
  m[i,11] = 'gDNA'
  m[i,12] = 'pH'
}

m.concat = rbind(m.concat, m)



# REPEAT FOR O HORIZON, temps < 50C:
m = data.frame('intercept'=0, 'mu.incub.trtmt'=0, 'mu.pH'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
               'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- 0
dT <- dt.WvW.O.lt50
for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'pb.wet.vs.SI.wet'
  m[i,7] = 'fast growth'
  m[i,8] = 'Mid thermo < 50C'
  m[i,9] = 'null'
  m[i,10] = 'O'
  m[i,11] = 'gDNA'
  m[i,12] = 'pH'
}

m.concat = rbind(m.concat, m)

# 
# 
# # REPEAT FOR A HORIZON, temps < 50C:
m = data.frame('intercept'=0, 'mu.incub.trtmt'=0, 'mu.pH'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
               'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- 0
dT <- dt.WvW.A.lt50
for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'pb.wet.vs.SI.wet'
  m[i,7] = 'fast growth'
  m[i,8] = 'Mid thermo < 50C'
  m[i,9] = 'null'
  m[i,10] = 'A'
  m[i,11] = 'gDNA'
  m[i,12] = 'pH'
}

m.concat = rbind(m.concat, m)




### Step 5. Clean up dataframe -----

# Let's bring back in our taxonomy from the tax table
SigOTUs = levels(as.factor(m.concat$OTU))
pruned = prune_taxa(SigOTUs,ps.norm.full)
taxtab = data.frame(tax_table(pruned))
taxtab$OTU = c(taxa_names(pruned))
df.joined = merge(m.concat,taxtab,by=c("OTU"))
head(df.joined)
dim(df.joined)



### Step 6. Save output: -----
#write.csv(df.joined,"data/sequence-data/LibCombined/corncob-output/Manuscript/Ex2_signOTUs_MidThermoCutoff_50C.csv", row.names = FALSE)
#   


