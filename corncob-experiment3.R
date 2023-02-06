# --------------------
# Title: FireSim2019 Experiment 3 corncob analysis
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

# Exclude all but 6 month incubation samples and remove taxa with 0 abundance:
ps.raw.prune <-  prune_samples(sample_data(ps.raw.full)$incub.trtmt %in% c('LIwA') &
                                 sample_data(ps.raw.full)$DNA.type == 'gDNA', ps.raw.full) 
ps.raw.prune <- prune_taxa(taxa_sums(ps.raw.prune)>0, ps.raw.prune) 


ps.norm.prune <-  prune_samples(sample_data(ps.norm.full)$incub.trtmt %in% c('LIwA') &
                                  sample_data(ps.norm.full)$DNA.type == 'gDNA', ps.norm.full) 
ps.norm.prune <- prune_taxa(taxa_sums(ps.norm.prune)>0, ps.norm.prune)

# Create phylum dataset:
#ps.raw.phylum <- tax_glom(ps.raw.prune, taxrank = 'Phylum')
#ps.raw.prune <- ps.raw.phylum





### Step 2. Create temp categories -----

# Subset dry burns
ps.LIwA.dry <-  prune_samples(sample_data(ps.raw.prune)$burn.trtmt != 'wet', ps.raw.prune)

# Want to pull out pairs (dry burn and control) of cores in which the dry burn
#   core reached temps greater than some T cutoff. To do this, we need to assign
#   a "max temp" of the core pair to both the dry burn and the control core. 
df.dry <- data.frame(sample_data(ps.LIwA.dry)) %>%
  group_by(site) %>%
  mutate(Dry.burn.max = max(Thermo.mid.max)) 

# Now combine the new sample data with the original phyloseq object:
ps.df <- sample_data(df.dry)

sample_names(ps.df) = df.dry$Full.id

sample_data(ps.LIwA.dry) <- ps.df



# # Subset wet burns
# ps.LIwA.wet <-  prune_samples(sample_data(ps.raw.prune)$burn.trtmt != 'dry', ps.raw.prune)
# 
# # Create dataframe from ps.LIwA.wet sample_data.
# df.wet <- data.frame(sample_data(ps.LIwA.wet)) %>%
#   group_by(site) %>%
#   mutate(Wet.burn.max = max(Thermo.mid.max)) 
# 
# # Now combine the new sample data with the original phyloseq object:
# ps.df <- sample_data(df.wet)
# 
# sample_names(ps.df) = df.wet$Full.id
# 
# sample_data(ps.LIwA.wet) <- ps.df




### Step 3. Subset phyloseq and run corncob ----- 

## Dry Vs. Control; O horizon; Mid thermo > 50 C
ps.corncob.dry.O.gt50 <-  prune_samples(sample_data(ps.LIwA.dry)$horizon == 'O' & 
                                           sample_data(ps.LIwA.dry)$Dry.burn.max>50, ps.LIwA.dry)
ps.corncob.dry.O.gt50 <- prune_taxa(taxa_sums(ps.corncob.dry.O.gt50)>0, ps.corncob.dry.O.gt50)
ps.corncob.dry.O.gt50

# Run corncob:
dt.LIwA.DryBurn.O.gt50 <- differentialTest(formula = ~burn.trtmt+cntl.pH,
                                           phi.formula = ~ burn.trtmt+cntl.pH,
                                           formula_null = ~ cntl.pH,
                                           phi.formula_null = ~burn.trtmt+ cntl.pH,
                                           test = 'Wald',boot=FALSE,
                                           data=ps.corncob.dry.O.gt50,
                                           fdr_cutoff = 0.05)



## Dry Vs. Control; A horizon; Mid thermo > 50 C
ps.corncob.dry.A.gt50 <-  prune_samples(sample_data(ps.LIwA.dry)$horizon == 'A' & 
                                           sample_data(ps.LIwA.dry)$Dry.burn.max>50, ps.LIwA.dry)

ps.corncob.dry.A.gt50 <- prune_taxa(taxa_sums(ps.corncob.dry.A.gt50)>0, ps.corncob.dry.A.gt50)
ps.corncob.dry.A.gt50

sample_names(ps.corncob.dry.A.gt50)
# Result = 

# Run corncob
dt.LIwA.DryBurn.A.gt50 <- differentialTest(formula = ~burn.trtmt+cntl.pH,
                                           phi.formula = ~ burn.trtmt+cntl.pH,
                                           formula_null = ~ cntl.pH,
                                           phi.formula_null = ~burn.trtmt+ cntl.pH,
                                           test = 'Wald',boot=FALSE,
                                           data=ps.corncob.dry.A.gt50,
                                           fdr_cutoff = 0.05)



# ##Dry Vs. Control; O horizon; Mid thermo < 50 C
# ps.corncob.dry.O.lt50 <-  prune_samples(sample_data(ps.LIwA.dry)$horizon == 'O' & 
#                                            sample_data(ps.LIwA.dry)$Dry.burn.max<50, ps.LIwA.dry)
# 
# ps.corncob.dry.O.lt50 <- prune_taxa(taxa_sums(ps.corncob.dry.O.lt50)>0, ps.corncob.dry.O.lt50)
# ps.corncob.dry.O.lt50
# 
# # Run corncob
# dt.LIwA.DryBurn.O.lt50 <- differentialTest(formula = ~burn.trtmt+cntl.pH,
#                                            phi.formula = ~ burn.trtmt+cntl.pH,
#                                            formula_null = ~ cntl.pH,
#                                            phi.formula_null = ~burn.trtmt+ cntl.pH,
#                                            test = 'Wald',boot=FALSE,
#                                            data=ps.corncob.dry.O.lt50,
#                                            fdr_cutoff = 0.05)
# 
#                                            
#                                            
# ## Dry Vs. Control; A horizon; Mid thermo < 50 C
# ps.corncob.dry.A.lt50 <-  prune_samples(sample_data(ps.LIwA.dry)$horizon == 'A' & 
#                                            sample_data(ps.LIwA.dry)$Dry.burn.max<50, ps.LIwA.dry)
# 
# ps.corncob.dry.A.lt50 <- prune_taxa(taxa_sums(ps.corncob.dry.A.lt50) >0, ps.corncob.dry.A.lt50)
# ps.corncob.dry.A.lt50
# 
# # If less than two samples make it through filters, don't run corncob.
# if (nsamples(ps.corncob.dry.A.lt50) <2) {
#   ps.corncob.dry.A.lt50 = 0
# }
# 
# # Run corncob
# dt.LIwA.DryBurn.A.lt50 <- differentialTest(formula = ~burn.trtmt+cntl.pH,
#                                            phi.formula = ~ burn.trtmt+cntl.pH,
#                                            formula_null = ~ cntl.pH,
#                                            phi.formula_null = ~burn.trtmt+ cntl.pH,
#                                            test = 'Wald',boot=FALSE,
#                                            data=ps.corncob.dry.A.lt50,
#                                            fdr_cutoff = 0.05)
# 
# 
# 
# ## Wet Vs. Control; O horizon; Mid thermo > 50 C
# ps.corncob.wet.O.gt50 <-  prune_samples(sample_data(ps.LIwA.wet)$horizon == 'O' & 
#                                            sample_data(ps.LIwA.wet)$Wet.burn.max>50, ps.LIwA.wet)
# 
# ps.corncob.wet.O.gt50 <- prune_taxa(taxa_sums(ps.corncob.wet.O.gt50)>0, ps.corncob.wet.O.gt50)
# ps.corncob.wet.O.gt50
# 
# # Run corncob:
# dt.LIwA.wetBurn.O.gt50 <- differentialTest(formula = ~burn.trtmt+cntl.pH,
#                                            phi.formula = ~ burn.trtmt+cntl.pH,
#                                            formula_null = ~ cntl.pH,
#                                            phi.formula_null = ~burn.trtmt+ cntl.pH,
#                                            test = 'Wald',boot=FALSE,
#                                            data=ps.corncob.wet.O.gt50,
#                                            fdr_cutoff = 0.05)
# 
# 
# 
# ## Wet Vs. Control; A horizon; Mid thermo > 50 C
# ps.corncob.wet.A.gt50 <-  prune_samples(sample_data(ps.LIwA.wet)$horizon == 'A' & 
#                                            sample_data(ps.LIwA.wet)$Wet.burn.max>50, ps.LIwA.wet)
# 
# ps.corncob.wet.A.gt50 <- prune_taxa(taxa_sums(ps.corncob.wet.A.gt50)>0, ps.corncob.wet.A.gt50)
# ps.corncob.wet.A.gt50
# 
# # Run corncob
# dt.LIwA.wetBurn.A.gt50 <- differentialTest(formula = ~burn.trtmt+cntl.pH,
#                                            phi.formula = ~ burn.trtmt+cntl.pH,
#                                            formula_null = ~ cntl.pH,
#                                            phi.formula_null = ~burn.trtmt+ cntl.pH,
#                                            test = 'Wald',boot=FALSE,
#                                            data=ps.corncob.wet.A.gt50,
#                                            fdr_cutoff = 0.05)
# 
# 
# 
# ## Wet Vs. Control; O horizon; Mid thermo < 50 C
# ps.corncob.wet.O.lt50 <-  prune_samples(sample_data(ps.LIwA.wet)$horizon == 'O' & 
#                                            sample_data(ps.LIwA.wet)$Wet.burn.max<50, ps.LIwA.wet)
# 
# ps.corncob.wet.O.lt50 <- prune_taxa(taxa_sums(ps.corncob.wet.O.lt50)>0, ps.corncob.wet.O.lt50)
# ps.corncob.wet.O.lt50
# 
# # Run corncob
# dt.LIwA.wetBurn.O.lt50 <- differentialTest(formula = ~burn.trtmt+cntl.pH,
#                                            phi.formula = ~ burn.trtmt+cntl.pH,
#                                            formula_null = ~ cntl.pH,
#                                            phi.formula_null = ~burn.trtmt+ cntl.pH,
#                                            test = 'Wald',boot=FALSE,
#                                            data=ps.corncob.wet.O.lt50,
#                                            fdr_cutoff = 0.05)
# 
# 
# 
# ##  Wet Vs. Control; A horizon; Mid thermo < 50 C
# ps.corncob.wet.A.lt50 <-  prune_samples(sample_data(ps.LIwA.wet)$horizon == 'A' & 
#                                            sample_data(ps.LIwA.wet)$Wet.burn.max<50, ps.LIwA.wet)
# 
# ps.corncob.wet.A.lt50 <- prune_taxa(taxa_sums(ps.corncob.wet.A.lt50)>0, ps.corncob.wet.A.lt50)
# ps.corncob.wet.A.lt50
# 
# # Run corncob
# dt.LIwA.wetBurn.A.lt50 <- differentialTest(formula = ~burn.trtmt+cntl.pH,
#                                            phi.formula = ~ burn.trtmt+cntl.pH,
#                                            formula_null = ~ cntl.pH,
#                                            phi.formula_null = ~burn.trtmt+ cntl.pH,
#                                            test = 'Wald',boot=FALSE,
#                                            data=ps.corncob.wet.A.lt50,
#                                            fdr_cutoff = 0.05)





### Step 4. Clean up corncob output -----

# PULL OUT COEFS FOR DRY BURNS, O HORIZON:

# Create empty dataframe to fill
m = data.frame('intercept'=0, 'mu.burn.trtmt'=0, 'mu.cntl.pH'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
               'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- dt.LIwA.DryBurn.O.gt50
for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'LIwA.dry.v.LIwA.control'
  m[i,7] = 'affinity'
  m[i,8] = 'Mid thermo > 50C'
  m[i,9] = 'null'
  m[i,10] = 'O'
  m[i,11] = 'gDNA'
  m[i,12] = 'pH'
}
m.concat = m



# Repeat with A horizon
m = data.frame('intercept'=0, 'mu.burn.trtmt'=0, 'mu.cntl.pH'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
               'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- dt.LIwA.DryBurn.A.gt50
for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'LIwA.dry.v.LIwA.control'
  m[i,7] = 'affinity'
  m[i,8] = 'Mid thermo > 50C'
  m[i,9] = 'null'
  m[i,10] = 'A'
  m[i,11] = 'gDNA'
  m[i,12] = 'pH'
}

m.concat = rbind(m.concat, m)



# # Temp < 50C
# m = data.frame('intercept'=0, 'mu.burn.trtmt'=0, 'mu.cntl.pH'=0, 'OTU'='taxa','p_fdr'=0,
#                'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
#                'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')
# 
# dT <- dt.LIwA.DryBurn.O.lt50
# for (i in 1:length(dT$significant_taxa)) {
#   m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
#   m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
#   # Also grab the p_fdr estimate for that taxon's model
#   m[i,5] = dT$p_fdr[dT$significant_taxa][i]
#   m[i,6] = 'LIwA.dry.v.LIwA.control'
#   m[i,7] = 'affinity'
#   m[i,8] = 'Mid thermo < 50C'
#   m[i,9] = 'null'
#   m[i,10] = 'O'
#   m[i,11] = 'gDNA'
#   m[i,12] = 'pH'
# }
# 
# m.concat = rbind(m.concat, m)
# 
# 
# 
# # Repeat with A horizon
# m = data.frame('intercept'=0, 'mu.burn.trtmt'=0, 'mu.cntl.pH'=0, 'OTU'='taxa','p_fdr'=0,
#                'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
#                'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')
# 
# dT <- dt.LIwA.DryBurn.A.lt50
# for (i in 1:length(dT$significant_taxa)) {
#   m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
#   m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
#   # Also grab the p_fdr estimate for that taxon's model
#   m[i,5] = dT$p_fdr[dT$significant_taxa][i]
#   m[i,6] = 'LIwA.dry.v.LIwA.control'
#   m[i,7] = 'affinity'
#   m[i,8] = 'Mid thermo < 50C'
#   m[i,9] = 'null'
#   m[i,10] = 'A'
#   m[i,11] = 'gDNA'
#   m[i,12] = 'pH'
# }
# 
# m.concat = rbind(m.concat, m)
# 
# 
# 
# 
# # PULL OUT COEFS FOR WET BURNS, O HORIZON:
# 
# # Create empty dataframe to fill
# m = data.frame('intercept'=0, 'mu.burn.trtmt'=0, 'mu.cntl.pH'=0, 'OTU'='taxa','p_fdr'=0,
#                'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
#                'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')
# 
# dT <- dt.LIwA.WetBurn.O.gt50
# for (i in 1:length(dT$significant_taxa)) {
#   m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
#   m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
#   # Also grab the p_fdr estimate for that taxon's model
#   m[i,5] = dT$p_fdr[dT$significant_taxa][i]
#   m[i,6] = 'LIwA.wet.v.LIwA.control'
#   m[i,7] = 'affinity'
#   m[i,8] = 'Mid thermo > 50C'
#   m[i,9] = 'null'
#   m[i,10] = 'O'
#   m[i,11] = 'gDNA'
#   m[i,12] = 'pH'
# }
# m.concat = rbind(m.concat, m)
# 
# 
# 
# Repeat with A horizon
# m = data.frame('intercept'=0, 'mu.burn.trtmt'=0, 'mu.cntl.pH'=0, 'OTU'='taxa','p_fdr'=0,
#                'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
#                'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')
# 
# dT <- dt.LIwA.WetBurn.A.gt50
# for (i in 1:length(dT$significant_taxa)) {
#   m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
#   m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
#   # Also grab the p_fdr estimate for that taxon's model
#   m[i,5] = dT$p_fdr[dT$significant_taxa][i]
#   m[i,6] = 'LIwA.wet.v.LIwA.control'
#   m[i,7] = 'affinity'
#   m[i,8] = 'Mid thermo > 50C'
#   m[i,9] = 'null'
#   m[i,10] = 'A'
#   m[i,11] = 'gDNA'
#   m[i,12] = 'pH'
# }
# 
# m.concat = rbind(m.concat, m)
# 
# 
# 
# # Temp < 50C
# m = data.frame('intercept'=0, 'mu.burn.trtmt'=0, 'mu.cntl.pH'=0, 'OTU'='taxa','p_fdr'=0,
#                'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
#                'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')
# 
# dT <- dt.LIwA.WetBurn.O.lt50
# for (i in 1:length(dT$significant_taxa)) {
#   m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
#   m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
#   # Also grab the p_fdr estimate for that taxon's model
#   m[i,5] = dT$p_fdr[dT$significant_taxa][i]
#   m[i,6] = 'LIwA.wet.v.LIwA.control'
#   m[i,7] = 'affinity'
#   m[i,8] = 'Mid thermo < 50C'
#   m[i,9] = 'null'
#   m[i,10] = 'O'
#   m[i,11] = 'gDNA'
#   m[i,12] = 'pH'
# }
# 
# m.concat = rbind(m.concat, m)
# 
# 
# 
# # Repeat with A horizon
# m = data.frame('intercept'=0, 'mu.burn.trtmt'=0, 'mu.cntl.pH'=0, 'OTU'='taxa','p_fdr'=0,
#                'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
#                'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')
# 
# dT <- dt.LIwA.WetBurn.A.lt50
# for (i in 1:length(dT$significant_taxa)) {
#   m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
#   m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
#   # Also grab the p_fdr estimate for that taxon's model
#   m[i,5] = dT$p_fdr[dT$significant_taxa][i]
#   m[i,6] = 'LIwA.wet.v.LIwA.control'
#   m[i,7] = 'affinity'
#   m[i,8] = 'Mid thermo < 50C'
#   m[i,9] = 'null'
#   m[i,10] = 'A'
#   m[i,11] = 'gDNA'
#   m[i,12] = 'pH'
# }
# 
# m.concat = rbind(m.concat, m)



### Step 5. Clean up output -----

# Let's bring back in our taxonomy from the tax table
SigOTUs = levels(as.factor(m.concat$OTU))
pruned = prune_taxa(SigOTUs,ps.norm.prune)
taxtab = data.frame(tax_table(pruned))
taxtab$OTU = c(taxa_names(pruned))
df.joined = merge(m.concat,taxtab,by=c("OTU"))
head(df.joined)

### Step 6. Save output -----
# May want to save at this point
write.csv(df.joined,"data/sequence-data/LibCombined/corncob-output/Manuscript/Ex3_signOTUs_MidThermoCutoff_50C.csv", row.names = FALSE)

