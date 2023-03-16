# --------------------
# Title: FireSim2019 Experiment 1 corncob analysis
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

# NOT RELEVANT TO 24 HR POST-BURN SAMPLES:
# ps.raw.full <- prune_samples(sample_data(ps.raw.full)$Full.id != '19UW-WB-06-08-A-SI-duplicate' &
#                                sample_data(ps.raw.full)$Full.id != '19UW-WB-07-02-O-SI-duplicate' &
#                                sample_data(ps.raw.full)$Full.id != '19UW-WB-11-03-O-SI-duplicate'&
#                                sample_data(ps.raw.full)$Full.id != '19UW-WB-19-07-A-SI-duplicate' &
#                                sample_data(ps.raw.full)$Full.id != '19UW-WB-08-10-O-SI', ps.raw.full)
# ps.raw.full <- prune_taxa(taxa_sums(ps.raw.full)>0, ps.raw.full)
# 
# # normalize abundance seq data
# ps.norm.full <- prune_samples(sample_data(ps.norm.full)$Full.id != '19UW-WB-06-08-A-SI-duplicate' &
#                                 sample_data(ps.norm.full)$Full.id != '19UW-WB-07-02-O-SI-duplicate' &
#                                 sample_data(ps.norm.full)$Full.id != '19UW-WB-11-03-O-SI-duplicate'&
#                                 sample_data(ps.norm.full)$Full.id != '19UW-WB-19-07-A-SI-duplicate' &
#                                 sample_data(ps.norm.full)$Full.id != '19UW-WB-08-10-O-SI', ps.norm.full)
# ps.norm.full <- prune_taxa(taxa_sums(ps.norm.full)>0, ps.norm.full)

# Exclude all but 24 hr post-burn
# Create separate genomic DNA and an RNA phyloseq objects
ps.raw.RNA <- ps.raw.full %>%
  phyloseq::subset_samples(incub.trtmt %in% c('pb')) %>%
  subset_samples(DNA.type %in% c('cDNA')) 
ps.raw.RNA <- prune_taxa(taxa_sums(ps.raw.RNA)>0, ps.raw.RNA)


ps.raw.DNA <- ps.raw.full %>%
  phyloseq::subset_samples(incub.trtmt %in% c('pb')) %>%
  subset_samples(DNA.type %in% c('gDNA')) 
ps.raw.DNA <- prune_taxa(taxa_sums(ps.raw.DNA)>0, ps.raw.DNA)


# 
# ### Step 2. Subset for the most abundance taxa for use in testing -----
# 
# # Create list of taxa with a mean abundance > cut-off
# MeanAbunTaxa = taxa_names(filter_taxa(ps.norm.full, function(x) mean(x) > 0.0001, TRUE))
# length(MeanAbunTaxa)
# 
# # Create list of taxa with max abundance > cutoff
# MaxAbunTaxa <-  taxa_names(filter_taxa(ps.norm.full, function(x) max(x) > 0.005, TRUE))
# length(MaxAbunTaxa)
# 
# ntaxa(ps.norm.full)
# 
# x <- c(MeanAbunTaxa, MaxAbunTaxa)
# length(unique(x))
# 
# # Prune taxa to include only OTUs with a maximum rel. abun > Cutoff
# ps.prune.raw.DNA = prune_taxa(MeanAbunTaxa, ps.raw.DNA)
# ps.prune.raw.DNA


# Create phylum dataset:
#ps.raw.phylum <- tax_glom(ps.raw.DNA, taxrank = 'Phylum')
#ps.prune.raw.DNA <- ps.raw.phylum




### Step 3. Set temperature cutoff for burns -----
# Focus on dry soil burns: 
ps.DNA.dry <-  prune_samples(sample_data(ps.raw.DNA)$burn.trtmt != 'wet', ps.raw.DNA)

# Create dataframe ps.DNA.dry sample_data.
df.dry <- data.frame(sample_data(ps.DNA.dry))

# Want to pull out pairs (dry burn and control) of cores in which the dry burn
#   core reached temps greater than some T cutoff. To do this, we need to assign
#   a "max temp" of the core pair to both the dry burn and the control core. 
df.dry <- df.dry %>%
  group_by(site) %>%
  mutate(Dry.burn.max = max(Thermo.mid.max)) 

# Now combine the new sample data with the original phyloseq object:
ps.df <- sample_data(df.dry)

sample_names(ps.df) = df.dry$Full.id

sample_data(ps.DNA.dry) <- ps.df



# # Focus on wet soil burns: 
# ps.DNA.wet <-  prune_samples(sample_data(ps.raw.DNA)$burn.trtmt != 'dry', ps.raw.DNA)
# 
# # Create dataframe ps.DNA.dry sample_data.
# df.wet <- data.frame(sample_data(ps.DNA.wet))
# 
# # Assign maximum temp at bottom thermocouple as the Max temp. for the core:
# df.wet <- df.wet %>%
#   group_by(site) %>%
#   mutate(Wet.burn.max = max(Thermo.low.max)) 
# 
# # Now combine the new sample data with the original phyloseq object:
# ps.df <- sample_data(df.wet)
# 
# sample_names(ps.df) = df.wet$Full.id
# 
# sample_data(ps.DNA.wet) <- ps.df




### Step 4. Run corncob on O and A horizons -----

# Create O horizon phyloseq object
ps.corncob.dry.O.gt50 <-  prune_samples(sample_data(ps.DNA.dry)$horizon == 'O' & 
                                          sample_data(ps.DNA.dry)$Dry.burn.max>50, ps.DNA.dry)
# Clean up ps to reduce size
ps.corncob.dry.O.gt50 <- prune_taxa(taxa_sums(ps.corncob.dry.O.gt50)>0, ps.corncob.dry.O.gt50)
ps.corncob.dry.O.gt50
# This is an odd number of samples because not all the samples are paired dry
#   burn+controls due to <100% success in extracting and sequencing DNA post-burn. 


# How many samples make it through?
sample_names(ps.corncob.dry.O.gt50) # 31 samples 

# Run corncob controlling for pre-burn pH (cntl.pH)
dt.DNA.DryBurn.O.gt50 <- differentialTest(formula = ~burn.trtmt+cntl.pH,
                                          phi.formula = ~ burn.trtmt+cntl.pH,
                                          formula_null = ~ cntl.pH,
                                          phi.formula_null = ~burn.trtmt+ cntl.pH,
                                          test = 'Wald',boot=FALSE,
                                          data=ps.corncob.dry.O.gt50,
                                          fdr_cutoff = 0.05)



# Repeat for A horizon:
ps.corncob.dry.A.gt50 <-  prune_samples(sample_data(ps.DNA.dry)$horizon == 'A' & 
                                          sample_data(ps.DNA.dry)$Dry.burn.max>50, ps.DNA.dry)
# Clean up ps to reduce size
ps.corncob.dry.A.gt50 <- prune_taxa(taxa_sums(ps.corncob.dry.A.gt50)>0, ps.corncob.dry.A.gt50)
ps.corncob.dry.A.gt50 # 23 samples

# Run corncob controlling for pre-burn pH (cntl.pH)
dt.DNA.DryBurn.A.gt50 <- differentialTest(formula = ~burn.trtmt+cntl.pH,
                                          phi.formula = ~ burn.trtmt+cntl.pH,
                                          formula_null = ~ cntl.pH,
                                          phi.formula_null = ~burn.trtmt+ cntl.pH,
                                          test = 'Wald',boot=FALSE,
                                          data=ps.corncob.dry.A.gt50,
                                          fdr_cutoff = 0.05)



# ## Dry Vs. Control; O horizon; Low thermo < 50 C 
# ps.corncob.dry.O.lt50 <-  prune_samples(sample_data(ps.DNA.dry)$horizon == 'O' &  sample_data(ps.DNA.dry)$Dry.burn.max<50, ps.DNA.dry)
# ps.corncob.dry.O.lt50 <- prune_taxa(taxa_sums(ps.corncob.dry.O.lt50)>0, ps.corncob.dry.O.lt50)
# ps.corncob.dry.O.lt50
# 
# # Run corncob
# dt.DNA.DryBurn.O.lt50 <- differentialTest(formula = ~burn.trtmt+cntl.pH,
#                                           phi.formula = ~ burn.trtmt+cntl.pH,
#                                           formula_null = ~ cntl.pH,
#                                           phi.formula_null = ~burn.trtmt+ cntl.pH,
#                                           test = 'Wald',boot=FALSE,
#                                           data=ps.corncob.dry.O.lt50,
#                                           fdr_cutoff = 0.05)
# 
# 
# 
# ## Dry Vs. Control; A horizon; Low thermo < 50 C 
# ps.corncob.dry.A.lt50 <-  prune_samples(sample_data(ps.DNA.dry)$horizon == 'A' & sample_data(ps.DNA.dry)$Dry.burn.max<50, ps.DNA.dry)
# ps.corncob.dry.A.lt50 <- prune_taxa(taxa_sums(ps.corncob.dry.A.lt50)>0, ps.corncob.dry.A.lt50)
# ps.corncob.dry.A.lt50
# 
# 
# # Run corncob
# dt.DNA.DryBurn.A.lt50 <- differentialTest(formula = ~burn.trtmt+cntl.pH,
#                                           phi.formula = ~ burn.trtmt+cntl.pH,
#                                           formula_null = ~ cntl.pH,
#                                           phi.formula_null = ~burn.trtmt+ cntl.pH,
#                                           test = 'Wald',boot=FALSE,
#                                           data=ps.corncob.dry.A.lt50,
#                                           fdr_cutoff = 0.05)



# ## Wet Vs. Control; O horizon; Low thermo > 50 C - CONDITIONS NOT MET
# ps.corncob.wet.O.gt50 <-  prune_samples(sample_data(ps.DNA.wet)$horizon == 'O' &  sample_data(ps.DNA.wet)$Wet.burn.max>50, ps.DNA.wet)
# ps.corncob.wet.O.gt50 <- prune_taxa(taxa_sums(ps.corncob.wet.O.gt50)>0, ps.corncob.wet.O.gt50)
# ps.corncob.wet.O.gt50
# 
# # Run corncob:
# dt.DNA.wetBurn.O.gt50 <- differentialTest(formula = ~burn.trtmt+cntl.pH,
#                                           phi.formula = ~ burn.trtmt+cntl.pH,
#                                           formula_null = ~ cntl.pH,
#                                           phi.formula_null = ~burn.trtmt+ cntl.pH,
#                                           test = 'Wald',boot=FALSE,
#                                           data=ps.corncob.wet.O.gt50,
#                                           fdr_cutoff = 0.05)
# 
# 
# 
# ## Wet Vs. Control; A horizon; Low thermo > 50 C - CONDITIONS NOT MET
# ps.corncob.wet.A.gt50 <-  prune_samples(sample_data(ps.DNA.wet)$horizon == 'A' & sample_data(ps.DNA.wet)$Wet.burn.max>50, ps.DNA.wet)
# ps.corncob.wet.A.gt50 <- prune_taxa(taxa_sums(ps.corncob.wet.A.gt50)>0, ps.corncob.wet.A.gt50)
# ps.corncob.wet.A.gt50
# 
# # Run corncob
# dt.DNA.wetBurn.A.gt50 <- differentialTest(formula = ~burn.trtmt+cntl.pH,
#                                           phi.formula = ~ burn.trtmt+cntl.pH,
#                                           formula_null = ~ cntl.pH,
#                                           phi.formula_null = ~burn.trtmt+ cntl.pH,
#                                           test = 'Wald',boot=FALSE,
#                                           data=ps.corncob.wet.A.gt50,
#                                           fdr_cutoff = 0.05)


# 
# ## Wet Vs. Control; O horizon; Low thermo < 50 C
# ps.corncob.wet.O.lt50 <-  prune_samples(sample_data(ps.DNA.wet)$horizon == 'O' & sample_data(ps.DNA.wet)$Wet.burn.max<50, ps.DNA.wet)
# ps.corncob.wet.O.lt50 <- prune_taxa(taxa_sums(ps.corncob.wet.O.lt50)>0, ps.corncob.wet.O.lt50)
# ps.corncob.wet.O.lt50
# 
# dt.DNA.wetBurn.O.lt50 <- differentialTest(formula = ~burn.trtmt+cntl.pH,
#                                           phi.formula = ~ burn.trtmt+cntl.pH,
#                                           formula_null = ~ cntl.pH,
#                                           phi.formula_null = ~burn.trtmt+ cntl.pH,
#                                           test = 'Wald',boot=FALSE,
#                                           data=ps.corncob.wet.O.lt50,
#                                           fdr_cutoff = 0.05)
# 
# 
# 
# ## Wet Vs. Control; A horizon; Low thermo < 50 C
# ps.corncob.wet.A.lt50 <-  prune_samples(sample_data(ps.DNA.wet)$horizon == 'A' &sample_data(ps.DNA.wet)$Wet.burn.max<50, ps.DNA.wet)
# ps.corncob.wet.A.lt50 <- prune_taxa(taxa_sums(ps.corncob.wet.A.lt50)>0, ps.corncob.wet.A.lt50)
# ps.corncob.wet.A.lt50
# 
# dt.DNA.wetBurn.A.lt50 <- differentialTest(formula = ~burn.trtmt+cntl.pH,
#                                           phi.formula = ~ burn.trtmt+cntl.pH,
#                                           formula_null = ~ cntl.pH,
#                                           phi.formula_null = ~burn.trtmt+ cntl.pH,
#                                           test = 'Wald',boot=FALSE,
#                                           data=ps.corncob.wet.A.lt50,
#                                           fdr_cutoff = 0.05)



### Step 5. Clean up corncob output -----
# Create empty dataframes to fill
m = data.frame('intercept'=0, 'mu.burn.trtmt'=0, 'mu.cntl.pH'=0, 'OTU'='taxa','p_fdr'=0,
      'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
      'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- dt.DNA.DryBurn.O.gt50
for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'pb.dry.vs.pb.control'
  m[i,7] = 'survival'
  m[i,8] = 'Mid thermo > 50C'
  m[i,9] = 'null'
  m[i,10] = 'O'
  m[i,11] = 'gDNA'
  m[i,12] = 'pre-burn pH'
}

m.concat = m



# REPEAT FOR A HORIZON:
m = data.frame('intercept'=0, 'mu.burn.trtmt'=0, 'mu.cntl.pH'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
               'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- dt.DNA.DryBurn.A.gt50
for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'pb.dry.vs.pb.control'
  m[i,7] = 'survival'
  m[i,8] = 'Mid thermo > 50C'
  m[i,9] = 'null'
  m[i,10] = 'A'
  m[i,11] = 'gDNA'
  m[i,12] = 'pre-burn pH'
}

m.concat = rbind(m.concat, m)



# # Create empty dataframes to fill
# m = data.frame('intercept'=0, 'mu.burn.trtmt'=0, 'mu.cntl.pH'=0, 'OTU'='taxa','p_fdr'=0,
#                'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
#                'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')
# 
# dT <- dt.DNA.DryBurn.O.lt50
# for (i in 1:length(dT$significant_taxa)) {
#   m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
#   m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
#   # Also grab the p_fdr estimate for that taxon's model
#   m[i,5] = dT$p_fdr[dT$significant_taxa][i]
#   m[i,6] = 'pb.dry.vs.pb.control'
#   m[i,7] = 'survival'
#   m[i,8] = 'Low thermo < 50C'
#   m[i,9] = 'null'
#   m[i,10] = 'O'
#   m[i,11] = 'gDNA'
#   m[i,12] = 'pre-burn pH'
# }
# 
# m.concat = rbind(m.concat, m)
# 
# 
# 
# # REPEAT FOR A HORIZON:
# m = data.frame('intercept'=0, 'mu.burn.trtmt'=0, 'mu.cntl.pH'=0, 'OTU'='taxa','p_fdr'=0,
#                'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
#                'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')
# 
# dT <- dt.DNA.DryBurn.A.lt50
# for (i in 1:length(dT$significant_taxa)) {
#   m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
#   m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
#   # Also grab the p_fdr estimate for that taxon's model
#   m[i,5] = dT$p_fdr[dT$significant_taxa][i]
#   m[i,6] = 'pb.dry.vs.pb.control'
#   m[i,7] = 'survival'
#   m[i,8] = 'Low thermo < 50C'
#   m[i,9] = 'null'
#   m[i,10] = 'A'
#   m[i,11] = 'gDNA'
#   m[i,12] = 'pre-burn pH'
# }
# 
# m.concat = rbind(m.concat, m)
# 
# 
# 
# # Create empty dataframes to fill
# m = data.frame('intercept'=0, 'mu.burn.trtmt'=0, 'mu.cntl.pH'=0, 'OTU'='taxa','p_fdr'=0,
#                'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
#                'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')
# 
# dT <- dt.DNA.wetBurn.O.gt50
# for (i in 1:length(dT$significant_taxa)) {
#   m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
#   m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
#   # Also grab the p_fdr estimate for that taxon's model
#   m[i,5] = dT$p_fdr[dT$significant_taxa][i]
#   m[i,6] = 'pb.wet.vs.pb.control'
#   m[i,7] = 'survival'
#   m[i,8] = 'Low thermo > 50C'
#   m[i,9] = 'null'
#   m[i,10] = 'O'
#   m[i,11] = 'gDNA'
#   m[i,12] = 'pre-burn pH'
# }
# 
# m.concat = rbind(m.concat, m)
# 
# 
# 
# # REPEAT FOR A HORIZON:
# m = data.frame('intercept'=0, 'mu.burn.trtmt'=0, 'mu.cntl.pH'=0, 'OTU'='taxa','p_fdr'=0,
#                'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
#                'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')
# 
# dT <- dt.DNA.wetBurn.A.gt50
# for (i in 1:length(dT$significant_taxa)) {
#   m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
#   m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
#   # Also grab the p_fdr estimate for that taxon's model
#   m[i,5] = dT$p_fdr[dT$significant_taxa][i]
#   m[i,6] = 'pb.wet.vs.pb.control'
#   m[i,7] = 'survival'
#   m[i,8] = 'Low thermo > 50C'
#   m[i,9] = 'null'
#   m[i,10] = 'A'
#   m[i,11] = 'gDNA'
#   m[i,12] = 'pre-burn pH'
# }
# 
# m.concat = rbind(m.concat, m)
# 
# 
# 
# # Create empty dataframes to fill
# m = data.frame('intercept'=0, 'mu.burn.trtmt'=0, 'mu.cntl.pH'=0, 'OTU'='taxa','p_fdr'=0,
#                'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
#                'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')
# 
# dT <- dt.DNA.wetBurn.O.lt50
# for (i in 1:length(dT$significant_taxa)) {
#   m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
#   m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
#   # Also grab the p_fdr estimate for that taxon's model
#   m[i,5] = dT$p_fdr[dT$significant_taxa][i]
#   m[i,6] = 'pb.wet.vs.pb.control'
#   m[i,7] = 'survival'
#   m[i,8] = 'Low thermo < 50C'
#   m[i,9] = 'null'
#   m[i,10] = 'O'
#   m[i,11] = 'gDNA'
#   m[i,12] = 'pre-burn pH'
# }
# 
# m.concat = rbind(m.concat, m)
# 
# 
# 
# # REPEAT FOR A HORIZON:
# m = data.frame('intercept'=0, 'mu.burn.trtmt'=0, 'mu.cntl.pH'=0, 'OTU'='taxa','p_fdr'=0,
#                'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
#                'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')
# 
# dT <- dt.DNA.wetBurn.A.lt50
# for (i in 1:length(dT$significant_taxa)) {
#   m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
#   m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
#   # Also grab the p_fdr estimate for that taxon's model
#   m[i,5] = dT$p_fdr[dT$significant_taxa][i]
#   m[i,6] = 'pb.wet.vs.pb.control'
#   m[i,7] = 'survival'
#   m[i,8] = 'Low thermo < 50C'
#   m[i,9] = 'null'
#   m[i,10] = 'A'
#   m[i,11] = 'gDNA'
#   m[i,12] = 'pre-burn pH'
# }
# 
# m.concat = rbind(m.concat, m)




# Let's bring back in our taxonomy from the tax table
SigOTUs = levels(as.factor(m.concat$OTU))
pruned = prune_taxa(SigOTUs,ps.norm.full)
taxtab = data.frame(tax_table(pruned))
taxtab$OTU = c(taxa_names(pruned))
df.joined = merge(m.concat,taxtab,by=c("OTU"))
head(df.joined)






### Step 7. Save results -----

write.csv(df.joined,"data/sequence-data/LibCombined/corncob-output/Manuscript/Ex1_signOTUs_gDNA.csv", row.names = FALSE)

#otu_to_taxonomy(OTU = dt.RNA.wetBurn.A.lt50$significant_taxa, data = ps.RNA.wet)
#plot(dt.RNA.wetBurn.A.lt50, level = c('Phylum'))

