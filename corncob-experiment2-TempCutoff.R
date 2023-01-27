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
#ps.raw.phylum <- tax_glom(ps.raw.full, taxrank = 'Phylum')
#ps.raw.phylum


# EXPERIMENT 2 - pb Dry vs. 5-week dry -----

ps.DvD <- ps.raw.prune %>%
  phyloseq::subset_samples(burn.trtmt %in% c('dry')) 

# Set burn temp cutoff:
df.dry <- data.frame(sample_data(ps.DvD))

df.dry <- df.dry %>%
  group_by(site) %>%
  mutate(Dry.burn.max = max(Thermo.low.max)) 

ps.df <- sample_data(df.dry)

sample_names(ps.df) = df.dry$Full.id

sample_data(ps.DvD) <- ps.df


# 1.1. Set up corncob O horizon, dry burn, Max temp > 50 -----
ps.corncob.DvD.O.gt50 <-  prune_samples(sample_data(ps.DvD)$horizon == 'O' & 
                                          sample_data(ps.DvD)$Dry.burn.max>50, ps.DvD)

ps.corncob.DvD.O.gt50 <- prune_taxa(taxa_sums(ps.corncob.DvD.O.gt50)>0, ps.corncob.DvD.O.gt50)
ps.corncob.DvD.O.gt50

# pb vs SI, both horizon, control (no burn)
dt.DvD.O.gt50 <- differentialTest(formula = ~incub.trtmt+pH,
                                         phi.formula = ~ incub.trtmt+pH,
                                         formula_null = ~ pH,
                                         phi.formula_null = ~incub.trtmt+ pH,
                                         test = 'Wald',boot=FALSE,
                                         data=ps.corncob.DvD.O.gt50,
                                         fdr_cutoff = 0.05)

df.DvD.O.gt50 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.DvD.O.gt50$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.DvD.O.gt50$significant_models[[j]]
  # Pull out the coefficients as above
  # NOTE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_dry = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.DvD.O.gt50$p_fdr[dt.DvD.O.gt50$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_dry$p_fdr = p_fdr
  # Create a column with the OTU ID
  mu_dry$OTU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.DvD.O.gt50 = rbind(df.DvD.O.gt50,mu_dry)
}

df.DvD.O.gt50$Comparison = 'pb.dry.v.SI.dry'
df.DvD.O.gt50$Experiment = 'fast growth'
df.DvD.O.gt50$Temp.Cutoff = 'Low thermo > 50 C'
df.DvD.O.gt50$pH.cutoff = 'null'
df.DvD.O.gt50$horizon = 'O'
df.DvD.O.gt50$DNA.type = 'gDNA'
df.DvD.O.gt50$Controlling.for = 'pH'

colnames(df.DvD.O.gt50) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                              'Experiment', 'Temp.cutoff','pH.cutoff','horizon', 
                              'DNA.type', 'Controlling.for')






# 2.1. Set up corncob A horizon, dry burn, Max temp > 50 -----
ps.corncob.DvD.A.gt50 <-  prune_samples(sample_data(ps.DvD)$horizon == 'A' & 
                                          sample_data(ps.DvD)$Dry.burn.max>50, ps.DvD)


ps.corncob.DvD.A.gt50 <- prune_taxa(taxa_sums(ps.corncob.DvD.A.gt50)>0, ps.corncob.DvD.A.gt50)
ps.corncob.DvD.A.gt50

sample_names(ps.corncob.DvD.A.gt50)
# Results: 

# pb vs SI, both horizon, control (no burn)
dt.DvD.A.gt50 <- differentialTest(formula = ~incub.trtmt+pH,
                                  phi.formula = ~ incub.trtmt+pH,
                                  formula_null = ~ pH,
                                  phi.formula_null = ~incub.trtmt+ pH,
                                  test = 'Wald',boot=FALSE,
                                  data=ps.corncob.DvD.A.gt50,
                                  fdr_cutoff = 0.05)

df.DvD.A.gt50 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.DvD.A.gt50$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.DvD.A.gt50$significant_models[[j]]
  # Pull out the coefficients as above
  # NOTE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_dry = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.DvD.A.gt50$p_fdr[dt.DvD.A.gt50$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_dry$p_fdr = p_fdr
  # Create a column with the OTU ID
  mu_dry$OTU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.DvD.A.gt50 = rbind(df.DvD.A.gt50,mu_dry)
}

df.DvD.A.gt50$Comparison = 'pb.dry.v.SI.dry'
df.DvD.A.gt50$Experiment = 'fast growth'
df.DvD.A.gt50$Temp.Cutoff = 'Low thermo >50 C'
df.DvD.A.gt50$pH.cutoff = 'null'
df.DvD.A.gt50$horizon = 'A'
df.DvD.A.gt50$DNA.type = 'gDNA'
df.DvD.A.gt50$Controlling.for = 'pH'

colnames(df.DvD.A.gt50) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                            'Experiment', 'Temp.cutoff','pH.cutoff','horizon', 
                            'DNA.type', 'Controlling.for')







# 3.1. Set up corncob O horizon, dry burn, Max temp < 50 -----
ps.corncob.DvD.O.lt50 <-  prune_samples(sample_data(ps.DvD)$horizon == 'O' & 
                                          sample_data(ps.DvD)$Dry.burn.max<50, ps.DvD)


ps.corncob.DvD.O.lt50 <- prune_taxa(taxa_sums(ps.corncob.DvD.O.lt50)>0, ps.corncob.DvD.O.lt50)
ps.corncob.DvD.O.lt50

sample_names(ps.corncob.DvD.O.lt50)
# Results: 3

# pb vs SI, both horizon, control (no burn)
dt.DvD.O.lt50 <- differentialTest(formula = ~incub.trtmt+pH,
                                  phi.formula = ~ incub.trtmt+pH,
                                  formula_null = ~ pH,
                                  phi.formula_null = ~incub.trtmt+ pH,
                                  test = 'Wald',boot=FALSE,
                                  data=ps.corncob.DvD.O.lt50,
                                  fdr_cutoff = 0.05)

df.DvD.O.lt50 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.DvD.O.lt50$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.DvD.O.lt50$significant_models[[j]]
  # Pull out the coefficients as above
  # NOTE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_dry = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.DvD.O.lt50$p_fdr[dt.DvD.O.lt50$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_dry$p_fdr = p_fdr
  # Create a column with the OTU ID
  mu_dry$OTU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.DvD.O.lt50 = rbind(df.DvD.O.lt50,mu_dry)
}

df.DvD.O.lt50$Comparison = 'pb.dry.v.SI.dry'
df.DvD.O.lt50$Experiment = 'fast growth'
df.DvD.O.lt50$Temp.Cutoff = 'Low thermo less than 50 C'
df.DvD.O.lt50$pH.cutoff = 'null'
df.DvD.O.lt50$horizon = 'O'
df.DvD.O.lt50$DNA.type = 'gDNA'
df.DvD.O.lt50$Controlling.for = 'pH'

colnames(df.DvD.O.lt50) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                            'Experiment', 'Temp.cutoff','pH.cutoff','horizon', 
                            'DNA.type', 'Controlling.for')







# 4.1. Set up corncob A horizon, dry burn, Max temp < 50 -----
ps.corncob.DvD.A.lt50 <-  prune_samples(sample_data(ps.DvD)$horizon == 'A' & 
                                          sample_data(ps.DvD)$Dry.burn.max<50, ps.DvD)


ps.corncob.DvD.A.lt50 <- prune_taxa(taxa_sums(ps.corncob.DvD.A.lt50)>0, ps.corncob.DvD.A.lt50)
ps.corncob.DvD.A.lt50

sample_names(ps.corncob.DvD.A.lt50)
# Results: 

# pb vs SI, both horizon, control (no burn)
dt.DvD.A.lt50 <- differentialTest(formula = ~incub.trtmt+pH,
                                  phi.formula = ~ incub.trtmt+pH,
                                  formula_null = ~ pH,
                                  phi.formula_null = ~incub.trtmt+ pH,
                                  test = 'Wald',boot=FALSE,
                                  data=ps.corncob.DvD.A.lt50,
                                  fdr_cutoff = 0.05)

df.DvD.A.lt50 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.DvD.A.lt50$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.DvD.A.lt50$significant_models[[j]]
  # Pull out the coefficients as above
  # NOTE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_dry = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.DvD.A.lt50$p_fdr[dt.DvD.A.lt50$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_dry$p_fdr = p_fdr
  # Create a column with the OTU ID
  mu_dry$OTU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.DvD.A.lt50 = rbind(df.DvD.A.lt50,mu_dry)
}

df.DvD.A.lt50$Comparison = 'pb.dry.v.SI.dry'
df.DvD.A.lt50$Experiment = 'fast growth'
df.DvD.A.lt50$Temp.Cutoff = 'Low thermo less than 50 C'
df.DvD.A.lt50$pH.cutoff = 'null'
df.DvD.A.lt50$horizon = 'A'
df.DvD.A.lt50$DNA.type = 'gDNA'
df.DvD.A.lt50$Controlling.for = 'pH'

colnames(df.DvD.A.lt50) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                            'Experiment', 'Temp.cutoff','pH.cutoff','horizon', 
                            'DNA.type', 'Controlling.for')








# EXPERIMENT 2 - pb wet vs. 5-week wet -----

ps.WvW <- ps.raw.prune %>%
  phyloseq::subset_samples(burn.trtmt %in% c('wet')) 

# Set burn temp cutoff:
df.wet <- data.frame(sample_data(ps.WvW))

df.wet <- df.wet %>%
  group_by(site) %>%
  mutate(wet.burn.max = max(Thermo.low.max)) 

ps.df <- sample_data(df.wet)

sample_names(ps.df) = df.wet$Full.id

sample_data(ps.WvW) <- ps.df


# 1.2. Set up corncob O horizon, wet burn, Max temp > 50 -----
ps.corncob.WvW.O.gt50 <-  prune_samples(sample_data(ps.WvW)$horizon == 'O' & 
                                          sample_data(ps.WvW)$wet.burn.max>50, ps.WvW)

ps.corncob.WvW.O.gt50 <- prune_taxa(taxa_sums(ps.corncob.WvW.O.gt50)>0, ps.corncob.WvW.O.gt50)
ps.corncob.WvW.O.gt50

# pb vs SI, both horizon, control (no burn)
dt.WvW.O.gt50 <- differentialTest(formula = ~incub.trtmt+pH,
                                  phi.formula = ~ incub.trtmt+pH,
                                  formula_null = ~ pH,
                                  phi.formula_null = ~incub.trtmt+ pH,
                                  test = 'Wald',boot=FALSE,
                                  data=ps.corncob.WvW.O.gt50,
                                  fdr_cutoff = 0.05)

df.WvW.O.gt50 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.WvW.O.gt50$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.WvW.O.gt50$significant_models[[j]]
  # Pull out the coefficients as above
  # NOTE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_wet = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.WvW.O.gt50$p_fdr[dt.WvW.O.gt50$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_wet$p_fdr = p_fdr
  # Create a column with the OTU ID
  mu_wet$OTU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.WvW.O.gt50 = rbind(df.WvW.O.gt50,mu_wet)
}

df.WvW.O.gt50$Comparison = 'pb.wet.v.SI.wet'
df.WvW.O.gt50$Experiment = 'fast growth'
df.WvW.O.gt50$Temp.Cutoff = 'Low thermo >50 C'
df.WvW.O.gt50$pH.cutoff = 'null'
df.WvW.O.gt50$horizon = 'O'
df.WvW.O.gt50$DNA.type = 'gDNA'
df.WvW.O.gt50$Controlling.for = 'pH'

colnames(df.WvW.O.gt50) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                            'Experiment', 'Temp.cutoff','pH.cutoff','horizon', 
                            'DNA.type', 'Controlling.for')




# 2.2. Set up corncob A horizon, wet burn, Max temp > 50 -----
ps.corncob.WvW.A.gt50 <-  prune_samples(sample_data(ps.WvW)$horizon == 'A' & 
                                          sample_data(ps.WvW)$wet.burn.max>50, ps.WvW)

ps.corncob.WvW.A.gt50 <- prune_taxa(taxa_sums(ps.corncob.WvW.A.gt50)>0, ps.corncob.WvW.A.gt50)
ps.corncob.WvW.A.gt50

sample_names(ps.corncob.WvW.A.gt50)
# Results: 

# pb vs SI, both horizon, control (no burn)
dt.WvW.A.gt50 <- differentialTest(formula = ~incub.trtmt+pH,
                                  phi.formula = ~ incub.trtmt+pH,
                                  formula_null = ~ pH,
                                  phi.formula_null = ~incub.trtmt+ pH,
                                  test = 'Wald',boot=FALSE,
                                  data=ps.corncob.WvW.A.gt50,
                                  fdr_cutoff = 0.05)

df.WvW.A.gt50 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.WvW.A.gt50$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.WvW.A.gt50$significant_models[[j]]
  # Pull out the coefficients as above
  # NOTE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_wet = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.WvW.A.gt50$p_fdr[dt.WvW.A.gt50$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_wet$p_fdr = p_fdr
  # Create a column with the OTU ID
  mu_wet$OTU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.WvW.A.gt50 = rbind(df.WvW.A.gt50,mu_wet)
}

df.WvW.A.gt50$Comparison = 'pb.wet.v.SI.wet'
df.WvW.A.gt50$Experiment = 'fast growth'
df.WvW.A.gt50$Temp.Cutoff = 'Low thermo >50 C'
df.WvW.A.gt50$pH.cutoff = 'null'
df.WvW.A.gt50$horizon = 'A'
df.WvW.A.gt50$DNA.type = 'gDNA'
df.WvW.A.gt50$Controlling.for = 'pH'

colnames(df.WvW.A.gt50) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                            'Experiment', 'Temp.cutoff','pH.cutoff','horizon', 
                            'DNA.type', 'Controlling.for')







# 3.2. Set up corncob O horizon, wet burn, Max temp < 50 -----
ps.corncob.WvW.O.lt50 <-  prune_samples(sample_data(ps.WvW)$horizon == 'O' & 
                                          sample_data(ps.WvW)$wet.burn.max<50, ps.WvW)


ps.corncob.WvW.O.lt50 <- prune_taxa(taxa_sums(ps.corncob.WvW.O.lt50)>0, ps.corncob.WvW.O.lt50)
ps.corncob.WvW.O.lt50

sample_names(ps.corncob.WvW.O.lt50)
# Results: 


# pb vs SI, both horizon, control (no burn)
dt.WvW.O.lt50 <- differentialTest(formula = ~incub.trtmt+pH,
                                  phi.formula = ~ incub.trtmt+pH,
                                  formula_null = ~ pH,
                                  phi.formula_null = ~incub.trtmt+ pH,
                                  test = 'Wald',boot=FALSE,
                                  data=ps.corncob.WvW.O.lt50,
                                  fdr_cutoff = 0.05)

df.WvW.O.lt50 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.WvW.O.lt50$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.WvW.O.lt50$significant_models[[j]]
  # Pull out the coefficients as above
  # NOTE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_wet = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.WvW.O.lt50$p_fdr[dt.WvW.O.lt50$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_wet$p_fdr = p_fdr
  # Create a column with the OTU ID
  mu_wet$OTU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.WvW.O.lt50 = rbind(df.WvW.O.lt50,mu_wet)
}

df.WvW.O.lt50$Comparison = 'pb.wet.v.SI.wet'
df.WvW.O.lt50$Experiment = 'fast growth'
df.WvW.O.lt50$Temp.Cutoff = 'Low thermo less than 50 C'
df.WvW.O.lt50$pH.cutoff = 'null'
df.WvW.O.lt50$horizon = 'O'
df.WvW.O.lt50$DNA.type = 'gDNA'
df.WvW.O.lt50$Controlling.for = 'pH'

colnames(df.WvW.O.lt50) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                            'Experiment', 'Temp.cutoff','pH.cutoff','horizon', 
                            'DNA.type', 'Controlling.for')








# 4.2. Set up corncob A horizon, wet burn, Max temp < 50 -----
ps.corncob.WvW.A.lt50 <-  prune_samples(sample_data(ps.WvW)$horizon == 'A' & 
                                          sample_data(ps.WvW)$wet.burn.max<50, ps.WvW)


ps.corncob.WvW.A.lt50 <- prune_taxa(taxa_sums(ps.corncob.WvW.A.lt50)>0, ps.corncob.WvW.A.lt50)
ps.corncob.WvW.A.lt50

sample_names(ps.corncob.WvW.A.lt50)
# Results: 


# pb vs SI, both horizon, control (no burn)
dt.WvW.A.lt50 <- differentialTest(formula = ~incub.trtmt+pH,
                                  phi.formula = ~ incub.trtmt+pH,
                                  formula_null = ~ pH,
                                  phi.formula_null = ~incub.trtmt+ pH,
                                  test = 'Wald',boot=FALSE,
                                  data=ps.corncob.WvW.A.lt50,
                                  fdr_cutoff = 0.05)

df.WvW.A.lt50 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.WvW.A.lt50$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.WvW.A.lt50$significant_models[[j]]
  # Pull out the coefficients as above
  # NOTE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_wet = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.WvW.A.lt50$p_fdr[dt.WvW.A.lt50$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_wet$p_fdr = p_fdr
  # Create a column with the OTU ID
  mu_wet$OTU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.WvW.A.lt50 = rbind(df.WvW.A.lt50,mu_wet)
}

df.WvW.A.lt50$Comparison = 'pb.wet.v.SI.wet'
df.WvW.A.lt50$Experiment = 'fast growth'
df.WvW.A.lt50$Temp.Cutoff = 'Low thermo less than 50 C'
df.WvW.A.lt50$pH.cutoff = 'null'
df.WvW.A.lt50$horizon = 'A'
df.WvW.A.lt50$DNA.type = 'gDNA'
df.WvW.A.lt50$Controlling.for = 'pH'

colnames(df.WvW.A.lt50) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                            'Experiment', 'Temp.cutoff','pH.cutoff','horizon', 
                            'DNA.type', 'Controlling.for')







# EXPERIMENT 2 - pb control vs. 5-week control -----

ps.CvC <- ps.raw.prune %>%
  phyloseq::subset_samples(burn.trtmt %in% c('control')) 

# Set burn temp cutoff:
df.control <- data.frame(sample_data(ps.CvC))

df.control <- df.control %>%
  group_by(site) %>%
  mutate(control.burn.max = max(Thermo.low.max)) 

ps.df <- sample_data(df.control)

sample_names(ps.df) = df.control$Full.id

sample_data(ps.CvC) <- ps.df


# 1.3. Set up corncob O horizon, control burn, Max temp > 50 -----
ps.corncob.CvC.O.gt50 <-  prune_samples(sample_data(ps.CvC)$horizon == 'O' & 
                                          sample_data(ps.CvC)$control.burn.max>50, ps.CvC)

ps.corncob.CvC.O.gt50 <- prune_taxa(taxa_sums(ps.corncob.CvC.O.gt50)>0, ps.corncob.CvC.O.gt50)
ps.corncob.CvC.O.gt50

# pb vs SI, both horizon, control (no burn)
dt.CvC.O.gt50 <- differentialTest(formula = ~incub.trtmt+pH,
                                  phi.formula = ~ incub.trtmt+pH,
                                  formula_null = ~ pH,
                                  phi.formula_null = ~incub.trtmt+ pH,
                                  test = 'Wald',boot=FALSE,
                                  data=ps.corncob.CvC.O.gt50,
                                  fdr_cutoff = 0.05)

df.CvC.O.gt50 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.CvC.O.gt50$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.CvC.O.gt50$significant_models[[j]]
  # Pull out the coefficients as above
  # NOTE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_control = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.CvC.O.gt50$p_fdr[dt.CvC.O.gt50$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_control$p_fdr = p_fdr
  # Create a column with the OTU ID
  mu_control$OTU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.CvC.O.gt50 = rbind(df.CvC.O.gt50,mu_control)
}

df.CvC.O.gt50$Comparison = 'pb.control.v.SI.control'
df.CvC.O.gt50$Experiment = 'fast growth'
df.CvC.O.gt50$Temp.Cutoff = 'Low thermo >50 C'
df.CvC.O.gt50$pH.cutoff = 'null'
df.CvC.O.gt50$horizon = 'O'
df.CvC.O.gt50$DNA.type = 'gDNA'
df.CvC.O.gt50$Controlling.for = 'pH'

colnames(df.CvC.O.gt50) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                            'Experiment', 'Temp.cutoff','pH.cutoff','horizon', 
                            'DNA.type', 'Controlling.for')





# 2.3. Set up corncob A horizon, control burn, Max temp > 50 -----
ps.corncob.CvC.A.gt50 <-  prune_samples(sample_data(ps.CvC)$horizon == 'A' & 
                                          sample_data(ps.CvC)$control.burn.max>50, ps.CvC)


ps.corncob.CvC.A.gt50 <- prune_taxa(taxa_sums(ps.corncob.CvC.A.gt50)>0, ps.corncob.CvC.A.gt50)
ps.corncob.CvC.A.gt50

sample_names(ps.corncob.CvC.A.gt50)
# Results: 

# pb vs SI, both horizon, control (no burn)
dt.CvC.A.gt50 <- differentialTest(formula = ~incub.trtmt+pH,
                                  phi.formula = ~ incub.trtmt+pH,
                                  formula_null = ~ pH,
                                  phi.formula_null = ~incub.trtmt+ pH,
                                  test = 'Wald',boot=FALSE,
                                  data=ps.corncob.CvC.A.gt50,
                                  fdr_cutoff = 0.05)

df.CvC.A.gt50 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.CvC.A.gt50$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.CvC.A.gt50$significant_models[[j]]
  # Pull out the coefficients as above
  # NOTE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_control = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.CvC.A.gt50$p_fdr[dt.CvC.A.gt50$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_control$p_fdr = p_fdr
  # Create a column with the OTU ID
  mu_control$OTU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.CvC.A.gt50 = rbind(df.CvC.A.gt50,mu_control)
}

df.CvC.A.gt50$Comparison = 'pb.control.v.SI.control'
df.CvC.A.gt50$Experiment = 'fast growth'
df.CvC.A.gt50$Temp.Cutoff = 'Low thermo >50 C'
df.CvC.A.gt50$pH.cutoff = 'null'
df.CvC.A.gt50$horizon = 'A'
df.CvC.A.gt50$DNA.type = 'gDNA'
df.CvC.A.gt50$Controlling.for = 'pH'

colnames(df.CvC.A.gt50) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                            'Experiment', 'Temp.cutoff','pH.cutoff','horizon', 
                            'DNA.type', 'Controlling.for')







# 3.3. Set up corncob O horizon, control burn, Max temp < 50 -----
ps.corncob.CvC.O.lt50 <-  prune_samples(sample_data(ps.CvC)$horizon == 'O' & 
                                          sample_data(ps.CvC)$control.burn.max<50, ps.CvC)


ps.corncob.CvC.O.lt50 <- prune_taxa(taxa_sums(ps.corncob.CvC.O.lt50)>0, ps.corncob.CvC.O.lt50)
ps.corncob.CvC.O.lt50

sample_names(ps.corncob.CvC.O.lt50)
# Results: 


# pb vs SI, both horizon, control (no burn)
dt.CvC.O.lt50 <- differentialTest(formula = ~incub.trtmt+pH,
                                  phi.formula = ~ incub.trtmt+pH,
                                  formula_null = ~ pH,
                                  phi.formula_null = ~incub.trtmt+ pH,
                                  test = 'Wald',boot=FALSE,
                                  data=ps.corncob.CvC.O.lt50,
                                  fdr_cutoff = 0.05)

df.CvC.O.lt50 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.CvC.O.lt50$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.CvC.O.lt50$significant_models[[j]]
  # Pull out the coefficients as above
  # NOTE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_control = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.CvC.O.lt50$p_fdr[dt.CvC.O.lt50$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_control$p_fdr = p_fdr
  # Create a column with the OTU ID
  mu_control$OTU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.CvC.O.lt50 = rbind(df.CvC.O.lt50,mu_control)
}

df.CvC.O.lt50$Comparison = 'pb.control.v.SI.control'
df.CvC.O.lt50$Experiment = 'fast growth'
df.CvC.O.lt50$Temp.Cutoff = 'Low thermo less than 50 C'
df.CvC.O.lt50$pH.cutoff = 'null'
df.CvC.O.lt50$horizon = 'O'
df.CvC.O.lt50$DNA.type = 'gDNA'
df.CvC.O.lt50$Controlling.for = 'pH'

colnames(df.CvC.O.lt50) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                            'Experiment', 'Temp.cutoff','pH.cutoff','horizon', 
                            'DNA.type', 'Controlling.for')








# 4.3. Set up corncob A horizon, control burn, Max temp < 50 -----
ps.corncob.CvC.A.lt50 <-  prune_samples(sample_data(ps.CvC)$horizon == 'A' & 
                                          sample_data(ps.CvC)$control.burn.max<50, ps.CvC)


ps.corncob.CvC.A.lt50 <- prune_taxa(taxa_sums(ps.corncob.CvC.A.lt50)>0, ps.corncob.CvC.A.lt50)
ps.corncob.CvC.A.lt50

sample_names(ps.corncob.CvC.A.lt50)
# Results: 


# pb vs SI, both horizon, control (no burn)
dt.CvC.A.lt50 <- differentialTest(formula = ~incub.trtmt+pH,
                                  phi.formula = ~ incub.trtmt+pH,
                                  formula_null = ~ pH,
                                  phi.formula_null = ~incub.trtmt+ pH,
                                  test = 'Wald',boot=FALSE,
                                  data=ps.corncob.CvC.A.lt50,
                                  fdr_cutoff = 0.05)

df.CvC.A.lt50 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.CvC.A.lt50$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.CvC.A.lt50$significant_models[[j]]
  # Pull out the coefficients as above
  # NOTE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_control = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.CvC.A.lt50$p_fdr[dt.CvC.A.lt50$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_control$p_fdr = p_fdr
  # Create a column with the OTU ID
  mu_control$OTU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.CvC.A.lt50 = rbind(df.CvC.A.lt50,mu_control)
}

df.CvC.A.lt50$Comparison = 'pb.control.v.SI.control'
df.CvC.A.lt50$Experiment = 'fast growth'
df.CvC.A.lt50$Temp.Cutoff = 'Low thermo less than 50 C'
df.CvC.A.lt50$pH.cutoff = 'null'
df.CvC.A.lt50$horizon = 'A'
df.CvC.A.lt50$DNA.type = 'gDNA'
df.CvC.A.lt50$Controlling.for = 'pH'

colnames(df.CvC.A.lt50) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                            'Experiment', 'Temp.cutoff','pH.cutoff','horizon', 
                            'DNA.type', 'Controlling.for')









# CvC: Combined 1-4 -----
df.joined.CvC <- rbind(df.CvC.O.gt50, 
                       df.CvC.A.gt50,
                       df.CvC.O.lt50,
                       df.CvC.A.lt50)

# Let's bring back in our taxonomy from the tax table
SigOTUs = levels(as.factor(df.joined.CvC$OTU))
pruned = prune_taxa(SigOTUs,ps.norm.full)
taxtab = data.frame(tax_table(pruned))
taxtab$OTU = c(taxa_names(pruned))
df.joined.CvC = merge(df.joined.CvC,taxtab,by=c("OTU"))
head(df.joined.CvC)




# DvD: Combined 1-4 -----
df.joined.DvD <- rbind(df.DvD.O.gt50, 
                       df.DvD.A.gt50,
                       df.DvD.O.lt50,
                       df.DvD.A.lt50)

# Let's bring back in our taxonomy from the tax table
SigOTUs = levels(as.factor(df.joined.DvD$OTU))
pruned = prune_taxa(SigOTUs,ps.norm.full)
taxtab = data.frame(tax_table(pruned))
taxtab$OTU = c(taxa_names(pruned))
df.joined.DvD = merge(df.joined.DvD,taxtab,by=c("OTU"))
head(df.joined.DvD)

  

# WvW: Combined 1-4 -----
df.joined.WvW <- rbind(df.WvW.O.gt50, 
                       df.WvW.A.gt50,
                       df.WvW.O.lt50,
                       df.WvW.A.lt50)

# Let's bring back in our taxonomy from the tax table
SigOTUs = levels(as.factor(df.joined.WvW$OTU))
pruned = prune_taxa(SigOTUs,ps.norm.full)
taxtab = data.frame(tax_table(pruned))
taxtab$OTU = c(taxa_names(pruned))
df.joined.WvW = merge(df.joined.WvW,taxtab,by=c("OTU"))
head(df.joined.WvW)



# Save output: -----
df.joined <- rbind(df.joined.DvD, df.joined.WvW, df.joined.CvC)
write.csv(df.joined,"data/sequence-data/LibCombined/corncob-output/OTU level/All_taxa/Ex2 Fast growth/low-thermo-cutoff-50C.csv", row.names = FALSE)
#   


