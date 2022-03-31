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


# Corncob analysis for Manuscript: 
# Input data -----
#      Input to corncob is raw count data!

ps.raw.full <- readRDS('data/sequence-data/LibCombined/phyloseq-objects/ps.raw.full')
ps.raw.full

ps.norm.full <- readRDS('data/sequence-data/LibCombined/phyloseq-objects/ps.norm.full')
ps.norm.full


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


ps.raw.full <- prune_taxa(taxa_sums(ps.raw.full)>0, ps.raw.full)
ps.norm.full <- prune_taxa(taxa_sums(ps.norm.full)>0, ps.norm.full)

# Abundance cut off and Phylum level prep -----
# subset raw data to include OTUs with a normalized relative abundance greater than some cutoff

# Mean cut-off
AbunTaxa = taxa_names(filter_taxa(ps.norm.full, function(x) mean(x) > 0.0001, TRUE))

# Cutoff for rel. abundance MAXIMUM
AbunTest <-  taxa_names(filter_taxa(ps.norm.full, function(x) max(x) > 0.005, TRUE))

length(AbunTaxa)

length(AbunTest)

ntaxa(ps.norm.full)

x <- c(AbunTaxa, AbunTest)
length(x)

x <- unique(x)
length(x)

# Prune taxa to include only OTUs with a maximum rel. abun > Cutoff
#ps.prune.raw = prune_taxa(AbunTest, ps.raw.full)



### EXPERIMENT 1 -----

#     post-burn RNA data
ps.raw.RNA <- ps.raw.full %>%
  phyloseq::subset_samples(incub.trtmt %in% c('pb')) %>%
  subset_samples(DNA.type %in% c('cDNA')) 


# DRY BURNS -----
#     Create Temp cutoff 
ps.RNA.dry <-  prune_samples(sample_data(ps.raw.RNA)$burn.trtmt != 'wet', ps.raw.RNA)

#    Create dataframe from ps.RNA.dry sample_data.
df.dry <- data.frame(sample_data(ps.RNA.dry))

#    Create temperature categories:
df.dry <- df.dry %>%
  group_by(site) %>%
  mutate(Dry.burn.max = max(Thermo.low.max)) 

#     Now combine the new sample data with the original phyloseq object:
ps.df <- sample_data(df.dry)

sample_names(ps.df) = df.dry$Full.id

sample_data(ps.RNA.dry) <- ps.df





# 1. Dry Vs. Control; O horizon; Low thermo > 50 C -----
ps.corncob.dry.O.gt50 <-  prune_samples(sample_data(ps.RNA.dry)$horizon == 'O' & 
                                          sample_data(ps.RNA.dry)$Dry.burn.max>50, ps.RNA.dry)

ps.corncob.dry.O.gt50 <- prune_taxa(taxa_sums(ps.corncob.dry.O.gt50)>0, ps.corncob.dry.O.gt50)

ps.corncob.dry.O.gt50

sample_names(ps.corncob.dry.O.gt50)
# Result 

# Run corncob:
dt.RNA.DryBurn.O.gt50 <- differentialTest(formula = ~burn.trtmt+cntl.pH,
                                           phi.formula = ~ burn.trtmt+cntl.pH,
                                           formula_null = ~ cntl.pH,
                                           phi.formula_null = ~burn.trtmt+ cntl.pH,
                                           test = 'Wald',boot=FALSE,
                                           data=ps.corncob.dry.O.gt50,
                                           fdr_cutoff = 0.05)

df.Dry.O.gt50 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.RNA.DryBurn.O.gt50$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.RNA.DryBurn.O.gt50$significant_models[[j]]
  # Pull out the coefficients as above
  # NOTE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_dry = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.RNA.DryBurn.O.gt50$p_fdr[dt.RNA.DryBurn.O.gt50$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_dry$p_fdr = p_fdr
  # Create a column with the OTU ID
  mu_dry$OTU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.Dry.O.gt50 = rbind(df.Dry.O.gt50,mu_dry)
}

df.Dry.O.gt50$Comparison = 'pb.dry.v.pb.control'
df.Dry.O.gt50$Experiment = 'survival'
df.Dry.O.gt50$Temp.Cutoff = 'Low thermo greater than 50 C'
df.Dry.O.gt50$pH.cutoff = 'null'
df.Dry.O.gt50$horizon = 'O'

colnames(df.Dry.O.gt50) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                            'Experiment', 'Temp.cutoff','pH.cutoff','horizon')




# 2. Dry Vs. Control; A horizon; Low thermo > 50 C -----
ps.corncob.dry.A.gt50 <-  prune_samples(sample_data(ps.RNA.dry)$horizon == 'A' & 
                                         sample_data(ps.RNA.dry)$Dry.burn.max>50, ps.RNA.dry)

ps.corncob.dry.A.gt50 <- prune_taxa(taxa_sums(ps.corncob.dry.A.gt50)>0, ps.corncob.dry.A.gt50)
ps.corncob.dry.A.gt50

sample_names(ps.corncob.dry.A.gt50)
# Result = 23 samples, 13 sites.
#      11-10-A = dry, no corresponding control A horizon
#      14-07-A = dry, no corresponding control A horizon
#      15-09-A = control, no corresponding dry A horizon 

# Run corncob
dt.RNA.DryBurn.A.gt50 <- differentialTest(formula = ~burn.trtmt+cntl.pH,
                                          phi.formula = ~ burn.trtmt+cntl.pH,
                                          formula_null = ~ cntl.pH,
                                          phi.formula_null = ~burn.trtmt+ cntl.pH,
                                          test = 'Wald',boot=FALSE,
                                          data=ps.corncob.dry.A.gt50,
                                          fdr_cutoff = 0.05)

df.Dry.A.gt50 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.RNA.DryBurn.A.gt50$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.RNA.DryBurn.A.gt50$significant_models[[j]]
  # Pull out the coefficients as above
  # NOTE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_dry = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.RNA.DryBurn.A.gt50$p_fdr[dt.RNA.DryBurn.A.gt50$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_dry$p_fdr = p_fdr
  # Create a column with the OTU ID
  mu_dry$OTU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.Dry.A.gt50 = rbind(df.Dry.A.gt50,mu_dry)
}

df.Dry.A.gt50$Comparison = 'pb.dry.v.pb.control'
df.Dry.A.gt50$Experiment = 'survival'
df.Dry.A.gt50$Temp.Cutoff = 'Low thermo greater than 50 C'
df.Dry.A.gt50$pH.cutoff = 'null'
df.Dry.A.gt50$horizon = 'A'

colnames(df.Dry.A.gt50) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                            'Experiment', 'Temp.cutoff','pH.cutoff','horizon')



# 3. Dry Vs. Control; O horizon; Low thermo < 50 C -----
ps.corncob.dry.O.lt50 <-  prune_samples(sample_data(ps.RNA.dry)$horizon == 'O' & 
                                          sample_data(ps.RNA.dry)$Dry.burn.max<50, ps.RNA.dry)

ps.corncob.dry.O.lt50 <- prune_taxa(taxa_sums(ps.corncob.dry.O.lt50)>0, ps.corncob.dry.O.lt50)
ps.corncob.dry.O.lt50

sample_names(ps.corncob.dry.O.lt50)
# Result = 5 samples, 3 sites 
#      08-05-O is a control core with no corresponding dry burn (failed RNA extraction)

# Run corncob
dt.RNA.DryBurn.O.lt50 <- differentialTest(formula = ~burn.trtmt+cntl.pH,
                                          phi.formula = ~ burn.trtmt+cntl.pH,
                                          formula_null = ~ cntl.pH,
                                          phi.formula_null = ~burn.trtmt+ cntl.pH,
                                          test = 'Wald',boot=FALSE,
                                          data=ps.corncob.dry.O.lt50,
                                          fdr_cutoff = 0.05)
#     All models failed to converge at temp cutoff of <50C


df.Dry.O.lt50 = data.frame()


# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.RNA.DryBurn.O.lt50$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.RNA.DryBurn.O.lt50$significant_models[[j]]
  # Pull out the coefficients as above
  # NOTE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_dry = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.RNA.DryBurn.O.lt50$p_fdr[dt.RNA.DryBurn.O.lt50$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_dry$p_fdr = p_fdr
  # Create a column with the OTU ID
  mu_dry$OTU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.Dry.O.lt50 = rbind(df.Dry.O.lt50,mu_dry)
}

df.Dry.O.lt50$Comparison = 'pb.dry.v.pb.control'
df.Dry.O.lt50$Experiment = 'survival'
df.Dry.O.lt50$Temp.Cutoff = 'Low thermo less than 50 C'
df.Dry.O.lt50$pH.cutoff = 'null'
df.Dry.O.lt50$horizon = 'O'

colnames(df.Dry.O.lt50) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                            'Experiment', 'Temp.cutoff','pH.cutoff','horizon')




# 4. Dry Vs. Control; A horizon; Low thermo < 50 C -----

ps.corncob.dry.A.lt50 <-  prune_samples(sample_data(ps.RNA.dry)$horizon == 'A' & 
                                          sample_data(ps.RNA.dry)$Dry.burn.max<50, ps.RNA.dry)

ps.corncob.dry.A.lt50 <- prune_taxa(taxa_sums(ps.corncob.dry.A.lt50)>0, ps.corncob.dry.A.lt50)

(ps.corncob.dry.A.lt50)
# Result = 1 sample, 1 site
#      04-03-A = dry, no corresponding control A horizon

# Run corncob
dt.RNA.DryBurn.A.lt50 <- differentialTest(formula = ~burn.trtmt+cntl.pH,
                                          phi.formula = ~ burn.trtmt+cntl.pH,
                                          formula_null = ~ cntl.pH,
                                          phi.formula_null = ~burn.trtmt+ cntl.pH,
                                          test = 'Wald',boot=FALSE,
                                          data=ps.corncob.dry.A.lt50,
                                          fdr_cutoff = 0.05)
#     All models failed to converge at temp cutoff of <50C

df.Dry.A.lt50 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.RNA.DryBurn.A.lt50$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.RNA.DryBurn.A.lt50$significant_models[[j]]
  # Pull out the coefficients as above
  # NOTE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_dry = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.RNA.DryBurn.A.lt50$p_fdr[dt.RNA.DryBurn.A.lt50$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_dry$p_fdr = p_fdr
  # Create a column with the OTU ID
  mu_dry$OTU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.Dry.A.lt50 = rbind(df.Dry.A.lt50,mu_dry)
}

df.Dry.A.lt50$Comparison = 'pb.dry.v.pb.control'
df.Dry.A.lt50$Experiment = 'survival'
df.Dry.A.lt50$Temp.Cutoff = 'Low thermo less than 50 C'
df.Dry.A.lt50$pH.cutoff = 'null'
df.Dry.A.lt50$horizon = 'A'

colnames(df.Dry.A.lt50) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                            'Experiment', 'Temp.cutoff','pH.cutoff','horizon')




# Combine 1-4 -----
df.joined.dry <- rbind(df.Dry.O.gt50, df.Dry.A.gt50,
                       df.Dry.O.lt50, df.Dry.A.lt50)


### WET BURNS: -----

ps.RNA.wet <-  prune_samples(sample_data(ps.raw.RNA)$burn.trtmt != 'dry', ps.raw.RNA)

#    Create dataframe from ps.RNA.wet sample_data.
df.wet <- data.frame(sample_data(ps.RNA.wet))

#    Create temperature categories:
df.wet <- df.wet %>%
  group_by(site) %>%
  mutate(Wet.burn.max = max(Thermo.low.max)) 

#     Now combine the new sample data with the original phyloseq object:
ps.df <- sample_data(df.wet)

sample_names(ps.df) = df.wet$Full.id

sample_data(ps.RNA.wet) <- ps.df


# 1. Wet Vs. Control; O horizon; Low thermo > 50 C -----

ps.corncob.wet.O.gt50 <-  prune_samples(sample_data(ps.RNA.wet)$horizon == 'O' & 
                                          sample_data(ps.RNA.wet)$Wet.burn.max>50, ps.RNA.wet)

ps.corncob.wet.O.gt50 <- prune_taxa(taxa_sums(ps.corncob.wet.O.gt50)>0, ps.corncob.wet.O.gt50)
ps.corncob.wet.O.gt50

sample_names(ps.corncob.wet.O.gt50)
# Result = 0 samples

# Run corncob:
dt.RNA.wetBurn.O.gt50 <- differentialTest(formula = ~burn.trtmt+cntl.pH,
                                          phi.formula = ~ burn.trtmt+cntl.pH,
                                          formula_null = ~ cntl.pH,
                                          phi.formula_null = ~burn.trtmt+ cntl.pH,
                                          test = 'Wald',boot=FALSE,
                                          data=ps.corncob.wet.O.gt50,
                                          fdr_cutoff = 0.05)
# Notes: 

df.wet.O.gt50 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.RNA.wetBurn.O.gt50$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.RNA.wetBurn.O.gt50$significant_models[[j]]
  # Pull out the coefficients as above
  # NOTE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_dry = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.RNA.wetBurn.O.gt50$p_fdr[dt.RNA.wetBurn.O.gt50$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_dry$p_fdr = p_fdr
  # Create a column with the OTU ID
  mu_dry$OTU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.wet.O.gt50 = rbind(df.wet.O.gt50,mu_dry)
}

df.wet.O.gt50$Comparison = 'pb.wet.v.pb.control'
df.wet.O.gt50$Experiment = 'survival'
df.wet.O.gt50$Temp.Cutoff = 'Low thermo greater than 50 C'
df.wet.O.gt50$pH.cutoff = 'null'
df.wet.O.gt50$horizon = 'O'

colnames(df.wet.O.gt50) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                            'Experiment', 'Temp.cutoff','pH.cutoff','horizon')




# 2.  Wet Vs. Control; A horizon; Low thermo > 50 C -----
ps.corncob.wet.A.gt50 <-  prune_samples(sample_data(ps.RNA.wet)$horizon == 'A' & 
                                          sample_data(ps.RNA.wet)$Wet.burn.max>50, ps.RNA.wet)

ps.corncob.wet.A.gt50 <- prune_taxa(taxa_sums(ps.corncob.wet.A.gt50)>0, ps.corncob.wet.A.gt50)
ps.corncob.wet.A.gt50

sample_names(ps.corncob.wet.A.gt50)
# Result = 0 samples


# Run corncob
dt.RNA.wetBurn.A.gt50 <- differentialTest(formula = ~burn.trtmt+cntl.pH,
                                          phi.formula = ~ burn.trtmt+cntl.pH,
                                          formula_null = ~ cntl.pH,
                                          phi.formula_null = ~burn.trtmt+ cntl.pH,
                                          test = 'Wald',boot=FALSE,
                                          data=ps.corncob.wet.A.gt50,
                                          fdr_cutoff = 0.05)
# Notes: 
df.wet.A.gt50 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.RNA.wetBurn.A.gt50$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.RNA.wetBurn.A.gt50$significant_models[[j]]
  # Pull out the coefficients as above
  # NOTE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_dry = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.RNA.wetBurn.A.gt50$p_fdr[dt.RNA.wetBurn.A.gt50$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_dry$p_fdr = p_fdr
  # Create a column with the OTU ID
  mu_dry$OTU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.wet.A.gt50 = rbind(df.wet.A.gt50,mu_dry)
}

df.wet.A.gt50$Comparison = 'pb.wet.v.pb.control'
df.wet.A.gt50$Experiment = 'survival'
df.wet.A.gt50$Temp.Cutoff = 'Low thermo greater than 50 C'
df.wet.A.gt50$pH.cutoff = 'null'
df.wet.A.gt50$horizon = 'A'

colnames(df.wet.A.gt50) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                            'Experiment', 'Temp.cutoff','pH.cutoff','horizon')



# 3.  Wet Vs. Control; O horizon; Low thermo < 50 C -----
ps.corncob.wet.O.lt50 <-  prune_samples(sample_data(ps.RNA.wet)$horizon == 'O' & 
                                          sample_data(ps.RNA.wet)$Wet.burn.max<50, ps.RNA.wet)

ps.corncob.wet.O.lt50 <- prune_taxa(taxa_sums(ps.corncob.wet.O.lt50)>0, ps.corncob.wet.O.lt50)
ps.corncob.wet.O.lt50

sample_names(ps.corncob.wet.O.lt50)
# Result = 38 samples, 19 sites

# Run corncob
dt.RNA.wetBurn.O.lt50 <- differentialTest(formula = ~burn.trtmt+cntl.pH,
                                          phi.formula = ~ burn.trtmt+cntl.pH,
                                          formula_null = ~ cntl.pH,
                                          phi.formula_null = ~burn.trtmt+ cntl.pH,
                                          test = 'Wald',boot=FALSE,
                                          data=ps.corncob.wet.O.lt50,
                                          fdr_cutoff = 0.05)
# Notes:  


df.wet.O.lt50 = data.frame()


# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.RNA.wetBurn.O.lt50$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.RNA.wetBurn.O.lt50$significant_models[[j]]
  # Pull out the coefficients as above
  # NOTE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_dry = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.RNA.wetBurn.O.lt50$p_fdr[dt.RNA.wetBurn.O.lt50$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_dry$p_fdr = p_fdr
  # Create a column with the OTU ID
  mu_dry$OTU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.wet.O.lt50 = rbind(df.wet.O.lt50,mu_dry)
}

df.wet.O.lt50$Comparison = 'pb.wet.v.pb.control'
df.wet.O.lt50$Experiment = 'survival'
df.wet.O.lt50$Temp.Cutoff = 'Low thermo less than 50 C'
df.wet.O.lt50$pH.cutoff = 'null'
df.wet.O.lt50$horizon = 'O'

colnames(df.wet.O.lt50) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                            'Experiment', 'Temp.cutoff','pH.cutoff','horizon')




# 4.  Wet Vs. Control; A horizon; Low thermo < 50 C -----
ps.corncob.wet.A.lt50 <-  prune_samples(sample_data(ps.RNA.wet)$horizon == 'A' & 
                                          sample_data(ps.RNA.wet)$Wet.burn.max<50, ps.RNA.wet)

ps.corncob.wet.A.lt50 <- prune_taxa(taxa_sums(ps.corncob.wet.A.lt50)>0, ps.corncob.wet.A.lt50)
ps.corncob.wet.A.lt50

sample_names(ps.corncob.wet.A.lt50)
# Result = 23 samples, thirteen sites:
# 11-03-A 
# 14-04-A
# 15-09-A

# Run corncob
dt.RNA.wetBurn.A.lt50 <- differentialTest(formula = ~burn.trtmt+cntl.pH,
                                          phi.formula = ~ burn.trtmt+cntl.pH,
                                          formula_null = ~ cntl.pH,
                                          phi.formula_null = ~burn.trtmt+ cntl.pH,
                                          test = 'Wald',boot=FALSE,
                                          data=ps.corncob.wet.A.lt50,
                                          fdr_cutoff = 0.05)
# Notes:  


df.wet.A.lt50 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.RNA.wetBurn.A.lt50$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.RNA.wetBurn.A.lt50$significant_models[[j]]
  # Pull out the coefficients as above
  # NOTE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_dry = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.RNA.wetBurn.A.lt50$p_fdr[dt.RNA.wetBurn.A.lt50$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_dry$p_fdr = p_fdr
  # Create a column with the OTU ID
  mu_dry$OTU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.wet.A.lt50 = rbind(df.wet.A.lt50,mu_dry)
}

df.wet.A.lt50$Comparison = 'pb.wet.v.pb.control'
df.wet.A.lt50$Experiment = 'survival'
df.wet.A.lt50$Temp.Cutoff = 'Low thermo less than 50 C'
df.wet.A.lt50$pH.cutoff = 'null'
df.wet.A.lt50$horizon = 'A'

colnames(df.wet.A.lt50) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                            'Experiment', 'Temp.cutoff','pH.cutoff','horizon')





# Combined 1-4 -----
df.joined.wet <- rbind(df.wet.O.gt50, df.wet.A.gt50,
                       df.wet.O.lt50, df.wet.A.lt50)

# Combine dry and wet -----

head(df.joined.dry)
head(df.joined.wet)

df.joined <- rbind(df.joined.dry, df.joined.wet)



# Let's bring back in our taxonomy from the tax table
SigOTUs = levels(as.factor(df.joined$OTU))
pruned = prune_taxa(SigOTUs,ps.norm.full)
taxtab = data.frame(tax_table(pruned))
taxtab$OTU = c(taxa_names(pruned))
df.joined = merge(df.joined,taxtab,by=c("OTU"))
head(df.joined)


# May want to save at this point
write.csv(df.joined,"data/sequence-data/LibCombined/corncob-output/OTU level/All_taxa/Ex1 Fire survival/low-thermo-cutoff-50C.csv", row.names = FALSE)

#otu_to_taxonomy(OTU = dt.RNA.wetBurn.A.lt50$significant_taxa, data = ps.RNA.wet)
#plot(dt.RNA.wetBurn.A.lt50, level = c('Phylum'))


