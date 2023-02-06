# --------------------
# Title: FireSim2019 Experiment 2: sliding pH categories x 3
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

# Remove taxa with 0 abundance and filter for only genomic DNA:
ps.raw.prune <-  prune_samples(sample_data(ps.raw.full)$incub.trtmt %in% c('pb','SI') &
                                sample_data(ps.raw.full)$DNA.type == 'gDNA', ps.raw.full) 
ps.raw.prune <- prune_taxa(taxa_sums(ps.raw.prune)>0, ps.raw.prune) 


ps.norm.prune <-  prune_samples(sample_data(ps.norm.full)$incub.trtmt %in% c('pb','SI') &
                                 sample_data(ps.norm.full)$DNA.type == 'gDNA', ps.norm.full) 
ps.norm.prune <- prune_taxa(taxa_sums(ps.norm.prune)>0, ps.norm.prune)



### Step 2. Subset dataset for initial analyses -----
# Create phylum dataset:
#ps.raw.phylum <- tax_glom(ps.raw.prune, taxrank = 'Phylum')
#ps.raw.prune <- ps.raw.phylum


### Step 3. Create sliding pH categories -----

df <- data.frame(sample_data(ps.raw.prune))

df$pH.3to7 = ''
df$pH.4to8 = ''
df$pH.5to9 = ''

for (i in 1:nrow(df)) {
  if (df$pH[i] > 3 & df$pH[i] < 4) {
    df$pH.3to7[i] = 'y'
    df$pH.4to8[i] = 'n'
    df$pH.5to9[i] = 'n'
  } else if (df$pH[i] >= 4 & df$pH[i] < 5) {
    df$pH.3to7[i] = 'y'
    df$pH.4to8[i] = 'y'
    df$pH.5to9[i] = 'n'
  } else if (df$pH[i] >= 5 & df$pH[i] < 6) {
    df$pH.3to7[i] = 'y'
    df$pH.4to8[i] = 'y'
    df$pH.5to9[i] = 'y'
  } else if (df$pH[i] >= 6 & df$pH[i] < 7) {
    df$pH.3to7[i] = 'y'
    df$pH.4to8[i] = 'y'
    df$pH.5to9[i] = 'y'
  } else if (df$pH[i] >= 7 & df$pH[i] < 8) {
    df$pH.3to7[i] = 'n'
    df$pH.4to8[i] = 'y'
    df$pH.5to9[i] = 'y'
  } else if (df$pH[i] >= 8 & df$pH[i] < 9) {
    df$pH.3to7[i] = 'n'
    df$pH.4to8[i] = 'n'
    df$pH.5to9[i] = 'y'
  } 
}

# Merge new categories with phyloseq object
ps.df <- sample_data(df)

sample_names(ps.df) = df$Full.id

sample_data(ps.raw.prune) <- ps.df




### OPTIONAL. Set temperature cutoff for burns -----
# # Focus on dry soil burns: 
# ps.DNA.dry <-  prune_samples(sample_data(ps.raw.prune)$burn.trtmt != 'wet', ps.raw.prune)
# 
# # Create dataframe ps.DNA.dry sample_data.
# df.dry <- data.frame(sample_data(ps.DNA.dry))
# 
# # Want to pull out pairs (dry burn and control) of cores in which the dry burn
# #   core reached temps greater than some T cutoff. To do this, we need to assign
# #   a "max temp" of the core pair to both the dry burn and the control core. 
# df.dry <- df.dry %>%
#   group_by(site) %>%
#   mutate(Dry.burn.max = max(Thermo.mid.max)) 
# 
# # Now combine the new sample data with the original phyloseq object:
# ps.df <- sample_data(df.dry)
# 
# sample_names(ps.df) = df.dry$Full.id
# 
# sample_data(ps.DNA.dry) <- ps.df
# 
# ps.DvD <- prune_samples(sample_data(ps.DNA.dry)$burn.trtmt == 'dry' &
#                           sample_data(ps.DNA.dry)$Dry.burn.max > 50, ps.DNA.dry)
# 


### Step 4: Subset phyloseq by burn.trtmt and run corncob ----
ps.DvD <- prune_samples(sample_data(ps.raw.prune)$burn.trtmt == 'dry', ps.raw.prune)
ps.WvW <- prune_samples(sample_data(ps.raw.prune)$burn.trtmt == 'wet', ps.raw.prune)
ps.CvC <- prune_samples(sample_data(ps.raw.prune)$burn.trtmt == 'control', ps.raw.prune)

### pH 3 to 7, Dry burns, O horizon
ps.corncob.DvD.O.pH3to7 <-  prune_samples(sample_data(ps.DvD)$horizon == 'O' &
                                            sample_data(ps.DvD)$pH.3to7 == 'y', ps.DvD)

ps.corncob.DvD.O.pH3to7 <- prune_taxa(taxa_sums(ps.corncob.DvD.O.pH3to7)>0, ps.corncob.DvD.O.pH3to7)
ps.corncob.DvD.O.pH3to7

# Run corncob
dt.DvD.O.pH3to7 <- differentialTest(formula = ~incub.trtmt+pH,
                                    phi.formula = ~ incub.trtmt+pH,
                                    formula_null = ~ pH,
                                    phi.formula_null = ~incub.trtmt+ pH,
                                    test = 'Wald',boot=FALSE,
                                    data=ps.corncob.DvD.O.pH3to7,
                                    fdr_cutoff = 0.05)


### pH 4 to 8, Dry burns, O horizon
ps.corncob.DvD.O.pH4to8 <-  prune_samples(sample_data(ps.DvD)$horizon == 'O' &
                                            sample_data(ps.DvD)$pH.4to8 == 'y', ps.DvD)

ps.corncob.DvD.O.pH4to8 <- prune_taxa(taxa_sums(ps.corncob.DvD.O.pH4to8)>0, ps.corncob.DvD.O.pH4to8)
ps.corncob.DvD.O.pH4to8

# Run corncob
dt.DvD.O.pH4to8 <- differentialTest(formula = ~incub.trtmt+pH,
                                    phi.formula = ~ incub.trtmt+pH,
                                    formula_null = ~ pH,
                                    phi.formula_null = ~incub.trtmt+ pH,
                                    test = 'Wald',boot=FALSE,
                                    data=ps.corncob.DvD.O.pH4to8,
                                    fdr_cutoff = 0.05)



### pH 5 to 9, Dry burns, O horizon
ps.corncob.DvD.O.pH5to9 <-  prune_samples(sample_data(ps.DvD)$horizon == 'O' &
                                            sample_data(ps.DvD)$pH.5to9 == 'y', ps.DvD)

ps.corncob.DvD.O.pH5to9 <- prune_taxa(taxa_sums(ps.corncob.DvD.O.pH5to9)>0, ps.corncob.DvD.O.pH5to9)
sample_names(ps.corncob.DvD.O.pH5to9)
ps.corncob.DvD.O.pH5to9


# Run corncob
dt.DvD.O.pH5to9 <- differentialTest(formula = ~incub.trtmt+pH,
                                    phi.formula = ~ incub.trtmt+pH,
                                    formula_null = ~ pH,
                                    phi.formula_null = ~incub.trtmt+ pH,
                                    test = 'Wald',boot=FALSE,
                                    data=ps.corncob.DvD.O.pH5to9,
                                    fdr_cutoff = 0.05)


### pH 3 to 7, Dry burns, A horizon
ps.corncob.DvD.A.pH3to7 <-  prune_samples(sample_data(ps.DvD)$horizon == 'A' &
                                            sample_data(ps.DvD)$pH.3to7 == 'y', ps.DvD)

ps.corncob.DvD.A.pH3to7 <- prune_taxa(taxa_sums(ps.corncob.DvD.A.pH3to7)>0, ps.corncob.DvD.A.pH3to7)
ps.corncob.DvD.A.pH3to7

# Run corncob
dt.DvD.A.pH3to7 <- differentialTest(formula = ~incub.trtmt+pH,
                                    phi.formula = ~ incub.trtmt+pH,
                                    formula_null = ~ pH,
                                    phi.formula_null = ~incub.trtmt+ pH,
                                    test = 'Wald',boot=FALSE,
                                    data=ps.corncob.DvD.A.pH3to7,
                                    fdr_cutoff = 0.05)


### pH 4 to 8, Dry burns, A horizon
ps.corncob.DvD.A.pH4to8 <-  prune_samples(sample_data(ps.DvD)$horizon == 'A' &
                                            sample_data(ps.DvD)$pH.4to8 == 'y', ps.DvD)

ps.corncob.DvD.A.pH4to8 <- prune_taxa(taxa_sums(ps.corncob.DvD.A.pH4to8)>0, ps.corncob.DvD.A.pH4to8)
ps.corncob.DvD.A.pH4to8

# Run corncob
dt.DvD.A.pH4to8 <- differentialTest(formula = ~incub.trtmt+pH,
                                    phi.formula = ~ incub.trtmt+pH,
                                    formula_null = ~ pH,
                                    phi.formula_null = ~incub.trtmt+ pH,
                                    test = 'Wald',boot=FALSE,
                                    data=ps.corncob.DvD.A.pH4to8,
                                    fdr_cutoff = 0.05)


### pH 5 to 9, Dry burns, A horizon
ps.corncob.DvD.A.pH5to9 <-  prune_samples(sample_data(ps.DvD)$horizon == 'A' &
                                            sample_data(ps.DvD)$pH.5to9 == 'y', ps.DvD)

ps.corncob.DvD.A.pH5to9 <- prune_taxa(taxa_sums(ps.corncob.DvD.A.pH5to9)>0, ps.corncob.DvD.A.pH5to9)
ps.corncob.DvD.A.pH5to9


# Run corncob
dt.DvD.A.pH5to9 <- differentialTest(formula = ~incub.trtmt+pH,
                                    phi.formula = ~ incub.trtmt+pH,
                                    formula_null = ~ pH,
                                    phi.formula_null = ~incub.trtmt+ pH,
                                    test = 'Wald',boot=FALSE,
                                    data=ps.corncob.DvD.A.pH5to9,
                                    fdr_cutoff = 0.05)



### pH 3 to 7, wet burns, O horizon
ps.corncob.WvW.O.pH3to7 <-  prune_samples(sample_data(ps.WvW)$horizon == 'O' &
                                            sample_data(ps.WvW)$pH.3to7 == 'y', ps.WvW)

ps.corncob.WvW.O.pH3to7 <- prune_taxa(taxa_sums(ps.corncob.WvW.O.pH3to7)>0, ps.corncob.WvW.O.pH3to7)
ps.corncob.WvW.O.pH3to7

# Run corncob
dt.WvW.O.pH3to7 <- differentialTest(formula = ~incub.trtmt+pH,
                                    phi.formula = ~ incub.trtmt+pH,
                                    formula_null = ~ pH,
                                    phi.formula_null = ~incub.trtmt+ pH,
                                    test = 'Wald',boot=FALSE,
                                    data=ps.corncob.WvW.O.pH3to7,
                                    fdr_cutoff = 0.05)



### pH 4 to 8, wet burns, O horizon
ps.corncob.WvW.O.pH4to8 <-  prune_samples(sample_data(ps.WvW)$horizon == 'O' &
                                            sample_data(ps.WvW)$pH.4to8 == 'y', ps.WvW)

ps.corncob.WvW.O.pH4to8 <- prune_taxa(taxa_sums(ps.corncob.WvW.O.pH4to8)>0, ps.corncob.WvW.O.pH4to8)
ps.corncob.WvW.O.pH4to8

# Run corncob
dt.WvW.O.pH4to8 <- differentialTest(formula = ~incub.trtmt+pH,
                                    phi.formula = ~ incub.trtmt+pH,
                                    formula_null = ~ pH,
                                    phi.formula_null = ~incub.trtmt+ pH,
                                    test = 'Wald',boot=FALSE,
                                    data=ps.corncob.WvW.O.pH4to8,
                                    fdr_cutoff = 0.05)



### pH 5 to 9, wet burns, O horizon
ps.corncob.WvW.O.pH5to9 <-  prune_samples(sample_data(ps.WvW)$horizon == 'O' &
                                            sample_data(ps.WvW)$pH.5to9 == 'y', ps.WvW)

ps.corncob.WvW.O.pH5to9 <- prune_taxa(taxa_sums(ps.corncob.WvW.O.pH5to9)>0, ps.corncob.WvW.O.pH5to9)
ps.corncob.WvW.O.pH5to9


# Run corncob
dt.WvW.O.pH5to9 <- differentialTest(formula = ~incub.trtmt+pH,
                                    phi.formula = ~ incub.trtmt+pH,
                                    formula_null = ~ pH,
                                    phi.formula_null = ~incub.trtmt+ pH,
                                    test = 'Wald',boot=FALSE,
                                    data=ps.corncob.WvW.O.pH5to9,
                                    fdr_cutoff = 0.05)



### pH 3 to 7, wet burns, A horizon
ps.corncob.WvW.A.pH3to7 <-  prune_samples(sample_data(ps.WvW)$horizon == 'A' &
                                            sample_data(ps.WvW)$pH.3to7 == 'y', ps.WvW)

ps.corncob.WvW.A.pH3to7 <- prune_taxa(taxa_sums(ps.corncob.WvW.A.pH3to7)>0, ps.corncob.WvW.A.pH3to7)
ps.corncob.WvW.A.pH3to7

# Run corncob
dt.WvW.A.pH3to7 <- differentialTest(formula = ~incub.trtmt+pH,
                                    phi.formula = ~ incub.trtmt+pH,
                                    formula_null = ~ pH,
                                    phi.formula_null = ~incub.trtmt+ pH,
                                    test = 'Wald',boot=FALSE,
                                    data=ps.corncob.WvW.A.pH3to7,
                                    fdr_cutoff = 0.05)



### pH 4 to 8, wet burns, A horizon
ps.corncob.WvW.A.pH4to8 <-  prune_samples(sample_data(ps.WvW)$horizon == 'A' &
                                            sample_data(ps.WvW)$pH.4to8 == 'y', ps.WvW)

ps.corncob.WvW.A.pH4to8 <- prune_taxa(taxa_sums(ps.corncob.WvW.A.pH4to8)>0, ps.corncob.WvW.A.pH4to8)
ps.corncob.WvW.A.pH4to8

# Run corncob
dt.WvW.A.pH4to8 <- differentialTest(formula = ~incub.trtmt+pH,
                                    phi.formula = ~ incub.trtmt+pH,
                                    formula_null = ~ pH,
                                    phi.formula_null = ~incub.trtmt+ pH,
                                    test = 'Wald',boot=FALSE,
                                    data=ps.corncob.WvW.A.pH4to8,
                                    fdr_cutoff = 0.05)



### pH 5 to 9, wet burns, A horizon
ps.corncob.WvW.A.pH5to9 <-  prune_samples(sample_data(ps.WvW)$horizon == 'A' &
                                            sample_data(ps.WvW)$pH.5to9 == 'y', ps.WvW)

ps.corncob.WvW.A.pH5to9 <- prune_taxa(taxa_sums(ps.corncob.WvW.A.pH5to9)>0, ps.corncob.WvW.A.pH5to9)
ps.corncob.WvW.A.pH5to9

# Run corncob
dt.WvW.A.pH5to9 <- differentialTest(formula = ~incub.trtmt+pH,
                                    phi.formula = ~ incub.trtmt+pH,
                                    formula_null = ~ pH,
                                    phi.formula_null = ~incub.trtmt+ pH,
                                    test = 'Wald',boot=FALSE,
                                    data=ps.corncob.WvW.A.pH5to9,
                                    fdr_cutoff = 0.05)



### pH 3 to 7, control burns, O horizon
ps.corncob.CvC.O.pH3to7 <-  prune_samples(sample_data(ps.CvC)$horizon == 'O' &
                                            sample_data(ps.CvC)$pH.3to7 == 'y', ps.CvC)

ps.corncob.CvC.O.pH3to7 <- prune_taxa(taxa_sums(ps.corncob.CvC.O.pH3to7)>0, ps.corncob.CvC.O.pH3to7)
ps.corncob.CvC.O.pH3to7

# Run corncob
dt.CvC.O.pH3to7 <- differentialTest(formula = ~incub.trtmt+pH,
                                    phi.formula = ~ incub.trtmt+pH,
                                    formula_null = ~ pH,
                                    phi.formula_null = ~incub.trtmt+ pH,
                                    test = 'Wald',boot=FALSE,
                                    data=ps.corncob.CvC.O.pH3to7,
                                    fdr_cutoff = 0.05)



### pH 4 to 8, control burns, O horizon
ps.corncob.CvC.O.pH4to8 <-  prune_samples(sample_data(ps.CvC)$horizon == 'O' &
                                            sample_data(ps.CvC)$pH.4to8 == 'y', ps.CvC)

ps.corncob.CvC.O.pH4to8 <- prune_taxa(taxa_sums(ps.corncob.CvC.O.pH4to8)>0, ps.corncob.CvC.O.pH4to8)
ps.corncob.CvC.O.pH4to8

# Run corncob
dt.CvC.O.pH4to8 <- differentialTest(formula = ~incub.trtmt+pH,
                                    phi.formula = ~ incub.trtmt+pH,
                                    formula_null = ~ pH,
                                    phi.formula_null = ~incub.trtmt+ pH,
                                    test = 'Wald',boot=FALSE,
                                    data=ps.corncob.CvC.O.pH4to8,
                                    fdr_cutoff = 0.05)



### pH 5 to 9, control burns, O horizon
ps.corncob.CvC.O.pH5to9 <-  prune_samples(sample_data(ps.CvC)$horizon == 'O' &
                                            sample_data(ps.CvC)$pH.5to9 == 'y', ps.CvC)

ps.corncob.CvC.O.pH5to9 <- prune_taxa(taxa_sums(ps.corncob.CvC.O.pH5to9)>0, ps.corncob.CvC.O.pH5to9)
ps.corncob.CvC.O.pH5to9

# Run corncob
dt.CvC.O.pH5to9 <- differentialTest(formula = ~incub.trtmt+pH,
                                    phi.formula = ~ incub.trtmt+pH,
                                    formula_null = ~ pH,
                                    phi.formula_null = ~incub.trtmt+ pH,
                                    test = 'Wald',boot=FALSE,
                                    data=ps.corncob.CvC.O.pH5to9,
                                    fdr_cutoff = 0.05)



### pH 3 to 7, control burns, A horizon
ps.corncob.CvC.A.pH3to7 <-  prune_samples(sample_data(ps.CvC)$horizon == 'A' &
                                            sample_data(ps.CvC)$pH.3to7 == 'y', ps.CvC)

ps.corncob.CvC.A.pH3to7 <- prune_taxa(taxa_sums(ps.corncob.CvC.A.pH3to7)>0, ps.corncob.CvC.A.pH3to7)
ps.corncob.CvC.A.pH3to7

# Run corncob
dt.CvC.A.pH3to7 <- differentialTest(formula = ~incub.trtmt+pH,
                                    phi.formula = ~ incub.trtmt+pH,
                                    formula_null = ~ pH,
                                    phi.formula_null = ~incub.trtmt+ pH,
                                    test = 'Wald',boot=FALSE,
                                    data=ps.corncob.CvC.A.pH3to7,
                                    fdr_cutoff = 0.05)



### pH 4 to 8, control burns, A horizon
ps.corncob.CvC.A.pH4to8 <-  prune_samples(sample_data(ps.CvC)$horizon == 'A' &
                                            sample_data(ps.CvC)$pH.4to8 == 'y', ps.CvC)

ps.corncob.CvC.A.pH4to8 <- prune_taxa(taxa_sums(ps.corncob.CvC.A.pH4to8)>0, ps.corncob.CvC.A.pH4to8)
ps.corncob.CvC.A.pH4to8

# Run corncob
dt.CvC.A.pH4to8 <- differentialTest(formula = ~incub.trtmt+pH,
                                    phi.formula = ~ incub.trtmt+pH,
                                    formula_null = ~ pH,
                                    phi.formula_null = ~incub.trtmt+ pH,
                                    test = 'Wald',boot=FALSE,
                                    data=ps.corncob.CvC.A.pH4to8,
                                    fdr_cutoff = 0.05)



### pH 5 to 9, control burns, A horizon
ps.corncob.CvC.A.pH5to9 <-  prune_samples(sample_data(ps.CvC)$horizon == 'A' &
                                            sample_data(ps.CvC)$pH.5to9 == 'y', ps.CvC)

ps.corncob.CvC.A.pH5to9 <- prune_taxa(taxa_sums(ps.corncob.CvC.A.pH5to9)>0, ps.corncob.CvC.A.pH5to9)
ps.corncob.CvC.A.pH5to9

# Run corncob
dt.CvC.A.pH5to9 <- differentialTest(formula = ~incub.trtmt+pH,
                                    phi.formula = ~ incub.trtmt+pH,
                                    formula_null = ~ pH,
                                    phi.formula_null = ~incub.trtmt+ pH,
                                    test = 'Wald',boot=FALSE,
                                    data=ps.corncob.CvC.A.pH5to9,
                                    fdr_cutoff = 0.05)




### Step 5. Clean up corncob output -----

# PULL OUT COEFS FOR DRY BURNS, O HORIZON:

# pH 3 to 7

# Create empty dataframes to fill
m = data.frame('intercept'=0, 'mu.incub.trtmt'=0, 'mu.pH'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
               'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- dt.DvD.O.pH3to7
for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'pb.dry.vs.SI.dry'
  m[i,7] = 'fast growth'
  m[i,8] = 'null'
  m[i,9] = '3 to 7'
  m[i,10] = 'O'
  m[i,11] = 'gDNA'
  m[i,12] = 'pH'
}
m.concat = m

# pH 4 TO 8
m = data.frame('intercept'=0, 'mu.incub.trtmt'=0, 'mu.pH'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
               'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- dt.DvD.O.pH4to8
for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'pb.dry.vs.SI.dry'
  m[i,7] = 'fast growth'
  m[i,8] = 'null'
  m[i,9] = '4 to 8'
  m[i,10] = 'O'
  m[i,11] = 'gDNA'
  m[i,12] = 'pH'
}

m.concat = rbind(m.concat, m)



# pH 5 to 9
m = data.frame('intercept'=0, 'mu.incub.trtmt'=0, 'mu.pH'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
               'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- dt.DvD.O.pH5to9
for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'pb.dry.vs.SI.dry'
  m[i,7] = 'fast growth'
  m[i,8] = 'null'
  m[i,9] = '5 to 9'
  m[i,10] = 'O'
  m[i,11] = 'gDNA'
  m[i,12] = 'pH'
}

m.concat = rbind(m.concat, m)


# REPEAT FOR A HORIZON
# pH 3 to 7
m = data.frame('intercept'=0, 'mu.incub.trtmt'=0, 'mu.pH'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
               'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- dt.DvD.A.pH3to7
for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'pb.dry.vs.SI.dry'
  m[i,7] = 'fast growth'
  m[i,8] = 'null'
  m[i,9] = '3 to 7'
  m[i,10] = 'A'
  m[i,11] = 'gDNA'
  m[i,12] = 'pH'
}

m.concat = rbind(m.concat, m)



# pH 4 to 8
m = data.frame('intercept'=0, 'mu.incub.trtmt'=0, 'mu.pH'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
               'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- dt.DvD.A.pH4to8
for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'pb.dry.vs.SI.dry'
  m[i,7] = 'fast growth'
  m[i,8] = 'null'
  m[i,9] = '4 to 8'
  m[i,10] = 'A'
  m[i,11] = 'gDNA'
  m[i,12] = 'pH'
}

m.concat = rbind(m.concat, m)



# pH 5 to 9
m = data.frame('intercept'=0, 'mu.incub.trtmt'=0, 'mu.pH'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
               'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- dt.DvD.A.pH5to9
for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'pb.dry.vs.SI.dry'
  m[i,7] = 'fast growth'
  m[i,8] = 'null'
  m[i,9] = '5 to 9'
  m[i,10] = 'A'
  m[i,11] = 'gDNA'
  m[i,12] = 'pH'
}

m.concat = rbind(m.concat, m)




# PULL OUT COEFS FOR WET BURNS, O HORIZON:

# pH 3 to 7

# Create empty dataframes to fill
m = data.frame('intercept'=0, 'mu.incub.trtmt'=0, 'mu.pH'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
               'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- dt.WvW.O.pH3to7
for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'pb.wet.vs.SI.wet'
  m[i,7] = 'fast growth'
  m[i,8] = 'null'
  m[i,9] = '3 to 7'
  m[i,10] = 'O'
  m[i,11] = 'gDNA'
  m[i,12] = 'pH'
}
m.concat = rbind(m.concat, m)


# pH 4 TO 8
m = data.frame('intercept'=0, 'mu.incub.trtmt'=0, 'mu.pH'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
               'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- dt.WvW.O.pH4to8
for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'pb.wet.vs.SI.wet'
  m[i,7] = 'fast growth'
  m[i,8] = 'null'
  m[i,9] = '4 to 8'
  m[i,10] = 'O'
  m[i,11] = 'gDNA'
  m[i,12] = 'pH'
}

m.concat = rbind(m.concat, m)



# pH 5 to 9
m = data.frame('intercept'=0, 'mu.incub.trtmt'=0, 'mu.pH'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
               'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- dt.WvW.O.pH5to9
for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'pb.wet.vs.SI.wet'
  m[i,7] = 'fast growth'
  m[i,8] = 'null'
  m[i,9] = '5 to 9'
  m[i,10] = 'O'
  m[i,11] = 'gDNA'
  m[i,12] = 'pH'
}

m.concat = rbind(m.concat, m)


# REPEAT FOR A HORIZON
# pH 3 to 7
m = data.frame('intercept'=0, 'mu.incub.trtmt'=0, 'mu.pH'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
               'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- dt.WvW.A.pH3to7
for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'pb.wet.vs.SI.wet'
  m[i,7] = 'fast growth'
  m[i,8] = 'null'
  m[i,9] = '3 to 7'
  m[i,10] = 'A'
  m[i,11] = 'gDNA'
  m[i,12] = 'pH'
}

m.concat = rbind(m.concat, m)



# pH 4 to 8
m = data.frame('intercept'=0, 'mu.incub.trtmt'=0, 'mu.pH'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
               'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- dt.WvW.A.pH4to8
for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'pb.wet.vs.SI.wet'
  m[i,7] = 'fast growth'
  m[i,8] = 'null'
  m[i,9] = '4 to 8'
  m[i,10] = 'A'
  m[i,11] = 'gDNA'
  m[i,12] = 'pH'
}

m.concat = rbind(m.concat, m)



# pH 5 to 9
m = data.frame('intercept'=0, 'mu.incub.trtmt'=0, 'mu.pH'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
               'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- dt.WvW.A.pH5to9
for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'pb.wet.vs.SI.wet'
  m[i,7] = 'fast growth'
  m[i,8] = 'null'
  m[i,9] = '5 to 9'
  m[i,10] = 'A'
  m[i,11] = 'gDNA'
  m[i,12] = 'pH'
}

m.concat = rbind(m.concat, m)



# PULL OUT COEFS FOR CNTL BURNS, O HORIZON:

# pH 3 to 7

# Create empty dataframes to fill
m = data.frame('intercept'=0, 'mu.incub.trtmt'=0, 'mu.pH'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
               'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- dt.CvC.O.pH3to7
for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'pb.cntl.vs.SI.cntl'
  m[i,7] = 'fast growth'
  m[i,8] = 'null'
  m[i,9] = '3 to 7'
  m[i,10] = 'O'
  m[i,11] = 'gDNA'
  m[i,12] = 'pH'
}
m.concat = rbind(m.concat, m)



# pH 4 TO 8
m = data.frame('intercept'=0, 'mu.incub.trtmt'=0, 'mu.pH'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
               'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- dt.CvC.O.pH4to8
for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'pb.cntl.vs.SI.cntl'
  m[i,7] = 'fast growth'
  m[i,8] = 'null'
  m[i,9] = '4 to 8'
  m[i,10] = 'O'
  m[i,11] = 'gDNA'
  m[i,12] = 'pH'
}

m.concat = rbind(m.concat, m)



# pH 5 to 9
m = data.frame('intercept'=0, 'mu.incub.trtmt'=0, 'mu.pH'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
               'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- dt.CvC.O.pH5to9
for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'pb.cntl.vs.SI.cntl'
  m[i,7] = 'fast growth'
  m[i,8] = 'null'
  m[i,9] = '5 to 9'
  m[i,10] = 'O'
  m[i,11] = 'gDNA'
  m[i,12] = 'pH'
}

m.concat = rbind(m.concat, m)


# REPEAT FOR A HORIZON
# pH 3 to 7
m = data.frame('intercept'=0, 'mu.incub.trtmt'=0, 'mu.pH'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
               'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- dt.CvC.A.pH3to7
for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'pb.cntl.vs.SI.cntl'
  m[i,7] = 'fast growth'
  m[i,8] = 'null'
  m[i,9] = '3 to 7'
  m[i,10] = 'A'
  m[i,11] = 'gDNA'
  m[i,12] = 'pH'
}

m.concat = rbind(m.concat, m)



# pH 4 to 8
m = data.frame('intercept'=0, 'mu.incub.trtmt'=0, 'mu.pH'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
               'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- dt.CvC.A.pH4to8
for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'pb.cntl.vs.SI.cntl'
  m[i,7] = 'fast growth'
  m[i,8] = 'null'
  m[i,9] = '4 to 8'
  m[i,10] = 'A'
  m[i,11] = 'gDNA'
  m[i,12] = 'pH'
}

m.concat = rbind(m.concat, m)



# pH 5 to 9
m = data.frame('intercept'=0, 'mu.incub.trtmt'=0, 'mu.pH'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
               'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- dt.CvC.A.pH5to9
for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'pb.cntl.vs.SI.cntl'
  m[i,7] = 'fast growth'
  m[i,8] = 'null'
  m[i,9] = '5 to 9'
  m[i,10] = 'A'
  m[i,11] = 'gDNA'
  m[i,12] = 'pH'
}

m.concat = rbind(m.concat, m)




### Alternative step 4: Don't separate by burn.trtmt; run corncob ----
#
# ### pH 3 to 7, O horizon
# ps.corncob.O.pH3to7 <-  prune_samples(sample_data(ps.raw.prune)$horizon == 'O' & 
#                                             sample_data(ps.raw.prune)$pH.3to7 == 'y', ps.raw.prune)
# 
# ps.corncob.O.pH3to7 <- prune_taxa(taxa_sums(ps.corncob.O.pH3to7)>0, ps.corncob.O.pH3to7)
# ps.corncob.O.pH3to7
# 
# # Run corncob
# dt.O.pH3to7 <- differentialTest(formula = ~incub.trtmt+pH,
#                                     phi.formula = ~ incub.trtmt+pH,
#                                     formula_null = ~ pH,
#                                     phi.formula_null = ~incub.trtmt+ pH,
#                                     test = 'Wald',boot=FALSE,
#                                     data=ps.corncob.O.pH3to7,
#                                     fdr_cutoff = 0.05)
# rm(ps.corncob.O.pH3to7)
# 
# 
# ### pH 4 to 8, O horizon
# ps.corncob.O.pH4to8 <-  prune_samples(sample_data(ps.raw.prune)$horizon == 'O' & 
#                                             sample_data(ps.raw.prune)$pH.4to8 == 'y', ps.raw.prune)
# 
# ps.corncob.O.pH4to8 <- prune_taxa(taxa_sums(ps.corncob.O.pH4to8)>0, ps.corncob.O.pH4to8)
# ps.corncob.O.pH4to8
# 
# # Run corncob
# dt.O.pH4to8 <- differentialTest(formula = ~incub.trtmt+pH,
#                                     phi.formula = ~ incub.trtmt+pH,
#                                     formula_null = ~ pH,
#                                     phi.formula_null = ~incub.trtmt+ pH,
#                                     test = 'Wald',boot=FALSE,
#                                     data=ps.corncob.O.pH4to8,
#                                     fdr_cutoff = 0.05)
# rm(ps.corncob.O.pH4to8)
# 
# 
# 
# ### pH 5 to 9, Dry burns, O horizon
# ps.corncob.O.pH5to9 <-  prune_samples(sample_data(ps.raw.prune)$horizon == 'O' & 
#                                             sample_data(ps.raw.prune)$pH.5to9 == 'y', ps.raw.prune)
# 
# ps.corncob.O.pH5to9 <- prune_taxa(taxa_sums(ps.corncob.O.pH5to9)>0, ps.corncob.O.pH5to9)
# sample_names(ps.corncob.O.pH5to9)
# ps.corncob.O.pH5to9
# 
# 
# # Run corncob
# dt.O.pH5to9 <- differentialTest(formula = ~incub.trtmt+pH,
#                                     phi.formula = ~ incub.trtmt+pH,
#                                     formula_null = ~ pH,
#                                     phi.formula_null = ~incub.trtmt+ pH,
#                                     test = 'Wald',boot=FALSE,
#                                     data=ps.corncob.O.pH5to9,
#                                     fdr_cutoff = 0.05)
# rm(ps.corncob.O.pH5to9)
# 
# 
# 
# ### pH 3 to 7, DA horizon
# ps.corncob.A.pH3to7 <-  prune_samples(sample_data(ps.raw.prune)$horizon == 'A' &
#                                             sample_data(ps.raw.prune)$pH.3to7 == 'y', ps.raw.prune)
# 
# ps.corncob.A.pH3to7 <- prune_taxa(taxa_sums(ps.corncob.A.pH3to7)>0, ps.corncob.A.pH3to7)
# ps.corncob.A.pH3to7
# 
# # Run corncob
# dt.A.pH3to7 <- differentialTest(formula = ~incub.trtmt+pH,
#                                     phi.formula = ~ incub.trtmt+pH,
#                                     formula_null = ~ pH,
#                                     phi.formula_null = ~incub.trtmt+ pH,
#                                     test = 'Wald',boot=FALSE,
#                                     data=ps.corncob.A.pH3to7,
#                                     fdr_cutoff = 0.05)
# rm(ps.corncob.A.pH3to7)
# 
# 
# 
# ### pH 4 to 8, A horizon
# ps.corncob.A.pH4to8 <-  prune_samples(sample_data(ps.raw.prune)$horizon == 'A' & 
#                                             sample_data(ps.raw.prune)$pH.4to8 == 'y', ps.raw.prune)
# 
# ps.corncob.A.pH4to8 <- prune_taxa(taxa_sums(ps.corncob.A.pH4to8)>0, ps.corncob.A.pH4to8)
# ps.corncob.A.pH4to8
# 
# # Run corncob
# dt.A.pH4to8 <- differentialTest(formula = ~incub.trtmt+pH,
#                                     phi.formula = ~ incub.trtmt+pH,
#                                     formula_null = ~ pH,
#                                     phi.formula_null = ~incub.trtmt+ pH,
#                                     test = 'Wald',boot=FALSE,
#                                     data=ps.corncob.A.pH4to8,
#                                     fdr_cutoff = 0.05)
# rm(ps.corncob.A.pH4to8)
# 
# 
# 
# ### pH 5 to 9, A horizon
# ps.corncob.A.pH5to9 <-  prune_samples(sample_data(ps.raw.prune)$horizon == 'A' & 
#                                             sample_data(ps.raw.prune)$pH.5to9 == 'y', ps.raw.prune)
# 
# ps.corncob.A.pH5to9 <- prune_taxa(taxa_sums(ps.corncob.A.pH5to9)>0, ps.corncob.A.pH5to9)
# ps.corncob.A.pH5to9
# 
# 
# # Run corncob
# dt.A.pH5to9 <- differentialTest(formula = ~incub.trtmt+pH,
#                                     phi.formula = ~ incub.trtmt+pH,
#                                     formula_null = ~ pH,
#                                     phi.formula_null = ~incub.trtmt+ pH,
#                                     test = 'Wald',boot=FALSE,
#                                     data=ps.corncob.A.pH5to9,
#                                     fdr_cutoff = 0.05)
# rm(ps.corncob.A.pH5to9)





### Alterantive Step 5. Clean up corncob output -----

# PULL OUT COEFS FOR O HORIZON:

# pH 3 to 7
# Create empty dataframes to fill
m = data.frame('intercept'=0, 'mu.incub.trtmt'=0, 'mu.pH'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
               'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- dt.O.pH3to7
for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'pb.vs.SI'
  m[i,7] = 'fast growth'
  m[i,8] = 'null'
  m[i,9] = '3 to 7'
  m[i,10] = 'O'
  m[i,11] = 'gDNA'
  m[i,12] = 'pH'
}
m.concat = m

# pH 4 TO 8
m = data.frame('intercept'=0, 'mu.incub.trtmt'=0, 'mu.pH'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
               'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- dt.O.pH4to8
for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'pb.vs.SI'
  m[i,7] = 'fast growth'
  m[i,8] = 'null'
  m[i,9] = '4 to 8'
  m[i,10] = 'O'
  m[i,11] = 'gDNA'
  m[i,12] = 'pH'
}

m.concat = rbind(m.concat, m)



# pH 5 to 9
m = data.frame('intercept'=0, 'mu.incub.trtmt'=0, 'mu.pH'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
               'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- dt.O.pH5to9
for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'pb.vs.SI'
  m[i,7] = 'fast growth'
  m[i,8] = 'null'
  m[i,9] = '5 to 9'
  m[i,10] = 'O'
  m[i,11] = 'gDNA'
  m[i,12] = 'pH'
}

m.concat = rbind(m.concat, m)


# REPEAT FOR A HORIZON
# pH 3 to 7
m = data.frame('intercept'=0, 'mu.incub.trtmt'=0, 'mu.pH'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
               'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- dt.A.pH3to7
for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'pb.vs.SI'
  m[i,7] = 'fast growth'
  m[i,8] = 'null'
  m[i,9] = '3 to 7'
  m[i,10] = 'A'
  m[i,11] = 'gDNA'
  m[i,12] = 'pH'
}

m.concat = rbind(m.concat, m)



# pH 4 to 8
m = data.frame('intercept'=0, 'mu.incub.trtmt'=0, 'mu.pH'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
               'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- dt.A.pH4to8
for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'pb.vs.SI'
  m[i,7] = 'fast growth'
  m[i,8] = 'null'
  m[i,9] = '4 to 8'
  m[i,10] = 'A'
  m[i,11] = 'gDNA'
  m[i,12] = 'pH'
}

m.concat = rbind(m.concat, m)



# pH 5 to 9
m = data.frame('intercept'=0, 'mu.incub.trtmt'=0, 'mu.pH'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','experiment'='experiment','T.cutoff'='cutoff',
               'pH.cutoff'=0,'horizon'='horizon','DNA.type'='type','Controlling.for'='text')

dT <- dt.A.pH5to9
for (i in 1:length(dT$significant_taxa)) {
  m[i, 1:3] = dT$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  m[i, 4] = dT$significant_taxa[[i]] # pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  m[i,5] = dT$p_fdr[dT$significant_taxa][i]
  m[i,6] = 'pb.vs.SI'
  m[i,7] = 'fast growth'
  m[i,8] = 'null'
  m[i,9] = '5 to 9'
  m[i,10] = 'A'
  m[i,11] = 'gDNA'
  m[i,12] = 'pH'
}

m.concat = rbind(m.concat, m)


### Step 6. Clean up output -----
rm(df.joined)

# Let's bring back in our taxonomy from the tax table
SigOTUs = levels(as.factor(m.concat$OTU))
pruned = prune_taxa(SigOTUs,ps.norm.prune)
taxtab = data.frame(tax_table(pruned))
taxtab$OTU = c(taxa_names(pruned))
df.joined = merge(m.concat,taxtab,by=c("OTU"))
head(df.joined)


### Step 7. Save output -----
# Save output: 

write.csv(df.joined,"data/sequence-data/LibCombined/corncob-output/Manuscript/Ex2_signOTUs_3sliding_pH_categories_byBurnTrtmt.csv", row.names = FALSE)







