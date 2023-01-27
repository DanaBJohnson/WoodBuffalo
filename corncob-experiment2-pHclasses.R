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

# Remove taxa with 0 abundance:
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


### Step 4. Create sliding pH categories -----

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


### Step 4 XXX, dry burns ----

ps.DvD <- prune_samples(sample_data(ps.raw.prune)$burn.trtmt == 'dry', ps.raw.prune)

### pH 3 to 7, Dry burns, O horizon -----
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

length(dt.DvD.O.pH3to7$significant_taxa)
# Result = 104 significant taxa


df.DvD.O.pH3to7 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.DvD.O.pH3to7$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.DvD.O.pH3to7$significant_models[[j]]
  # Pull out the coefficients as above
  # NOTE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_dry = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.DvD.O.pH3to7$p_fdr[dt.DvD.O.pH3to7$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_dry$p_fdr = p_fdr
  # Create a column with the OTU ID
  mu_dry$OTU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.DvD.O.pH3to7 = rbind(df.DvD.O.pH3to7,mu_dry)
}

df.DvD.O.pH3to7$Comparison = 'pb.dry.v.SI.dry'
df.DvD.O.pH3to7$Experiment = 'fast growth'
df.DvD.O.pH3to7$Temp.Cutoff = 'null'
df.DvD.O.pH3to7$pH.cutoff = '3 to 7'
df.DvD.O.pH3to7$horizon = 'O'
df.DvD.O.pH3to7$DNA.type = 'gDNA'
df.DvD.O.pH3to7$Controlling.for = 'pH'

colnames(df.DvD.O.pH3to7) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                              'Experiment', 'Temp.cutoff','pH.cutoff','horizon', 
                              'DNA.type', 'Controlling.for')



### pH 4 to 8, Dry burns, O horizon -----

ps.corncob.DvD.O.pH4to8 <-  prune_samples(sample_data(ps.DvD)$horizon == 'O' & 
                                            sample_data(ps.DvD)$pH.4to8 == 'y', ps.DvD)

ps.corncob.DvD.O.pH4to8 <- prune_taxa(taxa_sums(ps.corncob.DvD.O.pH4to8)>0, ps.corncob.DvD.O.pH4to8)
sample_names(ps.corncob.DvD.O.pH4to8)
ps.corncob.DvD.O.pH4to8


# Run corncob
dt.DvD.O.pH4to8 <- differentialTest(formula = ~incub.trtmt+pH,
                                    phi.formula = ~ incub.trtmt+pH,
                                    formula_null = ~ pH,
                                    phi.formula_null = ~incub.trtmt+ pH,
                                    test = 'Wald',boot=FALSE,
                                    data=ps.corncob.DvD.O.pH4to8,
                                    fdr_cutoff = 0.05)
length(dt.DvD.O.pH4to8$significant_taxa)
# Result = 140 significant taxa


df.DvD.O.pH4to8 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.DvD.O.pH4to8$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.DvD.O.pH4to8$significant_models[[j]]
  # Pull out the coefficients as above
  # NOTE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_dry = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.DvD.O.pH4to8$p_fdr[dt.DvD.O.pH4to8$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_dry$p_fdr = p_fdr
  # Create a column with the OTU ID
  mu_dry$OTU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.DvD.O.pH4to8 = rbind(df.DvD.O.pH4to8,mu_dry)
}

df.DvD.O.pH4to8$Comparison = 'pb.dry.v.SI.dry'
df.DvD.O.pH4to8$Experiment = 'fast growth'
df.DvD.O.pH4to8$Temp.Cutoff = 'null'
df.DvD.O.pH4to8$pH.cutoff = '4 to 8'
df.DvD.O.pH4to8$horizon = 'O'
df.DvD.O.pH4to8$DNA.type = 'gDNA'
df.DvD.O.pH4to8$Controlling.for = 'pH'

colnames(df.DvD.O.pH4to8) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                              'Experiment', 'Temp.cutoff','pH.cutoff','horizon', 
                              'DNA.type', 'Controlling.for')

### pH 5 to 9, Dry burns, O horizon -----
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
length(dt.DvD.O.pH5to9$significant_taxa)
# Result = 89 significant taxa

df.DvD.O.pH5to9 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.DvD.O.pH5to9$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.DvD.O.pH5to9$significant_models[[j]]
  # Pull out the coefficients as above
  # NOTE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_dry = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.DvD.O.pH5to9$p_fdr[dt.DvD.O.pH5to9$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_dry$p_fdr = p_fdr
  # Create a column with the OTU ID
  mu_dry$OTU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.DvD.O.pH5to9 = rbind(df.DvD.O.pH5to9,mu_dry)
}

df.DvD.O.pH5to9$Comparison = 'pb.dry.v.SI.dry'
df.DvD.O.pH5to9$Experiment = 'fast growth'
df.DvD.O.pH5to9$Temp.Cutoff = 'null'
df.DvD.O.pH5to9$pH.cutoff = '5 to 9'
df.DvD.O.pH5to9$horizon = 'O'
df.DvD.O.pH5to9$DNA.type = 'gDNA'
df.DvD.O.pH5to9$Controlling.for = 'pH'

colnames(df.DvD.O.pH5to9) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                              'Experiment', 'Temp.cutoff','pH.cutoff','horizon', 
                              'DNA.type', 'Controlling.for')

### pH 3 to 7, Dry burns, A horizon -----
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
length(dt.DvD.A.pH3to7$significant_taxa)
# Result = 30 significant taxa



df.DvD.A.pH3to7 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.DvD.A.pH3to7$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.DvD.A.pH3to7$significant_models[[j]]
  # Pull out the coefficients as above
  # NATE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_dry = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.DvD.A.pH3to7$p_fdr[dt.DvD.A.pH3to7$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_dry$p_fdr = p_fdr
  # Create a column with the ATU ID
  mu_dry$ATU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.DvD.A.pH3to7 = rbind(df.DvD.A.pH3to7,mu_dry)
}

df.DvD.A.pH3to7$Comparison = 'pb.dry.v.SI.dry'
df.DvD.A.pH3to7$Experiment = 'fast growth'
df.DvD.A.pH3to7$Temp.Cutoff = 'null'
df.DvD.A.pH3to7$pH.cutoff = '3 to 7'
df.DvD.A.pH3to7$horizon = 'A'
df.DvD.A.pH3to7$DNA.type = 'gDNA'
df.DvD.A.pH3to7$Controlling.for = 'pH'

colnames(df.DvD.A.pH3to7) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                              'Experiment', 'Temp.cutoff','pH.cutoff','horizon', 
                              'DNA.type', 'Controlling.for')

### pH 4 to 8, Dry burns, A horizon -----

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
length(dt.DvD.A.pH4to8$significant_taxa)
# Result = 54 significant taxa


df.DvD.A.pH4to8 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.DvD.A.pH4to8$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.DvD.A.pH4to8$significant_models[[j]]
  # Pull out the coefficients as above
  # NATE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_dry = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.DvD.A.pH4to8$p_fdr[dt.DvD.A.pH4to8$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_dry$p_fdr = p_fdr
  # Create a column with the ATU ID
  mu_dry$ATU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.DvD.A.pH4to8 = rbind(df.DvD.A.pH4to8,mu_dry)
}

df.DvD.A.pH4to8$Comparison = 'pb.dry.v.SI.dry'
df.DvD.A.pH4to8$Experiment = 'fast growth'
df.DvD.A.pH4to8$Temp.Cutoff = 'null'
df.DvD.A.pH4to8$pH.cutoff = '4 to 8'
df.DvD.A.pH4to8$horizon = 'A'
df.DvD.A.pH4to8$DNA.type = 'gDNA'
df.DvD.A.pH4to8$Controlling.for = 'pH'

colnames(df.DvD.A.pH4to8) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                              'Experiment', 'Temp.cutoff','pH.cutoff','horizon', 
                              'DNA.type', 'Controlling.for')

### pH 5 to 9, Dry burns, A horizon -----
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
length(dt.DvD.A.pH5to9$significant_taxa)
# Result = 96 significant taxa


df.DvD.A.pH5to9 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.DvD.A.pH5to9$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.DvD.A.pH5to9$significant_models[[j]]
  # Pull out the coefficients as above
  # NATE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_dry = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.DvD.A.pH5to9$p_fdr[dt.DvD.A.pH5to9$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_dry$p_fdr = p_fdr
  # Create a column with the ATU ID
  mu_dry$ATU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.DvD.A.pH5to9 = rbind(df.DvD.A.pH5to9,mu_dry)
}

df.DvD.A.pH5to9$Comparison = 'pb.dry.v.SI.dry'
df.DvD.A.pH5to9$Experiment = 'fast growth'
df.DvD.A.pH5to9$Temp.Cutoff = 'null'
df.DvD.A.pH5to9$pH.cutoff = '5 to 9'
df.DvD.A.pH5to9$horizon = 'A'
df.DvD.A.pH5to9$DNA.type = 'gDNA'
df.DvD.A.pH5to9$Controlling.for = 'pH'

colnames(df.DvD.A.pH5to9) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                              'Experiment', 'Temp.cutoff','pH.cutoff','horizon', 
                              'DNA.type', 'Controlling.for')


### Combine pH category output -----
df.DvD.pH.groups <- rbind(df.DvD.O.pH3to7,df.DvD.O.pH4to8,df.DvD.O.pH5to9,
                          df.DvD.A.pH3to7,df.DvD.A.pH4to8,df.DvD.A.pH5to9)



dim(df.DvD.O.pH3to7)
dim(df.DvD.O.pH4to8)
dim(df.DvD.O.pH5to9)
dim(df.DvD.A.pH3to7)
dim(df.DvD.A.pH4to8)
dim(df.DvD.A.pH5to9)



# Let's bring back in our taxonomy from the tax table
SigOTUs = levels(as.factor(df.DvD.pH.groups$OTU))
pruned = prune_taxa(SigOTUs,ps.norm.full)
taxtab = data.frame(tax_table(pruned))
taxtab$OTU = c(taxa_names(pruned))
df.joined.DvD.pH = merge(df.DvD.pH.groups,taxtab,by=c("OTU"))
head(df.joined.DvD.pH)




# EXPERIMENT 2 - pH categories, wet -----

ps.WvW <- prune_samples(sample_data(ps.raw.prune)$burn.trtmt == 'wet', ps.raw.prune)

### pH 3 to 7, wet burns, O horizon -----
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

df.WvW.O.pH3to7 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.WvW.O.pH3to7$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.WvW.O.pH3to7$significant_models[[j]]
  # Pull out the coefficients as above
  # NOTE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_wet = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.WvW.O.pH3to7$p_fdr[dt.WvW.O.pH3to7$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_wet$p_fdr = p_fdr
  # Create a column with the OTU ID
  mu_wet$OTU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.WvW.O.pH3to7 = rbind(df.WvW.O.pH3to7,mu_wet)
}

df.WvW.O.pH3to7$Comparison = 'pb.wet.v.SI.wet'
df.WvW.O.pH3to7$Experiment = 'fast growth'
df.WvW.O.pH3to7$Temp.Cutoff = 'null'
df.WvW.O.pH3to7$pH.cutoff = '3 to 7'
df.WvW.O.pH3to7$horizon = 'O'
df.WvW.O.pH3to7$DNA.type = 'gDNA'
df.WvW.O.pH3to7$Controlling.for = 'pH'

colnames(df.WvW.O.pH3to7) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                              'Experiment', 'Temp.cutoff','pH.cutoff','horizon', 
                              'DNA.type', 'Controlling.for')


### pH 4 to 8, wet burns, O horizon -----

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

df.WvW.O.pH4to8 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.WvW.O.pH4to8$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.WvW.O.pH4to8$significant_models[[j]]
  # Pull out the coefficients as above
  # NOTE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_wet = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.WvW.O.pH4to8$p_fdr[dt.WvW.O.pH4to8$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_wet$p_fdr = p_fdr
  # Create a column with the OTU ID
  mu_wet$OTU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.WvW.O.pH4to8 = rbind(df.WvW.O.pH4to8,mu_wet)
}

df.WvW.O.pH4to8$Comparison = 'pb.wet.v.SI.wet'
df.WvW.O.pH4to8$Experiment = 'fast growth'
df.WvW.O.pH4to8$Temp.Cutoff = 'null'
df.WvW.O.pH4to8$pH.cutoff = '4 to 8'
df.WvW.O.pH4to8$horizon = 'O'
df.WvW.O.pH4to8$DNA.type = 'gDNA'
df.WvW.O.pH4to8$Controlling.for = 'pH'

colnames(df.WvW.O.pH4to8) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                              'Experiment', 'Temp.cutoff','pH.cutoff','horizon', 
                              'DNA.type', 'Controlling.for')




### pH 5 to 9, wet burns, O horizon -----
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

df.WvW.O.pH5to9 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.WvW.O.pH5to9$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.WvW.O.pH5to9$significant_models[[j]]
  # Pull out the coefficients as above
  # NOTE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_wet = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.WvW.O.pH5to9$p_fdr[dt.WvW.O.pH5to9$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_wet$p_fdr = p_fdr
  # Create a column with the OTU ID
  mu_wet$OTU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.WvW.O.pH5to9 = rbind(df.WvW.O.pH5to9,mu_wet)
}

df.WvW.O.pH5to9$Comparison = 'pb.wet.v.SI.wet'
df.WvW.O.pH5to9$Experiment = 'fast growth'
df.WvW.O.pH5to9$Temp.Cutoff = 'null'
df.WvW.O.pH5to9$pH.cutoff = '5 to 9'
df.WvW.O.pH5to9$horizon = 'O'
df.WvW.O.pH5to9$DNA.type = 'gDNA'
df.WvW.O.pH5to9$Controlling.for = 'pH'

colnames(df.WvW.O.pH5to9) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                              'Experiment', 'Temp.cutoff','pH.cutoff','horizon', 
                              'DNA.type', 'Controlling.for')


### pH 3 to 7, wet burns, A horizon -----
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

df.WvW.A.pH3to7 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.WvW.A.pH3to7$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.WvW.A.pH3to7$significant_models[[j]]
  # Pull out the coefficients as above
  # NATE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_wet = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.WvW.A.pH3to7$p_fdr[dt.WvW.A.pH3to7$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_wet$p_fdr = p_fdr
  # Create a column with the ATU ID
  mu_wet$ATU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.WvW.A.pH3to7 = rbind(df.WvW.A.pH3to7,mu_wet)
}

df.WvW.A.pH3to7$Comparison = 'pb.wet.v.SI.wet'
df.WvW.A.pH3to7$Experiment = 'fast growth'
df.WvW.A.pH3to7$Temp.Cutoff = 'null'
df.WvW.A.pH3to7$pH.cutoff = '3 to 7'
df.WvW.A.pH3to7$horizon = 'A'
df.WvW.A.pH3to7$DNA.type = 'gDNA'
df.WvW.A.pH3to7$Controlling.for = 'pH'

colnames(df.WvW.A.pH3to7) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                              'Experiment', 'Temp.cutoff','pH.cutoff','horizon', 
                              'DNA.type', 'Controlling.for')

### pH 4 to 8, wet burns, A horizon -----

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

df.WvW.A.pH4to8 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.WvW.A.pH4to8$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.WvW.A.pH4to8$significant_models[[j]]
  # Pull out the coefficients as above
  # NATE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_wet = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.WvW.A.pH4to8$p_fdr[dt.WvW.A.pH4to8$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_wet$p_fdr = p_fdr
  # Create a column with the ATU ID
  mu_wet$ATU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.WvW.A.pH4to8 = rbind(df.WvW.A.pH4to8,mu_wet)
}

df.WvW.A.pH4to8$Comparison = 'pb.wet.v.SI.wet'
df.WvW.A.pH4to8$Experiment = 'fast growth'
df.WvW.A.pH4to8$Temp.Cutoff = 'null'
df.WvW.A.pH4to8$pH.cutoff = '4 to 8'
df.WvW.A.pH4to8$horizon = 'A'
df.WvW.A.pH4to8$DNA.type = 'gDNA'
df.WvW.A.pH4to8$Controlling.for = 'pH'

colnames(df.WvW.A.pH4to8) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                              'Experiment', 'Temp.cutoff','pH.cutoff','horizon', 
                              'DNA.type', 'Controlling.for')

### pH 5 to 9, wet burns, A horizon -----
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

df.WvW.A.pH5to9 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.WvW.A.pH5to9$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.WvW.A.pH5to9$significant_models[[j]]
  # Pull out the coefficients as above
  # NATE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_wet = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.WvW.A.pH5to9$p_fdr[dt.WvW.A.pH5to9$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_wet$p_fdr = p_fdr
  # Create a column with the ATU ID
  mu_wet$ATU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.WvW.A.pH5to9 = rbind(df.WvW.A.pH5to9,mu_wet)
}

df.WvW.A.pH5to9$Comparison = 'pb.wet.v.SI.wet'
df.WvW.A.pH5to9$Experiment = 'fast growth'
df.WvW.A.pH5to9$Temp.Cutoff = 'null'
df.WvW.A.pH5to9$pH.cutoff = '5 to 9'
df.WvW.A.pH5to9$horizon = 'A'
df.WvW.A.pH5to9$DNA.type = 'gDNA'
df.WvW.A.pH5to9$Controlling.for = 'pH'

colnames(df.WvW.A.pH5to9) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                              'Experiment', 'Temp.cutoff','pH.cutoff','horizon', 
                              'DNA.type', 'Controlling.for')

### Combine pH category output -----
ps.corncob.WvW.O.pH3to7
ps.corncob.WvW.O.pH4to8
ps.corncob.WvW.O.pH5to9
ps.corncob.WvW.A.pH3to7
ps.corncob.WvW.A.pH4to8
ps.corncob.WvW.A.pH5to9

length(dt.WvW.O.pH3to7$significant_taxa)
length(dt.WvW.O.pH4to8$significant_taxa)
length(dt.WvW.O.pH5to9$significant_taxa)
length(dt.WvW.A.pH3to7$significant_taxa)
length(dt.WvW.A.pH4to8$significant_taxa)
length(dt.WvW.A.pH5to9$significant_taxa)

df.WvW.pH.groups <- rbind(df.WvW.O.pH3to7,
                          df.WvW.O.pH4to8,
                          df.WvW.O.pH5to9,
                          df.WvW.A.pH3to7,
                          df.WvW.A.pH4to8,
                          df.WvW.A.pH5to9)


dim(df.WvW.O.pH3to7)
dim(df.WvW.O.pH4to8)
dim(df.WvW.O.pH5to9)
dim(df.WvW.A.pH3to7)
dim(df.WvW.A.pH4to8)
dim(df.WvW.A.pH5to9)

# Let's bring back in our taxonomy from the tax table
SigOTUs = levels(as.factor(df.WvW.pH.groups$OTU))
pruned = prune_taxa(SigOTUs,ps.norm.prune)
taxtab = data.frame(tax_table(pruned))
taxtab$OTU = c(taxa_names(pruned))
df.joined.WvW.pH = merge(df.WvW.pH.groups,taxtab,by=c("OTU"))
head(df.joined.WvW.pH)


# EXPERIMENT 2 - pH categories, control -----

ps.CvC <- prune_samples(sample_data(ps.raw.prune)$burn.trtmt == 'control', ps.raw.prune)

### pH 3 to 7, control burns, O horizon -----
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

df.CvC.O.pH3to7 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.CvC.O.pH3to7$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.CvC.O.pH3to7$significant_models[[j]]
  # Pull out the coefficients as above
  # NOTE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_control = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.CvC.O.pH3to7$p_fdr[dt.CvC.O.pH3to7$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_control$p_fdr = p_fdr
  # Create a column with the OTU ID
  mu_control$OTU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.CvC.O.pH3to7 = rbind(df.CvC.O.pH3to7,mu_control)
}

df.CvC.O.pH3to7$Comparison = 'pb.control.v.SI.control'
df.CvC.O.pH3to7$Experiment = 'fast growth'
df.CvC.O.pH3to7$Temp.Cutoff = 'null'
df.CvC.O.pH3to7$pH.cutoff = '3 to 7'
df.CvC.O.pH3to7$horizon = 'O'
df.CvC.O.pH3to7$DNA.type = 'gDNA'
df.CvC.O.pH3to7$Controlling.for = 'pH'

colnames(df.CvC.O.pH3to7) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                              'Experiment', 'Temp.cutoff','pH.cutoff','horizon', 
                              'DNA.type', 'Controlling.for')

### pH 4 to 8, control burns, O horizon -----

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

df.CvC.O.pH4to8 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.CvC.O.pH4to8$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.CvC.O.pH4to8$significant_models[[j]]
  # Pull out the coefficients as above
  # NOTE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_control = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.CvC.O.pH4to8$p_fdr[dt.CvC.O.pH4to8$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_control$p_fdr = p_fdr
  # Create a column with the OTU ID
  mu_control$OTU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.CvC.O.pH4to8 = rbind(df.CvC.O.pH4to8,mu_control)
}

df.CvC.O.pH4to8$Comparison = 'pb.control.v.SI.control'
df.CvC.O.pH4to8$Experiment = 'fast growth'
df.CvC.O.pH4to8$Temp.Cutoff = 'null'
df.CvC.O.pH4to8$pH.cutoff = '4 to 8'
df.CvC.O.pH4to8$horizon = 'O'
df.CvC.O.pH4to8$DNA.type = 'gDNA'
df.CvC.O.pH4to8$Controlling.for = 'pH'

colnames(df.CvC.O.pH4to8) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                              'Experiment', 'Temp.cutoff','pH.cutoff','horizon', 
                              'DNA.type', 'Controlling.for')

### pH 5 to 9, control burns, O horizon -----
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

df.CvC.O.pH5to9 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.CvC.O.pH5to9$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.CvC.O.pH5to9$significant_models[[j]]
  # Pull out the coefficients as above
  # NOTE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_control = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.CvC.O.pH5to9$p_fdr[dt.CvC.O.pH5to9$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_control$p_fdr = p_fdr
  # Create a column with the OTU ID
  mu_control$OTU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.CvC.O.pH5to9 = rbind(df.CvC.O.pH5to9,mu_control)
}

df.CvC.O.pH5to9$Comparison = 'pb.control.v.SI.control'
df.CvC.O.pH5to9$Experiment = 'fast growth'
df.CvC.O.pH5to9$Temp.Cutoff = 'null'
df.CvC.O.pH5to9$pH.cutoff = '5 to 9'
df.CvC.O.pH5to9$horizon = 'O'
df.CvC.O.pH5to9$DNA.type = 'gDNA'
df.CvC.O.pH5to9$Controlling.for = 'pH'

colnames(df.CvC.O.pH5to9) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                              'Experiment', 'Temp.cutoff','pH.cutoff','horizon', 
                              'DNA.type', 'Controlling.for')

### pH 3 to 7, control burns, A horizon -----
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

df.CvC.A.pH3to7 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.CvC.A.pH3to7$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.CvC.A.pH3to7$significant_models[[j]]
  # Pull out the coefficients as above
  # NATE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_control = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.CvC.A.pH3to7$p_fdr[dt.CvC.A.pH3to7$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_control$p_fdr = p_fdr
  # Create a column with the ATU ID
  mu_control$ATU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.CvC.A.pH3to7 = rbind(df.CvC.A.pH3to7,mu_control)
}

df.CvC.A.pH3to7$Comparison = 'pb.control.v.SI.control'
df.CvC.A.pH3to7$Experiment = 'fast growth'
df.CvC.A.pH3to7$Temp.Cutoff = 'null'
df.CvC.A.pH3to7$pH.cutoff = '3 to 7'
df.CvC.A.pH3to7$horizon = 'A'
df.CvC.A.pH3to7$DNA.type = 'gDNA'
df.CvC.A.pH3to7$Controlling.for = 'pH'

colnames(df.CvC.A.pH3to7) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                              'Experiment', 'Temp.cutoff','pH.cutoff','horizon', 
                              'DNA.type', 'Controlling.for')

### pH 4 to 8, control burns, A horizon -----

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

df.CvC.A.pH4to8 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.CvC.A.pH4to8$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.CvC.A.pH4to8$significant_models[[j]]
  # Pull out the coefficients as above
  # NATE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_control = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.CvC.A.pH4to8$p_fdr[dt.CvC.A.pH4to8$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_control$p_fdr = p_fdr
  # Create a column with the ATU ID
  mu_control$ATU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.CvC.A.pH4to8 = rbind(df.CvC.A.pH4to8,mu_control)
}

df.CvC.A.pH4to8$Comparison = 'pb.control.v.SI.control'
df.CvC.A.pH4to8$Experiment = 'fast growth'
df.CvC.A.pH4to8$Temp.Cutoff = 'null'
df.CvC.A.pH4to8$pH.cutoff = '4 to 8'
df.CvC.A.pH4to8$horizon = 'A'
df.CvC.A.pH4to8$DNA.type = 'gDNA'
df.CvC.A.pH4to8$Controlling.for = 'pH'

colnames(df.CvC.A.pH4to8) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                              'Experiment', 'Temp.cutoff','pH.cutoff','horizon', 
                              'DNA.type', 'Controlling.for')

### pH 5 to 9, control burns, A horizon -----
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

df.CvC.A.pH5to9 = data.frame()

# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.CvC.A.pH5to9$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.CvC.A.pH5to9$significant_models[[j]]
  # Pull out the coefficients as above
  # NATE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_control = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.CvC.A.pH5to9$p_fdr[dt.CvC.A.pH5to9$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_control$p_fdr = p_fdr
  # Create a column with the ATU ID
  mu_control$ATU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.CvC.A.pH5to9 = rbind(df.CvC.A.pH5to9,mu_control)
}

df.CvC.A.pH5to9$Comparison = 'pb.control.v.SI.control'
df.CvC.A.pH5to9$Experiment = 'fast growth'
df.CvC.A.pH5to9$Temp.Cutoff = 'null'
df.CvC.A.pH5to9$pH.cutoff = '5 to 9'
df.CvC.A.pH5to9$horizon = 'A'
df.CvC.A.pH5to9$DNA.type = 'gDNA'
df.CvC.A.pH5to9$Controlling.for = 'pH'

colnames(df.CvC.A.pH5to9) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                              'Experiment', 'Temp.cutoff','pH.cutoff','horizon', 
                              'DNA.type', 'Controlling.for')


### Combine pH category output -----
ps.corncob.CvC.O.pH3to7
ps.corncob.CvC.O.pH4to8
ps.corncob.CvC.O.pH5to9
ps.corncob.CvC.A.pH3to7
ps.corncob.CvC.A.pH4to8
ps.corncob.CvC.A.pH5to9

length(dt.CvC.O.pH3to7$significant_taxa)
length(dt.CvC.O.pH4to8$significant_taxa)
length(dt.CvC.O.pH5to9$significant_taxa)
length(dt.CvC.A.pH3to7$significant_taxa)
length(dt.CvC.A.pH4to8$significant_taxa)
length(dt.CvC.A.pH5to9$significant_taxa)



df.CvC.pH.groups <- rbind(df.CvC.O.pH3to7,
                          df.CvC.O.pH4to8,
                          df.CvC.O.pH5to9,
                          df.CvC.A.pH3to7,
                          df.CvC.A.pH4to8,
                          df.CvC.A.pH5to9)

# Let's bring back in our taxonomy from the tax table
SigOTUs = levels(as.factor(df.CvC.pH.groups$OTU))
pruned = prune_taxa(SigOTUs,ps.norm.prune)
taxtab = data.frame(tax_table(pruned))
taxtab$OTU = c(taxa_names(pruned))
df.joined.CvC.pH = merge(df.CvC.pH.groups,taxtab,by=c("OTU"))
head(df.joined.CvC.pH)




df.DvD.pH.groups
df.WvW.pH.groups
df.CvC.pH.groups

# Save output: 
write.csv(df.joined.DvD.pH,"data/sequence-data/LibCombined/corncob-output/OTU level/All_taxa/Ex2 Fast growth/DvD-3-sliding-pH-categories.csv", row.names = FALSE)
#    

# Save output: 
write.csv(df.joined.CvC.pH,"data/sequence-data/LibCombined/corncob-output/OTU level/All_taxa/Ex2 Fast growth/CvC-3-sliding-pH-categories.csv", row.names = FALSE)
#    

# Save output: 
write.csv(df.joined.WvW.pH,"data/sequence-data/LibCombined/corncob-output/OTU level/All_taxa/Ex2 Fast growth/WvW-3-sliding-pH-categories.csv", row.names = FALSE)
#  


