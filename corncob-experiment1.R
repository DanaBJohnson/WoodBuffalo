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


ps.raw.full <- prune_samples(sample_data(ps.raw.full)$Full.id != '19UW-WB-06-08-A-SI-duplicate' &
                               sample_data(ps.raw.full)$Full.id != '19UW-WB-07-02-O-SI-duplicate' &
                               sample_data(ps.raw.full)$Full.id != '19UW-WB-11-03-O-SI-duplicate'&
                               sample_data(ps.raw.full)$Full.id != '19UW-WB-19-07-A-SI-duplicate' &
                               sample_data(ps.raw.full)$Full.id != '19UW-WB-08-10-O-SI', ps.raw.full)
ps.raw.full <- prune_taxa(taxa_sums(ps.raw.full)>0, ps.raw.full)

# normalize abundance seq data
ps.norm.full <- prune_samples(sample_data(ps.norm.full)$Full.id != '19UW-WB-06-08-A-SI-duplicate' &
                                sample_data(ps.norm.full)$Full.id != '19UW-WB-07-02-O-SI-duplicate' &
                                sample_data(ps.norm.full)$Full.id != '19UW-WB-11-03-O-SI-duplicate'&
                                sample_data(ps.norm.full)$Full.id != '19UW-WB-19-07-A-SI-duplicate' &
                                sample_data(ps.norm.full)$Full.id != '19UW-WB-08-10-O-SI', ps.norm.full)
ps.norm.full <- prune_taxa(taxa_sums(ps.norm.full)>0, ps.norm.full)

# Exclude all but 24 hr post-burn
# Create separate genomic DNA and an RNA phyloseq objects
ps.raw.RNA <- ps.raw.full %>%
  phyloseq::subset_samples(incub.trtmt %in% c('pb')) %>%
  subset_samples(DNA.type %in% c('cDNA')) 

ps.raw.DNA <- ps.raw.full %>%
  phyloseq::subset_samples(incub.trtmt %in% c('pb')) %>%
  subset_samples(DNA.type %in% c('gDNA')) 



### Step 2. Subset for the most abundance taxa for use in testing -----

# Create list of taxa with a mean abundance greater than cut-off
MeanAbunTaxa = taxa_names(filter_taxa(ps.norm.full, function(x) mean(x) > 0.0001, TRUE))
length(MeanAbunTaxa)

# Create list of taxa with max abundance greater than cutoff
MaxAbunTaxa <-  taxa_names(filter_taxa(ps.norm.full, function(x) max(x) > 0.005, TRUE))
length(MaxAbunTaxa)

ntaxa(ps.norm.full)

x <- c(MeanAbunTaxa, MaxAbunTaxa)
length(unique(x))

# Prune taxa to include only OTUs with a maximum rel. abun > Cutoff
ps.prune.raw.DNA = prune_taxa(MaxAbunTaxa, ps.raw.DNA)



### Step 3. Set temperature cutoff for dry soil burns -----
# Focus on dry soil burns: 
ps.DNA.dry <-  prune_samples(sample_data(ps.prune.raw.DNA)$burn.trtmt != 'wet', ps.prune.raw.DNA)

# Create dataframe ps.DNA.dry sample_data.
df.dry <- data.frame(sample_data(ps.DNA.dry))

# Assign maximum temp at bottom thermocouple as the Max temp. for the core:
df.dry <- df.dry %>%
  group_by(site) %>%
  mutate(Dry.burn.max = max(Thermo.low.max)) 

# Now combine the new sample data with the original phyloseq object:
ps.df <- sample_data(df.dry)

sample_names(ps.df) = df.dry$Full.id

sample_data(ps.DNA.dry) <- ps.df




### Step 4. Run corncob on O and A horizons -----

# Create O horizon phyloseq object
ps.corncob.dry.O.gt50 <-  prune_samples(sample_data(ps.DNA.dry)$horizon == 'O' & 
                                          sample_data(ps.DNA.dry)$Dry.burn.max>50, ps.DNA.dry)
# Clean up ps to reduce size
ps.corncob.dry.O.gt50 <- prune_taxa(taxa_sums(ps.corncob.dry.O.gt50)>0, ps.corncob.dry.O.gt50)
ps.corncob.dry.O.gt50

# How many samples make it through?
sample_names(ps.corncob.dry.O.gt50)
# Results in 25 samples (some of the dry soils burns don't make the 50C cutoff)

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
ps.corncob.dry.A.gt50

# How many samples make it through?
sample_names(ps.corncob.dry.A.gt50)
# Results in 17 samples (not all cores have A horizons AND some of the dry soils 
#     burns don't make the 50C cutoff)

# Run corncob controlling for pre-burn pH (cntl.pH)
dt.DNA.DryBurn.A.gt50 <- differentialTest(formula = ~burn.trtmt+cntl.pH,
                                          phi.formula = ~ burn.trtmt+cntl.pH,
                                          formula_null = ~ cntl.pH,
                                          phi.formula_null = ~burn.trtmt+ cntl.pH,
                                          test = 'Wald',boot=FALSE,
                                          data=ps.corncob.dry.A.gt50,
                                          fdr_cutoff = 0.05)




### Step 5. Clean up corncob output -----
# Create empty dataframes to fill
df.Dry.O.gt50 = data.frame()
df.Dry.A.gt50 = data.frame()


# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.DNA.DryBurn.O.gt50$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.DNA.DryBurn.O.gt50$significant_models[[j]]
  # Pull out the coefficients as above
  # NOTE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_dry = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.DNA.DryBurn.O.gt50$p_fdr[dt.DNA.DryBurn.O.gt50$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_dry$p_fdr = p_fdr
  # Create a column with the OTU ID
  mu_dry$OTU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.Dry.O.gt50 = rbind(df.Dry.O.gt50,mu_dry)
}

# Fill in relavent information for later analyses:
df.Dry.O.gt50$Comparison = 'pb.dry.v.pb.control'
df.Dry.O.gt50$Experiment = 'survival'
df.Dry.O.gt50$Temp.Cutoff = 'Low thermo greater than 50 C'
df.Dry.O.gt50$pH.cutoff = 'null'
df.Dry.O.gt50$horizon = 'O'
df.Dry.O.gt50$DNA.type = 'gDNA'
df.Dry.O.gt50$Controlling.for = 'pre-burn pH'

colnames(df.Dry.O.gt50) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                            'Experiment', 'Temp.cutoff','pH.cutoff','horizon', 
                            'DNA.type', 'Controlling.for')


# REPEAT FOR A HORIZON:
# Loop to pull out coefficients for each taxon
for (j in 1:length(dt.DNA.DryBurn.A.gt50$significant_taxa)){
  # Get the significant model for that taxon
  sig_models = dt.DNA.DryBurn.A.gt50$significant_models[[j]]
  # Pull out the coefficients as above
  # NOTE: you have to make sure the row reference is correct, as determined above
  # Here we have it as 7.
  mu_dry = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  # Also grab the p_fdr estimate for that taxon's model
  p_fdr = dt.DNA.DryBurn.A.gt50$p_fdr[dt.DNA.DryBurn.A.gt50$significant_taxa][j]
  # Add that estimate onto our coefficient data frame
  mu_dry$p_fdr = p_fdr
  # Create a column with the OTU ID
  mu_dry$OTU= paste(row.names(data.frame(p_fdr)))
  # Add this row onto the df dataframe, which will collect the results
  # for all taxa as it iterates through this loop.
  df.Dry.A.gt50 = rbind(df.Dry.A.gt50,mu_dry)
}

# Fill in relavent information for later analyses:
df.Dry.A.gt50$Comparison = 'pb.dry.v.pb.control'
df.Dry.A.gt50$Experiment = 'survival'
df.Dry.A.gt50$Temp.Cutoff = 'Low thermo greater than 50 C'
df.Dry.A.gt50$pH.cutoff = 'null'
df.Dry.A.gt50$horizon = 'A'
df.Dry.A.gt50$DNA.type = 'gDNA'
df.Dry.A.gt50$Controlling.for = 'pre-burn pH'

colnames(df.Dry.A.gt50) = c("Estimate","SE","t","p","p_fdr","OTU", 'Comparison', 
                            'Experiment', 'Temp.cutoff','pH.cutoff','horizon', 
                            'DNA.type', 'Controlling.for')


df.joined <- rbind(df.Dry.O.gt50, df.Dry.A.gt50)

# Let's bring back in our taxonomy from the tax table
SigOTUs = levels(as.factor(df.joined$OTU))
pruned = prune_taxa(SigOTUs,ps.norm.full)
taxtab = data.frame(tax_table(pruned))
taxtab$OTU = c(taxa_names(pruned))
df.joined = merge(df.joined,taxtab,by=c("OTU"))
head(df.joined)



### Step 6. Pull out "Survivors" with an RNA:DNA ratio > [cutoff] -----

# Subset ps.norm by list of survivor OTUs
ps.norm.prune <- prune_taxa(taxa_names(ps.norm.full) %in% df.joined$OTU, ps.norm.full)

# Create separate dataframes for RNA  (cDNA) and genomic DNA 
df.norm.prune <- psmelt(ps.norm.prune)

df.DNA.prune <- df.norm.prune %>%
  subset(DNA.type == 'gDNA') %>%
  subset(select = c(OTU, core.id.hor.incub, Abundance))

df.RNA.prune <- df.norm.prune %>%
  subset(DNA.type == 'cDNA') %>%
  subset(select = c(OTU, core.id.hor.incub, Abundance))

# Change abundance column heading 
colnames(df.DNA.prune)[3] = 'Abun.DNA'
colnames(df.RNA.prune)[3] = 'Abun.RNA'

# Merge dataframes
df.merge <- merge(df.DNA.prune, df.RNA.prune, by = c('OTU','core.id.hor.incub'))

# Calculate RNA:DNA ratio
df.ratio <- df.merge %>%
  mutate(RNA.DNA.ratio = Abun.RNA/Abun.DNA) %>%
  # SUbset survivor list to include only OTUs with ratio > [cutoff]
  subset(RNA.DNA.ratio > 1)


# Prune survivors to include only OTUs meeting RNA:DNA ratio [cutoff]
length(unique(df.ratio$OTU))

df.joined.prune <- df.joined %>%
  subset(OTU %in% df.ratio$OTU)



### Step 7. Save results -----

write.csv(df.joined.prune,"data/sequence-data/LibCombined/corncob-output/OTU level/All_taxa/Ex1 Fire survival/gDNA-Temp-cutoff-50C.csv", row.names = FALSE)

otu_to_taxonomy(OTU = dt.RNA.wetBurn.A.lt50$significant_taxa, data = ps.RNA.wet)
#plot(dt.RNA.wetBurn.A.lt50, level = c('Phylum'))

