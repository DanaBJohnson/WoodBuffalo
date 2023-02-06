# Today's goal (21.Oct.2021) - for each experiment, create a dataframe with the
#    positive responders and a composite mu value. I can calculate this mu value
#    in different ways - max, min, mean - and create a weighted mu value.
#    So the resulting dataframe will have a mu.survival.max (or mu.mean or mu.min) and a 
#    mu.survival.weighted

library(ggplot2)
library(dplyr)
library(phyloseq)
library(tidyr)
library(phyloseq)

#setwd('../Box/WhitmanLab/Projects/WoodBuffalo/FireSim2019')


### STEP 1. Import list of enriched/depleted taxa generated from each experiment -----

# Experiment 1
df.exp1 <- read.csv('data/sequence-data/LibCombined/corncob-output/Manuscript/Ex1_signOTUs_gDNA.csv')
# Fix colnames for merging later:
colnames(df.exp1)[6] <- 'comparison'
colnames(df.exp1)[7] <- 'experiment'


# Experiment 2
df.exp2 <- read.csv('data/sequence-data/LibCombined/corncob-output/Manuscript/Ex2_signOTUs_3sliding_pH_categories_byBurnTrtmt.csv')

# Experiment 3
df.exp3 <- read.csv('data/sequence-data/LibCombined/corncob-output/Manuscript/Ex3_signOTUs_MidThermoCutoff_50C.csv')


# Combined into full dataset
colnames(df.exp1)
colnames(df.exp2)
colnames(df.exp3)

df.exp13 <- rbind(df.exp1,df.exp3) %>%
  mutate(mu.incub.trtmt = NA,
         mu.pH = NA)

df.exp2 <- df.exp2 %>%
  mutate(mu.burn.trtmt = NA,
         mu.cntl.pH = NA)

df.raw.combined <- rbind(df.exp13,df.exp2)
write.csv(df.raw.combined, 'data/sequence-data/LibCombined/corncob-output/Manuscript/All-responders.csv')





### STEP 2: For experiment 1, filter for "Survivors" with an RNA:DNA ratio > [cutoff] -----
# Import the full phyloseq object, filter for samples 24 hrs post-burn, with OTU
#   relative abundance > 0.
ps.norm.full <- readRDS('data/sequence-data/LibCombined/phyloseq-objects/ps.norm.full') 
ps.norm.prune <- prune_samples(sample_data(ps.norm.full)$incub.trtmt == 'pb', ps.norm.full)
ps.norm.prune <- prune_taxa(taxa_sums(ps.norm.prune)>0, ps.norm.prune)

# Create subset of ps.norm that only includes OTUs identified as differentially
#   enriched/depleted. This will include both positive and negative responders.
ps.responsive <- prune_taxa(taxa_names(ps.norm.prune) %in% df.exp1$OTU, ps.norm.prune)

# Create separate dataframes for RNA  (cDNA) and genomic DNA 
df.responsive <- psmelt(ps.responsive)

df.DNA.prune <- df.responsive %>%
  subset(DNA.type == 'gDNA') %>%
  subset(select = c(OTU, core.id.hor.incub, Abundance))

df.RNA.prune <- df.responsive %>%
  subset(DNA.type == 'cDNA') %>%
  subset(select = c(OTU, core.id.hor.incub, Abundance))

# Change abundance column heading 
colnames(df.DNA.prune)[3] = 'Abun.DNA'
colnames(df.RNA.prune)[3] = 'Abun.RNA'

# Merge dataframes and calculate RNA:DNA ratios
df.ratio <- merge(df.DNA.prune, df.RNA.prune, by = c('OTU','core.id.hor.incub')) %>%
  subset(Abun.DNA > 0 & Abun.RNA > 0) %>%
  # Calculate RNA:DNA for OTUs within each sample
  mutate(RNA.DNA.ratio = Abun.RNA/Abun.DNA) %>%
  # But then group by OTU and calculate the mean ratio for each OTU
  group_by(OTU) %>%
  mutate(max.RNA.DNA.ratio = max(RNA.DNA.ratio)) %>%
  subset(select = c('OTU','max.RNA.DNA.ratio')) %>%
  unique()


df.survivors <- merge(df.exp1, df.ratio, by = "OTU") %>%
  # Pull out positive responders
  subset(mu.burn.trtmt > 0) %>%
  #subset(mu.burn.trtmt > mean(mu.burn.trtmt)-1*sd(mu.burn.trtmt)) %>%
  # Remove taxa with an RNA:DNA ratio before [cutoff]
  subset(max.RNA.DNA.ratio > 0) %>%
  # Calculate max mu for each OTU
  group_by(OTU, horizon) %>%
  mutate(mu.survival.max = max(mu.burn.trtmt),
         mu.fast.max = NA,
         mu.affinity.max = NA) %>%
  subset(select = c(OTU, comparison, experiment, horizon,
                    Domain, Phylum, Class, Order, Family, Genus, Species, 
                    mu.survival.max, mu.fast.max, mu.affinity.max,max.RNA.DNA.ratio)) %>% 
  unique()

# Check for OTUs that are both positive and negative responders
# Check how number of responders changes with mu > mean(mu)-1*sd(mu)


# Number of fire surviving taxa
dim(df.survivors)
# Number of unique fire surviving taxa
length(unique(df.survivors)$OTU)




### Step 3. -----
# Calculate max mu value for each OTU
df.fast <- df.exp2 %>%
  # Pull out positive responders
  subset(mu.incub.trtmt > 0) %>%
  #subset(comparison == 'pb.dry.vs.SI.dry') %>%
  #subset(mu.incub.trtmt > mean(mu.incub.trtmt)-1*sd(mu.incub.trtmt)) %>%
  group_by(OTU,horizon) %>%
  mutate(mu.fast.max = max(mu.incub.trtmt),
         mu.survival.max = NA,
         mu.affinity.max = NA,
         max.RNA.DNA.ratio = NA) %>%
  subset(select = c(OTU, comparison, experiment, horizon,
                    Domain, Phylum, Class, Order, Family, Genus, Species, 
                    mu.survival.max, mu.fast.max, mu.affinity.max,max.RNA.DNA.ratio)) %>% 
  unique()
dim(df.fast)

# repeat for third experiment
df.affinity <- df.exp3 %>%
  # Pull out positive responders
  subset(mu.burn.trtmt > 0) %>%
  #subset(mu.burn.trtmt > mean(mu.burn.trtmt)-1*sd(mu.burn.trtmt)) %>%
  group_by(OTU, horizon) %>%
  mutate(mu.affinity.max = max(mu.burn.trtmt),
         mu.fast.max = NA,
         mu.survival.max = NA,
         max.RNA.DNA.ratio = NA) %>%
  subset(select = c(OTU, comparison, experiment, horizon,
                    Domain, Phylum, Class, Order, Family, Genus, Species, 
                    mu.survival.max, mu.fast.max, mu.affinity.max, max.RNA.DNA.ratio)) %>% 
  unique()
dim(df.affinity)

# Check that colnames match (and fix if they don't)
df.combined <- rbind(df.survivors, df.fast, df.affinity)



### Step 4. Save output?
#write.csv(df.combined, 'data/sequence-data/LibCombined/corncob-output/Manuscript/Trait-responder-OTUs.csv')























# #setwd('../Box/WhitmanLab/Projects/WoodBuffalo/FireSim2019/data')
# # Experiment 1: -----
# #     Import data
# 
# df.survival <- read.csv('data/sequence-data/LibCombined/corncob-output/Ex1 Fire survival/Ex1_signOTUs_gDNA.csv')
# 
# # OLD CODE
# # df.1.low <- read.csv('sequence-data/LibCombined/corncob-output/Ex1 Fire survival/low-thermo-cutoff-50C.csv')
# # df.1.mid <- read.csv('sequence-data/LibCombined/corncob-output/Ex1 Fire survival/Mid-thermo-cutoff.csv')
# # 
# # df.1 <- rbind(df.1.low,df.1.mid)
# # 
# 
# head(df.survival,2)
# dim(df.survival) # total of 532 responders (compared to 951 with RNA data)
# length(unique(df.survival$OTU)) # 379 unique OTU responders (compared to 388 with RNA data)
# 
# # Create dataframe with unique OTU identifies
# y <- data.frame(unique(df.survival$OTU))
# colnames(y)[1] = 'OTU'
# y$Responder.in.n.analyses = 0
# y$Pos.count = 0
# y$Neg.count = 0
# 
# # Calculate the number of times that each OTU was identified as a pos/neg fire responder
# for (i in 1:nrow(y)) {
#   for (j in 1:nrow(df.survival)) {
#     if (df.survival$OTU[j] == y$OTU[i]) {
#       # Add 1 every time a specific OTU appears
#       y$Responder.in.n.analyses[i] = y$Responder.in.n.analyses[i]+1
#     } else {
#       y$Responder.in.n.analyses[i] = y$Responder.in.n.analyses[i]
#     }
#   }
# }
# 
# for (i in 1:nrow(y)) {
#   for(j in 1:nrow(df.survival)) {
#     # If the intercept > 0, then this is INCREASED diff. abundance
#     if (df.survival$OTU[j] == y$OTU[i] & df.survival$intercept[j] > 0) {
#       # So add one to positive count
#       y$Pos.count[i] = y$Pos.count[i] + 1
#       # If intercept < 0, this is DECREASED diff. abun
#     } else if (df.survival$OTU[j] == y$OTU[i] & df.survival$intercept[j] < 0) {
#       # So add one to negative count
#       y$Neg.count[i] = y$Neg.count[i] + 1
#     } else {
#       y$Pos.count[i] = y$Pos.count[i]
#       y$Neg.count[i] = y$Neg.count[i]
#     }
#   }
# }
# 
# # Combine this data with the full corncob output:
# z <- merge(df.survival, y)
# hist(z$Responder.in.n.analyses)
# 
# # Create a list of the OTUs that were identified as both pos. and neg. responders
# Problem.OTUs <- subset(y, Pos.count > 0 & Neg.count > 0) 
# 
# df.problem <- subset(df.survival, OTU %in% Problem.OTUs$OTU) %>%
#   merge(Problem.OTUs)
# 
# 
# # Use only the dry v. control data, and the max temperature cutoff.
# df.fire.survivors.dry <- df.survival %>%
#   # Pull out only DRY and CONTROL samples, 24 HRS POST-BURN
#   subset(comparison == 'pb.dry.vs.pb.control') %>%
#   # Pull out only comparisons run at temperatures greater than cutoff
#   subset(T.cutoff %in% c('Low thermo > 50C')) %>%
#   # Only use positive responders
#   subset(intercept > 0) %>%
#   group_by(OTU, horizon) %>%
#   mutate(mu.survival.max = max(intercept),
#          mu.affinity.max = 'null',
#          mu.fast.max = 'null') %>%
#   subset(select = c(OTU, comparison, DNA.type, experiment, horizon, 
#                     Domain, Phylum, Class, Order, Family, Genus, Species, 
#                     mu.fast.max, mu.survival.max, mu.affinity.max)) %>%
#   # Make sure there are no duplicate rows
#   unique() %>%
#   # Remove dual responders (pos. & neg.)
#   subset(!OTU %in% Problem.OTUs$OTU)
# 
# 
# # Repeat for negative responders
# df.fire.susceptible.dry <- df.survival %>%
#   # Pull out only DRY and CONTROL samples, 24 HRS POST-BURN
#   subset(comparison == 'pb.dry.vs.pb.control') %>%
#   # Pull out only comparisons run at temperatures greater than cutoff
#   subset(T.cutoff %in% c('Low thermo > 50C')) %>%
#   # Only use positive responders
#   subset(intercept < 0) %>%
#   group_by(OTU, horizon) %>%
#   mutate(mu.survival.min = min(intercept)) %>%
#   subset(select = c(OTU, comparison, DNA.type, experiment, horizon, 
#                     Domain, Phylum, Class, Order, Family, Genus, Species, 
#                     mu.survival.min)) %>%
#   # Make sure there are no duplicate rows - there will be duplicate OTUs when 
#   #   an OTU was identified as a responder in both O and A horizons
#   unique() %>%
#   # Remove dual responders (pos. & neg.)
#   subset(!OTU %in% Problem.OTUs$OTU)
# 
# dim(df.fire.survivors.dry)
# dim(df.fire.susceptible.dry)
# 
# 
# 
# # Also keep a list of all the survivors:
# df.fire.survivors <- df.survival %>%
#   subset(intercept>0) %>%
#   group_by(OTU, horizon) %>%
#   mutate(mu.survival.max = max(intercept)) %>%
#   subset(select = c(OTU, comparison, experiment, horizon,
#                     Domain, Phylum, Class, Order, Family, Genus, Species, 
#                     mu.survival.max)) %>%
#   unique() %>%
#   subset(!OTU %in% Problem.OTUs$OTU)
# 
# df.fire.susceptible <- df.survival %>%
#   subset(intercept<0) %>%
#   group_by(OTU, horizon) %>%
#   mutate(mu.survival.min = min(intercept)) %>%
#   subset(select = c(OTU, comparison, experiment, horizon,
#                     Domain, Phylum, Class, Order, Family, Genus, Species, 
#                     mu.survival.min)) %>%
#   unique() %>%
#   subset(!OTU %in% Problem.OTUs$OTU)
# 
# 
# 
# ggplot(df.fire.survivors.dry, aes(x=mu.survival.max, y = Phylum, color = Phylum), alpha = 0.8) + 
#   geom_jitter() + 
#   facet_grid(horizon~comparison) + 
#   geom_vline(xintercept = 0) + 
#   theme_bw()
# 
# 
# # write.csv(df.fire.survivors, 'data/sequence-data/LibCombined/corncob-output/Manuscript/fire-survivors-all-burns.csv', row.names = FALSE)
# # write.csv(df.fire.susceptible, 'data/sequence-data/LibCombined/corncob-output/Manuscript/fire-susceptible-all-burns.csv', row.names = FALSE)
# # write.csv(df.fire.survivors.dry, 'data/sequence-data/LibCombined/corncob-output/Manuscript/fire-survivors-dry-burns.csv', row.names = FALSE)
# # write.csv(df.fire.susceptible.dry, 'data/sequence-data/LibCombined/corncob-output/Manuscript/fire-susceptible-dry-burns.csv', row.names = FALSE)
# 
# 
# 
# # Experiment 2: -----
# #    Import data.
# 
# df.fast.pH <- read.csv('data/sequence-data/LibCombined/corncob-output/Ex2 Fast growth/Ex2_signOTUs_3sliding_pH_categories.csv')
# df.fast.temp<-read.csv('data/sequence-data/LibCombined/corncob-output/Ex2 Fast growth/Ex2_signOTUs_LowThermoCutoff_50C.csv')
# df.fast <- rbind(df.fast.pH, df.fast.temp)
# 
# # Old code
# # df.2.pH <- read.csv('sequence-data/LibCombined/corncob-output/Ex2 Fast growth/3-sliding-pH-categories.csv')
# # df.2.low <- read.csv('sequence-data/LibCombined/corncob-output/Ex2 Fast growth/low-thermo-cutoff-50C.csv')
# # df.2.mid <- read.csv('sequence-data/LibCombined/corncob-output/Ex2 Fast growth/Mid-thermo-cutoff-50.csv')
# # df.2.mid100<- read.csv('sequence-data/LibCombined/corncob-output/Ex2 Fast growth/Mid-thermo-cutoff-100.csv')
# # df.fast <- rbind(df.fast.pH, df.fast.low, df.fast.mid, df.fast.mid100)
# 
# dim(df.fast) # 2462 OTUs
# length(unique(df.fast$OTU)) # 948 unique OTUs
# 
# y <- data.frame(unique(df.fast$OTU))
# 
# # Create dataframe to hold OTU counts
# head(df.fast)
# colnames(y)[1] = 'OTU'
# y$Responder.in.n.analyses = 0
# y$Pos.count = 0
# y$Neg.count = 0
# 
# # Count the number of times each OTU appears in dataframe
# for (i in 1:nrow(y)) {
#   for (j in 1:nrow(df.fast)) {
#     if (df.fast$OTU[j] == y$OTU[i]) {
#       y$Responder.in.n.analyses[i] = y$Responder.in.n.analyses[i]+1
#     } else {
#       y$Responder.in.n.analyses[i] = y$Responder.in.n.analyses[i]
#     }
#   }
# }
# 
# # Count the number of times each OTU is identified as pos/neg responder
# for (i in 1:nrow(y)) {
#   for(j in 1:nrow(df.fast)) {
#     if (df.fast$OTU[j] == y$OTU[i] & df.fast$mu.incub.trtmt[j] > 0) {
#       y$Pos.count[i] = y$Pos.count[i] + 1
#     } else if (df.fast$OTU[j] == y$OTU[i] & df.fast$mu.incub.trtmt[j] < 0) {
#       y$Neg.count[i] = y$Neg.count[i] + 1
#     } else {
#       y$Pos.count[i] = y$Pos.count[i]
#       y$Neg.count[i] = y$Neg.count[i]
#     }
#   }
# }
# 
# # Combine dataframe and count data
# z <- merge(df.fast, y)
# hist(z$Responder.in.n.analyses)
# 
# # Make a list of dual responders (ID as both pos. and neg. response)
# Dual.Responders <- subset(y, Pos.count > 0 & Neg.count > 0) 
# dim(Dual.Responders)
# 
# df.dual.responders <- subset(df.fast, OTU %in% Dual.Responders$OTU) %>%
#   merge(Dual.Responders)
# 
# # PLOT 
# ggplot(df.dual.responders, aes(x=intercept, y=Phylum, color = Phylum)) +
#   geom_point(size=2) + 
#   theme_bw()+
#   facet_grid(horizon~comparison) + 
#   geom_vline(xintercept = 0)
# 
# length(unique(df.dual.responders$OTU)) # 148 dual responders
# 
# 
# # Create final dataset of fast/slow growing OTUs:
# # Subset dataframe based on treatment used:
# df.dry <- subset(df.fast, comparison == 'pb.dry.vs.SI.dry') # Dry burns
# df.control <- subset(df.fast, comparison == 'pb.cntl.vs.SI.cntl') # No burning
# df.wet <- subset(df.fast, comparison == 'pb.wet.vs.SI.wet') # Wet (moist) burns
# 
# # Create dataframe of fast growing OTUs (pos. responders) after dry burns
# df.fast.growers.after.fire <- df.dry %>%
#   # Intercept > 0 indicates increased diff. abundance
#   subset(intercept > 0) %>%
#   # Exclude OTUs identified in comparisons of wet and control cores
#   subset(!OTU %in% subset(df.wet, intercept > 1)$OTU) %>%
#   subset(!OTU %in% subset(df.control, intercept > 1)$OTU) %>%
#   group_by(OTU, horizon) %>%
#   # Take the max intercept value for each OTU
#   mutate(mu.fast.max = max(intercept),
#          mu.survival.max = 'null',
#          mu.affinity.max = 'null') %>%
#   subset(select = c(OTU, comparison, DNA.type, experiment, horizon, 
#                     Domain, Phylum, Class, Order, Family, Genus, Species, 
#                     mu.fast.max, mu.survival.max, mu.affinity.max)) %>%
#   unique() 
# 
# length(unique(df.fast.growers.after.fire$OTU)) # 172 taxa
# 
# # Create dataframe with all fast growers (from all burn treatments)
# df.all.fast.growers <- subset(df.fast, intercept > 0) %>%
#   group_by(OTU, horizon) %>%
#   mutate(mu.fast.max = max(intercept)) %>%
#   subset(select = c(OTU, experiment, comparison, horizon, Domain, Phylum,Class,
#                     Order, Family, Genus, Species, mu.fast.max)) %>%
#   unique()
# 
# dim(df.all.fast.growers)
# length(unique(df.all.fast.growers$OTU)) #435 taxa
# 
# 
# ggplot(df.fast.growers.after.fire) + 
#   geom_jitter(aes(x=mu.fast.max, y = Phylum, color = Phylum), alpha = 0.8) + 
#   facet_grid(horizon~comparison) + 
#   geom_vline(xintercept = 0) + 
#   theme_bw()
# 
# # Save data
# # write.csv(df.fast.growers.after.fire, 'data/sequence-data/LibCombined/corncob-output/Manuscript/fast-growers-after-dry-burns.csv', row.names = FALSE)
# # write.csv(df.all.fast.growers, 'data/sequence-data/LibCombined/corncob-output/Manuscript/all-fast-growers.csv', row.names = FALSE)
# 
# 
# 
# 
# 
# 
# 
# # Experiment 3: -----
# #    Import data.
# df.affinity <- read.csv('data/sequence-data/LibCombined/corncob-output/Ex3 Affinity/Ex3_signOTUs_MidThermoCutoff_50C.csv')
# 
# # OLD CODE
# # df.3.low <- read.csv('data/sequence-data/LibCombined/corncob-output/Ex3 Affinity/low-thermo-cutoff-50C.csv')
# # df.3.mid <- read.csv('data/sequence-data/LibCombined/corncob-output/Ex3 Affinity/Mid-thermo-cutoff-50C.csv')
# # 
# # df.3 <- rbind(df.3.low,df.3.mid)
# # 
# # for (i in 1:nrow(df.affinity)) {
# #   if (df.affinity$T.cutoff[i] == 'low thermo less than 50 C') {
# #     df.affinity$T.cutoff[i] = 'Low thermo less than 50 C'
# #   } else if (df.affinity$T.cutoff[i] == 'low thermo greater than 50 C') {
# #     df.affinity$T.cutoff[i] = 'Low thermo greater than 50 C'
# #   } else if (df.affinity$T.cutoff[i] == 'mid thermo less than 50 C') {
# #     df.affinity$T.cutoff[i] = 'Mid thermo less than 50 C'
# #   } else if (df.affinity$T.cutoff[i] == 'mid thermo greater than 50 C') {
# #     df.affinity$T.cutoff[i] = 'Mid thermo greater than 50 C'
# #   }
# # } 
# 
# # Focus on samples that saw a meaningful increase in temperature
# df.affinity <- df.affinity %>%
#   subset(T.cutoff %in% c('Mid thermo > 50'))
# 
# 
# ggplot(df.affinity, aes(x=intercept, y = Class, color = Phylum), alpha = 0.8) + 
#   geom_jitter() + 
#   facet_grid(horizon~experiment) + 
#   geom_vline(xintercept = 0) + 
#   theme_bw()
# 
# length(unique(subset(df.affinity)$OTU)) # 51 taxa
# dim(df.affinity)
# 
# 
# 
# 
# 
# 
# 
# y <- data.frame(unique(df.affinity$OTU))
# 
# # Create dataframe to hold OTU counts
# head(df.affinity)
# colnames(y)[1] = 'OTU'
# y$Responder.in.n.analyses = 0
# y$Pos.count = 0
# y$Neg.count = 0
# 
# # Count the number of times each OTU appears in dataframe
# for (i in 1:nrow(y)) {
#   for (j in 1:nrow(df.affinity)) {
#     if (df.affinity$OTU[j] == y$OTU[i]) {
#       y$Responder.in.n.analyses[i] = y$Responder.in.n.analyses[i]+1
#     } else {
#       y$Responder.in.n.analyses[i] = y$Responder.in.n.analyses[i]
#     }
#   }
# }
# 
# # Count the number of times each OTU is identified as pos/neg responder
# for (i in 1:nrow(y)) {
#   for(j in 1:nrow(df.affinity)) {
#     if (df.affinity$OTU[j] == y$OTU[i] & df.affinity$intercept[j] > 0) {
#       y$Pos.count[i] = y$Pos.count[i] + 1
#     } else if (df.affinity$OTU[j] == y$OTU[i] & df.affinity$intercept[j] < 0) {
#       y$Neg.count[i] = y$Neg.count[i] + 1
#     } else {
#       y$Pos.count[i] = y$Pos.count[i]
#       y$Neg.count[i] = y$Neg.count[i]
#     }
#   }
# }
# 
# # Combine dataframe and count data
# z <- merge(df.affinity, y)
# hist(z$Responder.in.n.analyses)
# 
# # Combine this data with the full corncob output:
# Problem.OTUs <- subset(y, Pos.count > 0 & Neg.count > 0) 
# dim(Problem.OTUs)
# 
# # Remove the OTUs identified as being both pos. and neg. responders:
# df.affinity <- subset(df.affinity, !OTU %in% Problem.OTUs$OTU)
# 
# 
# # Create final dataset of fire affinity OTUs:
# df.fire.affinity <- df.affinity %>%
#   subset(intercept > 0) %>%
#   group_by(OTU, horizon) %>%
#   # Set intercept at max intercept for repeated OTUs
#   mutate(mu.affinity.max = max(intercept),
#          mu.survival.max = 'null',
#          mu.fast.max = 'null') %>%
#   subset(select = c(OTU, comparison, DNA.type, experiment, horizon, 
#                     Domain, Phylum, Class, Order, Family, Genus, Species, 
#                     mu.fast.max, mu.survival.max, mu.affinity.max)) %>%
#   unique() 
# 
# dim(df.fire.affinity) #52 taxa
# 
# # Save output: 
# # write.csv(df.fire.affinity, 'data/sequence-data/LibCombined/corncob-output/Manuscript/affinity-dry-burns.csv', row.names = FALSE)
# 
# 
# 
# # df.3.LIn.low <- read.csv('sequence-data/LibCombined/corncob-output/Ex3 Affinity/LIn/low-thermo-cutoff-50C.csv')
# # df.3.LIn.mid <- read.csv('sequence-data/LibCombined/corncob-output/Ex3 Affinity/LIn/mid-thermo-cutoff-50C.csv')
# # 
# # df.3.LIn <- rbind(df.3.LIn.low, df.3.LIn.mid, df.3)
# # 
# # ggplot(subset(df.3.LIn, comparison %in% c('LIn.dry.v.LIn.control','LIwA.dry.v.LIwA.control')), aes(x = intercept, y = Phylum, color = Phylum)) + 
# #   geom_jitter(size=2) + 
# #   theme_bw() + 
# #   geom_vline(xintercept = 0) + 
# #   facet_grid(~comparison)
# 
# 
# 
# 
# 
# 
# # Combining the three experiments: -----
# # Create matching columns in all dataframes:
# df.survival <- df.survival %>% 
#   mutate(mu.incub.trtmt = 'null',
#          mu.pH = 'null')
# 
# df.fast <- rbind(df.fast.pH, df.fast.temp) %>% 
#   mutate(mu.burn.trtmt = 'null', 
#          mu.cntl.pH = 'null')
# 
# df.affinity <- df.affinity %>%
#   mutate(mu.incub.trtmt = 'null',
#          mu.pH = 'null')
# 
# df.full <- rbind(df.survival, df.fast, df.affinity) 
# 
# # Save output: 
# #write.csv(df.full, 'data/sequence-data/LibCombined/corncob-output/Manuscript/all-responders.csv', row.names = FALSE)
# 
# df.full <- read.csv('data/sequence-data/LibCombined/corncob-output/all-resopnders.csv')
# 
# ggplot(df.full, aes(x=intercept, y = Phylum, color = comparison)) + 
#   geom_point(width = 0.1, alpha = 0.5) + 
#   theme_bw() + 
#   geom_vline(xintercept = 0) + 
#   facet_grid(horizon~experiment)
# 
# 
# 
# 
# 
# 
# 
# 
# ### STOPPED HERE _ Jan 31, 2023
# ps.norm.pb <- readRDS('sequence-data/LibCombined/phyloseq-objects/ps.norm.pb')
# ps.norm.SI <- readRDS('sequence-data/LibCombined/phyloseq-objects/ps.norm.SI')
# ps.norm.LI <- readRDS('sequence-data/LibCombined/phyloseq-objects/ps.norm.LI')
# 
# ps.norm.pb.cDNA <- prune_samples(sample_data(ps.norm.pb)$DNA.type == 'cDNA', ps.norm.pb)
# ps.norm.pb.cDNA <- prune_taxa(taxa_sums(ps.norm.pb.cDNA)>0, ps.norm.pb.cDNA)
# 
# ps.norm.pb.gDNA <- prune_samples(sample_data(ps.norm.pb)$DNA.type == 'gDNA', ps.norm.pb)
# ps.norm.pb.gDNA <- prune_taxa(taxa_sums(ps.norm.pb.gDNA)>0, ps.norm.pb.gDNA)
# 
# ps.norm.SI <- prune_taxa(taxa_sums(ps.norm.SI)>0, ps.norm.SI)
# ps.norm.LIwA <- prune_samples(sample_data(ps.norm.LI)$incub.trtmt == 'LIwA', ps.norm.LI)
# ps.norm.LIwA <- prune_taxa(taxa_sums(ps.norm.LIwA)>0, ps.norm.LIwA)
# 
# 
# df.pb.cDNA <- psmelt(ps.norm.pb.cDNA)
# df.pb.gDNA <- psmelt(ps.norm.pb.gDNA)
# 
# df.SI <- psmelt(ps.norm.SI)
# df.LIwA <- psmelt(ps.norm.LIwA)
# 
# 
# # PLOTS -----
# # Objective: Create boxplot with relative abundance of fire surviving and fast-growing
# #     taxa in the LIwA data.
# 
# # Survival
# df.OTU.survival.dry <- df.full %>%
#   subset(experiment == 'survival') %>%
#   subset(intercept > 0 & 
#            comparison %in% c('pb.dry.v.pb.control') &
#            T.cutoff %in% c('Low thermo greater than 50 C',
#                               'Mid thermo greater than 100 C',
#                               'Mid thermo greater than 50 C')) %>%
#   mutate(Survival = 'Fire survival dry') %>%
#   subset(select = c(OTU, Survival)) %>%
#   unique()
# 
# df.OTU.survival.wet <- df.full %>%
#   subset(experiment == 'survival') %>%
#   subset(intercept > 0 & 
#            comparison %in% c('pb.wet.v.pb.control')) %>%
#   mutate(Survival = 'Fire survival - wet') %>%
#   subset(select = c(OTU, Survival)) %>%
#   subset(!OTU %in% df.OTU.survival.dry$OTU) %>%
#   unique()
# 
# # Fast growth
# df.OTU.fast.dry <- df.full %>%
#   subset(experiment == 'fast growth') %>%
#   subset(intercept > 0 & 
#            comparison %in% c('pb.dry.v.SI.dry')) %>%
#   mutate(Fast = 'Fast growth - dry') %>%
#   subset(select = c(OTU, Fast)) %>%
#   unique()
# 
# df.OTU.fast.wet <- df.full %>%
#   subset(experiment == 'fast growth') %>%
#   subset(intercept > 0 & 
#            comparison %in% c('pb.wet.v.SI.wet')) %>%
#   mutate(Fast = 'Fast growth - wet') %>%
#   subset(select = c(OTU, Fast)) %>%
#   subset(!OTU %in% df.OTU.fast.dry$OTU) %>%
#   unique()
# 
# # Affinity
# df.OTU.affinity.dry <- df.full %>%
#   subset(experiment == 'affinity') %>%
#   subset(intercept > 0 & 
#            comparison %in% c('LIwA.dry.v.LIwA.control') &
#            T.cutoff %in% c('Low thermo greater than 50 C', "Mid thermo greater than 50 C")) %>%
#   mutate(Affinity = 'affinity for post-burn environ. - dry') %>%
#   subset(select = c(OTU, Affinity)) %>%
#   unique()
# 
# df.OTU.affinity.wet <- df.full %>%
#   subset(experiment == 'affinity') %>%
#   subset(intercept > 0 & 
#            comparison %in% c('LIwA.wet.v.LIwA.control')) %>%
#   mutate(Affinity = 'affinity for post-burn environ. - wet') %>%
#   subset(select = c(OTU, Affinity)) %>%
#   subset(!OTU %in% df.OTU.affinity.dry$OTU) %>%
#   unique()
# 
# 
# 
# 
# df.trait <- merge(rbind(df.OTU.survival.dry, df.OTU.survival.wet), 
#                   rbind(df.OTU.fast.dry, df.OTU.fast.wet), by = 'OTU', all = TRUE)
# df.trait <- merge(df.trait, rbind(df.OTU.affinity.dry,df.OTU.affinity.wet), by = 'OTU', all = TRUE)
# 
# 
# df.trait$Survival <- replace_na(df.trait$Survival, 'no')
# df.trait$Fast <- replace_na(df.trait$Fast, 'no')
# df.trait$Affinity <- replace_na(df.trait$Affinity, 'no')
# 
# head(df.trait)
# 
# for (i in 1:nrow(df.trait)) {
#   if (df.trait$Survival[i] != 'no' &
#       df.trait$Fast[i] != 'no' &
#       df.trait$Affinity[i] != 'no') {
#     df.trait$Trait[i]= 'Fire survivor, fast grower, and affinity'
#   } else if (df.trait$Survival[i] != 'no' &
#              df.trait$Fast[i] != 'no' &
#              df.trait$Affinity[i] == 'no') {
#     df.trait$Trait[i]= 'Fire survivor and fast grower'
#   }  else if (df.trait$Survival[i] == 'no' &
#               df.trait$Fast[i] != 'no' &
#               df.trait$Affinity[i] != 'no') {
#     df.trait$Trait[i]= 'Fast grower and affinity'
#   }  else if (df.trait$Survival[i] != 'no' &
#               df.trait$Fast[i] == 'no' &
#               df.trait$Affinity[i] != 'no') {
#     df.trait$Trait[i]= 'Fire survivor and affinity'
#   }  else if (df.trait$Survival[i] != 'no' &
#               df.trait$Fast[i] == 'no' &
#               df.trait$Affinity[i] == 'no') {
#     df.trait$Trait[i]= 'Fire survivor'
#   }  else if (df.trait$Survival[i] == 'no' &
#               df.trait$Fast[i] != 'no' &
#               df.trait$Affinity[i] == 'no') {
#     df.trait$Trait[i]= 'Fast grower'
#   }  else if (df.trait$Survival[i] == 'no' &
#               df.trait$Fast[i] == 'no' &
#               df.trait$Affinity[i] != 'no') {
#     df.trait$Trait[i]= 'Affinity'
#   } else {
#     df.trait$Trait[i]= 'no trait'
#   }
# }
# 
# 
# #df.trait <- subset(df.trait, select = c(OTU, Trait))
# 
# # pb
# df.pb.dry  <- subset(df.pb.cDNA, burn.trtmt %in% c('control','wet','dry') & Abundance > 0) %>%
#   merge(df.trait, by = c('OTU'), all = TRUE) %>%
#   subset(Abundance > 0)
# 
# # SI
# df.SI.dry <- subset(df.SI, burn.trtmt %in% c('control','wet','dry') & Abundance > 0) %>%
#   merge(df.trait, by = c('OTU'), all = TRUE) %>%
#   subset(Abundance > 0)
# 
# # LIwA
# df.LI.dry <- subset(df.LIwA, burn.trtmt %in% c('control','wet','dry') & Abundance > 0) %>%
#   merge(df.trait, by = c('OTU'), all = TRUE) %>%
#   subset(Abundance > 0)
# 
# 
# 
# df.dry <- rbind(df.pb.dry, df.SI.dry, df.LI.dry)
# 
# df.dry$Survival <- replace_na(df.dry$Survival, 'no')
# df.dry$Fast <- replace_na(df.dry$Fast, 'no')
# df.dry$Affinity <- replace_na(df.dry$Affinity, 'no')
# 
# 
# 
# df.dry$burn.trtmt <- factor(df.dry$burn.trtmt, levels = c('control','wet','dry'))
# df.dry$horizon <- factor(df.dry$horizon, levels = c('O','A'))
# df.dry$incub.trtmt <- factor(df.dry$incub.trtmt, levels = c('pb','SI','LIwA'))
# df.dry$Trait <- factor(df.dry$Trait, levels = c('Fire survivor',
#                                                 'Fast grower',
#                                                 'Fire survivor and fast grower',
#                                                 'Affinity',
#                                                 'Fire survivor and affinity',
#                                                 'Fast grower and affinity',
#                                                 'NA'))
# 
# 
# Burn.labels = c('Dry soil burn','Wet soil burn','Control, no burn')
# names(Burn.labels) = c("dry",'wet','control')
# 
# Incubation.labels = c('24-hours post-burn','5 weeks', 'Autoclave + 6 months')
# names(Incubation.labels) = c('pb','SI','LIwA')
# 
# 
# p = ggplot(subset(df.dry, horizon == 'O'), aes(x=site, y=Abundance, fill=Order))+
#   geom_bar(stat = 'identity') +
#   facet_grid(incub.trtmt~burn.trtmt, scales = 'free', 
#              labeller = labeller(incub.trtmt = Incubation.labels, 
#                                  burn.trtmt = Burn.labels)) +
#   labs(title = 'Traits yes/no, pb = cDNA, Horizon = O',
#        y = 'Relative abundance')+
#   theme(legend.position = 'none')
# 
# p
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # Other stuff -----
# 
# df.pb.full <- merge(df.pb, subset(df.full, select = c(OTU, intercept, comparison, experiment, T.cutoff, pH.cutoff, horizon)))
# df.pb.full <- subset(df.pb.full, Abundance > 0) %>%
#   subset(DNA.type == 'cDNA')
# 
# ggplot(df.pb.full, aes(x = burn.trtmt, y =intercept, color = Phylum)) + 
#   geom_jitter(size=2, alpha = 0.5) + 
#   facet_grid(horizon~experiment)
# 
# 
# ggplot(df.pb.full, aes(y = Abundance)) + 
#   geom_bar(aes(fill=Phylum), position = 'stack')
# 
# 
# df.affinity <- read.csv('sequence-data/LibCombined/corncob-output/Ex3 Affinity/affinity-dry-burns.csv')
# # OTUs with positive mu value from dry burns with low thermo cutoff > 50 C
# 
# df.fire.survivors  <- read.csv('sequence-data/LibCombined/corncob-output/Ex1 Fire survival/fire-survivors-all-burns.csv')
# # OTUs with positive mu value from dry burns exluding wet & control taxa
# 
# df.fast.growers.after.fire  <- read.csv('sequence-data/LibCombined/corncob-output/Ex2 Fast growth/fast-growers-after-dry-burns.csv')
# # OTUs with positive mu value in dry burns
# 
# colnames(df.affinity)
# colnames(df.fire.survivors)
# colnames(df.fast.growers.after.fire)
# 
# x <- subset(df.affinity, select = c(OTU,Domain,Phylum,Class,Order,Family,Genus, Species))
# y <- subset(df.fire.survivors, select = c(OTU,Domain,Phylum,Class,Order,Family,Genus, Species))
# z <- subset(df.fast.growers.after.fire, select = c(OTU,Domain,Phylum,Class,Order,Family,Genus, Species))
# 
# 
# df.OTU.info <- unique(rbind(x,y,z))
# 
# 
# df.full <- df.fire.survivors %>%
#   subset(select = c(OTU, mu.survival.max)) %>%
#   merge(subset(df.fast.growers.after.fire, select = c(OTU, mu.fastgrowth.max)), all = TRUE) %>%
#   merge(subset(df.affinity, select = c(OTU, mu.affinity.max)), all = TRUE) %>%
#   merge(df.OTU.info, all = TRUE)
# 
# # How to normalize response trait, mu? Normalizee like relative abundance?
# 
# ggplot(df.full) +
#   geom_point(aes(y=Phylum, x = mu.fastgrowth.max), color='green', size=2) +
#   geom_point(aes(y=Phylum, x = mu.affinity.max), color = 'blue', size=2) + 
#   geom_point(aes(y=Phylum, x = mu.survival.max), color = 'orange',size=2) + 
#   theme_bw() 
# 
# df.subset <- subset(df.full, mu.survival.max != 0 & mu.fastgrowth.max != 0) %>%
#   merge(subset(df.full, mu.fastgrowth.max !=0 & mu.affinity.max != 0), all = TRUE) %>%
#   merge(subset(df.full, mu.affinity.max!= 0 & mu.survival.max !=0), all = TRUE)
# 
# 
# 
# 
# 
# 
# 
# df.subset <- merge(df.survival, df.fast, by = 'OTU',all = FALSE)
# df.subset <- merge(df.subset, df.affinity, by = 'OTU', all = FALSE) 
# df.subset <- subset(df.subset, intercept > 0 & intercept.x > 0 & intercept.y > 0)
# 
# A = ggplot(df.subset) + 
#   geom_point(aes(x=intercept.x, y = intercept.y, color = Phylum), alpha=0.6,size=3)+
#   geom_vline(xintercept = 0) + 
#   geom_hline(yintercept = 0) + 
#   theme(legend.position = 'none')
# 
# C = ggplot(df.subset) + 
#   geom_point(aes(x=intercept.x, y = intercept, color = Phylum), alpha=0.6, size=3)+
#   geom_vline(xintercept = 0) + 
#   geom_hline(yintercept = 0) + 
#   theme(legend.position = 'none')
# 
# B = ggplot(df.subset) + 
#   geom_point(aes(x=intercept.y, y = intercept, color = Phylum), alpha=0.6, size=3) +
#   geom_vline(xintercept = 0) + 
#   geom_hline(yintercept = 0) 
#   theme(legend.position = 'none')
# 
# 
#   
# 
# cowplot::plot_grid(A,C,B, ncol=3, labels = c('A','C','B'))
# 
# 
# 
# 
# 
