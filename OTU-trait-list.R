# Today's goal (21.Oct.2021) - for each experiment, create a dataframe with the
#    positive responders and a composite mu value. I can calculate this mu value
#    in different ways - max, min, mean - and create a weighted mu value.
#    So the resulting dataframe will have a mu.survival.max (or mu.mean or mu.min) and a 
#    mu.survival.weighted

library(ggplot2)
library(dplyr)
library(phyloseq)
library(tidyr)

setwd('../Box/WhitmanLab/Projects/WoodBuffalo/FireSim2019/data')
# Experiment 1: -----
#     Import data
df.1.low <- read.csv('sequence-data/LibCombined/corncob-output/Ex1 Fire survival/low-thermo-cutoff-50C.csv')
df.1.mid <- read.csv('sequence-data/LibCombined/corncob-output/Ex1 Fire survival/Mid-thermo-cutoff.csv')

df.1 <- rbind(df.1.low,df.1.mid)

head(df.1)
dim(df.1) # total of 951 responders
length(unique(df.1$OTU)) # 388 unique OTU responders

# Create dataframe with unique OTU identifies
y <- data.frame(unique(df.1$OTU))
colnames(y)[1] = 'OTU'
y$Responder.in.n.analyses = 0
y$Pos.count = 0
y$Neg.count = 0

# Calculate the number of times that each OTU was identified as a pos/neg fire responder
for (i in 1:nrow(y)) {
  for (j in 1:nrow(df.1)) {
    if (df.1$OTU[j] == y$OTU[i]) {
      y$Responder.in.n.analyses[i] = y$Responder.in.n.analyses[i]+1
    } else {
      y$Responder.in.n.analyses[i] = y$Responder.in.n.analyses[i]
    }
  }
}

for (i in 1:nrow(y)) {
  for(j in 1:nrow(df.1)) {
    if (df.1$OTU[j] == y$OTU[i] & df.1$Estimate[j] > 0) {
      y$Pos.count[i] = y$Pos.count[i] + 1
    } else if (df.1$OTU[j] == y$OTU[i] & df.1$Estimate[j] < 0) {
      y$Neg.count[i] = y$Neg.count[i] + 1
    } else {
      y$Pos.count[i] = y$Pos.count[i]
      y$Neg.count[i] = y$Neg.count[i]
    }
  }
}

# Combine this data with the full corncob output:
z <- merge(df.1, y)
hist(z$Responder.in.n.analyses)

Problem.OTUs <- subset(y, Pos.count > 0 & Neg.count > 0) 

df.problem <- subset(df.1, OTU %in% Problem.OTUs$OTU) %>%
  merge(Problem.OTUs)


# Use only the dry v. control data, and the max temperature cutoff.
df.fire.survivors.dry <- df.1 %>%
  subset(Comparison == 'pb.dry.v.pb.control') %>%
  subset(Temp.cutoff %in% c('Mid thermo greater than 50 C', 'Low thermo greater than 50 C', 'Mid thermo greater than 100 C')) %>%
  subset(Estimate > 0) %>%
  group_by(OTU, horizon) %>%
  mutate(mu.survival.max = max(Estimate)) %>%
  subset(select = c(OTU, Comparison, Experiment, horizon,
                    Domain, Phylum, Class, Order, Family, Genus, Species, 
                    mu.survival.max)) %>%
  unique() %>%
  subset(!OTU %in% Problem.OTUs$OTU)


df.fire.susceptible.dry <- df.1 %>%
  subset(Comparison == 'pb.dry.v.pb.control') %>%
  subset(Temp.cutoff %in% c('Mid thermo greater than 50 C', 'Low thermo greater than 50 C', 'Mid thermo greater than 100 C')) %>%
  subset(Estimate < 0) %>%
  group_by(OTU, horizon) %>%
  mutate(mu.survival.min = min(Estimate)) %>%
  subset(select = c(OTU, Comparison, Experiment, horizon,
                    Domain, Phylum, Class, Order, Family, Genus, Species, 
                    mu.survival.min)) %>%
  unique() %>%
  subset(!OTU %in% Problem.OTUs$OTU)
dim(df.fire.susceptible.dry)

# Also keep a list of all the survivors:
df.fire.survivors <- df.1 %>%
  subset(Estimate>0) %>%
  group_by(OTU, horizon) %>%
  mutate(mu.survival.max = max(Estimate)) %>%
  subset(select = c(OTU, Comparison, Experiment, horizon,
                    Domain, Phylum, Class, Order, Family, Genus, Species, 
                    mu.survival.max)) %>%
  unique() %>%
  subset(!OTU %in% Problem.OTUs$OTU)

df.fire.susceptible <- df.1 %>%
  subset(Estimate<0) %>%
  group_by(OTU, horizon) %>%
  mutate(mu.survival.min = min(Estimate)) %>%
  subset(select = c(OTU, Comparison, Experiment, horizon,
                    Domain, Phylum, Class, Order, Family, Genus, Species, 
                    mu.survival.min)) %>%
  unique() %>%
  subset(!OTU %in% Problem.OTUs$OTU)



ggplot(df.fire.survivors.dry, aes(x=mu.survival.max, y = Phylum, color = Phylum), alpha = 0.8) + 
  geom_jitter() + 
  facet_grid(horizon~Comparison) + 
  geom_vline(xintercept = 0) + 
  theme_bw()


write.csv(df.fire.survivors, 'sequence-data/LibCombined/corncob-output/Ex1 Fire survival/fire-survivors-all-burns.csv', row.names = FALSE)
write.csv(df.fire.susceptible, 'sequence-data/LibCombined/corncob-output/Ex1 Fire survival/fire-susceptible-all-burns.csv', row.names = FALSE)
write.csv(df.fire.survivors.dry, 'sequence-data/LibCombined/corncob-output/Ex1 Fire survival/fire-survivors-dry-burns.csv', row.names = FALSE)
write.csv(df.fire.susceptible.dry, 'sequence-data/LibCombined/corncob-output/Ex1 Fire survival/fire-susceptible-dry-burns.csv', row.names = FALSE)



# Experiment 2: -----
#    Import data.
df.2.pH <- read.csv('sequence-data/LibCombined/corncob-output/Ex2 Fast growth/3-sliding-pH-categories.csv')
df.2.low <- read.csv('sequence-data/LibCombined/corncob-output/Ex2 Fast growth/low-thermo-cutoff-50C.csv')
df.2.mid <- read.csv('sequence-data/LibCombined/corncob-output/Ex2 Fast growth/Mid-thermo-cutoff-50.csv')
df.2.mid100<- read.csv('sequence-data/LibCombined/corncob-output/Ex2 Fast growth/Mid-thermo-cutoff-100.csv')

df.2 <- rbind(df.2.pH, df.2.low, df.2.mid, df.2.mid100)
dim(df.2)
length(unique(df.2$OTU))

y <- data.frame(unique(df.2$OTU))

head(df.2)
colnames(y)[1] = 'OTU'
y$Responder.in.n.analyses = 0
y$Pos.count = 0
y$Neg.count = 0

for (i in 1:nrow(y)) {
  for (j in 1:nrow(df.2)) {
    if (df.2$OTU[j] == y$OTU[i]) {
      y$Responder.in.n.analyses[i] = y$Responder.in.n.analyses[i]+1
    } else {
      y$Responder.in.n.analyses[i] = y$Responder.in.n.analyses[i]
    }
  }
}

for (i in 1:nrow(y)) {
  for(j in 1:nrow(df.2)) {
    if (df.2$OTU[j] == y$OTU[i] & df.2$Estimate[j] > 0) {
      y$Pos.count[i] = y$Pos.count[i] + 1
    } else if (df.2$OTU[j] == y$OTU[i] & df.2$Estimate[j] < 0) {
      y$Neg.count[i] = y$Neg.count[i] + 1
    } else {
      y$Pos.count[i] = y$Pos.count[i]
      y$Neg.count[i] = y$Neg.count[i]
    }
  }
}

z <- merge(df.2, y)
hist(z$Responder.in.n.analyses)

Dual.Responders <- subset(y, Pos.count > 0 & Neg.count > 0) 
dim(Dual.Responders)

df.dual.responders <- subset(df.2, OTU %in% Dual.Responders$OTU) %>%
  merge(Dual.Responders)

ggplot(df.dual.responders, aes(x=Estimate, y=Phylum, color = Phylum)) +
  geom_point(size=2) + 
  theme_bw()+
  facet_grid(horizon~Comparison) + 
  geom_vline(xintercept = 0)

length(unique(df.dual.responders$OTU))

# Create final dataset of fast/slow growing OTUs:
# Output:
df.dry <- subset(df.2, Comparison == 'pb.dry.v.SI.dry')
df.control <- subset(df.2, Comparison == 'pb.control.v.SI.control')
df.wet <- subset(df.2, Comparison == 'pb.wet.v.SI.wet')


df.fast.growers.after.fire <- df.dry %>%
  subset(Estimate > 0) %>%
  subset(!OTU %in% subset(df.wet, Estimate > 1)$OTU) %>%
  subset(!OTU %in% subset(df.control, Estimate > 1)$OTU) %>%
  group_by(OTU, horizon) %>%
  mutate(mu.fastgrowth.max = max(Estimate)) %>%
  subset(select = c(OTU, Comparison, Experiment, horizon,
                    Domain, Phylum, Class, Order, Family, Genus, Species, 
                    mu.fastgrowth.max)) %>%
  unique() 

length(unique(df.fast.growers.after.fire$OTU))


df.all.fast.growers <- subset(df.2, Estimate > 0) %>%
  group_by(OTU, horizon) %>%
  mutate(mu.fastgrowth.max = max(Estimate)) %>%
  subset(select = c(OTU, Experiment, Comparison, horizon, Domain, Phylum,Class,
                    Order, Family, Genus, Species, mu.fastgrowth.max)) %>%
  unique()

dim(df.all.fast.growers)
length(unique(df.all.fast.growers$OTU))




ggplot(df.fast.growers.after.fire) + 
  geom_jitter(aes(x=mu.fastgrowth.max, y = Phylum, color = Phylum), alpha = 0.8) + 
  facet_grid(horizon~Comparison) + 
  geom_vline(xintercept = 0) + 
  theme_bw()

write.csv(df.fast.growers.after.fire, 'sequence-data/LibCombined/corncob-output/Ex2 Fast growth/fast-growers-after-dry-burns.csv', row.names = FALSE)
write.csv(df.all.fast.growers, 'sequence-data/LibCombined/corncob-output/Ex2 Fast growth/all-fast-growers.csv', row.names = FALSE)







# Experiment 3: -----
#    Import data.
df.3.low <- read.csv('data/sequence-data/LibCombined/corncob-output/Ex3 Affinity/low-thermo-cutoff-50C.csv')
df.3.mid <- read.csv('data/sequence-data/LibCombined/corncob-output/Ex3 Affinity/Mid-thermo-cutoff-50C.csv')

df.3 <- rbind(df.3.low,df.3.mid)

for (i in 1:nrow(df.3)) {
  if (df.3$Temp.cutoff[i] == 'low thermo less than 50 C') {
    df.3$Temp.cutoff[i] = 'Low thermo less than 50 C'
  } else if (df.3$Temp.cutoff[i] == 'low thermo greater than 50 C') {
    df.3$Temp.cutoff[i] = 'Low thermo greater than 50 C'
  } else if (df.3$Temp.cutoff[i] == 'mid thermo less than 50 C') {
    df.3$Temp.cutoff[i] = 'Mid thermo less than 50 C'
  } else if (df.3$Temp.cutoff[i] == 'mid thermo greater than 50 C') {
    df.3$Temp.cutoff[i] = 'Mid thermo greater than 50 C'
  }
} 

# Focus on samples that saw a meaningful increase in temperature
df.3 <- df.3 %>%
  subset(Temp.cutoff %in% c('low thermo greater than 50 C', 'mid thermo greater than 50 C'))


ggplot(df.3, aes(x=Estimate, y = Class, color = Phylum), alpha = 0.8) + 
  geom_jitter() + 
  facet_grid(horizon~Experiment) + 
  geom_vline(xintercept = 0) + 
  theme_bw()

length(unique(subset(df.3)$OTU))
dim(df.3)


# Create dataframe with unique OTU identifies
df.unique.OTUs <- data.frame(unique(df.3$OTU))
colnames(df.unique.OTUs)[1] = 'OTU'
df.unique.OTUs$Responder.in.n.analyses = 0

# Calculate the number of times that each OTU was identified as a pos/neg fire responder
for (i in 1:nrow(df.unique.OTUs)) {
  for (j in 1:nrow(df.3)) {
    if (df.3$OTU[j] == df.unique.OTUs$OTU[i]) {
      df.unique.OTUs$Responder.in.n.analyses[i] = df.unique.OTUs$Responder.in.n.analyses[i]+1
    } else {
      df.unique.OTUs$Responder.in.n.analyses[i] = df.unique.OTUs$Responder.in.n.analyses[i]
    }
  }
}

hist(df.unique.OTUs$Responder.in.n.analyses)

df.unique.OTUs$Pos.count = 0
df.unique.OTUs$Neg.count = 0

for (i in 1:nrow(df.unique.OTUs)) {
  for(j in 1:nrow(df.3)) {
    if (df.3$OTU[j] == df.unique.OTUs$OTU[i] & df.3$Estimate[j] > 0) {
      df.unique.OTUs$Pos.count[i] = df.unique.OTUs$Pos.count[i] + 1
    } else if (df.3$OTU[j] == df.unique.OTUs$OTU[i] & df.3$Estimate[j] < 0) {
      df.unique.OTUs$Neg.count[i] = df.unique.OTUs$Neg.count[i] + 1
    } else {
      df.unique.OTUs$Pos.count[i] = df.unique.OTUs$Pos.count[i]
      df.unique.OTUs$Neg.count[i] = df.unique.OTUs$Neg.count[i]
    }
  }
}


# Combine this data with the full corncob output:
Problem.OTUs <- subset(df.unique.OTUs, Pos.count > 0 & Neg.count > 0) 
dim(Problem.OTUs)

# Remove the OTUs identified as being both pos. and neg. responders:
df.3 <- subset(df.3, !OTU %in% Problem.OTUs$OTU)


# Create final dataset of fire surviving/susceptible OTUs:
df.affinity <- df.3 %>%
  subset(Estimate > 0) %>%
  group_by(OTU, horizon) %>%
  mutate(mu.affinity.max = max(Estimate)) %>%
  subset(select = c(OTU, Comparison, Experiment, horizon,
                    Domain, Phylum, Class, Order, Family, Genus, Species, 
                    mu.affinity.max)) %>%
  unique() 

dim(df.affinity)

write.csv(df.affinity, 'sequence-data/LibCombined/corncob-output/Ex3 Affinity/affinity-dry-burns.csv', row.names = FALSE)



df.3.LIn.low <- read.csv('sequence-data/LibCombined/corncob-output/Ex3 Affinity/LIn/low-thermo-cutoff-50C.csv')
df.3.LIn.mid <- read.csv('sequence-data/LibCombined/corncob-output/Ex3 Affinity/LIn/mid-thermo-cutoff-50C.csv')

df.3.LIn <- rbind(df.3.LIn.low, df.3.LIn.mid, df.3)

ggplot(subset(df.3.LIn, Comparison %in% c('LIn.dry.v.LIn.control','LIwA.dry.v.LIwA.control')), aes(x = Estimate, y = Phylum, color = Phylum)) + 
  geom_jitter(size=2) + 
  theme_bw() + 
  geom_vline(xintercept = 0) + 
  facet_grid(~Comparison)






# Combining the three experiments: -----
df.full <- rbind(df.1, df.2, df.3) 

#write.csv(df.full, 'sequence-data/LibCombined/corncob-output/all-resopnders.csv', row.names = FALSE)

df.full <- read.csv('sequence-data/LibCombined/corncob-output/all-resopnders.csv')

ggplot(df.full, aes(x=Estimate, y = Phylum, color = Comparison)) + 
  geom_point(width = 0.1, alpha = 0.5) + 
  theme_bw() + 
  geom_vline(xintercept = 0) + 
  facet_grid(horizon~Experiment)


ps.norm.pb <- readRDS('sequence-data/LibCombined/phyloseq-objects/ps.norm.pb')
ps.norm.SI <- readRDS('sequence-data/LibCombined/phyloseq-objects/ps.norm.SI')
ps.norm.LI <- readRDS('sequence-data/LibCombined/phyloseq-objects/ps.norm.LI')

ps.norm.pb.cDNA <- prune_samples(sample_data(ps.norm.pb)$DNA.type == 'cDNA', ps.norm.pb)
ps.norm.pb.cDNA <- prune_taxa(taxa_sums(ps.norm.pb.cDNA)>0, ps.norm.pb.cDNA)

ps.norm.pb.gDNA <- prune_samples(sample_data(ps.norm.pb)$DNA.type == 'gDNA', ps.norm.pb)
ps.norm.pb.gDNA <- prune_taxa(taxa_sums(ps.norm.pb.gDNA)>0, ps.norm.pb.gDNA)

ps.norm.SI <- prune_taxa(taxa_sums(ps.norm.SI)>0, ps.norm.SI)
ps.norm.LIwA <- prune_samples(sample_data(ps.norm.LI)$incub.trtmt == 'LIwA', ps.norm.LI)
ps.norm.LIwA <- prune_taxa(taxa_sums(ps.norm.LIwA)>0, ps.norm.LIwA)


df.pb.cDNA <- psmelt(ps.norm.pb.cDNA)
df.pb.gDNA <- psmelt(ps.norm.pb.gDNA)

df.SI <- psmelt(ps.norm.SI)
df.LIwA <- psmelt(ps.norm.LIwA)


# PLOTS -----
# Objective: Create boxplot with relative abundance of fire surviving and fast-growing
#     taxa in the LIwA data.

# Survival
df.OTU.survival.dry <- df.full %>%
  subset(Experiment == 'survival') %>%
  subset(Estimate > 0 & 
           Comparison %in% c('pb.dry.v.pb.control') &
           Temp.cutoff %in% c('Low thermo greater than 50 C',
                              'Mid thermo greater than 100 C',
                              'Mid thermo greater than 50 C')) %>%
  mutate(Survival = 'Fire survival dry') %>%
  subset(select = c(OTU, Survival)) %>%
  unique()

df.OTU.survival.wet <- df.full %>%
  subset(Experiment == 'survival') %>%
  subset(Estimate > 0 & 
           Comparison %in% c('pb.wet.v.pb.control')) %>%
  mutate(Survival = 'Fire survival - wet') %>%
  subset(select = c(OTU, Survival)) %>%
  subset(!OTU %in% df.OTU.survival.dry$OTU) %>%
  unique()

# Fast growth
df.OTU.fast.dry <- df.full %>%
  subset(Experiment == 'fast growth') %>%
  subset(Estimate > 0 & 
           Comparison %in% c('pb.dry.v.SI.dry')) %>%
  mutate(Fast = 'Fast growth - dry') %>%
  subset(select = c(OTU, Fast)) %>%
  unique()

df.OTU.fast.wet <- df.full %>%
  subset(Experiment == 'fast growth') %>%
  subset(Estimate > 0 & 
           Comparison %in% c('pb.wet.v.SI.wet')) %>%
  mutate(Fast = 'Fast growth - wet') %>%
  subset(select = c(OTU, Fast)) %>%
  subset(!OTU %in% df.OTU.fast.dry$OTU) %>%
  unique()

# Affinity
df.OTU.affinity.dry <- df.full %>%
  subset(Experiment == 'affinity') %>%
  subset(Estimate > 0 & 
           Comparison %in% c('LIwA.dry.v.LIwA.control') &
           Temp.cutoff %in% c('Low thermo greater than 50 C', "Mid thermo greater than 50 C")) %>%
  mutate(Affinity = 'affinity for post-burn environ. - dry') %>%
  subset(select = c(OTU, Affinity)) %>%
  unique()

df.OTU.affinity.wet <- df.full %>%
  subset(Experiment == 'affinity') %>%
  subset(Estimate > 0 & 
           Comparison %in% c('LIwA.wet.v.LIwA.control')) %>%
  mutate(Affinity = 'affinity for post-burn environ. - wet') %>%
  subset(select = c(OTU, Affinity)) %>%
  subset(!OTU %in% df.OTU.affinity.dry$OTU) %>%
  unique()




df.trait <- merge(rbind(df.OTU.survival.dry, df.OTU.survival.wet), 
                  rbind(df.OTU.fast.dry, df.OTU.fast.wet), by = 'OTU', all = TRUE)
df.trait <- merge(df.trait, rbind(df.OTU.affinity.dry,df.OTU.affinity.wet), by = 'OTU', all = TRUE)


df.trait$Survival <- replace_na(df.trait$Survival, 'no')
df.trait$Fast <- replace_na(df.trait$Fast, 'no')
df.trait$Affinity <- replace_na(df.trait$Affinity, 'no')

head(df.trait)

for (i in 1:nrow(df.trait)) {
  if (df.trait$Survival[i] != 'no' &
      df.trait$Fast[i] != 'no' &
      df.trait$Affinity[i] != 'no') {
    df.trait$Trait[i]= 'Fire survivor, fast grower, and affinity'
  } else if (df.trait$Survival[i] != 'no' &
             df.trait$Fast[i] != 'no' &
             df.trait$Affinity[i] == 'no') {
    df.trait$Trait[i]= 'Fire survivor and fast grower'
  }  else if (df.trait$Survival[i] == 'no' &
              df.trait$Fast[i] != 'no' &
              df.trait$Affinity[i] != 'no') {
    df.trait$Trait[i]= 'Fast grower and affinity'
  }  else if (df.trait$Survival[i] != 'no' &
              df.trait$Fast[i] == 'no' &
              df.trait$Affinity[i] != 'no') {
    df.trait$Trait[i]= 'Fire survivor and affinity'
  }  else if (df.trait$Survival[i] != 'no' &
              df.trait$Fast[i] == 'no' &
              df.trait$Affinity[i] == 'no') {
    df.trait$Trait[i]= 'Fire survivor'
  }  else if (df.trait$Survival[i] == 'no' &
              df.trait$Fast[i] != 'no' &
              df.trait$Affinity[i] == 'no') {
    df.trait$Trait[i]= 'Fast grower'
  }  else if (df.trait$Survival[i] == 'no' &
              df.trait$Fast[i] == 'no' &
              df.trait$Affinity[i] != 'no') {
    df.trait$Trait[i]= 'Affinity'
  } else {
    df.trait$Trait[i]= 'no trait'
  }
}


#df.trait <- subset(df.trait, select = c(OTU, Trait))

# pb
df.pb.dry  <- subset(df.pb.cDNA, burn.trtmt %in% c('control','wet','dry') & Abundance > 0) %>%
  merge(df.trait, by = c('OTU'), all = TRUE) %>%
  subset(Abundance > 0)

# SI
df.SI.dry <- subset(df.SI, burn.trtmt %in% c('control','wet','dry') & Abundance > 0) %>%
  merge(df.trait, by = c('OTU'), all = TRUE) %>%
  subset(Abundance > 0)

# LIwA
df.LI.dry <- subset(df.LIwA, burn.trtmt %in% c('control','wet','dry') & Abundance > 0) %>%
  merge(df.trait, by = c('OTU'), all = TRUE) %>%
  subset(Abundance > 0)



df.dry <- rbind(df.pb.dry, df.SI.dry, df.LI.dry)

df.dry$Survival <- replace_na(df.dry$Survival, 'no')
df.dry$Fast <- replace_na(df.dry$Fast, 'no')
df.dry$Affinity <- replace_na(df.dry$Affinity, 'no')



df.dry$burn.trtmt <- factor(df.dry$burn.trtmt, levels = c('control','wet','dry'))
df.dry$horizon <- factor(df.dry$horizon, levels = c('O','A'))
df.dry$incub.trtmt <- factor(df.dry$incub.trtmt, levels = c('pb','SI','LIwA'))
df.dry$Trait <- factor(df.dry$Trait, levels = c('Fire survivor',
                                                'Fast grower',
                                                'Fire survivor and fast grower',
                                                'Affinity',
                                                'Fire survivor and affinity',
                                                'Fast grower and affinity',
                                                'NA'))


Burn.labels = c('Dry soil burn','Wet soil burn','Control, no burn')
names(Burn.labels) = c("dry",'wet','control')

Incubation.labels = c('24-hours post-burn','5 weeks', 'Autoclave + 6 months')
names(Incubation.labels) = c('pb','SI','LIwA')


p = ggplot(subset(df.dry, horizon == 'O'), aes(x=site, y=Abundance, fill=Order))+
  geom_bar(stat = 'identity') +
  facet_grid(incub.trtmt~burn.trtmt, scales = 'free', 
             labeller = labeller(incub.trtmt = Incubation.labels, 
                                 burn.trtmt = Burn.labels)) +
  labs(title = 'Traits yes/no, pb = cDNA, Horizon = O',
       y = 'Relative abundance')+
  theme(legend.position = 'none')

p












# Other stuff -----

df.pb.full <- merge(df.pb, subset(df.full, select = c(OTU, Estimate, Comparison, Experiment, Temp.cutoff, pH.cutoff, horizon)))
df.pb.full <- subset(df.pb.full, Abundance > 0) %>%
  subset(DNA.type == 'cDNA')

ggplot(df.pb.full, aes(x = burn.trtmt, y =Estimate, color = Phylum)) + 
  geom_jitter(size=2, alpha = 0.5) + 
  facet_grid(horizon~Experiment)


ggplot(df.pb.full, aes(y = Abundance)) + 
  geom_bar(aes(fill=Phylum), position = 'stack')


df.affinity <- read.csv('sequence-data/LibCombined/corncob-output/Ex3 Affinity/affinity-dry-burns.csv')
# OTUs with positive mu value from dry burns with low thermo cutoff > 50 C

df.fire.survivors  <- read.csv('sequence-data/LibCombined/corncob-output/Ex1 Fire survival/fire-survivors-all-burns.csv')
# OTUs with positive mu value from dry burns exluding wet & control taxa

df.fast.growers.after.fire  <- read.csv('sequence-data/LibCombined/corncob-output/Ex2 Fast growth/fast-growers-after-dry-burns.csv')
# OTUs with positive mu value in dry burns

colnames(df.affinity)
colnames(df.fire.survivors)
colnames(df.fast.growers.after.fire)

x <- subset(df.affinity, select = c(OTU,Domain,Phylum,Class,Order,Family,Genus, Species))
y <- subset(df.fire.survivors, select = c(OTU,Domain,Phylum,Class,Order,Family,Genus, Species))
z <- subset(df.fast.growers.after.fire, select = c(OTU,Domain,Phylum,Class,Order,Family,Genus, Species))


df.OTU.info <- unique(rbind(x,y,z))


df.full <- df.fire.survivors %>%
  subset(select = c(OTU, mu.survival.max)) %>%
  merge(subset(df.fast.growers.after.fire, select = c(OTU, mu.fastgrowth.max)), all = TRUE) %>%
  merge(subset(df.affinity, select = c(OTU, mu.affinity.max)), all = TRUE) %>%
  merge(df.OTU.info, all = TRUE)

# How to normalize response trait, mu? Normalizee like relative abundance?

ggplot(df.full) +
  geom_point(aes(y=Phylum, x = mu.fastgrowth.max), color='green', size=2) +
  geom_point(aes(y=Phylum, x = mu.affinity.max), color = 'blue', size=2) + 
  geom_point(aes(y=Phylum, x = mu.survival.max), color = 'orange',size=2) + 
  theme_bw() 

df.subset <- subset(df.full, mu.survival.max != 0 & mu.fastgrowth.max != 0) %>%
  merge(subset(df.full, mu.fastgrowth.max !=0 & mu.affinity.max != 0), all = TRUE) %>%
  merge(subset(df.full, mu.affinity.max!= 0 & mu.survival.max !=0), all = TRUE)







df.subset <- merge(df.1, df.2, by = 'OTU',all = FALSE)
df.subset <- merge(df.subset, df.3, by = 'OTU', all = FALSE) 
df.subset <- subset(df.subset, Estimate > 0 & Estimate.x > 0 & Estimate.y > 0)

A = ggplot(df.subset) + 
  geom_point(aes(x=Estimate.x, y = Estimate.y, color = Phylum), alpha=0.6,size=3)+
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  theme(legend.position = 'none')

C = ggplot(df.subset) + 
  geom_point(aes(x=Estimate.x, y = Estimate, color = Phylum), alpha=0.6, size=3)+
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  theme(legend.position = 'none')

B = ggplot(df.subset) + 
  geom_point(aes(x=Estimate.y, y = Estimate, color = Phylum), alpha=0.6, size=3) +
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) 
  theme(legend.position = 'none')


  

cowplot::plot_grid(A,C,B, ncol=3, labels = c('A','C','B'))





