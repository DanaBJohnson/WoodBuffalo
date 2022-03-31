# --------------------
# Title: Creating phyloseq object
# Author: Dana Johnson
# Date: 2021-April-6
#
# Purpose: Make figures from seq data 

# Libraries -----
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(corncob)


# OBJECTIVE: Combine reads from Lanes 1 and 2 ----- 

# Create master sample list:
### Lib1:
SamList <- read.csv('../data/sequence-data/Lib1-CombinedLanes/Seq-IDs-by-lane-lib1.csv',fileEncoding="UTF-8-BOM")
### Lib2:
SamList <- read.csv('../data/sequence-data/Lib2-CombinedLanes/Seq-IDs-by-lane-lib2.csv', fileEncoding = "UTF-8-BOM")
### Combined Libs:
SamList <- read.csv('../Box/WhitmanLab/Projects/WoodBuffalo/FireSim2019/data/sequence-data/LibCombined/Seq-IDs-combined-libs.csv', fileEncoding = "UTF-8-BOM")

SamList <- read.csv('data-exploration/data/LibCombined/Seq-IDs-combined-libs.csv')
head(SamList)

# End of Lane 1-specific code; skip to line 33

# Converge ID list to phyloseq object
SamList <- sample_data(SamList)

# Assign the sample names to match the samples' IDs
sample_names(SamList) = SamList$Lane.IDs

# Import necessary seq data
### Lib1:
ps <- import_biom('../data/sequence-data/Lib1-CombinedLanes/full-diversity_OTU_table/feature-table-metaD-tax.biom')
### Lib2: 
ps <- import_biom('../data/sequence-data/Lib2-CombinedLanes/full-diversity_OTU_table/feature-table-metaD-tax.biom')
ps
### Combined Libs:
ps <- import_biom('../Box/WhitmanLab/Projects/WoodBuffalo/FireSim2019/data/sequence-data/LibCombined/full-diversity_OTU_table/feature-table-metaD-tax.biom')

ps <- import_biom('data-exploration/data/LibCombined/full-diversity_OTU_table/feature-table-metaD-tax.biom')

# Create a new phyloseq object including the master sample list
ps.Lane = phyloseq(otu_table(ps),tax_table(ps),SamList)
ps.Lane

otu_table(ps)[1:10,1:10]
otu_table(ps.Lane)[1:10, 1:10]

# -------
# 1. Look at blanks
ps.blank <- prune_samples(sample_data(ps.Lane)$Blank == 'Yes', ps.Lane)

n_seq = sample_sums(ps.blank)
n_seq = data.frame(names(n_seq),n_seq)
colnames(n_seq)=c("Sample.ID","Total.seq")

# Look at n_seq
min(n_seq$Total.seq)
max(n_seq$Total.seq)
n_seq

#Let's look at this in histogram form  to confirm. 
p = qplot(n_seq$Total.seq, geom="histogram", binwidth=100)+ 
  scale_x_continuous(labels = function(x) sprintf("%g", x))
p

ps.blank.glom <- ps.blank %>%
  tax_glom('Rank2')

df.blank <- psmelt(ps.blank.glom)

p = ggplot(df.blank,aes(x=Sample.ID, y=Abundance, fill=Rank2)) + 
  geom_bar(stat='identity') + 
  coord_flip()+
  facet_wrap(Library~DNA.type, scales = 'free') + 
  theme(axis.text.x = element_text(angle = 90, vjust=0))
p

#write.csv(df.blank, 'data-exploration/data/LibCombined/phyloseq-objects/df.blank')

# -------
# 2. Look at DNeasy extraction. 
### Look at kit comparison: RNeasy vs. DNeasy
ps.Kit <- prune_samples(sample_data(ps.Lane)$Rerun == 'KitComparison', ps.Lane)

ps.Kit.glom <- ps.Kit %>%
  tax_glom('Rank2') %>%
  clean_taxa_names()
ps.Kit.rel <-  transform_sample_counts(ps.Kit.glom, function(x) x / sum(x))

df.Kit <- psmelt(ps.Kit.glom)
df.Kit.rel <- psmelt(ps.Kit.rel)
colnames(df.Kit)

p = ggplot(df.Kit.rel,aes(x=Sample.ID, y=Abundance, fill=Rank2)) + 
  geom_bar(stat='identity') + 
  coord_flip()+
  facet_wrap(~DNA.type, scales = 'free') + 
  theme(axis.text.x = element_text(angle = 90, vjust=0))


# -------
# 3. Compare SI using DNeasy from Lib1 and Lib2
### 06-08-A-SI
### 07-02-O-SI
### 11-03-O-SI
### 19-07-A-SI

# DNeasy run
ps.SI <- prune_samples(sample_data(ps.Lane)$incub.trtmt == 'SI', ps.Lane)

ps.SI.glom <- ps.SI %>%
  tax_glom('Phylum')

ps.SI.rel <- transform_sample_counts(ps.SI.glom, function(x) x / sum(x))

df.SI <- psmelt(ps.SI.rel)

df.SI.compare <- subset(df.SI, core.id.hor.incub %in% c('19UW-WB-19-07-A-SI','19UW-WB-11-03-O-SI',
                                                        '19UW-WB-07-02-O-SI',
                                                        '19UW-WB-06-08-A-SI', 'Blank-Ext-4-SI'))
dim(df.SI.compare)

p = ggplot(df.SI.compare,aes(x=core.id.hor.incub, y=Abundance, fill=Phylum)) + 
  geom_bar(stat='identity') + 
  coord_flip()+
  theme(axis.text.x = element_text(angle = 90, vjust=0))

# RNeasy run:
ps.SI.2 <- prune_samples(sample_data(ps.Lane)$Timepoint == 'SI',ps.Lane)
ps.SI.2 <- merge_samples(ps.SI.2, "Sample.ID")

ps.SI.2.glom <- ps.SI.2 %>%
  tax_glom('Rank2')
ps.SI.2.rel <- transform_sample_counts(ps.SI.2.glom, function(x) x / sum(x))

df.SI.2 <- psmelt(ps.SI.2.rel)

p = ggplot(df.SI.2,aes(x=Sample, y=Abundance, fill=Rank2)) + 
  geom_bar(stat='identity') + 
  coord_flip()+
  theme(axis.text.x = element_text(angle = 90, vjust=0))

# -------
# 4. Compare Redo
# Let's look at just site 16
ps.Redo <- prune_samples(sample_data(ps.Lane)$Rerun == 'REDO',  ps.Lane)

ps.Redo.rel <- transform_sample_counts(ps.Redo, function(x) x / sum(x))
 
ps.Redo.phy <- ps.Redo.rel %>%
  tax_glom('Rank2')  #group taxa at genus level -> I can group by Phylum when plotting later
 
# Normalize reads:
df.Redo <- psmelt(ps.Redo.phy)
colnames(df.Redo)
# Plot sequencing data (
p = ggplot(df.Redo,aes(x=Sample.ID, y=Abundance, fill=Rank2)) + 
  geom_bar(stat='identity') + 
  coord_flip()+
  facet_wrap(~DNA.type, scales = 'free') + 
  theme(axis.text.x = element_text(angle = 90, vjust=0))

# -------
# 5. Calculate num. seq
n_seq = sample_sums(ps.Lane)
n_seq = data.frame(names(n_seq),n_seq)
colnames(n_seq)=c("Sample.ID","Total.seq")

#write.csv(n_seq, 'data-exploration/data/LibCombined/phyloseq-objects/CombinedLanes_n-seq.csv')

# -------
# 6. Look at DNeasy samples (used to account for low abun in RNeasy gDNA)
ps.DNez <- prune_samples(sample_data(ps.Lane)$Rerun == 'DNez',  ps.Lane)

ps.DNez.rel <- transform_sample_counts(ps.DNez, function(x) x / sum(x))

ps.DNez.phy <- ps.DNez.rel %>%
  tax_glom('Rank2')  #group taxa at genus level -> I can group by Phylum when plotting later
# Normalize reads:
df.DNez <- psmelt(ps.DNez.phy)
colnames(df.DNez)
# Plot sequencing data (
p = ggplot(df.DNez,aes(x=Sample.ID, y=Abundance, fill=Rank2)) + 
  geom_bar(stat='identity') + 
  coord_flip()+
  facet_wrap(~DNA.type, scales = 'free') + 
  theme(axis.text.x = element_text(angle = 90, vjust=0))



# -------
# 7. Look at DNeasy + Redo's
ps.DNez.Redo <- prune_samples(sample_data(ps.Lane)$Rerun == 'Dnez-REDO',  ps.Lane)

ps.DNez.Redo.rel <- transform_sample_counts(ps.DNez.Redo, function(x) x / sum(x))

ps.DNez.Redo.phy <- ps.DNez.Redo.rel %>%
  tax_glom('Rank2')  #group taxa at genus level -> I can group by Phylum when plotting later
# Normalize reads:
df.DNez.Redo <- psmelt(ps.DNez.Redo.phy)
colnames(df.DNez.Redo)
# Plot sequencing data (
p = ggplot(df.DNez.Redo,aes(x=Sample.ID, y=Abundance, fill=Rank2)) + 
  geom_bar(stat='identity') + 
  coord_flip()+
  facet_wrap(~DNA.type, scales = 'free') + 
  theme(axis.text.x = element_text(angle = 90, vjust=0))

# -------
# 8. Look at DNeasy + Redo's
ps.2 <- prune_samples(sample_data(ps.Lane)$Rerun == '2',  ps.Lane)

ps.2.rel <- transform_sample_counts(ps.2, function(x) x / sum(x))

ps.2.phy <- ps.2.rel %>%
  tax_glom('Rank2')  #group taxa at genus level -> I can group by Phylum when plotting later
# Normalize reads:
df.2 <- psmelt(ps.2.phy)
colnames(df.2)
# Plot sequencing data (
p = ggplot(df.2,aes(x=Sample.ID, y=Abundance, fill=Rank2)) + 
  geom_bar(stat='identity') + 
  coord_flip()+
  facet_wrap(~DNA.type, scales = 'free') + 
  theme(axis.text.x = element_text(angle = 90, vjust=0))




# -------
# Clean up ps:
### Remove blanks
ps.prune <- prune_samples(sample_data(ps.Lane)$Blank == 'No', ps.Lane)

### Remove lab std
ps.prune <- prune_samples(sample_data(ps.prune)$Timepoint != 'Lab.std', ps.prune)

### Remove Redos/DNez extractions with <1000 seqs
ps.prune <- prune_samples(sample_data(ps.prune)$Status == 'Keep', ps.prune)

ps.prune <- prune_taxa(taxa_sums(ps.prune)>0, ps.prune)

# merge sample reads by sample.id (there are two samples per each Sample.ID)
ps.merge <- merge_samples(ps.prune, "Sample.ID")

# This switches the OTU table columns and rows. Undo this:
otu_table(ps.merge) <- t(otu_table(ps.merge))

# Take a look at OTU table:
otu_table(ps.merge)[1:3,1:3]
otu_table(ps.prune)[1:3,1:3]

ps.merge # 470 samples

# OBJECTIVE: Fixing taxonomy import -----
head(tax_table(ps.merge))
# Yup, looks ugly.

##Fixing the tax_table object to "clean-up" the column names

x = data.frame(tax_table(ps.merge))
# Making a dummy variable to store the taxonomy data

colnames(x) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# Assigning the proper column names instead of SILVA ranks

x$Domain = gsub("d__", "", as.character(x$Domain))
x$Phylum = gsub("p__", "", as.character(x$Phylum))
x$Class = gsub("c__", "", as.character(x$Class))
x$Order = gsub("o__", "", as.character(x$Order))
x$Family = gsub("f__", "", as.character(x$Family))
x$Genus = gsub("g__", "", as.character(x$Genus))
x$Species = gsub("s__", "", as.character(x$Species))
# Substituting the characters we don't want with nothing in the taxonomy

x=tax_table(as.matrix(x,dimnames=list(row.names(x),colnames(x))))
# Turning it into a taxonomy table, while saving the rownames and column names
tax_table(ps.merge)=x
# Reassigning the taxonomy table in ps_xxx to the new modified one

head(tax_table(ps.merge))
# Check for success
# Looks good.



# OBJECTIVE: Add metadata -----

# Read in the mapping file, which has sample data
df.samdat = read.csv("../Box/WhitmanLab/Projects/WoodBuffalo/FireSim2019/data/all-sample-metadata.csv", fileEncoding="UTF-8-BOM")

#Check out your sample data
colnames(df.samdat)
dim(df.samdat)
class(df.samdat$site)
df.samdat$site <- as.character(df.samdat$site)

# Reformat site id's
df.samdat$site[df.samdat$site == '1'] <- '01'
df.samdat$site[df.samdat$site == '2'] <- '02'
df.samdat$site[df.samdat$site == '3'] <- '03'
df.samdat$site[df.samdat$site == '4'] <- '04'
df.samdat$site[df.samdat$site == '5'] <- '05'
df.samdat$site[df.samdat$site == '6'] <- '06'
df.samdat$site[df.samdat$site == '7'] <- '07'
df.samdat$site[df.samdat$site == '8'] <- '08'
df.samdat$site[df.samdat$site == '9'] <- '09'

class(df.samdat$core)
df.samdat$core <- as.character(df.samdat$core)
df.samdat$core[df.samdat$core == '1'] <- '01'
df.samdat$core[df.samdat$core == '2'] <- '02'
df.samdat$core[df.samdat$core == '3'] <- '03'
df.samdat$core[df.samdat$core == '4'] <- '04'
df.samdat$core[df.samdat$core == '5'] <- '05'
df.samdat$core[df.samdat$core == '6'] <- '06'
df.samdat$core[df.samdat$core == '7'] <- '07'
df.samdat$core[df.samdat$core == '8'] <- '08'
df.samdat$core[df.samdat$core == '9'] <- '09'

# Assign degree hours based on horizon so mnrl soil degree hours = Low.degree.C.hours
#   and O horizon degree hours = Mid.degree.C.hours. For cores with all O horizon, 
#   avg degree hrs across mid and low thermocouple.
for (i in 1:nrow(df.samdat)) {
  if (df.samdat$PRE.O.hor.thickness.cm[i] == 10) {
    df.samdat$hrzn.degree.hrs[i] = (df.samdat$Low.degree.C.hours[i] + df.samdat$Mid.degree.C.hours[i])/2
  } else if (df.samdat$horizon[i] == 'A') {
    df.samdat$hrzn.degree.hrs[i] = df.samdat$Low.degree.C.hours[i]
  } else {
    df.samdat$hrzn.degree.hrs[i] = df.samdat$Mid.degree.C.hours[i]
  }
}

# Turn it into phyloseq sample data format
samdat = sample_data(df.samdat)

# Assign the sample names to match the samples' IDs
sample_names(samdat) = samdat$Full.id

# Create a new phyloseq object with all the data
ps.full = phyloseq(otu_table(ps.merge),tax_table(ps.merge),samdat)
ps.full

# Sample clean up: 
### All libraries: Remove core 11-02 (it fell on the ground, so we burned a 
###    fourth core from this site)
ps.sam <- prune_samples(sample_data(ps.full)$core.id != '19UW-WB-11-02', ps.full)
ps.sam
### Lib1 and Combined Libs:  Remove "19UW-WB-08-10-O-SI" because the reads are really low
ps.sam <- prune_samples(sample_data(ps.sam)$Full.id != '19UW-WB-08-10-O-SI', ps.sam)
ps.sam


# OBJECTIVE: Plotting phylogenetic tree -----

tree <- ape::read.tree('../Box/WhitmanLab/Projects/WoodBuffalo/FireSim2019/data/sequence-data/LibCombined/exported-tree/tree.nwk')

tree <- phy_tree(tree)
tax <- tax_table(ps.sam)
otu <- otu_table(ps.sam)
sam <- sample_data(ps.sam)

ps.tree.raw = phyloseq::phyloseq(otu,tax, sam,tree)
ps.tree.raw

# Plot tree
#plot_tree(ps.tree, color='incub.trtmt', label.tips='Rank1',
#          shape='burn.trtmt', size='pH')
#plot_tree(ps.tree, "treeonly", nodeplotblank)

# Clean taxa names now:

# Temporarily remove  samples that had low reads:
ps.tree.norm = transform_sample_counts(ps.tree.raw, function(x) x / sum(x))


# -------
## Lib 1: 
### Create separate -SI and -LI phyloseq objects
ps.raw.SI <- prune_samples(sample_data(ps.tree.raw)$incub.trtmt == 'SI', ps.tree.raw)
ps.raw.LI <- prune_samples(sample_data(ps.tree.raw)$incub.trtmt != 'SI', ps.tree.raw)
ps.norm.SI <- prune_samples(sample_data(ps.tree.norm)$incub.trtmt == 'SI', ps.tree.norm)
ps.norm.LI <- prune_samples(sample_data(ps.tree.norm)$incub.trtmt != 'SI', ps.tree.norm)

# SAVE phyloseq object:
#saveRDS(ps.raw.SI, '../data/sequence-data/Lib1-CombinedLanes/ps.raw.SI')
#saveRDS(ps.raw.LI, '../data/sequence-data/Lib1-CombinedLanes/ps.raw.LI')
#saveRDS(ps.norm.SI, '../data/sequence-data/Lib1-CombinedLanes/ps.norm.SI')
#saveRDS(ps.norm.LI, '../data/sequence-data/Lib1-CombinedLanes/ps.norm.LI')



# -------
## Lib 2: 
### Create separate ps for -pb and -SI samples:
ps.raw.pb <- prune_samples(sample_data(ps.tree.raw)$incub.trtmt == 'pb', ps.tree.raw)
ps.raw.SI.redo <- prune_samples(sample_data(ps.tree.raw)$core.id.hor.incub == '19UW-WB-08-10-O-SI', ps.tree.raw)

ps.norm.pb <- prune_samples(sample_data(ps.tree.norm)$incub.trtmt == 'pb', ps.tree.norm)
ps.norm.SI.redo <- prune_samples(sample_data(ps.tree.norm)$core.id.hor.incub == '19UW-WB-08-10-O-SI', ps.tree.norm)

# SAVE phyloseq object:
#saveRDS(ps.raw.pb, '../data/sequence-data/Lib2-CombinedLanes/ps.raw.pb')
#saveRDS(ps.raw.SI.redo, '../data/sequence-data/Lib2-CombinedLanes/ps.raw.SI.redo')
#saveRDS(ps.norm.pb, '../data/sequence-data/Lib2-CombinedLanes/ps.norm.pb')
#saveRDS(ps.norm.SI.redo, '../data/sequence-data/Lib2-CombinedLanes/ps.norm.SI.redo')


# -------
## LibCombined:
### Create separate ps for -pb, -LI, and -SI samples

#ps.tree.norm <- readRDS('../Box/WhitmanLab/Projects/WoodBuffalo/FireSim2019/data/sequence-data/LibCombined/phyloseq-objects/ps.norm.full')
ps.raw.pb <- prune_samples(sample_data(ps.tree.raw)$incub.trtmt == 'pb', ps.tree.raw)
ps.raw.SI <- prune_samples(sample_data(ps.tree.raw)$incub.trtmt == 'SI', ps.tree.raw)
ps.raw.LI <- prune_samples(sample_data(ps.tree.raw)$incub.trtmt != "SI" &
                             sample_data(ps.tree.raw)$incub.trtmt != 'pb', ps.tree.raw)
ps.raw.full <- ps.tree.raw

ps.norm.pb <- prune_samples(sample_data(ps.tree.norm)$incub.trtmt == 'pb', ps.tree.norm)
ps.norm.SI <- prune_samples(sample_data(ps.tree.norm)$incub.trtmt == 'SI', ps.tree.norm)
ps.norm.LI <- prune_samples(sample_data(ps.tree.norm)$incub.trtmt != "SI" &
                              sample_data(ps.tree.norm)$incub.trtmt != 'pb', ps.tree.norm)
ps.norm.full <- ps.tree.norm


ps.raw.SI <- prune_samples(sample_data(ps.raw.SI)$Full.id != '19UW-WB-06-08-A-SI-duplicate' &
                             sample_data(ps.raw.SI)$Full.id != '19UW-WB-07-02-O-SI-duplicate' &
                             sample_data(ps.raw.SI)$Full.id != '19UW-WB-11-03-O-SI-duplicate'&
                             sample_data(ps.raw.SI)$Full.id != '19UW-WB-19-07-A-SI-duplicate' &
                             sample_data(ps.raw.SI)$Full.id != '19UW-WB-08-10-O-SI', ps.raw.SI)

ps.norm.SI <- prune_samples(sample_data(ps.norm.SI)$Full.id != '19UW-WB-06-08-A-SI-duplicate' &
                              sample_data(ps.norm.SI)$Full.id != '19UW-WB-07-02-O-SI-duplicate' &
                              sample_data(ps.norm.SI)$Full.id != '19UW-WB-11-03-O-SI-duplicate'&
                              sample_data(ps.norm.SI)$Full.id != '19UW-WB-19-07-A-SI-duplicate' &
                              sample_data(ps.norm.SI)$Full.id != '19UW-WB-08-10-O-SI', ps.norm.SI)
# SAVE phyloseq object:
#saveRDS(ps.raw.pb,   '../../../WhitmanLab/Projects/WoodBuffalo/FireSim2019/data/sequence-data/LibCombined/phyloseq-objects/ps.raw.pb')
#saveRDS(ps.raw.SI,   '../Box/WhitmanLab/Projects/WoodBuffalo/FireSim2019/data/sequence-data/LibCombined/phyloseq-objects/ps.raw.SI')
#saveRDS(ps.raw.LI,   '../../../WhitmanLab/Projects/WoodBuffalo/FireSim2019/data/sequence-data/LibCombined/phyloseq-objects/ps.raw.LI')
#saveRDS(ps.raw.full, '../../../WhitmanLab/Projects/WoodBuffalo/FireSim2019/data/sequence-data/LibCombined/phyloseq-objects/ps.raw.full')
#saveRDS(ps.norm.pb,  '../../../WhitmanLab/Projects/WoodBuffalo/FireSim2019/data/sequence-data/LibCombined/phyloseq-objects/ps.norm.pb')
#saveRDS(ps.norm.SI,  '../Box/WhitmanLab/Projects/WoodBuffalo/FireSim2019/data/sequence-data/LibCombined/phyloseq-objects/ps.norm.SI')
#saveRDS(ps.norm.LI,  '../../../WhitmanLab/Projects/WoodBuffalo/FireSim2019/data/sequence-data/LibCombined/phyloseq-objects/ps.norm.LI')
#saveRDS(ps.norm.full,'../../../WhitmanLab/Projects/WoodBuffalo/FireSim2019/data/sequence-data/LibCombined/phyloseq-objects/ps.norm.full')

# Visualize -----
# Take a quick look at how many sequences are in each sample.
n_seq = sample_sums(ps.raw.pb)
n_seq = data.frame(names(n_seq),n_seq)
colnames(n_seq)=c("Sample.ID","Total.seq")

# Look at n_seq
min(n_seq$Total.seq)
tail(n_seq,100)

# Let's look at this in histogram form  to confirm. 
p = qplot(n_seq$Total.seq, geom="histogram", binwidth=100)+ 
  scale_x_continuous(labels = function(x) sprintf("%g", x))
p

# More files to save: 
# Separate gDNA and RNA sequences:
ps.cDNA.raw <- prune_samples(sample_data(ps.raw)$DNA.type == 'cDNA', ps.raw)
ps.gDNA.raw <- prune_samples(sample_data(ps.raw)$DNA.type == 'gDNA', ps.raw)
ps.cDNA.norm <- prune_samples(sample_data(ps.norm)$DNA.type == 'cDNA', ps.norm)
ps.gDNA.norm <- prune_samples(sample_data(ps.norm)$DNA.type == 'gDNA', ps.norm)

#saveRDS(ps.cDNA.raw, '../data/sequence-data/LibCombined/phyloseq-objects/ps.cDNA.raw.pb')
#saveRDS(ps.gDNA.raw, '../data/sequence-data/LibCombined/phyloseq-objects/ps.gDNA.raw.pb')
#saveRDS(ps.cDNA.norm, '../data/sequence-data/LibCombined/phyloseq-objects/ps.cDNA.norm.pb')
#saveRDS(ps.gDNA.norm, '../data/sequence-data/LibCombined/phyloseq-objects/ps.gDNA.norm.pb')

