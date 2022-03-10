##### Does mean predicted 16S copy number increase with degree hours?

# Adding RDP Classifier data
# In this case, I made a .fasta file from the OTU names, which were the sequences (rather than OTU1, OTU2, or whatever)
# One could just use the .fasta file from the sequence processing workflow
# There might be some additional steps required in that case to match the OTU names, but should work.

DNA = Biostrings::DNAStringSet(merged$OTU)
names(DNA)=merged$OTU
Biostrings::writeXStringSet(DNA,"mergedResponders.fasta")

# We already ran rrnDB classifier RDP Classifier version 2.12 https://rrndb.umms.med.umich.edu/estimate/run_classifier
# Using 0.8 cutoff on our .fasta file of OTUs
# Importing the results of that run

# Sequence counts
# This file has what they classified each read as
RDP = read.csv('../../../../Box/WhitmanLab/Projects/WoodBuffalo/FireSim2019/data/sequence-data/LibCombined/full-diversity_OTU_table/dna-sequences.tsv',header=FALSE,sep=";")
head(RDP)

# We want to extract the genus they assigned for each OTU.
# Create function to extract genus, if present
GenusGenerator <- function(taxon) {
  Genus = strsplit(gsub(".*family\\s*|genus*", "", taxon),split="\t")[[1]][2]
  return(Genus)
}

# Extract the genus
RDP$GenusRDP = sapply(RDP$V1,GenusGenerator)
head(RDP)

# Might as well pull out OTU ID to be certain
OTUGenerator <- function(taxon) {
  OTU = strsplit(paste(taxon),split="\t")[[1]][1]
  return(OTU)
}
RDP$OTU = sapply(RDP$V1,OTUGenerator)
head(RDP)

# Ok, we've got what we need.
# Trim it down
RDP = RDP[,c(2,3)]

# Can now pull data from RRNDB.
# Reading in the rrnDB v5.5 file
rrnDB = read.csv('../../../../Box/WhitmanLab/Projects/WoodBuffalo/FireSim2019/data/sequence-data/LibCombined/full-diversity_OTU_table/rrnDB-5.7_pantaxa_stats_RDP.tsv',sep="\t")

head(rrnDB)

# Creating a list of genera in the DB
rrnDBGenera = as.character(rrnDB[rrnDB$rank=="genus",]$name)

# Matching up genus name with mean predicted copy number
for (i in 1:length(RDP$GenusRDP)){
  GenusRDP = paste(RDP$GenusRDP[i])
  CopyNum = ifelse(GenusRDP %in% rrnDBGenera, rrnDB[rrnDB$name==GenusRDP,9],"")
  RDP$CopyNum[i] = CopyNum
}
head(RDP)

write.csv(RDP,'../Box/WhitmanLab/Projects/WoodBuffalo/FireSim2019/data/All-OTUs-CopyNum.csv', row.names = FALSE)

RDP <- read.csv('../../../../Box/WhitmanLab/Projects/WoodBuffalo/FireSim2019/data/All-OTUs-CopyNum.csv')
ps <- readRDS('../../../../Box/WhitmanLab/Projects/WoodBuffalo/FireSim2019/data/sequence-data/LibCombined/phyloseq-objects/ps.raw.full')

#ps.SI <- prune_samples(sample_data(ps)$incub.trtmt == 'SI', ps)
#ps.SI

ps.LIwA <- prune_samples(sample_data(ps)$incub.trtmt == 'LIwA', ps)

ps.LIn <- prune_samples(sample_data(ps)$incub.trtmt == 'LIn', ps)

ps.pb <- prune_samples(sample_data(ps)$incub.trtmt == 'pb' &
                         sample_data(ps)$DNA.type == 'gDNA', ps)
ps.pb

# How should I pick this cutoff? Just setting it at >0 removes over 10000 taxa.
ps.mini <- prune_taxa(taxa_sums(ps) > 0, ps.LIn)
ps.mini

ps.mini.norm <- transform_sample_counts(ps.mini, function(x) x / sum(x))


# Working with melted phyloseq object
mdf = psmelt(ps.mini.norm) 


# Add the rrnDB copy number data to the melted phyloseq object
mdf = plyr::join(mdf,RDP,by="OTU")
mdf$CopyNum = as.numeric(mdf$CopyNum)
mdf$Abundance = as.numeric(mdf$Abundance)

# From Nemergut et al. (2016) - OTU data were then normalized (standardized) for copy number 
# by dividing by copy number. For each sample, we calculated the community aggregated trait value 
# (weighted mean) by taking the product of the estimated operon copy number and the relative abundance 
# for each OTU, and summing this value across all OTUs in a sample. 

# So, first, we divide abundance by copy number
# Then, we re-calculate the relative abundnace, now adjusted for copy number
# The risk there, is, for any organisms without assigned copy numbers, they are excluded from this calculation.
# However, I think we have a pretty good fraction of the community with copy numbers
# To check:

d = mdf %>%
  dplyr::group_by(Sample)%>%
  dplyr::filter(is.na(CopyNum))%>%
  dplyr::summarize(NoCopyNum = sum(Abundance))
hist(d$NoCopyNum)
# Not too bad, but a wide range
# Could assume that unassigned taxa have mean copy number across dataset:
meanCopyNum = RDP%>%
  group_by(GenusRDP,CopyNum)%>%
  summarize(N=n())
meanCopyNum = mean(as.numeric(meanCopyNum$CopyNum),na.rm=TRUE)

# Calculating weighted mean copy numbers (with optional mean replacement; not used):
df = mdf %>%
  #dplyr::mutate(CopyNum = ifelse(is.na(CopyNum) | CopyNum =="NA",meanCopyNum,CopyNum))%>%
  dplyr::filter(!is.na(CopyNum))%>%
  dplyr::mutate(WtAbund = Abundance/CopyNum)%>%
  dplyr::group_by(Sample)%>%
  dplyr::mutate(AdjAbund = WtAbund/sum(WtAbund))%>%
  dplyr::mutate(WtCopyNum = CopyNum*AdjAbund)%>%
  dplyr::group_by(Full.id,core.id, site, 
                  burn.trtmt,
                  incub.trtmt, 
                  PRE.hor.thickness.cm,
                  POST.hor.thickness.cm,
                  PRE.O.hor.thickness.cm,
                  POST.O.hor.thickness.cm,
                  Veg.type,
                  texture,pH,horizon,
                  C.N.ratio,
                  hor.thickness.loss.cm,
                  Mid.degree.C.hours,
                  Low.degree.C.hours)%>%
  dplyr::summarize(WtMeanCopyNum = sum(WtCopyNum,na.rm=TRUE))
hist(df$WtMeanCopyNum)

#write.csv(df, '../../../../Box/WhitmanLab/Projects/WoodBuffalo/FireSim2019/data/WtMeanCopyNum_LIn.csv')
#write.csv(df, '../data/WtMeanCopyNum_SI.csv')

df.CopyNum.SI <- read.csv('../../../../Box/WhitmanLab/Projects/WoodBuffalo/FireSim2019/data/WtMeanCopyNum_SI.csv')
df.CopyNum.LIn <- read.csv('../../../../Box/WhitmanLab/Projects/WoodBuffalo/FireSim2019/data/WtMeanCopyNum_LIn.csv')
df.CopyNum.LIwA <- read.csv('../../../../Box/WhitmanLab/Projects/WoodBuffalo/FireSim2019/data/WtMeanCopyNum_LIwA.csv')
df.CopyNum.pb.gDNA <- read.csv('../../../../Box/WhitmanLab/Projects/WoodBuffalo/FireSim2019/data/WtMeanCopyNum_pb_gDNA.csv')

df.CopyNum.full <- rbind(df.CopyNum.SI, df.CopyNum.LIn, df.CopyNum.LIwA, df.CopyNum.pb.gDNA)

write.csv(df.CopyNum.full, '../../../../Box/WhitmanLab/Projects/WoodBuffalo/FireSim2019/data/WtMeanCopyNum_AllgDNAsamples.csv')
