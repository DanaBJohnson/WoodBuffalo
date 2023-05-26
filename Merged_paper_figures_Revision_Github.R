# Load packages
library(phyloseq)
library(ggplot2)
library(dplyr)

# Setting some plot parameters used throughout
year.labs = c("One year post-fire","Five years post-fire")
names(year.labs) = c(1,5)
horizon.labs = c("Organic","Mineral")

# Importing the phyloseq object with the two years of field data and the lab data
ps = readRDS("ps.1519FireSim.merged")
tree <- ape::read.tree('tree/merged-exported-tree/tree.nwk')

ps <- phyloseq(sample_data(ps), otu_table(ps), tax_table(ps), tree)
ps

# Add column with pairs info
sample_data(ps)$Pairs = paste(sample_data(ps)$Site_ID,sample_data(ps)$Org_or_Min,sep="_")
head(sample_data(ps)$Pairs)

# Normalize sample counts
ps.norm = transform_sample_counts(ps,function(x) x/sum(x))

# Add a column with which study the samples is from
sample_data(ps.norm)$Study = sample_data(ps.norm)$Years_Since_Fire
sample_data(ps.norm)$Study[is.na(sample_data(ps.norm)$Study)]="Lab"

# Drop non-autoclaved lab samples, cDNA, and treed wetland field samples for ordination
sample_data(ps.norm)$incub.trtmt = ifelse(is.na(sample_data(ps.norm)$incub.trtmt), "field",sample_data(ps.norm)$incub.trtmt)
sample_data(ps.norm)$DNA.type[is.na(sample_data(ps.norm)$DNA.type)]="gDNA"
sample_data(ps.norm)$Veg.type = ifelse(is.na(sample_data(ps.norm)$Veg.type), paste(sample_data(ps.norm)$Veg_Comm),sample_data(ps.norm)$Veg.type)

ps.norm = subset_samples(ps.norm,sample_data(ps.norm)$incub.trtmt %in% c("field","SI","LIwA","pb"))
ps.norm = subset_samples(ps.norm,sample_data(ps.norm)$DNA.type == "gDNA")
ps.norm = subset_samples(ps.norm,sample_data(ps.norm)$Veg.type != "Treed Wetland")
ps.norm

# Ordinate plots
ord.nmds = ordinate(ps.norm,method = "NMDS",distance="bray",k=3, trymax=200)
ord.nmds

# Making figure for paper
# Get ordination data
x = data.frame(ord.nmds$points)$MDS1
y = data.frame(ord.nmds$points)$MDS2
z = data.frame(ord.nmds$points)$MDS3

# Add to sample data to plot
full.plot = sample_data(ps.norm)
full.plot$x = x
full.plot$y = y
full.plot$z = z

# Add field ID to burn treatment
full.plot$burn.trtmt = ifelse(is.na(full.plot$burn.trtmt) & full.plot$Burned_Unburned == "Burned","field burned",
                              ifelse(is.na(full.plot$burn.trtmt) & full.plot$Burned_Unburned == "Unburned","field unburned",full.plot$burn.trtmt))

full.plot$burn.trtmt = as.factor(full.plot$burn.trtmt)
full.plot$burn.trtmt = ordered(full.plot$burn.trtmt,levels=c("control","wet","dry","field unburned","field burned"))

# Make horizon variable same for both sets
full.plot$horizon = ifelse(is.na(full.plot$horizon), paste(full.plot$Org_or_Min),full.plot$horizon)
full.plot$horizon[full.plot$horizon=="M"]="A"
full.plot$horizon[full.plot$horizon=="A"]="Mineral"
full.plot$horizon[full.plot$horizon=="O"]="Organic"

# Make new year-agnostic study variable
full.plot$StudyShape = full.plot$Study
full.plot$StudyShape[full.plot$StudyShape %in% c("1","5")]="Field"

# Factor horizon variable
full.plot$horizon <- factor(full.plot$horizon, levels = c('Organic','Mineral'))

# Reformat full plot
full.plot = data.frame(full.plot)

plot.data <- full.plot %>%
  tidyr::unite(burn.trtmt.hor, c(burn.trtmt,horizon), sep='-', remove = FALSE)

# Save statistical source data Figure 3
ssd.3 = plot.data[,c("x","y","z","burn.trtmt","burn.trtmt.hor")]
#write.csv(ssd.3,"../FireSim2019/paper/FinalRevision/Johnson_SourceDSata_Fig3.csv")

# Plot figure showing field and lab samples ordinated together
p = ggplot(plot.data) + theme_bw()
p = p + geom_point(aes(x=x,y=y,shape=burn.trtmt.hor,color=burn.trtmt),size=3, alpha=0.7)
p = p + scale_color_manual(values=c("grey50","orange","red3","grey50","red3"))
p = p + scale_shape_manual(values=c('field burned-Mineral'=0,
                                    'field burned-Organic'=1,
                                    'field unburned-Organic'=1,
                                    'field unburned-Mineral'=0,
                                    'wet-Mineral'=15,
                                    'wet-Organic'=16,
                                    'dry-Mineral'=15,
                                    'dry-Organic'=16,
                                    'control-Mineral'=15,
                                    'control-Organic'=16),
                           labels = c('Field, Mineral horizon',
                                      'Field, Organic horizon',
                                      'Field, Organic horizon',
                                      'Field, Mineral horizon',
                                      'Lab, Mineral horizon',
                                      'Lab, Organic horizon',
                                      'Lab, Mineral horizon',
                                      'Lab, Organic horizon',
                                      'Lab, Mineral horizon',
                                      'Lab, Organic horizon'))
p = p + xlab("NMDS1") + ylab("NMDS2")
p = p + guides(shape=guide_legend(title="Study"),color=guide_legend(title="Burn\nTreatment"))
p


# Our ultimate goal is to bring in trait assignments, match to these taxa, and trace abundance
# with burn severity and time
# Now that we have the OTU IDs, we can actually drop all the lab samples (and their unique taxa).
ps.norm.field = subset_samples(ps.norm,sample_data(ps.norm)$Study != "Lab")
ps.norm.field = prune_taxa(taxa_sums(ps.norm.field)>0,ps.norm.field)
ps.norm.field

mdf = psmelt(ps.norm.field)

# Import responder data, calculated from lab samples using corncob
Responders = read.csv("../FireSim2019/data/sequence-data/LibCombined/corncob-output/Manuscript/Trait-responder-OTUs.csv")

# Get list of survivors
Survivors = Responders %>%
  filter(!is.na(mu.survival.max))

# Get list of fire-survivor OTU IDs
Survivors = unique(Survivors$OTU)

# Add trait classification to full dataframe
mdf$Survivors = ifelse(mdf$taxa_OTU %in% Survivors, "Survivor", "None")

# Pull out data of interest and sum total abundance of survivor taxa for relevant veg types
mdf.Survivors = mdf %>%
  group_by(Sample,Burn_Severity_Index,Severity_Class,Org_or_Min,Veg_Comm,Burned_Unburned,Survivors,Years_Since_Fire,Land_Class)%>%
  summarize(TotalAbundance = sum(Abundance))%>%
  filter(Survivors == "Survivor")%>%
  filter(Veg_Comm %in% c("Jack Pine","Mixedwood","Black Spruce"))

# Get max and run linear models with and without interaction
max(mdf.Survivors$TotalAbundance)
summary(lm(data=mdf.Survivors,TotalAbundance~Years_Since_Fire*Burn_Severity_Index))
summary(lm(data=mdf.Survivors,TotalAbundance~Years_Since_Fire+Burn_Severity_Index))

# Re-plot figures with linear models
survivors.lm = summary(lm(data=mdf.Survivors,TotalAbundance~Years_Since_Fire+Burn_Severity_Index))

# Assign coefficients
Intercept = survivors.lm$coefficients[1,1]
YSF5 = survivors.lm$coefficients[2,1]
BSI = survivors.lm$coefficients[3,1]
# Store coeffs for 1 and 5 years post-fire
# Form = intercept, YSF, BSI
survivors.lm.coeffs = data.frame(t(c("Intercept"=Intercept,"YSF5"=YSF5,"BSI"=BSI)))

Int.1.s=survivors.lm.coeffs$Intercept+survivors.lm.coeffs$BSI*1
Int.5.s=survivors.lm.coeffs$Intercept+survivors.lm.coeffs$YSF5+survivors.lm.coeffs$BSI*1
Slope.s=survivors.lm.coeffs$BSI
Int.1.s
Int.5.s
Slope.s

# This should give us the data frame we need to plot
# Derive the sets of co-ordinates we want to plot and add to our dataset
mdf.Survivors$LinearFit = survivors.lm.coeffs$Intercept+ifelse(mdf.Survivors$Years_Since_Fire=="5",survivors.lm.coeffs$YSF5,0)+(mdf.Survivors$Burn_Severity_Index)*survivors.lm.coeffs$BSI

# Saving statistical source data
ssd.4.A = mdf.Survivors[,c("Burn_Severity_Index","TotalAbundance","Org_or_Min","Years_Since_Fire","LinearFit")]
ssd.4.A$Burn_Severity_Index=ssd.4.A$Burn_Severity_Index-1
#write.csv(ssd.4.A,"../FireSim2019/paper/FinalRevision/Johnson_SourceData_Fig4A.csv")

# Make plot with linear fits
p = ggplot(mdf.Survivors)
p = p + geom_point(aes(x=Burn_Severity_Index-1,y=TotalAbundance,color=Org_or_Min,shape=Org_or_Min),size=2)
p = p + facet_grid(~Years_Since_Fire, labeller = labeller(Years_Since_Fire = year.labs))
p = p + theme_bw()
p = p + ylab("Survivor total read \nrelative abundance")
p = p + xlab("Burn severity index")
p = p + scale_color_manual(values = c("slateblue4","palegreen3"),labels=horizon.labs)
p = p + scale_shape_manual(values = c(16,17),labels=horizon.labs)
p = p + labs(color="Soil horizon",shape="Soil horizon")
p = p + theme(strip.text = element_text(size=14),
              axis.text = element_text(size=14), 
              axis.title.y = element_text(size=14, margin = margin(r=10)),
              axis.title.x = element_text(size=14, margin = margin(t=10)),
              legend.text = element_text(size=12),
              legend.title = element_text(size=12))
p = p + geom_line(aes(x=Burn_Severity_Index-1,y=LinearFit))
p


# Fire Environment 

# Select fire-environment-associated taxa
FireEnv = Responders %>%
  filter(!is.na(mu.affinity.max))

# Get OTU IDs
FireEnv = unique(FireEnv$OTU)

# Add fire environment characteristic to full dataset
mdf$FireEnv = ifelse(mdf$taxa_OTU %in% FireEnv, "FireEnv", "None")

# Select parameters of interest and sum total abundance of fire env responders
mdf.FireEnv = mdf %>%
  group_by(Sample,Burn_Severity_Index,Org_or_Min,Severity_Class,pH,Veg_Comm,TSLF,Burned_Unburned,FireEnv,Years_Since_Fire,Land_Class)%>%
  summarize(TotalAbundance = sum(Abundance))%>%
  filter(FireEnv == "FireEnv")%>%
  filter(Veg_Comm %in% c("Jack Pine","Mixedwood","Black Spruce"))

# Look at max and run linear models with and without interaction
max(mdf.FireEnv$TotalAbundance)
summary(lm(data=mdf.FireEnv,TotalAbundance~Years_Since_Fire*Burn_Severity_Index))
summary(lm(data=mdf.FireEnv,TotalAbundance~Years_Since_Fire+Burn_Severity_Index))

# Save linear model
FireEnv.lm = summary(lm(data=mdf.FireEnv,TotalAbundance~Years_Since_Fire+Burn_Severity_Index))

# Assign coefficients
Intercept = FireEnv.lm$coefficients[1,1]
YSF5 = FireEnv.lm$coefficients[2,1]
BSI = FireEnv.lm$coefficients[3,1]

# Store coeffs for 1 and 5 years post-fire
# Form = intercept, YSF, BSI
FireEnv.lm.coeffs = data.frame(t(c("Intercept"=Intercept,"YSF5"=YSF5,"BSI"=BSI)))

Int.1.fe=FireEnv.lm.coeffs$Intercept+FireEnv.lm.coeffs$BSI*1
Int.5.fe=FireEnv.lm.coeffs$Intercept+FireEnv.lm.coeffs$YSF5+FireEnv.lm.coeffs$BSI*1
Slope.fe=FireEnv.lm.coeffs$BSI
Int.1.fe
Int.5.fe
Slope.fe

# This should give us the data frame we need to plot
# Derive the sets of co-ordinates we want to plot and add to our dataset
mdf.FireEnv$LinearFit = FireEnv.lm.coeffs$Intercept+ifelse(mdf.FireEnv$Years_Since_Fire=="5",FireEnv.lm.coeffs$YSF5,0)+(mdf.FireEnv$Burn_Severity_Index)*FireEnv.lm.coeffs$BSI

# Saving statistical source data
ssd.4.C = mdf.FireEnv[,c("Burn_Severity_Index","TotalAbundance","Org_or_Min","Years_Since_Fire","LinearFit")]
ssd.4.C$Burn_Severity_Index=ssd.4.C$Burn_Severity_Index-1
#write.csv(ssd.4.C,"../FireSim2019/paper/FinalRevision/Johnson_SourceData_Fig4C.csv")


# Plot relative abundance of affinity taxa vs burn severity with linear fits
p = ggplot(mdf.FireEnv)
p = p + geom_point(aes(x=Burn_Severity_Index-1,y=TotalAbundance,color=Org_or_Min,shape=Org_or_Min),size=2)
p = p + facet_grid(~Years_Since_Fire, labeller = labeller(Years_Since_Fire = year.labs))
p = p + theme_bw()
p = p + ylab("Affinity total read \nrelative abundance")
p = p + xlab("Burn severity index")
p = p + scale_color_manual(values = c("slateblue4","palegreen3"),labels=horizon.labs)
p = p + scale_shape_manual(values = c(16,17),labels=horizon.labs)
p = p + labs(color="Soil horizon",shape="Soil horizon")
p = p + theme(strip.text = element_text(size=14),
              axis.text = element_text(size=14), 
              axis.title.y = element_text(size=14, margin = margin(r=10)),
              axis.title.x = element_text(size=14, margin = margin(t=10)),
              legend.text = element_text(size=12),
              legend.title = element_text(size=12))
p = p + geom_line(aes(x=Burn_Severity_Index-1,y=LinearFit))
p = p + scale_y_continuous(breaks=c(0.00,0.01,0.02,0.03,0.04,0.05))
p


# Select fast growing taxa
FastGrowers = Responders %>%
    filter(!is.na(mu.fast.max))

# Get OTU IDs
FastGrowers = unique(FastGrowers$OTU)

# Add trait to main dataset
mdf$FastGrowers = ifelse(mdf$taxa_OTU %in% FastGrowers, "FastGrower", "None")

# Select parameters of interest and sum abundance
mdf.Fast = mdf %>%
  group_by(Sample,Burn_Severity_Index,Severity_Class,Veg_Comm,Org_or_Min,Burned_Unburned,FastGrowers,Years_Since_Fire,Land_Class)%>%
  summarize(TotalAbundance = sum(Abundance))%>%
  filter(FastGrowers == "FastGrower")%>%
  filter(Veg_Comm %in% c("Jack Pine","Mixedwood","Black Spruce"))

# Look at max responder fraction
max(mdf.Fast$TotalAbundance)

# Ensure years since fire is coded as categorical, not numeric
mdf.Fast$Years_Since_Fire <- as.character(mdf.Fast$Years_Since_Fire)

# Run linear model with interaction
summary(lm(data=mdf.Fast,TotalAbundance~Years_Since_Fire*Burn_Severity_Index))

# Plot figure of relative abundance of fast growers vs. BSI and linear models
# Note interaction is significant in this one!
Fast.lm = summary(lm(data=mdf.Fast,TotalAbundance~Years_Since_Fire*(Burn_Severity_Index)))

# Assign coefficients
Intercept = Fast.lm$coefficients[1,1]
YSF5 = Fast.lm$coefficients[2,1]
BSI = Fast.lm$coefficients[3,1]
Interaction = Fast.lm$coefficients[4,1]
# Store coeffs for 1 and 5 years post-fire
# Form = intercept, YSF, BSI
Fast.lm.coeffs = data.frame(t(c("Intercept"=Intercept,"YSF5"=YSF5,"BSI"=BSI,"Interaction"=Interaction)))

Int.1.fg=Fast.lm.coeffs$Intercept+Fast.lm.coeffs$BSI*1
Int.5.fg=Fast.lm.coeffs$Intercept+Fast.lm.coeffs$YSF5+(Fast.lm.coeffs$BSI+Fast.lm.coeffs$Interaction)*1
Slope.1.fg=Fast.lm.coeffs$BSI
Slope.5.fg=Fast.lm.coeffs$BSI+Fast.lm.coeffs$Interaction

Int.1.fg
Int.5.fg
Slope.1.fg
Slope.5.fg

# This should give us the data frame we need to plot
# Derive the sets of co-ordinates we want to plot and add to our dataset
mdf.Fast$LinearFit = Fast.lm.coeffs$Intercept + 
  ifelse(mdf.Fast$Years_Since_Fire=="5",Fast.lm.coeffs$YSF5,0) +
  mdf.Fast$Burn_Severity_Index*ifelse(mdf.Fast$Years_Since_Fire=="5",Fast.lm.coeffs$BSI+Fast.lm.coeffs$Interaction,Fast.lm.coeffs$BSI)

# Saving statistical source data
ssd.4.B = mdf.Fast[,c("Burn_Severity_Index","TotalAbundance","Org_or_Min","Years_Since_Fire","LinearFit")]
ssd.4.B$Burn_Severity_Index=ssd.4.B$Burn_Severity_Index-1
#write.csv(ssd.4.B,"../FireSim2019/paper/FinalRevision/Johnson_SourceData_Fig4B.csv")


# Plot relative abundance of fast growers vs BSI and linear fits
p = ggplot(mdf.Fast)
p = p + geom_point(aes(x=Burn_Severity_Index-1,y=TotalAbundance,color=Org_or_Min,shape=Org_or_Min),size=2)
p = p + facet_wrap(~Years_Since_Fire, labeller = labeller(Years_Since_Fire = year.labs))
p = p + theme_bw() + ylab("Fast growers total abundance")
p = p + ylab("Fast grower total read \nrelative abundance")
p = p + xlab("Burn severity index")
p = p + scale_color_manual(values = c("slateblue4","palegreen3"),labels=horizon.labs)
p = p + scale_shape_manual(values = c(16,17),labels=horizon.labs)
p = p + labs(color="Soil horizon",shape="Soil horizon")
p = p + theme(strip.text = element_text(size=14),
              axis.text = element_text(size=12), 
              axis.title.y = element_text(size=14, margin = margin(r=10)),
              axis.title.x = element_text(size=14),
              legend.text = element_text(size=12),
              legend.title = element_text(size=12), 
              legend.position = 'right')
p = p + geom_line(aes(x=Burn_Severity_Index-1,y=LinearFit))
p

# Calculate mean fraction of community for BSI >=3 (coded as 4)
MeanFractionFast.1 = mdf.Fast %>%
  group_by()%>%
  filter(Burn_Severity_Index>=4 & Years_Since_Fire=="1")%>%
  summarize(Average=mean(TotalAbundance),N=n(),sd=sd(TotalAbundance))
MeanFractionFast.1

MeanFractionFast.5 = mdf.Fast %>%
  group_by()%>%
  filter(Burn_Severity_Index>=4 & Years_Since_Fire=="5")%>%
  summarize(Average=mean(TotalAbundance),N=n(),sd=sd(TotalAbundance))
MeanFractionFast.5


# Bringing all responders together
mdf = mdf %>%
    mutate(Responders = paste(Survivors,FastGrowers,FireEnv))

# Collecting all responders of any kind for burned samples
ResponderSummarySamples = mdf %>%
  filter(Veg_Comm %in% c("Jack Pine","Mixedwood","Black Spruce"))%>%
  filter(Burned_Unburned == "Burned")%>%
  filter(Responders != "None None None")%>%
  group_by(Sample,Severity_Class,Veg_Comm,Burned_Unburned,Years_Since_Fire)%>%
  summarize(TotalAbundance=sum(Abundance))

# Getting stats for all responders summed together
max(ResponderSummarySamples$TotalAbundance)
mean(ResponderSummarySamples$TotalAbundance)
sd(ResponderSummarySamples$TotalAbundance)


ResponderSummarySamplesHighBSI = mdf %>%
  filter(Veg_Comm %in% c("Jack Pine","Mixedwood","Black Spruce"))%>%
  filter(Burned_Unburned == "Burned")%>%
  filter(Responders != "None None None")%>%
  filter(Burn_Severity_Index>=4)%>%
  group_by(Sample,Severity_Class,Veg_Comm,Burned_Unburned,Years_Since_Fire)%>%
  summarize(TotalAbundance=sum(Abundance))
head(ResponderSummarySamples)

# Mean and sd 1-year post-fire for high BSI
mean(ResponderSummarySamplesHighBSI[ResponderSummarySamplesHighBSI$Years_Since_Fire==1,]$TotalAbundance)
sd(ResponderSummarySamplesHighBSI[ResponderSummarySamplesHighBSI$Years_Since_Fire==1,]$TotalAbundance)

# Mean and sd 5-year post-firefor high BSI
mean(ResponderSummarySamplesHighBSI[ResponderSummarySamplesHighBSI$Years_Since_Fire==5,]$TotalAbundance)
sd(ResponderSummarySamplesHighBSI[ResponderSummarySamplesHighBSI$Years_Since_Fire==5,]$TotalAbundance)

# Mean and sd both years post-fire
mean(ResponderSummarySamplesHighBSI$TotalAbundance)
sd(ResponderSummarySamplesHighBSI$TotalAbundance)


### Looking at 16S copy number data ###
# Essentially making sure general patterns still hold after trying
# to adjust by copy number
# Focusing on strongest responders and most likely co-correlations - fast growers

# Bring in RDP estimated copy nubmer data
RDP = read.csv('../FireSim2019/data/All-OTUs-CopyNum.csv')

# Make sure they match our OTU IDs
colnames(RDP)[2] = "taxa_OTU"
  
# Add the rrnDB copy number data to the melted phyloseq object
mdf.RDP = plyr::join(mdf,RDP,by="taxa_OTU")
mdf.RDP$CopyNum = as.numeric(mdf.RDP$CopyNum)
mdf.RDP$Abundance = as.numeric(mdf.RDP$Abundance)

# From Nemergut et al. (2016) - OTU data were then normalized (standardized) for copy number 
# by dividing by copy number. For each sample, we calculated the community aggregated trait value 
# (weighted mean) by taking the product of the estimated operon copy number and the relative abundance 
# for each OTU, and summing this value across all OTUs in a sample. 

# So, first, we divide abundance by copy number
# Then, we re-calculate the relative abundnace, now adjusted for copy number
# The risk there, is, for any organisms without assigned copy numbers, they are excluded from this calculation.
# However, I think we have a pretty good fraction of the community with copy numbers
# To check:

d = mdf.RDP %>%
  dplyr::group_by(Sample)%>%
  dplyr::filter(is.na(CopyNum))%>%
  dplyr::summarize(NoCopyNum = sum(Abundance))
hist(d$NoCopyNum)

# Moderate range

# Calculating mean copy number for unassigned taxa:
meanCopyNum = RDP%>%
  group_by(GenusRDP,CopyNum)%>%
  summarize(N=n())
medianCopyNum = median(as.numeric(meanCopyNum$CopyNum),na.rm=TRUE)
meanCopyNum = mean(as.numeric(meanCopyNum$CopyNum),na.rm=TRUE)

# Calculating weighted mean copy numbers (with mean or median replacement):
df = mdf.RDP %>%
  dplyr::mutate(CopyNum = ifelse(is.na(CopyNum) | CopyNum =="NA",medianCopyNum,CopyNum))%>%
  dplyr::filter(!is.na(CopyNum))%>%
  dplyr::mutate(WtAbund = Abundance/CopyNum)%>%
  dplyr::group_by(Sample)%>%
  dplyr::mutate(AdjAbund = WtAbund/sum(WtAbund))%>%
  dplyr::mutate(WtCopyNum = CopyNum*AdjAbund)%>%
  dplyr::group_by(Sample,Years_Since_Fire,Burn_Severity_Index,Org_or_Min,Veg_Comm)%>%
  dplyr::summarize(WtMeanCopyNum = sum(WtCopyNum,na.rm=TRUE))
hist(df$WtMeanCopyNum)
sum(is.na(df$WtMeanCopyNum))


# Look at fast growers weighted by mean copy number
# Pull out data of interest and sum total abundance of survivor taxa for relevant veg types
mdf.Fast.RDP = mdf.RDP %>%
  mutate(CopyNum = ifelse(is.na(CopyNum) | CopyNum =="NA",meanCopyNum,CopyNum))%>%
  mutate(RDPAbundance = Abundance/CopyNum)%>%
  group_by(Sample,Burn_Severity_Index,Severity_Class,Org_or_Min,Veg_Comm,Burned_Unburned,FastGrowers,Years_Since_Fire,Land_Class)%>%
  summarize(TotalAbundanceRDP = sum(RDPAbundance))%>%
  filter(FastGrowers == "FastGrower")%>%
  filter(Veg_Comm %in% c("Jack Pine","Mixedwood","Black Spruce"))

# Plot burn severity index vs. total read abundance
p = ggplot(mdf.Fast.RDP)
p = p + geom_point(aes(x=Burn_Severity_Index-1,y=TotalAbundanceRDP,color=Org_or_Min,shape=Org_or_Min),size=2)
p = p + facet_grid(~Years_Since_Fire, labeller = labeller(Years_Since_Fire = year.labs))
p = p + theme_bw()
p = p + ylab("Fast grower total read relative abundance\nweighted by copy number")
p = p + xlab("Burn severity index")
p = p + scale_color_manual(values = c("black","grey40"),labels=horizon.labs)
p = p + scale_shape_manual(values = c(16,17),labels=horizon.labs)
p = p + labs(color="Soil horizon",shape="Soil horizon")
p

summary(lm(data=mdf.Fast.RDP[mdf.Fast.RDP$Years_Since_Fire=="1",],TotalAbundanceRDP~Burn_Severity_Index))

# Does the weighted mean copy number also correspond to the fraction of fast growers?
# Want to combine WtMeanCopyNum with fraction of fast growers at sample level
FractFast = mdf%>%
  filter(FastGrowers=="FastGrower")%>%
  group_by(Sample)%>%
  summarize(RelabundFastGrowers = sum(Abundance))

# Merge dataframes
df.plot = merge(df,FractFast,by="Sample")

# Focus on one year post-fire
df.plot = df.plot%>%
  filter(Years_Since_Fire==1)

### rrnDB Result Plot ###
p = ggplot(df.plot)
p = p + geom_point(aes(x=WtMeanCopyNum,y=RelabundFastGrowers,shape=Org_or_Min,colour=Org_or_Min))
p = p + xlab("Weighted Mean Predicted 16S Copy Number")
p = p + ylab("Relative Abundance of Lab-identified Fast Growers")
p = p + theme_bw()
p = p + scale_color_manual(values = c("slateblue4","palegreen3"),labels=horizon.labs)
p = p + scale_shape_manual(values = c(16,17),labels=horizon.labs)
p = p + labs(color="Soil horizon",
             shape="Soil horizon")
p = p + theme(axis.text.y = element_text(size = 14),
              legend.title = element_text(size=14),
              axis.text.x = element_text(size = 13),
              axis.title.x = element_text(size = 14, margin = margin(t=10, b=10)),
              axis.title.y = element_text(size = 14, margin = margin(r=10, l=10)),
              strip.text = element_text(size=14),
              legend.text = element_text(size=14),
              legend.position = '')
p

fast.vs.copynum.stats = summary(lm(data=df.plot,formula=RelabundFastGrowers~WtMeanCopyNum))

# Assign coefficients
Intercept = fast.vs.copynum.stats$coefficients[1,1]
Slope = fast.vs.copynum.stats$coefficients[2,1]

# Store coeffs for 1 and 5 years post-fire
# Form = intercept, Slope

# Add slope to plot
p = p + geom_abline(aes(intercept=Intercept,slope=Slope))
p

# However, one could argue of course these are correlated. What if we plot the
# rrnDB-adjusted lab fast growers relabund?
FractFast.RDP = mdf.Fast.RDP%>%
  group_by()%>%
    select(Sample,TotalAbundanceRDP)

# Merge dataframes
df.RDP.plot = merge(df,FractFast.RDP,by="Sample")

# Select one year post-fire
df.RDP.plot = df.RDP.plot%>%
  filter(Years_Since_Fire==1)%>%
  filter(Veg_Comm %in% c("Jack Pine","Mixedwood","Black Spruce"))

# Saving statistical source data
ssd.ED6 = df.RDP.plot[,c("WtMeanCopyNum","TotalAbundanceRDP","Org_or_Min")]
#write.csv(ssd.ED6,"../FireSim2019/paper/FinalRevision/Johnson_SourceData_FigED6.csv")

### rrnDB Result Plot, accounting for 16S copy number in relabund ###
p = ggplot(df.RDP.plot)
p = p + geom_point(aes(x=WtMeanCopyNum,y=TotalAbundanceRDP,shape=Org_or_Min,colour=Org_or_Min), size=2)
p = p + xlab("Weighted Mean Predicted 16S Copy Number")
p = p + ylab("Relative Abundance of Fast Growers\nAdjusted for Predicted Copy Number")
p = p + theme_bw()
p = p + scale_color_manual(values = c("slateblue4","palegreen3"),labels=horizon.labs)
p = p + scale_shape_manual(values = c(16,17),labels=horizon.labs)
p = p + labs(color="Soil horizon",
             shape="Soil horizon")
p = p + theme(axis.text.y = element_text(size = 12),
              legend.title = element_text(size=12),
              axis.text.x = element_text(size = 13),
              axis.title.x = element_text(size = 12, margin = margin(t=10, b=10)),
              axis.title.y = element_text(size = 12, margin = margin(r=10, l=10)),
              strip.text = element_text(size=12),
              legend.text = element_text(size=12),
              legend.position = '', 
              plot.margin = margin(1,1,1,1,'cm'))
p

fast.vs.copynum.stats = summary(lm(data=df.RDP.plot,formula=TotalAbundanceRDP~WtMeanCopyNum))

# Assign coefficients
Intercept = fast.vs.copynum.stats$coefficients[1,1]
Slope = fast.vs.copynum.stats$coefficients[2,1]

# Store coeffs for 1 and 5 years post-fire
# Form = intercept, Slope

# Add line to plot
p = p + geom_abline(aes(intercept=Intercept,slope=Slope))
p

## Final figure
fast.vs.copynum.stats 

##############################
# We also want to compare BC-dissim for burned vs. unburned (same veg)
# from the field sampels against total relabund of trait-responders

# Load additional packages
library(reshape2)
library(vegan)
library(reshape)

# We're probably interested in seeing how the same site in year 1 differs from that site in year 5
ps = ps.norm.field

# Establish a Paired variable
sample_data(ps)$Pairs = paste(sample_data(ps)$Site_ID,sample_data(ps)$Org_or_Min,sep="_")

# Calculate dissimilarities
Dist.mb = as.matrix(vegan::vegdist(otu_table(ps), method="bray", type="samples"))
Dist.mb[upper.tri(Dist.mb)] = NA
df = data.frame(melt(Dist.mb))

#Generates a dataframe with each contrast and the dissimilarity for the mb comm
colnames(df) = c("Sample_ID_1","Sample_ID_2","Mb_dist")
df = df[df$Mb_dist!=0,]
df = df[!is.na(df$Mb_dist),]

# Making a matrix with one data entry for each site
SamDat = sample_data(ps)
SamDat$Sample_ID = row.names(SamDat)
SamDat$VegSoilYear = paste(SamDat$Veg_Comm,SamDat$Org_or_Min,SamDat$Years_Since_Fire)
SamDat$Burned_Unburned

# Matching same site (different years)
for (i in df$Sample_ID_1){
  df$Pairs_1[df$Sample_ID_1==i] = paste(SamDat$Pairs[SamDat$Sample_ID==i])
}
for (i in df$Sample_ID_2){
  df$Pairs_2[df$Sample_ID_2==i] = paste(SamDat$Pairs[SamDat$Sample_ID==i])
}
df$Pairs = ifelse(df$Pairs_1==df$Pairs_2,df$Pairs_1,"Different")


# Matching Veg and Org/Min for same year
for (i in df$Sample_ID_1){
  df$VegSoilYear_1[df$Sample_ID_1==i] = paste(SamDat$VegSoilYear[SamDat$Sample_ID==i])
}
for (i in df$Sample_ID_2){
  df$VegSoilYear_2[df$Sample_ID_2==i] = paste(SamDat$VegSoilYear[SamDat$Sample_ID==i])
}
df$VegSoilYear = ifelse(df$VegSoilYear_1==df$VegSoilYear_2,df$VegSoilYear_1,"Different")

# Matching Burned / unburned
for (i in df$Sample_ID_1){
  df$BurnedUnburned_1[df$Sample_ID_1==i] = paste(SamDat$Burned_Unburned[SamDat$Sample_ID==i])
}
for (i in df$Sample_ID_2){
  df$BurnedUnburned_2[df$Sample_ID_2==i] = paste(SamDat$Burned_Unburned[SamDat$Sample_ID==i])
}
df$BurnedUnburned = ifelse(df$BurnedUnburned_1==df$BurnedUnburned_2,df$BurnedUnburned_1,"Different")


# Selecting only comparisons between unburned and burned sites, 
# same horizon, same veg comm, same year
df = df[df$VegSoilYear != "Different",]
df = df[df$BurnedUnburned == "Different",]

# Remove the treed wetlands
TW = c("Treed Wetland O 5","Treed Wetland M 5","Treed Wetland O 1","Treed Wetland M 1")
df = df[!(df$VegSoilYear %in% TW),]

# Now we have a matrix with dissim indices for ...
# sites that were burned vs an unburned site,
# within the same year
# matched by veg type and horizon
# for the veg types relevant to this study
# We want to pull out the sample ID for the burned site

df$BurnedSiteID = ifelse(df$BurnedUnburned_1=="Burned",paste(df$Sample_ID_1),paste(df$Sample_ID_2))

# Add the data and clean it up
df = merge(df,data.frame(SamDat),by.x="BurnedSiteID",by.y="Sample_ID")%>%
  select(BurnedSiteID,Burned_Unburned,Veg_Comm,Burn_Severity_Index,Org_or_Min,Mb_dist)

# Want to summarize mean dist to unburned sites for each site
df = df%>%
  group_by(BurnedSiteID,Org_or_Min)%>%
  summarize(MeanBCDissimToUnburned=mean(Mb_dist))

# Now want to match these with the fraction of comm that are responders
FinalMatch = merge(df,ResponderSummarySamples,by.x="BurnedSiteID",by.y="Sample")

# Plot
p = ggplot(FinalMatch)
p = p + geom_point(aes(x=MeanBCDissimToUnburned,y=TotalAbundance,color=Org_or_Min, shape = Org_or_Min))
p = p + facet_grid(~Years_Since_Fire, labeller = labeller(Years_Since_Fire = year.labs))
p = p + scale_color_manual(values = c("slateblue4","palegreen3"),labels=horizon.labs)
p = p + scale_shape_manual(values = c(16,17),labels=horizon.labs)
p = p + labs(color="Soil horizon",
             shape="Soil horizon",
             x='Mean Bray-Curtis dissimilarity',
             y='Total relative abundance of traits')
p = p + theme_bw()
p = p + theme(axis.title.x = element_text(margin = margin(t=10, b=10)))
p

# Linear model
summary(lm(data=FinalMatch,TotalAbundance~Years_Since_Fire*MeanBCDissimToUnburned))

# Checking N
Enns.FinalMatch = FinalMatch %>%
  group_by(Years_Since_Fire,Org_or_Min)%>%
  summarize(N=n())
Enns.FinalMatch

# Ensure years since fire is categorical
FinalMatch$Years_Since_Fire <- as.character(FinalMatch$Years_Since_Fire)

# Re-plot figures with linear models
# Note interaction is significant.
BC.lm = summary(lm(data=FinalMatch,TotalAbundance~Years_Since_Fire*MeanBCDissimToUnburned))

# Assign coefficients
Intercept = BC.lm$coefficients[1,1]
YSF5 = BC.lm$coefficients[2,1]
BC = BC.lm$coefficients[3,1]
Interaction = BC.lm$coefficients[4,1]

# Store coeffs for 1 and 5 years post-fire
# Form = intercept, YSF, BC
BC.lm.coeffs = data.frame(t(c("Intercept"=Intercept,"YSF5"=YSF5,"BC"=BC,"Interaction"=Interaction)))
BC.lm.coeffs

# This should give us the data frame we need to plot
# Derive the sets of co-ordinates we want to plot and add to our dataset
FinalMatch$LinearFit = BC.lm.coeffs$Intercept + 
  ifelse(FinalMatch$Years_Since_Fire=="5",BC.lm.coeffs$YSF5,0) +
  FinalMatch$MeanBCDissimToUnburned*ifelse(FinalMatch$Years_Since_Fire=="5",BC.lm.coeffs$BC+BC.lm.coeffs$Interaction,BC.lm.coeffs$BC)

# Saving statistical source data
ssd.5 = FinalMatch[,c("MeanBCDissimToUnburned","TotalAbundance","Org_or_Min", "Years_Since_Fire","LinearFit")]
#write.csv(ssd.5,"../FireSim2019/paper/FinalRevision/Johnson_SourceData_Fig5.csv")

# Add lines to plot
p = ggplot(FinalMatch,aes(x=MeanBCDissimToUnburned,y=TotalAbundance))
p = p + geom_point(aes(color=Org_or_Min, shape = Org_or_Min))
p = p + facet_grid(~Years_Since_Fire, labeller = labeller(Years_Since_Fire = year.labs))
p = p + scale_color_manual(values = c("slateblue4","palegreen3"),labels=horizon.labs)
p = p + scale_shape_manual(values = c(16,17),labels=horizon.labs)
p = p + labs(color="Soil horizon",
             shape="Soil horizon",
             x='Mean Bray-Curtis dissimilarity',
             y='Total relative abundance')
p = p + theme_bw()
p = p + theme(axis.text.y = element_text(size = 12),
      legend.title = element_text(size=12),
      axis.text.x = element_text(size = 13),
      axis.title.x = element_text(size = 12, margin = margin(t=10, b=10)),
      axis.title.y = element_text(size = 12, margin = margin(r=10, l=10)),
      strip.text = element_text(size=12),
      legend.text = element_text(size=12),
      legend.position = '')
p = p + geom_line(aes(x=MeanBCDissimToUnburned,y=LinearFit))
p
