# --------------------
# Title: incubations-microbial-respiration.R
# Author: Dana Johnson
# Date: 2020-Nov-21
#
# Purpose: This is code to plot microbial respiration data from short and long term
# incubation of FireSim2019 cores. Similar to ec-alldays.R 
#
# Output: p1 <- day-vs-respiration-rate-SI.png
#         p2 <- day-vs-respiration-rate-LI.png

# --------------------
# Load packages
library(dplyr)
library(ggplot2)
library(tidyr)
# --------------------
setwd("C:/Users/danab/Box/WhitmanLab/Projects/WoodBuffalo/FireSim2019/data")


# --------------------
# OBJECTIVE: Calculate initial total C for each core

# Import C data and clean up dataframe
dfCN <- read.csv("raw-data/CN-data.csv")
dfCN <- dfCN[1:94, 1:6]
head(dfCN)

colnames(dfCN)[1] <- 'core.id.full'
colnames(dfCN)[4] <- 'mass.mg'
colnames(dfCN)[5] <- 'percent.total.N'
colnames(dfCN)[6] <- 'percent.total.C'

dfCN <- dfCN %>%
  separate(core.id.full, c("year", "project", "site", "core", "horizon"), sep = "-", remove = FALSE) %>%
  unite(core.id, c(year, project, site, core), sep = "-", remove = FALSE) 
head(dfCN)


# I need to get the initial dry mass numbers for O and mnrl horizons. Import the incubation set-up data
df_initial <- read.csv("raw-data/incubation-set-up.csv")
head(df_initial)
colnames(df_initial)

# Fix naming error on first column
colnames(df_initial)[1] <- 'core.id.hor'

# Pull out relevant dry mass info
df_initial <- subset(df_initial, select = c(core.id.hor,
                                            Horizon,
                                            SI.mass.dry.soil.g,
                                            LIwA.mass.dry.soil.g,
                                            LIn.mass.dry.soil.g,
                                            LIwA.mass.dry.inoculum..to.add.g,
                                            LIn.mass.dry.inoculum.to.add.g))
head(df_initial)

df_initial <- df_initial %>%
  separate(core.id.hor, c('year','project','site','core','hor'), sep = '-', remove=FALSE) %>%
  unite(core.id, c(year,project,site,core), sep='-',remove=TRUE)
head(df_initial)

df_initial$SI.mass.dry.soil.g <- as.numeric(as.character(df_initial$SI.mass.dry.soil.g))
df_initial$LIwA.mass.dry.soil.g <- as.numeric(as.character(df_initial$LIwA.mass.dry.soil.g))
df_initial$LIwA.mass.dry.inoculum..to.add.g <- as.numeric(as.character(df_initial$LIwA.mass.dry.inoculum..to.add.g))
df_initial$LIn.mass.dry.soil.g <- as.numeric(as.character(df_initial$LIn.mass.dry.soil.g))
df_initial$LIn.mass.dry.inoculum.to.add.g <- as.numeric(as.character(df_initial$LIn.mass.dry.inoculum.to.add.g))



# Calculate total dry mass in each jar - need to separate SI from LIwA and LIn. LIwA and LIn had additional soil added (inoculum)
df_SI <- subset(df_initial, select = c(core.id.hor,
                                       core.id,
                                       Horizon,
                                       SI.mass.dry.soil.g))
df_LIwA <- subset(df_initial, select = c(core.id.hor,
                                         core.id,
                                         Horizon,
                                         LIwA.mass.dry.soil.g,
                                         LIwA.mass.dry.inoculum..to.add.g))
df_LIn <- subset(df_initial, select = c(core.id.hor,
                                        core.id,
                                        Horizon,
                                        LIn.mass.dry.soil.g,
                                        LIn.mass.dry.inoculum.to.add.g))


# Format these dataframes so each horizon gets its own row: core.id, horizon, initial.dry.mass.g
# SI
head(df_SI)
df_SI_O <- subset(df_SI, Horizon == 'O')
df_SI_A <- subset(df_SI, Horizon == 'A')

colnames(df_SI_O)[4] <- 'O.dry.mass.g'
df_SI_O$incub <- 'SI'
df_SI_O <- df_SI_O %>%
  unite(core.id.incub, c(core.id,incub), sep='-', remove=FALSE)
head(df_SI_O)

colnames(df_SI_A)[4] <- 'A.dry.mass.g'
df_SI_A$incub <- 'SI'
df_SI_A <- df_SI_A %>%
  unite(core.id.incub, c(core.id,incub), sep='-', remove=FALSE)
head(df_SI_A)

df_SI_O <- subset(df_SI_O, select = c(core.id,
                                      O.dry.mass.g))
df_SI_A <- subset(df_SI_A, select = c(core.id,
                                      A.dry.mass.g))

SI <-merge(df_SI_A, df_SI_O, by = 'core.id')
SI$incub <- 'SI'
head(SI)

#LIwA
LIwA_O <- subset(df_LIwA, Horizon =='O')
LIwA_A <- subset(df_LIwA, Horizon =='A')

head(LIwA_O)
class(LIwA_O$O.dry.mass.g)
colnames(LIwA_O)[4] <- 'O.dry.mass.g'
colnames(LIwA_O)[5] <- 'O.inoculum.dry.mass.g'
LIwA_O$incub <- "LIwA"
LIwA_O <- LIwA_O %>%
  mutate(total.O.dry.mass.g = O.dry.mass.g + O.inoculum.dry.mass.g, na.rm=TRUE) %>%
  unite(core.id.incub, c(core.id,incub), sep='-', remove=FALSE)
head(LIwA_O)

head(LIwA_A)
colnames(LIwA_A)[4] <- 'A.dry.mass.g'
colnames(LIwA_A)[5] <- 'A.inoculum.dry.mass.g'
LIwA_A <- LIwA_A %>%
  mutate(incub = 'LIwA') %>%
  mutate(total.A.dry.mass.g = A.dry.mass.g + A.inoculum.dry.mass.g, na.rm=TRUE) %>%
  unite(core.id.incub, c(core.id, incub), sep='-', remove=FALSE)
head(LIwA_A)

LIwA_O <- subset(LIwA_O, select = c(core.id,
                                    O.dry.mass.g,
                                    O.inoculum.dry.mass.g,
                                    total.O.dry.mass.g))
LIwA_A <- subset(LIwA_A, select = c(core.id,
                                    A.dry.mass.g,
                                    A.inoculum.dry.mass.g,
                                    total.A.dry.mass.g))
LIwA <- merge(LIwA_O, LIwA_A, by ='core.id')
LIwA$incub <- 'LIwA'
head(LIwA)



# LIn
LIn_O <- subset(df_LIn, Horizon =='O')
LIn_A <- subset(df_LIn, Horizon =='A')
head(LIn_O)

colnames(LIn_O)[4] <- 'O.dry.mass.g'
colnames(LIn_O)[5] <- 'O.inoculum.dry.mass.g'

LIn_O$incub <- "LIn"
LIn_O <- LIn_O %>%
  mutate(total.O.dry.mass.g = O.dry.mass.g + O.inoculum.dry.mass.g, na.rm=TRUE) %>%
  unite(core.id.incub, c(core.id,incub), sep='-', remove=FALSE)
head(LIn_O)

head(LIn_A)
colnames(LIn_A)[4] <- 'A.dry.mass.g'
colnames(LIn_A)[5] <- 'A.inoculum.dry.mass.g'
LIn_A <- LIn_A %>%
  mutate(incub = 'LIn') %>%
  mutate(total.A.dry.mass.g = A.dry.mass.g + A.inoculum.dry.mass.g, na.rm=TRUE) %>%
  unite(core.id.incub, c(core.id, incub), sep='-', remove=FALSE)
head(LIn_A)

LIn_O <- subset(LIn_O, select = c(core.id,
                                    O.dry.mass.g,
                                    O.inoculum.dry.mass.g,
                                    total.O.dry.mass.g))
LIn_A <- subset(LIn_A, select = c(core.id,
                                    A.dry.mass.g,
                                    A.inoculum.dry.mass.g,
                                    total.A.dry.mass.g))

LIn <- merge(LIn_O, LIn_A, by = 'core.id')
LIn$incub <- 'LIn'
head(LIn)



# Take a look at SI, LIwA, and LIn and bring in C:N data
head(SI)
head(LIwA)
head(LIn)

head(dfCN)


# Split CN dataframe by horizon
dfCN_O <- subset(dfCN, horizon == 'O')
dfCN_A <- subset(dfCN, horizon == 'A')

colnames(dfCN_O)[12] <- 'O.percent.total.C'
colnames(dfCN_A)[12] <- 'A.percent.total.C'

dfCN_O <- subset(dfCN_O, select = c(core.id, O.percent.total.C))
dfCN_A <- subset(dfCN_A, select = c(core.id, A.percent.total.C))

# Here's my problem: I'm merging the dataframes and losing data of O horizons
CN <- merge(dfCN_O, dfCN_A, by = 'core.id', all.x = TRUE)
head(CN)




# Now combine SI and C:N data. Calculate dry mass C by horizon.
SI <- merge(SI, CN, by = 'core.id')
head(SI)

SI <- SI %>% 
  mutate(O.dry.mass.C.g=O.dry.mass.g*O.percent.total.C/100) %>%
  mutate(A.dry.mass.C.g=A.dry.mass.g*A.percent.total.C/100) 

# Combine intitial dry mass of C from O and A horizons:
SI <- SI %>%
  group_by(core.id) %>%
  mutate(total.initial.dry.mass.C.g = sum(O.dry.mass.C.g,A.dry.mass.C.g,na.rm=TRUE))
head(SI)





# Make control dataframe
dfSite <- read.csv("raw-data/sampling-site-characteristics.csv")

names(dfSite)[1] <- 'core.id'

dfControl <- merge(CN, subset(dfSite, select = c(core.id, Site, Core.trtmt))) %>%
  subset(Core.trtmt == 'control') %>%
  mutate(O.percent.total.C.inoculum = O.percent.total.C,
         A.percent.total.C.inoculum = A.percent.total.C) %>%
  subset(select = c(Site, O.percent.total.C.inoculum, A.percent.total.C.inoculum))

LIwA <- merge(LIwA, subset(dfSite, select = c(core.id, Site)))
        
LIwA <- merge(LIwA, dfControl, by = 'Site', all =TRUE)

LIwA <- merge(LIwA, CN, by = 'core.id')
head(LIwA)

LIwA <- LIwA %>% 
  mutate(O.dry.mass.C.g = O.dry.mass.g * O.percent.total.C/100) %>%
  mutate(A.dry.mass.C.g = A.dry.mass.g * A.percent.total.C/100) %>%
  mutate(O.inoculum.C.g = O.inoculum.dry.mass.g * O.percent.total.C.inoculum/100) %>%
  mutate(A.inoculum.C.g = A.inoculum.dry.mass.g * A.percent.total.C.inoculum/100)

# Combine intitial dry mass of C from O and A horizons:
LIwA <- LIwA %>%
  group_by(core.id) %>%
  mutate(total.initial.dry.mass.C.g = sum(O.dry.mass.C.g, A.dry.mass.C.g, O.inoculum.C.g, A.inoculum.C.g, na.rm = TRUE))

head(LIwA)
colnames(LIwA)





LIn <- merge(LIn, subset(dfSite, select = c(core.id, Site)))

LIn <- merge(LIn, dfControl, by = 'Site', all =TRUE)

LIn <- merge(LIn, CN, by = 'core.id')
head(LIn)

LIn <- LIn %>% 
  mutate(O.dry.mass.C.g = O.dry.mass.g * O.percent.total.C/100) %>%
  mutate(A.dry.mass.C.g = A.dry.mass.g * A.percent.total.C/100) %>%
  mutate(O.inoculum.C.g = O.inoculum.dry.mass.g * O.percent.total.C.inoculum/100) %>%
  mutate(A.inoculum.C.g = A.inoculum.dry.mass.g * A.percent.total.C.inoculum/100)

# Combine intitial dry mass of C from O and A horizons:
LIn <- LIn %>%
  group_by(core.id) %>%
  mutate(total.initial.dry.mass.C.g = sum(O.dry.mass.C.g, A.dry.mass.C.g,O.inoculum.C.g,A.inoculum.C.g, na.rm = TRUE))

head(LIn)
colnames(LIn)



# Combine SI, LIwA, and LIn
SI <- subset(SI, select = c(core.id,
                            incub,
                            total.initial.dry.mass.C.g))
LIwA <- subset(LIwA, select = c(core.id,
                                incub,
                                total.initial.dry.mass.C.g))
LIn <- subset(LIn, select = c(core.id,
                              incub,
                              total.initial.dry.mass.C.g))


SI <- SI %>%
  unite(core.id.incub, c(core.id, incub),sep='-', remove=TRUE) 
LIwA <- LIwA %>%
  unite(core.id.incub, c(core.id,incub),sep='-', remove=TRUE) 
LIn <- LIn %>%
  unite(core.id.incub, c(core.id, incub),sep='-', remove=TRUE) 

df_initialC <- rbind(SI, LIwA)
df_initialC <- rbind(df_initialC, LIn)
# Great. I've got a dataframe with ID's and total initial dry C mass
head(df_initialC)







# Import raw conductivity data
df = read.csv("raw-data/compiled-conductivity-data.csv")
head(df)

colnames(df)[1] <- 'core.id'

df <- df %>%
  unite(core.id.incub, c(core.id, Incubation.trtmt), sep = "-", remove = FALSE)

df <- subset(df, select = c(core.id, core.id.incub, Burn.trtmt, Incubation.trtmt, Dry.mass.g,
                            Core.set, Respirer, Day, Base.trap.M, Initial.cond.mS.cm,
                            Final.cond.mS.cm, Half.M.mS.cm, KOH.trap.time.period, 
                            KOH.trap.time.period.exact, Exact.day))
head(df)

df <- df %>%
  separate(core.id, c('year','project','site','core'), sep = '-', remove = FALSE)

# Removing sample that fell on floor post burn treatment
df <- subset(df, core.id != "19UW-WB-11-02")
head(df)

# I believe that I incorrectly entered 62.20 instead of 162.20 for this data point.
df$Final.cond.mS.cm[df$Final.cond.mS.cm == 62.20] <- 162.20

df <- merge(df, df_initialC, by = 'core.id.incub', all = TRUE)

head(df)



# ----------------
# OBJECTIVE: Calculate CO2 absorbed by KOH trap. 
# Make new column that is difference between ref and measurement
# Group by measurement/setup date
# Make new column where we get the reference blank value (mean of all blanks for that date)
# Make new column where we subtract the relevant blank mean for that date


# Calculate CO2 absorbtion in mg using [CO2 absorbed = P*V*C*R*M]
# P = (Change in conductivity)/(Conductivity at starting concentration - conductivity at 1/2 of starting concentration)
# V = Volume (mL) of KOH trap
# C = Concentration of KOH solution
# R = Ratio
# M = Molecular weight of C
# KOH base trap capacity = V x C x R x M
V = 0.015 # L
C = 0.5   # mol/L
R = 0.5   # Ratio: 1 mol KOH can absorbed 1/2 mol CO2
M = 12.01 # g/mol

# Calculate C-CO2 absorbed by KOH trap.

df = df %>%
  mutate(P = (Initial.cond.mS.cm-Final.cond.mS.cm)/(Initial.cond.mS.cm-Half.M.mS.cm)) %>% # Half.M.mS.cm = molarity of KOH trap at saturation
  mutate(C.CO2.absorbed.g = P*V*C*R*M)

####### CALCULATING ALTERNATIVE C respired ######
# Smirnova et al. 2014 - Calibration of CO2 trapping in alkaline solutions during soil incubation
# Meyer et al. 2018 - Soil respiration and its temp. sensitivity 
#A = 173 # constant: A = 273.75 mg CO2 for 15 mL of 0.5 M KOH solution, 
#df = df %>%
#  mutate(alternative = A *(Initial.cond.mS.cm-Final.cond.mS.cm)/Initial.cond.mS.cm * (12.01)/(44.01*1000)) %>%
#  #mutate(C.CO2.absorbed.g = alternative) %>%
#  mutate(test = (Initial.cond.mS.cm-Final.cond.mS.cm)/Initial.cond.mS.cm)




colnames(df)

# Adjust for C-CO2 absorbed minus blank.
df = df %>%
  unite(msm_ID, c(Respirer,Day), remove = FALSE) %>%
  group_by(msm_ID) %>%
  arrange(msm_ID , .by_group = TRUE) %>%
  mutate(Blank.C.CO2.absorbed.g = mean(C.CO2.absorbed.g[Incubation.trtmt=="Blank"])) 

df = df %>%
  # Normalize C-CO2 absorbtion by gram of initial dry soil
  mutate(C.CO2.absorbed.minus.blank.g = C.CO2.absorbed.g - Blank.C.CO2.absorbed.g) %>%
  mutate(C.CO2.absorbed.per.g.dry.mass = C.CO2.absorbed.minus.blank.g/Dry.mass.g)

# Remove the rows for the Blanks (no longer needed)
df = df %>%
    filter(Incubation.trtmt != "Blank")

### Calculate cumulative C-CO2 respired over the course of the incubation
df <- df %>%
  group_by(core.id.incub) %>%
  arrange(Day, .by_group = TRUE) %>%
  mutate(cum.C.CO2.absorbed.g = cumsum(C.CO2.absorbed.minus.blank.g)) 
# cumsum - order of df is important!



# Normalize C-CO2 absorbed by day
df <- df %>%
  group_by(core.id.incub)%>%
  mutate(C.CO2.absorbed.per.g.dry.mass.per.day = C.CO2.absorbed.per.g.dry.mass/KOH.trap.time.period.exact) %>%
# Calcualte fractional C remaining
  mutate(mass.C.remaining.g = total.initial.dry.mass.C.g - cum.C.CO2.absorbed.g) %>%
  mutate(C.CO2.absorbed.per.g.C = C.CO2.absorbed.minus.blank.g/total.initial.dry.mass.C.g) %>%
  mutate(C.CO2.absorbed.per.g.C.per.day = C.CO2.absorbed.per.g.C/KOH.trap.time.period.exact)


df <- df %>%
  mutate(fractional.C.remaining = (mass.C.remaining.g)/total.initial.dry.mass.C.g) %>%
  mutate(frac = (total.initial.dry.mass.C.g - C.CO2.absorbed.per.g.C.per.day * total.initial.dry.mass.C.g * KOH.trap.time.period.exact)/total.initial.dry.mass.C.g)

# Quick check.  
ggplot(subset(df, site == 10), aes(x=Exact.day, y = mass.C.remaining.g, color = Burn.trtmt)) +
  geom_point() +
  facet_grid(~Incubation.trtmt) + 
  theme_bw() + 
  geom_hline(yintercept = 0)



df$incub.trtmt = df$Incubation.trtmt
df$burn.trtmt = df$Burn.trtmt
colnames(df)


summary(df$C.CO2.absorbed.per.g.dry.mass.per.day)


# Format for saving
X <- subset(df, select = c(core.id,
                           core.id.incub,
                           site,
                           incub.trtmt,
                           burn.trtmt,
                           Day,
                           mass.C.remaining.g,
                           C.CO2.absorbed.per.g.C,
                           fractional.C.remaining))


#write.csv(X, "fractional-c-remaining.csv")

head(df)

# ---------------
# PLOTS:


# Make Site numeric
colnames(df)
class(df$site)
df$site <- as.numeric(df$site)


# Make Site a factor instead of a continuous numeric value
class(df$site)
df$site = as.factor(df$site)
levels(df$site)

# Make Day a factor -> "DayFactor"
class(df$Day)
df$DayFactor = as.factor(df$Day)
levels(df$DayFactor)

# Display burn treatments in set order: Control, Wet, Dry
df$Burn.trtmt = ordered(df$Burn.trtmt,levels=c("control","wet burn","dry burn"))
df$burn.trtmt = ordered(df$burn.trtmt,levels=c("control","wet burn","dry burn"))

# Import site data
dfSite <- read.csv("raw-data/sampling-site-characteristics.csv")

names(dfSite)[1] <- 'core.id'

dfCOMP <- merge(df, dfSite, by = 'core.id')

# Making a vector with 3 colours
palette = c("black","orange","red3")

# Edit names of Leading Species
levels(dfCOMP$Leading.Species)
class(dfCOMP$Leading.Species)
Leading.Species.labs = c("P. glauca",
                         "P. mariana",
                         "P. banksiana",
                         "P. banksiana overstory",
                         "P. tremuloides")
names(Leading.Species.labs) <- c("Picea_glauca",
                                 "Picea_mariana",
                                 "Pinus_banksiana",
                                 "Pinus_banksiana overstory piceglauc understory",
                                 "Populus_tremula")

# Edit Site Labels
dfCOMP$Site.char <- as.factor(dfCOMP$Site)
levels(dfCOMP$Site.char)
class(dfCOMP$Site.char)


Site.labs = c("Site 01","Site 02","Site 03","Site 04","Site 05","Site 06",
               "Site 07","Site 08","Site 09","Site 10","Site 11","Site 12",
               "Site 13","Site 14","Site 15","Site 16","Site 17","Site 18",
               "Site 19")
names(Site.labs) <- c("1","2","3","4","5","6","7","8","9","10","11",
                      "12","13","14","15","16","17","18","19")


head(df_LI)
black_palette = c('black','black','black','black','black','black','black','black',
                  'black','black','black','black','black','black','black','black',
                  'black','black','black','black','black','black','black','black',
                  'black','black','black','black','black','black','black','black',
                  'black','black','black','black','black','black','black','black',
                  'black','black','black','black','black','black','black','black')



for (i in 1:nrow(dfCOMP)) {
  if (dfCOMP$Leading.Species[i] == 'Pinus_banksiana') {
    dfCOMP$veg[i] = 'Pinus banksiana'
  } else if (dfCOMP$Leading.Species[i] == 'Populus_tremula') {
    dfCOMP$veg[i] = 'P. tremuloides'
  } else if (dfCOMP$Leading.Species[i] == 'Picea_glauca') {
    dfCOMP$veg[i] = 'Picea spp.'
  }  else if (dfCOMP$Leading.Species[i] == 'Picea_mariana') {
    dfCOMP$veg[i] = 'Picea spp.'
  }  else if (dfCOMP$Leading.Species[i] == 'Pinus_banksiana overstory piceglauc understory') {
    dfCOMP$veg[i] = 'Pinus banksiana'
  }
}


dfCOMP$Leading.Species <-  factor(dfCOMP$Leading.Species,
                                    levels = c("Picea_glauca",
                                               "Picea_mariana",
                                               "Pinus_banksiana",
                                               "Pinus_banksiana overstory piceglauc understory",
                                               "Populus_tremula"),
                                    labels = c("italic(P.~glauca)", 
                                               "italic(P.~mariana)", 
                                               "italic(P.~banksiana)",
                                               "italic(P.~banksiana~overstory)",
                                               "italic(P.~tremuloides)"))

dfCOMP$site <- factor(dfCOMP$site,
                      levels = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19),
                        labels = c("Site~01","Site~02","Site~03","Site~04","Site~05","Site~06",
                                   "Site~07","Site~08","Site~09","Site~10","Site~11","Site~12",
                                   "Site~13","Site~14","Site~15","Site~16","Site~17","Site~18",
                                   "Site~19"))

# Create an SI dataframe
df_SI<- dfCOMP %>%
  subset(Incubation.trtmt =="SI")

df_LI <- dfCOMP %>%
  subset(Incubation.trtmt != "SI")


# Supplementary Figure S5
p1 = ggplot(subset(dfCOMP, incub.trtmt %in% c('LIwA') & site != 10), 
            aes(x = DayFactor, 
                #y = C.CO2.absorbed.per.g.C.per.day, 
                y = C.CO2.absorbed.per.g.C.per.day,
                fill = Burn.trtmt)) +
  #geom_point(size=2) +
  geom_boxplot(aes(x=DayFactor),alpha=0.6)+
  theme_bw() +
  scale_fill_manual(values = palette,
                    limits = c('control','wet burn','dry burn'),
                    labels = c('Control','Moist soil','Dry soil')) +
  #scale_shape_manual(values = c('circle','triangle'),
   #                  limits = c('LIn','SI'),
    #                 labels = c('Post-fire environment affinity','Fast growth'))+
  #facet_wrap(~veg+site, scale = "free",
   #          labeller = labeller(Leading.Species = Leading.Species.labs,site = Site.labs)) +
  labs(x = "Day",
       y = expression(Microbial~respiration~rate~(g~CO[2]-C~g~total~C^{"-1"}~day^{"-1"})),
       fill = "Burn treatment",
       shape='Incubation') +
  theme(strip.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        title = element_text(size = 12, margin = margin(b=20)),
        axis.title.y = element_text(size = 12, margin = margin(r=20, l=20)),
        axis.title.x = element_text(size = 12, margin = margin(t=20, b=20)))

p1


# Supplementary Figure S5
p1B = ggplot(subset(dfCOMP, incub.trtmt %in% c('LIwA')& site != 10), 
            aes(x = DayFactor, 
                #y = C.CO2.absorbed.per.g.C.per.day, 
                y = C.CO2.absorbed.per.g.dry.mass.per.day,
                fill = Burn.trtmt)) +
  #geom_point(size=2) +
  geom_boxplot(aes(x=DayFactor),alpha=0.6)+
  theme_bw() +
  scale_fill_manual(values = palette,
                    limits = c('control','wet burn','dry burn'),
                    labels = c('control','moist soil','dry soil')) +
  #scale_shape_manual(values = c('circle','triangle'),
  #                  limits = c('LIn','SI'),
  #                 labels = c('Post-fire environment affinity','Fast growth'))+
  #facet_wrap(~veg+site, scale = "free",
  #          labeller = labeller(Leading.Species = Leading.Species.labs,site = Site.labs)) +
  labs(x = "Day",
       y = expression(Microbial~respiration~rate~(g~CO[2]-C~g~dry~soil^{"-1"}~day^{"-1"})),
       fill = "Burn treatment",
       shape='Incubation') +
  theme(strip.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        title = element_text(size = 12, margin = margin(b=20)),
        axis.title.y = element_text(size = 12, margin = margin(r=20, l=20)),
        axis.title.x = element_text(size = 12, margin = margin(t=20, b=20)))

p1B

test = aov(C.CO2.absorbed.per.g.C.per.day~burn.trtmt, subset(dfCOMP, incub.trtmt %in% c('LIwA','LIn') & Day == 103))
summary(test)
TukeyHSD(test)




# Supplementary Figure S6
df.coefs <- read.csv('../data/2-pool-decay-model-output-LIwA-adjusted-for-total-C.csv')
df.coefs <- mutate(df.coefs, Day = t)


df_SI <- merge(subset(dfCOMP, Incubation.trtmt =="SI"), subset(df.coefs, select = c(core.id.incub, Mt.fit, Day)))
df_LIwA <- merge(subset(dfCOMP, incub.trtmt =="LIwA"), subset(df.coefs, select = c(core.id.incub, Mt.fit, Day)))


p2 = ggplot(subset(df_LIwA, site != 10), aes(x = Day, 
                      y = fractional.C.remaining, 
                      group = core.id.incub,
                      color = burn.trtmt,
                      shape = burn.trtmt)) +
  geom_point(size=3) +
  geom_line(aes(x=Day, y = Mt.fit))+
  theme_bw()+
  scale_color_manual(values = palette,
                     limits = c('control','wet burn','dry burn'),
                     labels = c('control','moist soil','dry soil'))+
  scale_shape_manual(values = c('circle','square','triangle'),
                     limits = c('control','wet burn','dry burn'),
                     labels = c('control','moist soil','dry soil'))+
  facet_wrap(~Leading.Species+site, scale = "fixed",ncol=4,
             labeller = label_parsed) +
  labs(x = "Day",
       y = "C remaining as a fraction of total C",
       #title = "Microbial respiration during long-term incubation",
       color = "Burn treatment",
       shape = 'Burn treatment') +
  theme(axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 12),
        strip.text.x = element_text(size = 10))
p2







# Supplementary Figure S7
df.coefs <- read.csv('../data/2-pool-decay-model-output-LIwA-adjusted-for-total-C.csv')
df.coefs <- mutate(df.coefs, Day = t)

df_LI <- merge(df_LI, subset(df.coefs, select = c(core.id.incub, Mt.fit, Day)))

p3 = ggplot(subset(df_LI, Incubation.trtmt == 'LIwA' & site != 10), aes(x = Day, 
                       y = fractional.C.remaining, 
                       group = core.id.incub,
                       color = Burn.trtmt,
                       shape = Burn.trtmt)) +
  geom_point(size=3, alpha=0.7) +
  geom_line(aes(x=Day, y = Mt.fit))+
  theme_bw()+
  scale_color_manual(values = palette) +
  #scale_shape_manual(values = c('triangle','circle'),
  #                  limits = c('SI', 'LIwA'),
  #                 labels = c('5 weeks','6 months')) +
  facet_wrap(~Leading.Species+site, scale = "fixed",
             labeller = labeller(Leading.Species = Leading.Species.labs,site = Site.labs)) +
  labs(x = "Day",
       y = "C remaining as a fraction of initial total C",
       #title = "Microbial respiration during long-term incubation",
       color = "Burn treatment",
       shape = 'Burn treatment') +
  theme(axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 12),
        strip.text.x = element_text(size = 10))

p3





# Supplementary Figure S4
p4 = ggplot(subset(dfCOMP, incub.trtmt %in% c('SI','LIwA')), aes(x = Day, 
                       y = C.CO2.absorbed.per.g.C.per.day, 
                       group = core.id.incub,
                       color = Burn.trtmt,
                       shape = incub.trtmt)) +
  geom_point(alpha=0.7,size=2) +
  #geom_line() +
  theme_bw() +
  scale_color_manual(values = palette) +
  #scale_color_manual(values = c('green4','black'),limits = c('SI','LIwA'),labels = c('5 weeks','6 months')) +
  scale_shape_manual(values = c('triangle','circle'),
                     limits = c('SI', 'LIwA'),
                     labels = c('5 weeks','6 months')) +
  #facet_grid(rows = vars(Leading.Species), 
   facet_wrap(~veg + site,          
              scales = "free",
              ncol=4,
              labeller = labeller(Leading.Species = Leading.Species.labs, 
                                 Site.char = Site.labs)) +
  labs(x = "Day",
       y = expression(Rate~of~CO[2]-C~loss~(g~C~respired~g^{"-1"}~initial~total~C~day^{"-1"})),
       color = "Incubation \ntreatment", 
       shape = 'Incubation \ntreatment') +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        title = element_text(size = 12, margin = margin(b=20)),
        axis.title.y = element_text(size = 14, margin = margin(r=20, l=20)),
        axis.title.x = element_text(size = 14, margin = margin(t=20, b=20)),
        legend.text = element_text(size=12))
p4





# Plot 5 - Exploring differences in respiration with autoclaving and inoculation:
p5 = ggplot(subset(dfCOMP, incub.trtmt %in% c('SI','LIwA', 'LIn') & 
                     Burn.trtmt == 'dry burn' & site != 10 &
                     Day < 40), aes(x = Day, 
                                    y = C.CO2.absorbed.per.g.C.per.day, 
                                    group = core.id.incub,
                                    color = incub.trtmt,
                                    shape = incub.trtmt)) +
  geom_point(alpha=0.8,size=3) +
  #geom_line() +
  theme_bw() +
  #scale_color_manual(values = palette) +
  scale_color_manual(values = c('blue3','grey60','black'),
                     limits = c('SI','LIwA', 'LIn'),
                     labels = c('5 weeks','6 months', '6 months: No autoclaving')) +
  scale_shape_manual(values = c('triangle','circle', 'square'),
                     limits = c('SI', 'LIwA', 'LIn'),
                     labels = c('5 weeks','6 months', '6 months: No autoclaving')) +
  #facet_grid(rows = vars(Leading.Species), 
  facet_wrap(~Leading.Species+site,          
             scales = "fixed",
             ncol=4,
             labeller = label_parsed) +
  labs(x = "Day",
       y = expression(Rate~of~CO[2]-C~loss~(g~C~respired~g^{"-1"}~total~C~day^{"-1"})),
       color = "Incubation treatment", 
       shape = 'Incubation treatment') +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        title = element_text(size = 12, margin = margin(b=20)),
        axis.title.y = element_text(size = 14, margin = margin(r=20, l=20)),
        axis.title.x = element_text(size = 14, margin = margin(t=20, b=20)),
        legend.text = element_text(size=12))
p5






p6 = ggplot(subset(dfCOMP, site==10 & incub.trtmt %in% c('LIwA', 'LIn')), aes(x = Day, 
                      y = mass.C.remaining.g, 
                      color = Day,
                      shape = burn.trtmt)) +
  geom_point(alpha=0.8,size=3) +
  geom_abline(intercept=0, slope=1)+
  #geom_line() +
  theme_bw() +
  #scale_color_manual(values = c('blue3','grey60','black'),limits = c('SI','LIwA', 'LIn'),labels = c('5 weeks','6 months', '6 months: No autoclaving')) +
  #scale_shape_manual(values = c('triangle','circle', 'square'),
   #                  limits = c('SI', 'LIwA', 'LIn'),
    #                 labels = c('5 weeks','6 months', '6 months: No autoclaving')) +
  facet_wrap(~Site)+
  #geom_line()+
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        title = element_text(size = 12, margin = margin(b=20)),
        axis.title.y = element_text(size = 14, margin = margin(r=20, l=20)),
        axis.title.x = element_text(size = 14, margin = margin(t=20, b=20)),
        legend.text = element_text(size=12)) 
p6




p6 = ggplot(y, aes(x = O.dry.mass.C.g, 
                                    y = A.dry.mass.C.g, 
                      shape = burn.trtmt,
                      color=burn.trtmt)) +
  geom_point(alpha=0.8,size=3) +
  #geom_line() +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1)+
  #scale_color_manual(values = c('blue3','grey60','black'),limits = c('SI','LIwA', 'LIn'),labels = c('5 weeks','6 months', '6 months: No autoclaving')) +
  #scale_shape_manual(values = c('triangle','circle', 'square'),
   #                  limits = c('SI', 'LIwA', 'LIn'),
    #                 labels = c('5 weeks','6 months', '6 months: No autoclaving')) +
  facet_wrap(~Site)+
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        title = element_text(size = 12, margin = margin(b=20)),
        axis.title.y = element_text(size = 14, margin = margin(r=20, l=20)),
        axis.title.x = element_text(size = 14, margin = margin(t=20, b=20)),
        legend.text = element_text(size=12))
p6





