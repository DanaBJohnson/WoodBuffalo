# --------------------
# 
# Title: creating-burn-temp-profiles
# Author: Dana Johnson
# Date: 2020-Feb-3
#
# Purpose: This is code to plot the temperature profiles during simulated burns.
# Can I add this to Git?
# --------------------
# Import packages
library(tidyverse)
library(stringr)
library(ggplot2)
# --------------------

# Read temperature data
df <-read.csv('../Temp-all-burns-clean.csv')


# Quick look at maximum temperatures across vegetation type and burn severity
X <- df %>%
  subset(Leading_species == 'Picea_mariana') %>%
  subset(Thermo_position == 'Mid') %>%
  subset(Burn_severity == 'High_Sev') %>%
  group_by(Site) %>%
  mutate(max.Temp = max(Temp))
range(X$max.Temp)


# Fix column name
names(df)[1] <- 'core.id'

# Remove core that fell on ground
df <- subset(df, core.id != '19UW-WB-11-02')

# Change time to numeric
class(df$time_s)
df$time_s = as.numeric(df$time_s)


# Integrate area under the curve for "cumulative temperature" (terminology?)
df.full <- read.csv('../Temp-all-burns.csv')

# Calculate cumulative temperature above a baseline of 21 degrees C and convert 
#   to "degrees * hours"
# 
df.full <- df.full %>%
  subset(Thermo_mid >0) %>%
  subset(Thermo_low > 0) %>%
  group_by(Site_ID) %>%
  mutate(cumul.T.s.mid = sum(Thermo_mid-18)) %>%
  mutate(cumul.T.s.low = sum(Thermo_low-18)) %>%
  mutate(cumul.T.hr.mid = cumul.T.s.mid /3600) %>%
  mutate(cumul.T.hr.low = cumul.T.s.low/3600) %>%
# Subset this dataframe to include just the site ID and the cumulative temp data.
  subset(time_s == 501) %>%
  subset(select = c(Site_ID,
                    cumul.T.hr.mid,
                    cumul.T.hr.low))
# Save as .csv file
# write.csv(df.area, '../Cumulative-temp-above-18-degrees.csv')



# Create names vector 
names.vec <- seq(1,19, by=1)

# Create output dataframe
output.df <- data.frame("core.id" = NULL, 
                            "Site" = NULL, 
                            "Core" = NULL,
                            "time_s" = NULL,
                            "Thermo_position" = NULL,
                            "Temp" = NULL,
                            "Burn_severity" = NULL,
                            "Leading_species" = NULL,
                            "ten_s" = NULL,
                            "time_minutes" = NULL,
                            "time_hr" = NULL)

# Subset the data by taking every tenth time point
for (i in seq_along(names.vec)) {
  df.x <- df %>%
    subset(Site == names.vec[[i]]) %>%
    subset(time_s < 21600) %>%
    mutate(ten_s = time_s/100) %>%
    filter(ten_s == as.integer(ten_s)) %>%
    mutate(time_minutes = time_s/60) %>%
    mutate(time_hr = time_minutes/60) 
  
  # Filter out negative temperature values as a quality control
  df.x <- df.x %>%
    subset(Temp >= 0)
  
  output.df <- rbind(output.df, df.x)
}

head(output.df)

# Import site property data and merge with temp data.
dfSite <- read.csv("../../raw-data/sampling-site-characteristics.csv")
colnames(dfSite)
names(dfSite)[1] <- 'core.id'

df <- merge(dfSite, output.df)

df$Site <- as.numeric(df$Site)
levels(df$Site)

# Import wet mass loss data and merge with temp data
df_wet <-read.csv('../../wet.mass.loss.g.during.burn.csv')
head(df_wet)

# Remove core that fell on ground
df_wet <- subset(df_wet, core.id != '19UW-WB-11-02')


df_wet <- df_wet[,2:3]

df <- merge(df, df_wet)
head(df)



df <- df %>%
  separate(Core.trtmt, c('burn.trtmt','burn','duplicate'), sep = ' ', remove = TRUE) %>%
  subset(select = c(core.id, Site,Core,burn.trtmt,Leading.Species,Temp, Thermo_position,
                    time_s, Burn_severity, time_minutes,time_hr,wet.mass.loss.g,
                    Bulk.mass.g,Bulk.height.g,Bulk.diameter.cm,Bulk.density.dry.g.cm.3,
                    Bulk.density.wet.g.cm.3,Mass.wet.soil.g,Mass.dry.soil.g,Percent.moisture))

head(df)
#write.csv(df, '../Temp-all-burns-clean.csv')





### --------------------
# PLOTS


# Create custom color palettes
palette2 = c("red3", "orange")
palette5 = c('green3','blue1','gold','purple','blue4')


# Create labels for Leading species and Sites
levels(df$Leading_species)
df$Leading_species <- as.factor(df$Leading_species)
Leading.Species.labs = c("Picea glauca",
                         "Picea mariana",
                         "Pinus banksiana",
                         "Pinus banksiana overstory",
                         "Populus tremuloides")
names(Leading.Species.labs) <- c("Picea_glauca",
                                 "Picea_mariana",
                                 "Pinus_banksiana",
                                 "Pinus_banksiana overstory piceglauc understory",
                                 "Populus_tremula")

df$Site <- as.factor(df$Site)
levels(df$Site)

Site.labs <- c("Site 01","Site 02","Site 03","Site 04","Site 05","Site 06",
               "Site 07","Site 08","Site 09","Site 10","Site 11","Site 12",
               "Site 13","Site 14","Site 15","Site 16","Site 17","Site 18",
               "Site 19")
names(Site.labs) <- c("1","2","3","4","5","6","7","8","9","10","11",
                      "12","13","14","15","16","17","18","19")

df$Burn_severity <- as.factor(df$Burn_severity)
levels(df$Burn_severity)
Burn.severity.labs = c('Wet burn','Dry burn')
names(Burn.severity.labs) <- c("Low_Sev",'High_Sev')


### PLOT: Cumulative mins above 100C by site 
# Calculate cumulative time above 100 degrees C
df_100 <- subset(df, Temp >= 100)

# Calculate time above 100 degrees C
df_100 <- df_100 %>%
  group_by(core.id) %>%
  mutate(range.time.min = max(time_minutes) - min(time_minutes)) %>%
  mutate(min.time = min(time_minutes))




# Quick look at time over 100 C
X <- df_100 %>%
  subset(Leading.Species == 'Populus_tremula')
range(X$range.time.min)




p1 = ggplot(df_100, aes(x = Site, 
                          y = range.time.min, 
                          group=core.id, 
                          color = Leading.Species)) +
  geom_point(size=4) +
  labs(x = "Site",
       y = expression(paste("Cumulative minutes above 100 ( ",degree ~ C, ")")),
       shape = "Position of Thermocouple",
       color = "Burn severity treatment") + 
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  scale_x_discrete(limits = c('1','2','3','4','5','6','7','8','9','10',
  '11','12','13','14','15','16','17','18','19')) +
  scale_color_manual(values = palette5)
p1
#ggsave("../../../figures/time-over-100-C.png", plot = p1, width = 12, height = 6)



# PLOT cumulative time over 100 degrees vs wet mass loss (g)
p2 <- ggplot(df_100, aes(x = wet.mass.loss.g,
                    y = range.time.min/60)) + 
  geom_point(size=3) + 
  geom_text(aes(label=Site),
            hjust= 1,
            vjust=1.7) +
  facet_grid(~Leading.Species,
             labeller = labeller(Leading.Species = Leading.Species.labs)) + 
  labs(x = 'Wet mass loss (g)',
       y = expression(paste('Total hours above 100 ', degree, 'C'))) + 
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(16),
        strip.text.y = element_text(16)) + 
  ggtitle("XXX") +
  theme(plot.title = element_text(color = "black", hjust = 0.5, size = 18)) + 
  geom_hline(yintercept=0,color ='black')
p2

#ggsave('../../../figures/wet-mass-loss-vs-time-above-100.png',plot=p, height = 6, width = 10)

colnames(df)



# PLOT: High and Low severity temp profiles
df_High <- subset(df, Burn_severity == 'High_Sev')
df_Low <- subset(df, Burn_severity == 'Low_Sev')


p3 <- ggplot (df_High, aes (x = time_hr,
                       y = Temp,
                       color = Thermo_position)) + 
  geom_point() +
  labs(title = 'Dry burn, all sites',
       x = "Time from start of burn (hours)",
       y = expression(paste("Temperature of core (",degree~C,")"))) + 
  facet_wrap(~Leading.Species + Site,
             labeller = labeller(Leading.Species = Leading.Species.labs,
                                 Site = NA)) +
  theme(strip.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_line(color = 'grey'),
        panel.border = element_rect(color = 'black',fill = NA)) + 
  theme(axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.text.x = element_text(size = 14)) +
  scale_color_manual(name = 'Thermocouple position',
                       breaks = c('Mid','Low'),
                       labels = c('O/Mnrl interface', 'Core base'),
                       values = palette2) +
  ylim(18,600) 
#  xlim(0,5)
p3
#ggsave("../../../figures/dry-burn-temp-profiles.png", plot = p3, width=14, height=10)




# PLOT: High and Low severity profiles side-by-side
df_Site <- subset(df, Site == '12')
p4 <- ggplot (df_Site, aes (x = time_hr,
                               y = Temp,
                               color = Thermo_position)) + 
  geom_point() +
  labs(title = 'Site 12, Picea mariana',
       x = "Time from start of burn (hours)",
       y = expression(paste("Temperature of core (",degree~C,")"))) + 
  facet_wrap(~Burn_severity,
             labeller = labeller(Burn_severity = Burn.severity.labs)) +
  theme(axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.text.x = element_text(size = 14)) +
  scale_color_manual(name = 'Thermocouple position',
                     breaks = c('Mid','Low'),
                     labels = c('O/Mnrl interface', 'Core base'),
                     values = palette2) +
  geom_hline(yintercept=0,color ='black') +
  ylim(18,600) 
#  xlim(0,5)
p4

#ggsave("../../../figures/site-12-temp-profiles.png", plot = p4, width=6, height=4)



# PLOT - temp profiles by vegetation type
df_Picea <- subset()
p5 <- ggplot (df_High, aes (x = time_hr,
                            y = Temp,
                            color = Thermo_position)) + 
  geom_point() +
  labs(x = "Time from start of burn (hours)",
       y = "Temperature of core (degrees Celsius)") + 
  facet_wrap(~Leading.Species + Site,
             labeller = labeller(Leading.Species = Leading.Species.labs,
                                 Site = element_blank())) +
  theme(strip.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_line(color = 'grey'),
        panel.border = element_rect(color = 'black',fill = NA)) + 
  theme(axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.text.x = element_text(size = 14)) +
  scale_color_manual(name = 'Thermocouple position',
                     breaks = c('Mid','Low'),
                     labels = c('O/Mnrl interface', 'Core base'),
                     values = palette2) +
  ylim(18,600) 
#  xlim(0,5)
p5