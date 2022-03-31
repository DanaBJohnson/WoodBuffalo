library(dplyr)
library(corncob)
library(ggplot2)
getwd()
setwd('../Box/WhitmanLab/Projects/WoodBuffalo/FireSim2019/code')

df <- read.csv('../data/sequence-data/LibCombined/corncob-output/raw-df_L2FC.csv')
df.fast <- read.csv('../data/sequence-data/LibCombined/corncob-output/fast-Mid-and-low-thermo-L2FC.csv')


head(df,2)
head(df.fast,2)
unique(df$test)
unique(df.fast$test)

summary(df$cntl.pH)


df <- rbind(df, df.fast) %>%
  subset(!(test %in% c("dt.DvD.O.gt50","dt.DvD.A.gt50","dt.DvD.O.lt50","dt.DvD.A.lt50",
                       "dt.WvW.O.lt50" ,"dt.WvW.A.lt50","dt.CvC.O.lt50","dt.CvC.A.lt50")))

unique(df$test)

head(df)
dim(df)
unique(df$test)

df <- df %>% subset(test != 'test')
dim(df)

df$experiment = 'x'
df$burn.comparison = 'x'
df$horizon = 'x'
df$Temp.cutoff = 'x'
df$pH.cutoff = 'x'
df$controlling.for = 'x'

for (i in 1:nrow(df)) {
# Experiment 1, low thermo. 
  if (df$test[i] == 'dt.RNA.dryBurn.O.gt50') {
    df$test[i] = 'dt.low.RNA.dryBurn.O.gt50'
    df$experiment[i] = 'survival'
    df$burn.comparison[i] ='dry.v.control'
    df$horizon[i] = 'O'
    df$Temp.cutoff[i] = 'low thermo greater than 50 C'
    df$pH.cutoff[i] = 'null'
    df$controlling.for[i] = 'cntl.pH'
  }   else if (df$test[i] == 'dt.RNA.dryBurn.A.gt50') {
    df$test[i] = 'dt.low.RNA.dryBurn.A.gt50'
    df$experiment[i] = 'survival'
    df$burn.comparison[i] ='dry.v.control'
    df$horizon[i] = 'A'
    df$Temp.cutoff[i] = 'low thermo greater than 50 C'
    df$pH.cutoff[i] = 'null'
    df$controlling.for[i] = 'cntl.pH'
  }   else if (df$test[i] == 'dt.RNA.dryBurn.O.lt50') {
    df$test[i] = 'dt.low.RNA.dryBurn.O.lt50'
    df$experiment[i] = 'survival'
    df$burn.comparison[i] ='dry.v.control'
    df$horizon[i] = 'O'
    df$Temp.cutoff[i] = 'low thermo less than 50 C'
    df$pH.cutoff[i] = 'null'
    df$controlling.for[i] = 'cntl.pH'
  }   else if (df$test[i] == 'dt.RNA.wetBurn.O.lt50') {
    df$test[i] = 'dt.low.RNA.wetBurn.O.lt50'
    df$experiment[i] = 'survival'
    df$burn.comparison[i] ='wet.v.control'
    df$horizon[i] = 'O'
    df$Temp.cutoff[i] = 'low thermo less than 50 C'
    df$pH.cutoff[i] = 'null'
    df$controlling.for[i] = 'cntl.pH'
  }   else if (df$test[i] == 'dt.RNA.wetBurn.A.lt50') {
    df$test[i] = 'dt.low.RNA.wetBurn.A.lt50'
    df$experiment[i] = 'survival'
    df$burn.comparison[i] ='wet.v.control'
    df$horizon[i] = 'A'
    df$Temp.cutoff[i] = 'low thermo less than 50 C'
    df$pH.cutoff[i] = 'null'
    df$controlling.for[i] = 'cntl.pH'
    
# Experiment 1, mid thermo  
    } else if (df$test[i] == 'dt.mid.RNA.dryBurn.O.gt50') {
      df$test[i] = 'dt.mid.RNA.dryBurn.O.gt50'
      df$experiment[i] = 'survival'
      df$burn.comparison[i] ='dry.v.control'
      df$horizon[i] = 'O'
      df$Temp.cutoff[i] = 'mid thermo greater than 50 C'
      df$pH.cutoff[i] = 'null'
      df$controlling.for[i] = 'cntl.pH'
    } else if (df$test[i] == 'dt.mid.RNA.dryBurn.A.gt50') {
      df$test[i] = 'dt.mid.RNA.dryBurn.A.gt50'
      df$experiment[i] = 'survival'
      df$burn.comparison[i] ='dry.v.control'
      df$horizon[i] = 'A'
      df$Temp.cutoff[i] = 'mid thermo greater than 50 C'
      df$pH.cutoff[i] = 'null'
      df$controlling.for[i] = 'cntl.pH'
    }   else if (df$test[i] == 'dt.mid.RNA.dryBurn.O.lt50') {
      df$test[i] = 'dt.mid.RNA.dryBurn.O.lt50'
      df$experiment[i] = 'survival'
      df$burn.comparison[i] ='dry.v.control'
      df$horizon[i] = 'O'
      df$Temp.cutoff[i] = 'mid thermo less than 50 C'
      df$pH.cutoff[i] = 'null'
      df$controlling.for[i] = 'cntl.pH'
    }   else if (df$test[i] == 'dt.mid.RNA.dryBurn.A.lt50') {
      df$test[i] = 'dt.mid.RNA.dryBurn.A.lt50'
      df$experiment[i] = 'survival'
      df$burn.comparison[i] ='dry.v.control'
      df$horizon[i] = 'A'
      df$Temp.cutoff[i] = 'mid thermo less than 50 C'
      df$pH.cutoff[i] = 'null'
      df$controlling.for[i] = 'cntl.pH'
    }   else if (df$test[i] == 'dt.mid.RNA.wetBurn.O.lt50') {
      df$test[i] = 'dt.mid.RNA.wetBurn.O.lt50'
      df$experiment[i] = 'survival'
      df$burn.comparison[i] ='wet.v.control'
      df$horizon[i] = 'O'
      df$Temp.cutoff[i] = 'mid thermo less than 50 C'
      df$pH.cutoff[i] = 'null'
      df$controlling.for[i] = 'cntl.pH'
    }   else if (df$test[i] == 'dt.mid.RNA.wetBurn.A.lt50') {
      df$test[i] = 'dt.mid.RNA.wetBurn.A.lt50'
      df$experiment[i] = 'survival'
      df$burn.comparison[i] ='wet.v.control'
      df$horizon[i] = 'A'
      df$Temp.cutoff[i] = 'mid thermo less than 50 C'
      df$pH.cutoff[i] = 'null'
      df$controlling.for[i] = 'cntl.pH'
    }  else if (df$test[i] == 'dt.mid.RNA.wetBurn.O.gt50') {
      df$test[i] = 'dt.mid.RNA.wetBurn.O.gt50'
      df$experiment[i] = 'survival'
      df$burn.comparison[i] ='wet.v.control'
      df$horizon[i] = 'O'
      df$Temp.cutoff[i] = 'mid thermo greater than 50 C'
      df$pH.cutoff[i] = 'null'
      df$controlling.for[i] = 'cntl.pH'
    }   else if (df$test[i] == 'dt.mid.RNA.wetBurn.A.gt50') {
      df$test[i] = 'dt.mid.RNA.wetBurn.A.gt50'
      df$experiment[i] = 'survival'
      df$burn.comparison[i] ='wet.v.control'
      df$horizon[i] = 'A'
      df$Temp.cutoff[i] = 'mid thermo greater than 50 C'
      df$pH.cutoff[i] = 'null'
      df$controlling.for[i] = 'cntl.pH'
      
# Experiment 2
  }   else if (df$test[i] == 'dt.WvW.O.pH3to7') {
    df$experiment[i] = 'fast growth'
    df$burn.comparison[i] ='wet.pb.v.wet.SI'
    df$horizon[i] = 'O'
    df$Temp.cutoff[i] = 'null'
    df$pH.cutoff[i] = '3to7'
    df$controlling.for[i] = 'pH'
  }  else if (df$test[i] == 'dt.WvW.O.pH4to8') {
    df$experiment[i] = 'fast growth'
    df$burn.comparison[i] ='wet.pb.v.wet.SI'
    df$horizon[i] = 'O'
    df$Temp.cutoff[i] = 'null'
    df$pH.cutoff[i] = '4to8'
    df$controlling.for[i] = 'pH'
  }   else if (df$test[i] == 'dt.WvW.O.pH5to9') {
    df$experiment[i] = 'fast growth'
    df$burn.comparison[i] ='wet.pb.v.wet.SI'
    df$horizon[i] = 'O'
    df$Temp.cutoff[i] = 'null'
    df$pH.cutoff[i] = '5to9'
    df$controlling.for[i] = 'pH'
  }   else if (df$test[i] == 'dt.WvW.A.pH3to7') {
    df$experiment[i] = 'fast growth'
    df$burn.comparison[i] ='wet.pb.v.wet.SI'
    df$horizon[i] = 'A'
    df$Temp.cutoff[i] = 'null'
    df$pH.cutoff[i] = '3to7'
    df$controlling.for[i] = 'pH'
  }   else if (df$test[i] == 'dt.WvW.A.pH4to8') {
    df$experiment[i] = 'fast growth'
    df$burn.comparison[i] ='wet.pb.v.wet.SI'
    df$horizon[i] = 'A'
    df$Temp.cutoff[i] = 'null'
    df$pH.cutoff[i] = '4to8'
    df$controlling.for[i] = 'pH'
  }   else if (df$test[i] == 'dt.DvD.O.pH3to7') {
    df$experiment[i] = 'fast growth'
    df$burn.comparison[i] ='dry.pb.v.dry.SI'
    df$horizon[i] = 'O'
    df$Temp.cutoff[i] = 'null'
    df$pH.cutoff[i] = '3to7'
    df$controlling.for[i] = 'pH'
  }   else if (df$test[i] == 'dt.DvD.O.pH4to8') {
    df$experiment[i] = 'fast growth'
    df$burn.comparison[i] ='dry.pb.v.dry.SI'
    df$horizon[i] = 'O'
    df$Temp.cutoff[i] = 'null'
    df$pH.cutoff[i] = '4to8'
    df$controlling.for[i] = 'pH'
  }   else if (df$test[i] == 'dt.DvD.A.pH3to7') {
    df$experiment[i] = 'fast growth'
    df$burn.comparison[i] ='dry.pb.v.dry.SI'
    df$horizon[i] = 'A'
    df$Temp.cutoff[i] = 'null'
    df$pH.cutoff[i] = '3to7'
    df$controlling.for[i] = 'pH'
  }   else if (df$test[i] == 'dt.DvD.A.pH4to8') {
    df$experiment[i] = 'fast growth'
    df$burn.comparison[i] ='dry.pb.v.dry.SI'
    df$horizon[i] = 'A'
    df$Temp.cutoff[i] = 'null'
    df$pH.cutoff[i] = '4to8'
    df$controlling.for[i] = 'pH'
  }   else if (df$test[i] == 'dt.DvD.A.pH5to9') {
    df$experiment[i] = 'fast growth'
    df$burn.comparison[i] ='dry.pb.v.dry.SI'
    df$horizon[i] = 'A'
    df$Temp.cutoff[i] = 'null'
    df$pH.cutoff[i] = '5to9'
    df$controlling.for[i] = 'pH'
  }   else if (df$test[i] == 'dt.DvD.O.pH5to9') {
    df$experiment[i] = 'fast growth'
    df$burn.comparison[i] ='dry.pb.v.dry.SI'
    df$horizon[i] = 'O'
    df$Temp.cutoff[i] = 'null'
    df$pH.cutoff[i] = '5to9'
    df$controlling.for[i] = 'pH'
  } else if (df$test[i] == 'dt.CvC.O.pH3to7') {
    df$experiment[i] = 'fast growth'
    df$burn.comparison[i] ='control.pb.v.control.SI'
    df$horizon[i] = 'O'
    df$Temp.cutoff[i] = 'null'
    df$pH.cutoff[i] = '3to7'
    df$controlling.for[i] = 'pH'
  } else if (df$test[i] == 'dt.CvC.O.pH4to8') {
    df$experiment[i] = 'fast growth'
    df$burn.comparison[i] ='control.pb.v.control.SI'
    df$horizon[i] = 'O'
    df$Temp.cutoff[i] = 'null'
    df$pH.cutoff[i] = '4to8'
    df$controlling.for[i] = 'pH'
  } else if (df$test[i] == 'dt.CvC.A.pH3to7') {
    df$experiment[i] = 'fast growth'
    df$burn.comparison[i] ='control.pb.v.control.SI'
    df$horizon[i] = 'A'
    df$Temp.cutoff[i] = 'null'
    df$pH.cutoff[i] = '3to7'
    df$controlling.for[i] = 'pH'
  } else if (df$test[i] == 'dt.CvC.A.pH4to8') {
    df$experiment[i] = 'fast growth'
    df$burn.comparison[i] ='control.pb.v.control.SI'
    df$horizon[i] = 'A'
    df$Temp.cutoff[i] = 'null'
    df$pH.cutoff[i] = '4to8'
    df$controlling.for[i] = 'pH'
    
    # Experiment 2 switch to temp.cutoff - mid thermo
  } else if (df$test[i] == 'dt.mid.DvD.O.gt50') {
    df$experiment[i] = 'fast growth'
    df$burn.comparison[i] ='dry.pb.v.dry.SI'
    df$horizon[i] = 'O'
    df$Temp.cutoff[i] = 'mid thermo greater than 50 C'
    df$pH.cutoff[i] = 'null'
    df$controlling.for[i] = 'pH'
  } else if (df$test[i] == 'dt.mid.DvD.A.gt50') {
    df$experiment[i] = 'fast growth'
    df$burn.comparison[i] ='dry.pb.v.dry.SI'
    df$horizon[i] = 'A'
    df$Temp.cutoff[i] = 'mid thermo greater than 50 C'
    df$pH.cutoff[i] = 'null'
    df$controlling.for[i] = 'pH'
  }else if (df$test[i] == 'dt.mid.DvD.O.lt50') {
    df$experiment[i] = 'fast growth'
    df$burn.comparison[i] ='dry.pb.v.dry.SI'
    df$horizon[i] = 'O'
    df$Temp.cutoff[i] = 'mid thermo less than 50 C'
    df$pH.cutoff[i] = 'null'
    df$controlling.for[i] = 'pH'
  }else if (df$test[i] == 'dt.mid.DvD.A.lt50') {
    df$experiment[i] = 'fast growth'
    df$burn.comparison[i] ='dry.pb.v.dry.SI'
    df$horizon[i] = 'A'
    df$Temp.cutoff[i] = 'mid thermo less than 50 C'
    df$pH.cutoff[i] = 'null'
    df$controlling.for[i] = 'pH'
  }else if (df$test[i] == 'dt.mid.WvW.O.lt50') {
    df$experiment[i] = 'fast growth'
    df$burn.comparison[i] ='wet.pb.v.wet.SI'
    df$horizon[i] = 'O'
    df$Temp.cutoff[i] = 'mid thermo less than 50 C'
    df$pH.cutoff[i] = 'null'
    df$controlling.for[i] = 'pH'
  }else if (df$test[i] == 'dt.mid.WvW.A.lt50') {
    df$experiment[i] = 'fast growth'
    df$burn.comparison[i] ='wet.pb.v.wet.SI'
    df$horizon[i] = 'A'
    df$Temp.cutoff[i] = 'mid thermo less than 50 C'
    df$pH.cutoff[i] = 'null'
    df$controlling.for[i] = 'pH'
  }else if (df$test[i] == 'dt.mid.CvC.O.lt50') {
    df$experiment[i] = 'fast growth'
    df$burn.comparison[i] ='control.pb.v.control.SI'
    df$horizon[i] = 'O'
    df$Temp.cutoff[i] = 'mid thermo less than 50 C'
    df$pH.cutoff[i] = 'null'
    df$controlling.for[i] = 'pH'
  }else if (df$test[i] == 'dt.mid.CvC.A.lt50') {
    df$experiment[i] = 'fast growth'
    df$burn.comparison[i] ='control.pb.v.control.SI'
    df$horizon[i] = 'A'
    df$Temp.cutoff[i] = 'mid thermo less than 50 C'
    df$pH.cutoff[i] = 'null'
    df$controlling.for[i] = 'pH'
    
    
    # Experiment 2 switch to temp.cutoff - low thermo
  } else if (df$test[i] == 'dt.low.DvD.O.gt50') {
    df$experiment[i] = 'fast growth'
    df$burn.comparison[i] ='dry.pb.v.dry.SI'
    df$horizon[i] = 'O'
    df$Temp.cutoff[i] = 'low thermo greater than 50 C'
    df$pH.cutoff[i] = 'null'
    df$controlling.for[i] = 'pH'
  } else if (df$test[i] == 'dt.low.DvD.A.gt50') {
    df$experiment[i] = 'fast growth'
    df$burn.comparison[i] ='dry.pb.v.dry.SI'
    df$horizon[i] = 'A'
    df$Temp.cutoff[i] = 'low thermo greater than 50 C'
    df$pH.cutoff[i] = 'null'
    df$controlling.for[i] = 'pH'
  }else if (df$test[i] == 'dt.low.DvD.O.lt50') {
    df$experiment[i] = 'fast growth'
    df$burn.comparison[i] ='dry.pb.v.dry.SI'
    df$horizon[i] = 'O'
    df$Temp.cutoff[i] = 'low thermo less than 50 C'
    df$pH.cutoff[i] = 'null'
    df$controlling.for[i] = 'pH'
  }else if (df$test[i] == 'dt.low.DvD.A.lt50') {
    df$experiment[i] = 'fast growth'
    df$burn.comparison[i] ='dry.pb.v.dry.SI'
    df$horizon[i] = 'A'
    df$Temp.cutoff[i] = 'low thermo less than 50 C'
    df$pH.cutoff[i] = 'null'
    df$controlling.for[i] = 'pH'
  }else if (df$test[i] == 'dt.low.WvW.O.lt50') {
    df$experiment[i] = 'fast growth'
    df$burn.comparison[i] ='wet.pb.v.wet.SI'
    df$horizon[i] = 'O'
    df$Temp.cutoff[i] = 'low thermo less than 50 C'
    df$pH.cutoff[i] = 'null'
    df$controlling.for[i] = 'pH'
  }else if (df$test[i] == 'dt.low.WvW.A.lt50') {
    df$experiment[i] = 'fast growth'
    df$burn.comparison[i] ='wet.pb.v.wet.SI'
    df$horizon[i] = 'A'
    df$Temp.cutoff[i] = 'low thermo less than 50 C'
    df$pH.cutoff[i] = 'null'
    df$controlling.for[i] = 'pH'
  }else if (df$test[i] == 'dt.low.CvC.O.lt50') {
    df$experiment[i] = 'fast growth'
    df$burn.comparison[i] ='control.pb.v.control.SI'
    df$horizon[i] = 'O'
    df$Temp.cutoff[i] = 'low thermo less than 50 C'
    df$pH.cutoff[i] = 'null'
    df$controlling.for[i] = 'pH'
  }else if (df$test[i] == 'dt.low.CvC.A.lt50') {
    df$experiment[i] = 'fast growth'
    df$burn.comparison[i] ='control.pb.v.control.SI'
    df$horizon[i] = 'A'
    df$Temp.cutoff[i] = 'low thermo less than 50 C'
    df$pH.cutoff[i] = 'null'
    df$controlling.for[i] = 'pH'
    
    
    
  # Experiment 3
  } else if (df$test[i] == 'dT.mid.LIwA.dryBurn.O.gt50') {
    df$experiment[i] = 'affinity'
    df$burn.comparison[i] ='dry.v.control'
    df$horizon[i] = 'O'
    df$Temp.cutoff[i] = 'mid thermo greater than 50 C'
    df$pH.cutoff[i] = 'null'
    df$controlling.for[i] = 'cntl.pH'
  }else if (df$test[i] == 'dT.mid.LIwA.dryBurn.A.gt50') {
    df$experiment[i] = 'affinity'
    df$burn.comparison[i] ='dry.v.control'
    df$horizon[i] = 'A'
    df$Temp.cutoff[i] = 'mid thermo greater than 50 C'
    df$pH.cutoff[i] = 'null'
    df$controlling.for[i] = 'cntl.pH'
  }else if (df$test[i] == 'dT.mid.LIwA.wetBurn.A.lt50') {
    df$experiment[i] = 'affinity'
    df$burn.comparison[i] ='wet.v.control'
    df$horizon[i] = 'A'
    df$Temp.cutoff[i] = 'mid thermo less than 50 C'
    df$pH.cutoff[i] = 'null'
    df$controlling.for[i] = 'cntl.pH'
  } else if (df$test[i] == 'dT.mid.LIwA.wetBurn.O.lt50') {
    df$experiment[i] = 'affinity'
    df$burn.comparison[i] ='wet.v.control'
    df$horizon[i] = 'O'
    df$Temp.cutoff[i] = 'mid thermo less than 50 C'
    df$pH.cutoff[i] = 'null'
    df$controlling.for[i] = 'cntl.pH'
  }else if (df$test[i] == 'dT.low.LIwA.dryBurn.O.gt50') {
    df$experiment[i] = 'affinity'
    df$burn.comparison[i] ='dry.v.control'
    df$horizon[i] = 'O'
    df$Temp.cutoff[i] = 'low thermo greater than 50 C'
    df$pH.cutoff[i] = 'null'
    df$controlling.for[i] = 'cntl.pH'
  }else if (df$test[i] == 'dT.low.LIwA.dryBurn.O.lt50') {
    df$experiment[i] = 'affinity'
    df$burn.comparison[i] ='dry.v.control'
    df$horizon[i] = 'O'
    df$Temp.cutoff[i] = 'low thermo less than 50 C'
    df$pH.cutoff[i] = 'null'
    df$controlling.for[i] = 'cntl.pH'
  } else if (df$test[i] == 'dT.low.LIwA.dryBurn.A.gt50') {
  df$experiment[i] = 'affinity'
  df$burn.comparison[i] ='dry.v.control'
  df$horizon[i] = 'A'
  df$Temp.cutoff[i] = 'low thermo greater than 50 C'
  df$pH.cutoff[i] = 'null'
  df$controlling.for[i] = 'cntl.pH'
} else if (df$test[i] == 'dT.low.LIwA.wetBurn.O.lt50') {
  df$experiment[i] = 'affinity'
  df$burn.comparison[i] ='wet.v.control'
  df$horizon[i] = 'O'
  df$Temp.cutoff[i] = 'low thermo less than 50 C'
  df$pH.cutoff[i] = 'null'
  df$controlling.for[i] = 'cntl.pH'
} else if (df$test[i] == 'dT.low.LIwA.wetBurn.A.lt50') {
  df$experiment[i] = 'affinity'
  df$burn.comparison[i] ='wet.v.control'
  df$horizon[i] = 'A'
  df$Temp.cutoff[i] = 'low thermo less than 50 C'
  df$pH.cutoff[i] = 'null'
  df$controlling.for[i] = 'cntl.pH'
 }
}

unique(df$test)

df.cntlpH <- subset(df, experiment %in% c('survival', 'affinity'))
df.fast <- subset(df, experiment == 'fast growth')

df.fast <- df.fast %>%
  mutate(mean.sample.pH = cntl.pH,
         cntl.pH = NA) 
df.cntlpH <- df.cntlpH %>%
  mutate(mean.sample.pH = NA)


df <- rbind(df.fast, df.cntlpH)



for (i in 1:nrow(df)) {
  if (df$experiment[i] == 'survival') {
    df$RA.burned[i] = invlogit(df$intercept[i] + df$burn[i] + df$cntl.coef[i] * df$cntl.pH[i])
    df$RA.base[i] = invlogit(df$intercept[i] + df$cntl.coef[i]*df$cntl.pH[i])
    df$fold.change[i] = df$RA.burned[i]/df$RA.base[i]
    df$L2FC[i] = log(df$fold.change[i], base=2)
  } else if (df$experiment[i] == 'fast growth') {
    df$RA.burned[i] = invlogit(df$intercept[i] + df$burn[i] + df$cntl.coef[i] * df$mean.sample.pH[i])
    df$RA.base[i] = invlogit(df$intercept[i] + df$cntl.coef[i]*df$mean.sample.pH[i])
    df$fold.change[i] = df$RA.burned[i]/df$RA.base[i]
    df$L2FC[i] = log(df$fold.change[i], base=2)
  } else if (df$experiment[i] == 'affinity') {
    df$RA.burned[i] = invlogit(df$intercept[i] + df$burn[i] + df$cntl.coef[i] * df$cntl.pH[i])
    df$RA.base[i] = invlogit(df$intercept[i] + df$cntl.coef[i]*df$cntl.pH[i])
    df$fold.change[i] = df$RA.burned[i]/df$RA.base[i]
    df$L2FC[i] = log(df$fold.change[i], base=2)
  }
}

ggplot(data = df, aes(x=burn, y = L2FC)) + 
  geom_point(size=2) + 
  facet_wrap(~experiment)


head(df)
colnames(df)[2] = 'mu'
colnames(df)[5] = 'OTU'
colnames(df)[8] = 'Experiment'
#write.csv(df, '../data/sequence-data/LibCombined/corncob-output/Full_L2FC.csv', row.names = FALSE)

pHmean=mean(sample_data(ps.corncob.dry.O.gt50)$cntl.pH)

RA.burned.31=invlogit(m.31.coef[1]+m.31.coef[2]+m.31.coef[3]*pHmean)
RA.base.31=invlogit(m.31.coef[1]+m.31.coef[3]*pHmean)
FC.31 = RA.burned.31/RA.base.31
L2FC.31 = log(FC.31,base=2)





