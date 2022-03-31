ps.cDNA.wet <- prune_samples(sample_data(ps.cDNA.raw)$burn.trtmt == 'wet', ps.cDNA.raw)

ps.cDNA.wet <- prune_taxa(taxa_sums(ps.cDNA.wet)>0, ps.cDNA.wet)

df.cDNA.wet <- psmelt(ps.cDNA.wet)


df.control.counts <- df.cDNA.control %>%
  subset(Abundance > 0) %>%
  mutate(Frequence = 1) %>%
  group_by(site) %>%
  summarize(OTUnum = sum(Frequence))

df.dry.site.counts <- df.cDNA.dry %>%
  subset(Abundance > 0) %>%
  subset(select = c(site, OTU))


df.wet.site.counts <- df.cDNA.wet %>%
  subset(Abundance > 0) %>%
  subset(select = c(site, OTU))


df.counts <- df.cDNA.control %>%
  subset(Abundance > 0) %>%
  merge(df.dry.site.counts) %>%
  mutate(Frequency = 1) %>%
  group_by(site) %>%
  summarize(OTUnum = sum(Frequency))



ps.gDNA.dry <- prune_samples(sample_data(ps.gDNA.raw)$burn.trtmt == 'dry', ps.gDNA.raw)

ps.gDNA.dry <- prune_taxa(taxa_sums(ps.gDNA.dry)>0, ps.gDNA.dry)

df.gDNA.dry <- psmelt(ps.gDNA.dry)

ps.gDNA.control <- prune_samples(sample_data(ps.gDNA.raw)$burn.trtmt == 'control', ps.gDNA.raw)

ps.gDNA.control <- prune_taxa(taxa_sums(ps.gDNA.control)>0, ps.gDNA.control)

df.gDNA.control <- psmelt(ps.gDNA.control)

ps.gDNA.wet <- prune_samples(sample_data(ps.gDNA.raw)$burn.trtmt == 'wet', ps.gDNA.raw)

ps.gDNA.wet <- prune_taxa(taxa_sums(ps.gDNA.wet)>0, ps.gDNA.wet)

df.gDNA.wet <- psmelt(ps.gDNA.wet)


df.control.counts <- df.gDNA.control %>%
  subset(Abundance > 0) %>%
  mutate(Frequence = 1) %>%
  group_by(site) %>%
  summarize(OTUnum = sum(Frequence))

df.dry.site.counts <- df.gDNA.dry %>%
  subset(Abundance > 0) %>%
  subset(select = c(site, OTU))


df.wet.site.counts <- df.gDNA.wet %>%
  subset(Abundance > 0) %>%
  subset(select = c(site, OTU))


df.counts <- df.gDNA.control %>%
  subset(Abundance > 0) %>%
  merge(df.wet.site.counts) %>%
  mutate(Frequency = 1) %>%
  group_by(site) %>%
  summarize(OTUnum = sum(Frequency))




ps.SI.dry <- prune_samples(sample_data(ps.raw.SI)$burn.trtmt == 'dry', ps.raw.SI)

ps.SI.dry <- prune_taxa(taxa_sums(ps.SI.dry)>0, ps.SI.dry)

df.SI.dry <- psmelt(ps.SI.dry)

ps.SI.control <- prune_samples(sample_data(ps.raw.SI)$burn.trtmt == 'control', ps.raw.SI)

ps.SI.control <- prune_taxa(taxa_sums(ps.SI.control)>0, ps.SI.control)

df.SI.control <- psmelt(ps.SI.control)

ps.SI.wet <- prune_samples(sample_data(ps.raw.SI)$burn.trtmt == 'wet', ps.raw.SI)

ps.SI.wet <- prune_taxa(taxa_sums(ps.SI.wet)>0, ps.SI.wet)

df.SI.wet <- psmelt(ps.SI.wet)


df.control.counts <- df.SI.control %>%
  subset(Abundance > 0) %>%
  mutate(Frequence = 1) %>%
  group_by(site) %>%
  summarize(OTUnum = sum(Frequence))

df.dry.site.counts <- df.SI.dry %>%
  subset(Abundance > 0) %>%
  subset(select = c(site, OTU))


df.wet.site.counts <- df.SI.wet %>%
  subset(Abundance > 0) %>%
  subset(select = c(site, OTU))


df.counts <- df.SI.control %>%
  subset(Abundance > 0) %>%
  merge(df.wet.site.counts) %>%
  mutate(Frequency = 1) %>%
  group_by(site) %>%
  summarize(OTUnum = sum(Frequency))



ps.raw.LIwA <- prune_samples(sample_data(ps.raw.LI)$incub.trtmt == 'LIwA', ps.raw.LI)

ps.LIwA.dry <- prune_samples(sample_data(ps.raw.LIwA)$burn.trtmt == 'dry', ps.raw.LIwA)

ps.LIwA.dry <- prune_taxa(taxa_sums(ps.LIwA.dry)>0, ps.LIwA.dry)

df.LIwA.dry <- psmelt(ps.LIwA.dry)

ps.LIwA.control <- prune_samples(sample_data(ps.raw.LIwA)$burn.trtmt == 'control', ps.raw.LIwA)

ps.LIwA.control <- prune_taxa(taxa_sums(ps.LIwA.control)>0, ps.LIwA.control)

df.LIwA.control <- psmelt(ps.LIwA.control)

ps.LIwA.wet <- prune_samples(sample_data(ps.raw.LIwA)$burn.trtmt == 'wet', ps.raw.LIwA)

ps.LIwA.wet <- prune_taxa(taxa_sums(ps.LIwA.wet)>0, ps.LIwA.wet)

df.LIwA.wet <- psmelt(ps.LIwA.wet)


df.control.counts <- df.LIwA.control %>%
  subset(Abundance > 0) %>%
  mutate(Frequence = 1) %>%
  group_by(site) %>%
  summarize(OTUnum = sum(Frequence))

df.dry.LIwAte.counts <- df.LIwA.dry %>%
  subset(Abundance > 0) %>%
  subset(select = c(site, OTU))


df.wet.site.counts <- df.LIwA.wet %>%
  subset(Abundance > 0) %>%
  subset(select = c(site, OTU))


df.counts <- df.LIwA.control %>%
  subset(Abundance > 0) %>%
  merge(df.dry.site.counts) %>%
  mutate(Frequency = 1) %>%
  group_by(site) %>%
  summarize(OTUnum = sum(Frequency))




ps.raw.LIn <- prune_samples(sample_data(ps.raw.LI)$incub.trtmt == 'LIn', ps.raw.LI)

ps.LIn.dry <- prune_samples(sample_data(ps.raw.LIn)$burn.trtmt == 'dry', ps.raw.LIn)

ps.LIn.dry <- prune_taxa(taxa_sums(ps.LIn.dry)>0, ps.LIn.dry)

df.LIn.dry <- psmelt(ps.LIn.dry)

ps.LIn.control <- prune_samples(sample_data(ps.raw.LIn)$burn.trtmt == 'control', ps.raw.LIn)

ps.LIn.control <- prune_taxa(taxa_sums(ps.LIn.control)>0, ps.LIn.control)

df.LIn.control <- psmelt(ps.LIn.control)

ps.LIn.wet <- prune_samples(sample_data(ps.raw.LIn)$burn.trtmt == 'wet', ps.raw.LIn)

ps.LIn.wet <- prune_taxa(taxa_sums(ps.LIn.wet)>0, ps.LIn.wet)

df.LIn.wet <- psmelt(ps.LIn.wet)


df.control.counts <- df.LIn.control %>%
  subset(Abundance > 0) %>%
  mutate(Frequence = 1) %>%
  group_by(site) %>%
  summarize(OTUnum = sum(Frequence))

df.dry.LInte.counts <- df.LIn.dry %>%
  subset(Abundance > 0) %>%
  subset(select = c(site, OTU))


df.wet.site.counts <- df.LIn.wet %>%
  subset(Abundance > 0) %>%
  subset(select = c(site, OTU))


df.counts <- df.LIn.control %>%
  subset(Abundance > 0) %>%
  merge(df.dry.site.counts) %>%
  mutate(Frequency = 1) %>%
  group_by(site) %>%
  summarize(OTUnum = sum(Frequency))
