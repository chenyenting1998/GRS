library(dplyr)
library(Polychrome)
library(GRSmacrofauna)

# density rare ----
rank_den <-
  size %>%
  group_by(Taxon) %>%
  summarize(Count = n()) %>%
  mutate(percent = round(Count/sum(Count) * 100, 2)) %>%
  arrange(desc(percent))

## Idensity dominant taxa-----
# label rare taxa as others
rank_den$Taxa <- if_else(rank_den$percent < 1, "Others", "Dominant")
# use the original label as donimants
rank_den[rank_den$Taxa == "Dominant",]$Taxa <- rank_den[rank_den$Taxa == "Dominant",]$Taxon

# factorize
rank_den_order <- rank_den$Taxa[length(unique(rank_den$Taxa)):1]
rank_den$Taxa <- factor(rank_den$Taxa, rank_den_order)

# biomass rare----
rank_bio <-
  size %>%
  group_by(Taxon) %>%
  summarize(Biomass = sum(Size)) %>%
  mutate(percent = round(Biomass/sum(Biomass) * 100, 2)) %>%
  arrange(desc(percent))

# label rare taxa as others
rank_bio$Taxa <- if_else(rank_bio$percent < 1, "Others", "Dominant")
# use the original label as donimants
rank_bio[rank_bio$Taxa == "Dominant",]$Taxa <- rank_bio[rank_bio$Taxa == "Dominant",]$Taxon

# factorize
rank_bio_order <- rank_bio$Taxa[length(unique(rank_bio$Taxa)):1]
rank_bio$Taxa <- factor(rank_bio$Taxa, rank_bio_order)

# all taxa list
all_taxa <- unique(c(rank_den_order, rank_bio_order))
# grab colors
taxa_color <- kelly.colors(length(all_taxa) + 1)[-1]
# assign colors to taxa
names(taxa_color) <- all_taxa

# den_taxa_color
taxa_den_color <- taxa_color[names(taxa_color) %in% rank_den$Taxa]

# bio_taxa_color
taxa_bio_color <- taxa_color[names(taxa_color) %in% rank_bio$Taxa]


save(rank_den, rank_bio,
     file = "data/taxa_rank.Rdata")
save(taxa_den_color, taxa_bio_color,
     file = "data/taxa_color.Rdata")
