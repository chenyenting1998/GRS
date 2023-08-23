# Variance partitioning ------------
# match env_spatial dataframe wtih count_wide
env_spatial_expand <- 
  left_join(count_wide, env) %>% 
  select(all_of(env_variables_spatial))
# RDA
count_rda_spatial <- rda(count_chord ~ scale(env_spatial_expand))

# r.square
RsquareAdj(count_rda_spatial)
# backward selection
count_rda_spatial_back <- ordistep(count_rda_spatial, method = "backward")
# anova.cca
set.seed(10)
count_rda_spatial_axis <- anova(count_rda_spatial, by = "axis", permutations = 9999)
count_rda_spatial_margin <- anova(count_rda_spatial, by = "margin", permutations = 9999)

# plot_RDA
count_varpart <- 
  varpart(count_chord,
          ~ scale(env_spatial_expand),
          ~ scale(env_selected_expand))

showvarparts(2, bg = c("hotpink","skyblue"))
plot(count_varpart, bg = c("hotpink","skyblue"))

