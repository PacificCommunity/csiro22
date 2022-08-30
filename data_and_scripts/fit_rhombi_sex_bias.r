## Cobalt Squarehead (Rhombichthys cyanorosea) with sex-biased sampling

# Use original dataset, but subsample it differentially (fewer juves, too)
# Took me *ages* to not have bugs in this...

samp_rhombi_sex_bias <- subset_samples( samp_rhombi_D, 
    runif( 4000) < ifelse (is.na( Sex),
      0.4, ifelse( Sex=='M', 
      0.1,
      0.5)),
    guess=F)
with( samp_rhombi_sex_bias$Samps, table( Sex, useNA='ifany'))

# We do not need to code a lglk in this case, because the estimates are obvious; but 

env2 <- boring_data_prep_rhombi_D( samp_rhombi_sex_bias, prev_env=.GlobalEnv)
env2$n_comp_POP_S
env2$n_POP_s


# The truth...
samp_rhombi_sex_bias@secret$Nad_S # trooo

# Sex-specific estimates
Nobvious2 <- with( env2, n_comp_POP_S / n_POP_S)
sum( Nobvious2) # both sexes

# How would aggregate version go?
Nsexlumpo2 <- with( env2, 2 * sum( n_comp_POP_S) / sum( n_POP_S))
Nsexlumpo2 # Oooerr missus...

