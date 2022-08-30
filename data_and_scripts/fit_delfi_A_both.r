## Pinocchio's Dolphin (Delfinus mendax) from Acme Archipelago

# Lethal sampling; adult age & juve age know exactly OK; "mammals"...
# ... ie knife-edge maturity; all adults equal WRTO reprod & survival; constant adult survival
# "Sampling season" is after "breeding season"
# Only female samples used; only female pop dyn considered; POPs are actually Parent(other)O(ffspring)Ps
# POPs only (HSPs not used yet)

library( mvbutils)  # various low-level utilities, eg cq()
library( offarray)  # arrays/vectors don't have to start at 1 any more! Freedom!
library( debug)     # so we can watch the fun
library( atease)
source("ckmr_funs.r")
source("fit_delfi_A_both_funs.r")

### ACTION STARTS HERE

## Load and set up the data
od <- options( digits=7)
on.exit( options( od))
print( load( 'samp_delfi_A.rda')) # object samp_notog1
print( names( samp_delfi_A))
print( head( samp_delfi_A$Samps))
scatn( 'Delfi_A #POPs = %i', nrow( samp_delfi_A$POPs))

## samp_delfi_A <- "READ THE DAMN SCRIPT AND EDIT APPROPRIATELY; DON'T JUST RUN IT"

# Summarize sampling & kin & other data for delfi_A, and teach the lglk function about it
lglk_delfi_A_both <- generic_lglk_POP_ideal_mammal_both # a copy that will be taught about delfi_A data
envb <- boring_data_prep_delfi_A_both( samp_delfi_A, prev_env=environment( lglk_delfi_A_both))
environment( lglk_delfi_A_both) <- envb # now, lglk_<blah> will always know eg what...
# ... catches, POPs, comparisons etc were, without having to pass that stuff in every time it's called

starto <- c( log( 100000), 0.01, logit(0.45)) # logNad and RoI
lglk_delfi_A_both( starto) # just check it doesn't crash!

fitto <- nlminb( starto, NEG( lglk_delfi_A_both))

# Lots of stuff now available in env, which is _identical_ to environment( lglk_delfi_A)

###############
# Howzit lookin?
Amat <- envb$Amat

tru_Nadf <- sumover( samp_delfi_A@secret$N_sya[, envb$years, Amat:40], c('A', 'SEX'))
pred_Nadf <- sumover(envb$N_YS, 'S')

yl <- max(c(tru_Nadf, pred_Nadf))
plot(envb$years, tru_Nadf, type="l", ylim=c(0, yl))
points(envb$years, pred_Nadf, pch=19, col="dodgerblue")
