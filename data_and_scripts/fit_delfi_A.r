## Pinocchio's Dolphin (Delfinus mendax) from Acme Archipelago

# Lethal sampling; adult age & juve age know exactly OK; "mammals"...
# ... ie knife-edge maturity; all adults equal WRTO reprod & survival; constant adult survival
# "Sampling season" is after "breeding season"
# Only female samples used; only female pop dyn considered; MOPs are actually M(other)D(aughter)Ps
# MOPs only (HSPs not used yet)

library( mvbutils)  # various low-level utilities, eg cq()
library( offarray)  # arrays/vectors don't have to start at 1 any more! Freedom!
library( debug)     # so we can watch the fun
library( atease)
source("ckmr_funs.r")

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
lglk_delfi_A <- generic_lglk_MOP_ideal_mammal # a copy that will be taught about delfi_A data
env <- boring_data_prep_delfi_A( samp_delfi_A, prev_env=environment( lglk_delfi_A))
environment( lglk_delfi_A) <- env # now, lglk_<blah> will always know eg what...
# ... catches, MOPs, comparisons etc were, without having to pass that stuff in every time it's called

starto <- c( log( 100000), 0.01) # logNad and RoI
lglk_delfi_A( starto) # just check it doesn't crash!

fitto <- nlminb( starto, NEG( lglk_delfi_A))

# Lots of stuff now available in env, which is _identical_ to environment( lglk_delfi_A)

###############
# Howzit lookin?
Amat <- env$Amat
Brange <- dimseq( env$N_Y) # what years get estimated?

tru_Nadf_Y <- sumover( samp_delfi_A@secret$N_sya[SLICE='F', Brange, Amat:40], 'A')

maxyplot <- max( c( tru_Nadf_Y, env$N_Y)) * 1.2 # bit of extra room for CI later
plot( Brange, tru_Nadf_Y, type='l', ylim=c( 0, maxyplot))
points( Brange, env$N_Y, col='green', pch=16)

#################
# Diagnostics (pretty feeble attempt)
# Expected MOPs (for diagnostic stuff)
e_MOP_BYA <- with( env,
  n_comp_MOP_BYA * Pr_MOP_BYA)
sum( e_MOP_BYA)
sum( env$n_MOP_BYA) # equal, phew

# Check whether first half and second half are OK
# What's the pattern in the first few years?
# Formal check would use chi-sq here... DoF=???
eee <- sumover( e_MOP_BYA[ 2011:2015,,], cq( Yad, Aad)) # expected
sum( eee)
ooo <- sumover( env$n_MOP_BYA[ 2011:2015,,], cq( Yad, Aad)) # observed
sum( ooo)
# ... very close

######################
# And you can do usual inference stuff. There's nothing CKMR-specific here
# Will need Hessian...
dlglk <- function( params) numderiv( lglk_delfi_A, params)
Hess <- numderiv( dlglk, fitto$par)
V_par <- solve( -Hess)

# Now use delta-method on stat-of-interest--- here "current absolute annual increase/decrease"
my_goal <- function( params) {
  lglk_delfi_A( params) # to set up all the Ns etc
return( env$N_Y[2020] * env$RoI)
}

scatn( "Point est of current annual delta-N: %5.0f animals", my_goal( fitto$par))
dmy_goal <- numderiv( my_goal, fitto$par)
V_my_goal <- dmy_goal %*% V_par %*% dmy_goal
scatn( "Standard error: %5.1f", sqrt( V_my_goal))

# So we could put some kinda CI on that plot, too:

my_goal_log_all_N <- function( params) {
  lglk_delfi_A( params) # to set up all the Ns etc
return( log( unclass( env$N_Y))) # whole lot; prolly don't need unclass()
}

lan <- log( env$N_Y)
d_mglan <- numderiv( my_goal_log_all_N, fitto$par)
V_mglan <- d_mglan %*% V_par %*% t( d_mglan)
SE_mglan <- sqrt( diag( V_mglan))

# CI bars: R base graphics is not very intuitive...
arrows( x0=Brange, y0=exp( lan - 1.645 * SE_mglan), y1=exp( lan + 1.645 * SE_mglan),
    code=3, angle=90, length=par( 'pin')[1] * 0.1 / length( Brange), col='grey')


#################
# Construct an "annual index" version, instead of fitting an explicit pop dyn model
# Only even-vaguely-defensible for "mammals" with MOP-only and no other params
# Need to count comparisons only where adult was alive & mature
# For some unknown reason, B-range is different in the CKMR data vs 'env' in general...
# ... also true for adult age... so just use actual dimensions of n_MOP & n_comp_MOP
extract.named( with( dimseq( env$n_comp_MOP_BYA),
  autoloop( B=Bju, Y=Yad, A=Aad, {
    is_plausible <- ( Y >= B) & (A-(Y-B) >= Amat)
    ncomp2 <- env$n_comp_MOP_BYA[ B, Y, A] * is_plausible
    npop2 <- env$n_MOP_BYA[ B, Y, A] * is_plausible # same as n_MOP unless I slipped up with >/>=...
  returnList( ncomp2, npop2)
  })))

n_comp_MOP_B <- sumover( ncomp2, cq( Y, A))
n_MOP_B <- sumover( npop2, cq( Y, A))
Nhat_B <- n_comp_MOP_B / n_MOP_B

# Point estimate for each year is then ncomp/nMOP--- setting up a lglk would give the same result...
# ... more slowly

maxyplot <- max( c( env$N_Y, tru_Nadf_Y, Nhat_B)) * 1.1
plot( Brange, tru_Nadf_Y, type='l', ylim=c( 0, maxyplot))
points( Brange, env$N_Y, col='green', pch=16)
points( dimseq( n_MOP_B), Nhat_B, col='orange', pch='x')
