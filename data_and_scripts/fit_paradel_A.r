# Lethal sampling; HSP s
# ... ie knife-edge maturity; all adults equal WRTO reprod & survival; constant adult survival
# "Sampling season" is after "breeding season"
# Only female samples used; only female pop dyn considered; MHSPs only (i hope)


### ACTION STARTS HERE
library( mvbutils)  # various low-level utilities, eg cq()
library( offarray)  # arrays/vectors don't have to start at 1 any more! Freedom!
library( atease)    # x@a <=> attr( x, "a")
library( debug)     # so we can watch the fun

## Load and set up the data
od <- options( digits=7)
on.exit( options( od))

print( load( 'samp_paradel_A.rda'))
print( names( samp_paradel_A))
print( head( samp_paradel_A$Samps))
print( ls( samp_paradel_A@public))

# Data summaries. Notice anything?
with( samp_paradel_A, with( Samps, table( 
    Sex[ POPs[,1]], 
    Y[ POPs[,1]] - (Y[POPs[,2]]-A[POPs[,2]]
  ))))
with( samp_paradel_A, with( Samps, table( 
    (Y-A)[MHSPs[,2]] - (Y-A)[MHSPs[,1]]
  )))

lglk_paradel_A <- generic_lglk_skippy # a copy that will be taught about paradel_A data
# Data is just like Delfinus mendax, Acme Archipelago--- "only" lglk might be different
env <- boring_data_prep_delfi_A( samp_paradel_A, prev_env=environment( lglk_paradel_A))
environment( lglk_paradel_A) <- env # now, lglk_<blah> will always know eg what...

