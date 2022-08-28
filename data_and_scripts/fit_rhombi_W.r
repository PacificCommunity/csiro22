## Cobalt Squarehad (Rhombichthys cyanorosea) from Whyalla

library( mvbutils)  # various low-level utilities, eg cq()
library( offarray)  # arrays/vectors don't have to start at 1 any more! Freedom!
library( debug)     # so we can watch the fun

### ACTION STARTS HERE

## Load and set up the data
od <- options( digits=7)
on.exit( options( od))
print( load( 'samp_rhombi_W.rda')) # object samp_notog1
print( names( samp_rhombi_W))
print( head( samp_rhombi_W$Samps))
scatn( 'Rhombi W #POPs = %i', nrow( samp_rhombi_W$POPs))

samp_rhombi_W <- "READ THE DAMN SCRIPT AND EDIT APPROPRIATELY; DON'T JUST RUN IT"

# Summarize sampling & kin & other data for rhombi_W, and teach the lglk function about it
lglk_rhombi_W <- generic_lglk_Covar # a copy that will be taught about rhombi_W data
env <- boring_data_prep_rhombi_W( samp_rhombi_W, prev_env=environment( lglk_rhombi_W))
environment( lglk_rhombi_W) <- env # now, lglk_<blah> will always know about data

starto <- c( rep( log( 100000), 2), 0.5) #
lglk_rhombi_W( starto) # just check it doesn't crash!

fitto <- nlminb( starto, NEG( lglk_rhombi_W))

# How did we go?
cat( 'Abund:\n')
rbind( Est=env$Nad_S, Tru=samp_rhombi_W@secret$Nad_S)

cat( 'Fec\n')
rbind( Est=with( env, fec_Co_m / fec_Co_m[3]),
    Tru=with( samp_rhombi_W@secret, fec_Co / fec_Co[3])
  )

# Well Done Us!!! Go Us!!!

# How did Bozo et al. get on?
# Bozo to student: "_You_ have to get that PhD chapter published this week,
# ... else _we_ don't receive the Government $ ..."
# "Never mind the different sampling, you can write that up another time, you'll get another paper"
# Analyse data as per Derwent (unselective) samples
lglk_Bozo <- generic_lglk_cartoon
env <- boring_data_prep_rhombi_D( samp_rhombi_W, prev_env=environment( lglk_Bozo))
environment( lglk_Bozo) <- env
starto <- rep( log( 100000), 2)
lglk_Bozo( starto)
fitto <- nlminb( starto, NEG( lglk_Bozo))
cat( 'Est:\n')
env$Nad_S # estimated
cat( 'Tru:\n')
samp_rhombi_W@secret$Nad_S

# Easy way to get the Hessian: use optim()
fitto <- optim( starto, NEG( lglk_wrongo), hessian=TRUE)
fitto$hessian # NB the NEG has made this POSitive
solve( fitto$hessian)
#         [,1]     [,2]
#[1,] 0.004311 0.000000
#[2,] 0.000000 0.003745

# Males were 2nd, so...
sqrt( 0.003745)
fitto$par
# [1] 9.205 7.860

# 95% UCI from Bozo's model
exp( 7.86 + 2*0.06)
# [1] 2922

# Truth (3000) is outside 95% CI.
# But, never mind: Bozo gets the $$$! Happy times!


