library( mvbutils)  # various low-level utilities, eg cq()
library( offarray)  # arrays/vectors don't have to start at 1 any more! Freedom!
library( debug)     # so we can watch the fun
library( atease)
source("ckmr_funs.r")

print( load( 'samp_delfi_A.rda')) # object samp_notog1


"boring_data_prep_delfi_A_both" <-
function( sampo, prev_env=parent.env()){
## Given "CKMR sample data", probably from 'prep_from_sim2'...
## Make arrays with n_comps and n_kin, per kinship and per covars
## Stick them into child env of prev_env--- which means you
## can subsequently associate them with a lglk function
## so it will "know" about its specific data

## usage:
# lglk_mydata <- generic_lglk_for_this_kinda_CKMR_data  # a function
# env <- boring_data_prep_mydata( my_samps, prev_env=environment( lglk_mydata))
# environment( lglk_mydata) <- env
# lglk_mydata( params) # lglk_mydata can now refer internally to 'n_MOP' etc

## Warning: this is boring. Did the name not tip you off? You almost
## certainly do NOT need to understand what's in it. And there's esoteRica
## which I am NOT going to explain. So you should probably stop reading RIGHT NOW...

  extract.named( sampo[ cq( Samps, POPs)]) # drop HSPs ...

  extract.named( Samps)
  extract.named( sampo@public) # Amat
    AMAX <- lastel( C_sya, 3) # maximum age of adults
    
  # Year of birth
  B <- Y - A

  parposs <- A >= Amat
  offposs <- A < Amat

  # Only AJ POPs
  POPs <- POPs[ parposs[ POPs[,1]] & offposs[ POPs[,2]],]

  rPOPs <- pairid( POPs[,1], POPs[,2]) # combine into single real number

  # Package up stuff, to be used as environment for lglk function
  y0 <- min( B[ offposs]) # SHOULDN'T really be data-driven
  years <- min( B[ offposs]) %upto% max( Y[ parposs])
  envo <- list2env( mget( cq( rPOPs, B, Y, A, Amat, y0, years, offposs, parposs)), parent=prev_env)

  ## Stuff for aggregated version:
  # NB things will be evaluate over entire range, even if gaps
  # Dodge-able at C level (ie gaps could be handled), but maybe not in R
  Bju_range <- min( B[ offposs]) %upto% max( B[ offposs])
  Yad_range <- min( Y[ parposs]) %upto% max( Y[ parposs])
    Aad_range <- min( A[ parposs]) %upto%  max( A[ parposs])
    Sad_range <- c('M', 'F')

  # m_... is samp size
  m_ad_YAS <- offarray( table( Y=Y[ parposs], A=A[ parposs], S=Sex[ parposs]))
  m_ju_B <- offarray( table( B=B[ offposs]))

  # Number of comparisons: product of sample sizes by category
  # See ?offarray::noloop  or code of lglk_aggregate() below
    n_comp_POP_BYAS <- autoloop( Bju=Bju_range, Yad=Yad_range,
                               Aad=Aad_range, Sad=Sad_range, {
      Bad <- Yad - Aad
      # Only do comps where ju is born after adult
      # ... which also avoids double-counting and self-comparisons
      m_ju_B[ Bju] * m_ad_YAS[ Yad, Aad, Sad] * (Bju > Bad)
    })

  n_POP_BYAS <- offarray( table(
      Bju=B[ POPs[,2]],
      Yad=Y[ POPs[,1]],
      Aad=A[ POPs[,1]],
      Sad=Sex[ POPs[,1]]),
      template=n_comp_POP_BYAS) # template ensures full ranges used, even if no entries for some values

  # That's what's needed for fitting...
  # ... next only for "pedagogical" purposes (data inspection)
  POP_df <- boring_dfize( n_POP_BYAS, n_comp_POP_BYAS)

  # copy useful vars into envo... R magic, just trust me on this one ;)
  list2env( mget( cq( n_comp_POP_BYAS, n_POP_BYAS,
      Bju_range, Yad_range, Aad_range, Sad_range, # in sample
      AMAX,
      POP_df)),
      envo)

return( envo)


}


"generic_lglk_POP_ideal_mammal_both" <-
function( params){
  assign( 'last_params', params, environment( sys.function())) # debugging etc

  ## Unpack parameters
  Nad_y0 <- exp( params[ 1])
  RoI <- params[ 2]
  srM <- inv.logit( params[ 3]) # sex ratio for males
  # ... that's all,  folks

  ## Population dynamics
  N_YS <- offarray( 0, dimseq=list( Y=years, S=Sad_range))
  N_YS[, 'M'] <- Nad_y0 * srM * exp( RoI * (years-y0))
  N_YS[, 'F'] <- Nad_y0 * (1-srM) * exp( RoI * (years-y0))  

  # Ideal kinship probs, ie if we knew all important covariates......
  Pr_POP_BYAS <- autoloop( Bju=Bju_range, Yad=Yad_range, Aad=Aad_range, Sad=Sad_range, {
    Bad <- Yad - Aad
    ( Bju <= Yad ) *          # was ad still alive at B[ju] ?
    ( Bju >= Amat + Bad ) *   # was ad mature at B[ju] ?
    ( 1/N_YS[ Bju, Sad])              # ad's ERO at B[ju] if alive & mature, divided by TRO that year
  })

  # ...... which in this case, we DO! so no need to corrupt this to observables

  # Housekeeping: useful to keep some stuff after function exits. Trust me, this incantation works...
  list2env( mget( cq( Nad_y0, RoI, srM, N_YS, Pr_POP_BYAS)), envir=environment( sys.function()))

  lglk <- MAKE_FINITE( # anti-infinity guard
      sum( dpois( n_POP_BYAS, lambda=n_comp_POP_BYAS * Pr_POP_BYAS, log=TRUE))
    )

return( lglk)
}

# a copy that will be taught about delfi_A data
lglk_delfi_A_both <- generic_lglk_POP_ideal_mammal_both 
envb <- boring_data_prep_delfi_A_both( samp_delfi_A, prev_env=environment( lglk_delfi_A_both))
environment( lglk_delfi_A_both) <- envb # now, lglk_<blah> will always know eg what...
# ... catches, MOPs, comparisons etc were, without having to pass that stuff in every time it's called

starto <- c( log( 100000), 0.01, logit(0.45)) # logNad and RoI
lglk_delfi_A_both( starto) # just check it doesn't crash!

fitto <- nlminb( starto, NEG( lglk_delfi_A_both))

Amat <- envb$Amat
tru_Nadf <- sumover( samp_delfi_A@secret$N_sya[ , envb$years, Amat:40], c('A', 'SEX'))
pred_Nadf <- sumover(envb$N_YS, 'S')

yl <- max(c(tru_Nadf, pred_Nadf))
plot(envb$years, tru_Nadf, type='l', ylim=c(0, yl))
points(envb$years, pred_Nadf, pch=19, col='dodgerblue')
