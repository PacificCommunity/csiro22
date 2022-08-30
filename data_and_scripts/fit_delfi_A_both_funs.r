"boring_data_prep_delfi_A" <-
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

  # For THIS version, only look at FEMALE samples, drop all HSPs, and restrict
  # ... to adult-juve comps at time of sampling

  sampo <- subset_samples( sampo, Sex=='F') # drops all POPs/HSPs with any non-Female member
  extract.named( sampo[ cq( Samps, POPs)]) # drop HSPs ...

  extract.named( Samps)
  extract.named( sampo@public) # Amat
  AMAX <- lastel( C_sya, 3)
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

  # m_... is samp size
  m_ad_YA <- offarray( table( Y=Y[ parposs], A=A[ parposs]))
  m_ju_B <- offarray( table( B=B[ offposs]))

  # Number of comparisons: product of sample sizes by category
  # See ?offarray::noloop  or code of lglk_aggregate() below
  n_comp_MOP_BYA <- autoloop( Bju=Bju_range, Yad=Yad_range, Aad=Aad_range, {
      Bad <- Yad - Aad
      # Only do comps where ju is born after adult
      # ... which also avoids double-counting and self-comparisons
      m_ju_B[ Bju] * m_ad_YA[ Yad, Aad] * (Bju > Bad)
    })

  n_MOP_BYA <- offarray( table(
      Bju=B[ POPs[,2]],
      Yad=Y[ POPs[,1]],
      Aad=A[ POPs[,1]]),
      template=n_comp_MOP_BYA) # template ensures full ranges used, even if no entries for some values

  # That's what's needed for fitting...
  # ... next only for "pedagogical" purposes (data inspection)
  POP_df <- boring_dfize( n_MOP_BYA, n_comp_MOP_BYA)

  # copy useful vars into envo... R magic, just trust me on this one ;)
  list2env( mget( cq( n_comp_MOP_BYA, n_MOP_BYA,
      Bju_range, Yad_range, Aad_range, # in sample
      AMAX,
      POP_df)),
      envo)

return( envo)


}




"generic_lglk_MOP_ideal_mammal" <-
function( params){
  assign( 'last_params', params, environment( sys.function())) # debugging etc

  ## Unpack parameters
  Nfad_y0 <- exp( params[ 1])
  RoI <- params[ 2]
  # ... that's all,  folks

  ## Population dynamics

  N_Y <- offarray( 0, dimseq=list( Y=years))
  N_Y[] <- Nfad_y0 * exp( RoI * (years-y0))

  # Ideal kinship probs, ie if we knew all important covariates......
  Pr_MOP_BYA <- autoloop( Bju=Bju_range, Yad=Yad_range, Aad=Aad_range, {
    Bad <- Yad - Aad
    ( Bju <= Yad ) *          # was ad still alive at B[ju] ?
    ( Bju >= Amat + Bad ) *   # was ad mature at B[ju] ?
    ( 1/N_Y[ Bju])              # ad's ERO at B[ju] if alive & mature, divided by TRO that year
  })

  # ...... which in this case, we DO! so no need to corrupt this to observables

  # Housekeeping: useful to keep some stuff after function exits. Trust me, this incantation works...
  list2env( mget( cq( Nfad_y0, RoI, N_Y, Pr_MOP_BYA)), envir=environment( sys.function()))

  lglk <- MAKE_FINITE( # anti-infinity guard
      sum( dpois( n_MOP_BYA, lambda=n_comp_MOP_BYA * Pr_MOP_BYA, log=TRUE))
    )

return( lglk)
}
