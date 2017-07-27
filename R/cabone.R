
loc <- "/Users/kyleb/git/metrumresearchgroup/cabone/inst"
camodel <- "cabone"

##' Get the locaton of model source code.
##' 
##' @export
cablib <- function() system.file(package="cabone")


##' Load calcium / bone model.
##' 
##' @param ... passed to update
##' 
##' @export
cabone <- function(...) {
  update(mread_cache(camodel,cablib()),...)
}



##' Convert teriparatide doses
##' 
##' @param x teriparatide dose in micrograms
##' 
##' @return teriparatide dose in 
##' @export
amt_teri <- function(x) x*1E6/4117.8

##' Convert denosumab doses
##' 
##' @param x denosumab dose in milligrams
##' 
##' @return denosumab dose in mmol
##' @export
amt_denos <- function(x) x*1


##' Simulate with teriparatide dosing
##' 
##' @param dose teriparatide dose in micrograms
##' @param ii dosing interval in hours
##' @param dur number of doses to simulate
##' @param delta simulation time grid
##' @param request outputs to request
##' 
##' @export
sim_teri <- function(dose=20, ii=24, dur=27, delta=0.1, request="PTHpm,CaC") {
  mod <- cabone()
  cmtn <- mrgsolve::cmtn(mod,"TERISC")
  data <- expand.ev(amt=amt_teri(dose), ii=ii, addl=dur, cmt=cmtn)
  mrgsim(mod, data=data, delta=delta, end=(dur+1)*ii, Req=request)
}

##' Simulate with denosumab dosing
##' 
##' @param dose denosumab dose in micrograms
##' @parma ii dosing interval in months
##' @param dur number of doses to simulate
##' @param delta simulation time grid in hours
##' @param request outputs to request
##' @param tscale factor for rescaling time in simulated output
##' 
##' @export
sim_denos <- function(dose=60, ii=6, dur=3, delta=4, 
                      request="DENCP,DENMOL,BMDlsDENchange", 
                      tscale=1/(24*28)) {
  mod <- cabone()
  cmtn <- mrgsolve::cmtn(mod,"DENSC")
  ii <- ii*28*24
  data <- expand.ev(amt=amt_denos(dose), ii=ii, addl=dur-1, cmt=cmtn)
  mrgsim(mod, data=data, delta=delta, end=(dur+1)*ii,
         tscale=tscale,
         Req=request)
}




