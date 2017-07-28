
camodel <- "cabone"

##' Get the locaton of model source code.
##' 
##' @export
cablib <- function() system.file(package="cabone")


##' Calcium / bone homeostatis model.
##' 
##' 
##' @param ... passed to update
##' @export
cabone <- function(...) {
  update(mread_cache(camodel,cablib()),...)
}


##' @rdname cabone
##' @export
cabone_export <- function(file=tempfile(fileext=".cpp"), overwrite=FALSE) {
  if(is.null(file)) stop("please provide a file name to write model code.", call.=FALSE)
  if(!grepl(".*\\.cpp$",file)) stop("file must end in '.cpp'",call.=FALSE)
  file <- normalizePath(file,mustWork=FALSE)
  if(file.exists(file) & !overwrite) stop("output file already exists.", call.=FALSE)
  mod <- mread(camodel,cablib(),compile=FALSE)
  message("Writing model code to file ", file)
  writeLines(mod@code,file)
  return(file)
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


##' Simulate with scler dosing
##' 
##' @param dose teriparatide dose in micrograms
##' @param ii dosing interval in hours
##' @param dur number of doses to simulate
##' @param delta simulation time grid
##' @param request outputs to request
##' 
##' @export
sim_scler <- function(dose=20, ii=24, dur=27, delta=0.1, request="PTHpm,CaC") {
  mod <- cabone()
  cmtn <- mrgsolve::cmtn(mod,"TERISC")
  data <- expand.ev(amt=amt_teri(dose), ii=ii, addl=dur, cmt=cmtn)
  mrgsim(mod, data=data, delta=delta, end=(dur+1)*ii, Req=request)
}

##' Simulate with denosumab dosing
##' 
##' @param dose denosumab dose in micrograms
##' @param ii dosing interval in months
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

##' Simulate secondary hyperparathyroidism
##' 
##' @param GFRdelta change in GFR from baseline value
##' @param GFRtau time interval in years over which GFR changes 
##' @param dur number of doses to simulate
##' @param delta simulation time grid in hours
##' @param request outputs to request
##' 
##' @export
sim_2h <- function(GFRdelta = 84, GFRtau = 10, delta=24, 
                   request="GFR,CaC,PTHpm,OC,ECCPhos",cfb=TRUE) {
  if(GFRtau <= 0 ) stop("GFRtau must be greater than zero.", call.=FALSE)
  mod <- cabone()
  ini <- as.list(init(mod))
  if(GFRdelta <= 0 | GFRdelta >= ini$GFR*16.66666667) {
    stop("GFRdelta out of bounds.", call.=FALSE) 
  }
  out <- 
    mod %>%
    param(GFRdelta = GFRdelta, GFRtau = GFRtau) %>%
    mrgsim(end=24*7*52*GFRtau, delta=delta, Req=request,tscale=1/(24*7*52))
  
  if(!cfb) return(out)
  
  base <- dplyr::slice(out,1)
  vars <- c(out@request,out@outnames)
  base <- as.list(base[,vars])
  data <- as.data.frame(out)
  for(i in seq_along(vars)) {
    data[,vars[i]] <- data[,vars[i]]/base[[vars[i]]]
  }
  out@data <- data
  out
}


