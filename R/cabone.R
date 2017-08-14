
camodel <- "cabone"
camodel_scler <- "cabone2"

##' Plotting helper functions.
##' @param ... passed to \code{ggplot::scale_color_brewer}
##' @export
.colSet1 <- function(...) ggplot2::scale_color_brewer(palette="Set1",...)
##' @export
##' @rdname .colSet1 
.colSet2 <- function(...) ggplot2::scale_color_brewer(palette="Set2",...)

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

##' @export
##' @rdname cabone
cabone_scler <- function(...) {
  update(mread_cache(camodel_scler,cablib()),...)
}

##' @rdname cabone
##' @export
cabone_export <- function(file=tempfile(fileext=".cpp"), overwrite=FALSE) {
  if(is.null(file)) {
    stop("please provide a file name to write model code.", call.=FALSE)
  }
  if(!grepl(".*\\.cpp$",file)) {
    stop("file must end in '.cpp'",call.=FALSE)
  }
  file <- normalizePath(file,mustWork=FALSE)
  if(file.exists(file) & !overwrite) {
    stop("output file already exists.", call.=FALSE)
  }
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

##' Convert sclerostin doses
##' 
##' @param x sclerostin dose in milligrams
##' 
##' @return sclerostin dose in nmol
##' @export
amt_scler <- function(x) x/0.145 

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

##' Simulate with sclerostin dosing
##' 
##' @param dose sclerostin dose in milligrams
##' @param ii dosing interval in months
##' @param dur number of doses to simulate
##' @param delta simulation time grid
##' @param request outputs to request
##' @param tscale factor for rescaling time in simulated output
##' 
##' @export
sim_scler <- function(dose=210, ii=1*24*28, dur=12, delta=24*28, 
                      request="SOSTCP, lsBMDsimSCLER, OBchange") {
  mod <- cabone_scler()
  cmtn <- mrgsolve::cmtn(mod,"SOSTSC")
  data <- expand.ev(amt=amt_scler(dose), ii=ii, addl=dur, cmt=cmtn)
  mrgsim(mod, data=data, delta=delta, end=(dur+1)*ii, Req=request)
}

##' Simulate with denosumab dosing
##' 
##' @param dose denosumab dose in milligrams
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

##' This is a function to re-create plots from Eudy, et al.(2015), CPT:PSP; Figure 3 ## 
##' 
##' 
##' @examples
##' 
##' sims <- sim_scler_data()
##' 
##' \dontrun{
##'   require(ggplot2)
##' 
##'   ggplot(data=sims, aes(time,P1NPsim)) + 
##'     geom_line() + facet_wrap(~label)
##' 
##'   ggplot(data=sims, aes(time,CTXsim)) + 
##'     geom_line() + facet_wrap(~label)
##'   
##'   ggplot(data=sims, aes(time,lsBMDsimSCLER)) + 
##'     geom_line() + facet_wrap(~label)
##'   
##'   ggplot(data=sims, aes(time,thBMDsimSCLER)) + 
##'     geom_line() + facet_wrap(~label)
##'    
##' }
##' 
##' @export
##' 
sim_scler_data <- function(){
  
  .Req <- c("P1NPsim","CTXsim","lsBMDsimSCLER","thBMDsimSCLER")
  .end <- 8787
  .delta <- 2*24
  
  mod <- cabone_scler()
  
  ev1 <- ev(ID=1, amt=0, cmt="SOSTSC", ii=24*14, addl=25)
  out1 <- mrgsim(mod, ev1, end=.end, delta=.delta, Req=.Req)
  dat1 <- as.data.frame(out1)  %>% dplyr::mutate(label= "PBO Q2W")
  
  ev2 <- ev(ID=2, amt=1244.555, cmt="SOSTSC", ii=24*28, addl=11)
  out2 <- mrgsim(mod, ev2, end=.end, delta=.delta, Req=.Req)
  dat2 <- as.data.frame(out2) %>% dplyr::mutate(label= "180 mg Q4W")
  
  ev3 <- ev(ID=3, amt=1244.555, cmt="SOSTSC", ii=24*14, addl=25)
  out3 <- mrgsim(mod, ev3, end=.end, delta=.delta, Req=.Req)
  dat3 <- as.data.frame(out3) %>% dplyr::mutate(label= "180 mg Q2W")
  
  ev4 <- ev(ID=4, amt=1866.83, cmt="SOSTSC", ii=24*14, addl=25)
  out4 <- mrgsim(mod, ev4, end=.end, delta=.delta, Req=.Req)
  dat4 <- as.data.frame(out4) %>% dplyr::mutate(label= "270 mg Q2W")
  
  dplyr::bind_rows(dat1,dat2,dat3,dat4)
  
}

##' This is a function to simulate combination arms (ACOP 2015 poster) ## 
##' 
##' @examples
##' 
##' \dontrun{
##'   require(ggplot2)
##'   
##'   sims <- sim_combo_arms()
##'   
##'   ggplot(sims, aes(time,lsBMD,col=regimen)) + 
##'     geom_line(lwd=1) + theme(legend.position="top")
##' }
##' 
##' @export
sim_combo_arms <- function() {
  
  mod <- cabone_scler()
  
  events1 <- ev(amt = 4856.962,cmt="TERISC",ii=24, addl= 671) 
  events2 <- ev(amt = 6E7,cmt="DENSC",ii=4032, addl= 3)
  events3 <-  events1 + events2
  
  .end <- 16128
  .delta <- 24*7
  
  message("Simulating teriparatide data ...")
  out1 <- 
    mod %>% 
    ev(events1) %>% 
    mrgsim(end=.end, Req="lsBMDsimTERI", delta=.delta) %>%
    dplyr::mutate(regimen = "teriparatide", 
                  lsBMD = lsBMDsimTERI,
                  output = "lsBMDsimTERI", 
                  lsBMDsimTERI = NULL) 
  
  message("Simulating denosumab data ...")
  out2 <- 
    mod %>% 
    ev(events2) %>% 
    mrgsim(end=.end,Req="lsBMDsimDEN", delta=.delta) %>%
    dplyr::mutate(regimen = "denosumab", 
                  lsBMD = lsBMDsimDEN, 
                  output = "lsBMDsimDEN", 
                  lsBMDsimDEN = NULL)
  
  message("Simulating combination data ...")
  out3 <- 
    mod %>% 
    ev(events3) %>% 
    param(DEN_TERI_COMBO=1) %>% 
    mrgsim(end=.end, Req="lsBMDsimCOMBO", delta=.delta) %>%
    dplyr::mutate(regimen = "teri+denos", 
                  lsBMD = lsBMDsimCOMBO,
                  output = "lsBMDsimCOMBO", 
                  lsBMDsimCOMBO = NULL)
  
  message("Done.")
  out <- dplyr::bind_rows(out1,out2,out3)
  
  dplyr::mutate(out, regimen = factor(regimen, levels=unique(regimen)))
}
