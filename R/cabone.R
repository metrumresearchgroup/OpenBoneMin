
camodel <- "cabone2"

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

##' Convert sclerostin doses
##' 
##' @param x sclerostin dose in milligrams
##' 
##' @return sclerostin dose in nmol
##' @export
amt_scler <- function(x) x/0.145 


##' Simulate with teriparatide dosing
##' 
##' @param dose teriparatide dose in milligrams
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
  mod <- cabone()
  cmtn <- mrgsolve::cmtn(mod,"SOSTSC")
  data <- expand.ev(amt=amt_scler(dose), ii=ii, addl=dur, cmt=cmtn)
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

##' This is a function to re-create plots from Eudy, et al.(2015), CPT:PSP; Figure 3 ## 
##' 
##' @param endpoint ## can be P1NP, CTX, lsBMD, or thBMD
##' @export
##' 
## Bloso PBO ##
create_scler_plots<- function (endpoint){
  if (!endpoint %in% c("P1NP","CTX","lsBMD","thBMD")) stop("Endpoint must be P1NP, CTX, lsBMD, or thBMD.", call.=FALSE)
  
  mod <- cabone()
  out <- mrgsim(mod, ev(ID=1, amt=0, cmt="SOSTSC", time=0, ii=24*14, addl=25), end=8784, delta=2*24,Req=c("P1NPsim","CTXsim","lsBMDsimSCLER","thBMDsimSCLER"))
  
  ## Bloso 180 mg Q4W SQ ##
  
  out2 <- mrgsim(mod, ev(ID=2, amt=1244.555, cmt="SOSTSC", time=0, ii=24*28, addl=11), end=8784, delta=2*24,Req=c("P1NPsim","CTXsim","lsBMDsimSCLER","thBMDsimSCLER"))
  
  ## Bloso 180 mg Q2W SQ ##
  
  out3 <- mrgsim(mod, ev(ID=3, amt=1244.555, cmt="SOSTSC", time=0, ii=24*14, addl=25), end=8784, delta=2*24,Req=c("P1NPsim","CTXsim","lsBMDsimSCLER","thBMDsimSCLER"))
  
  ## Bloso 270 mg Q2W SQ ##
  
  out4 <- mrgsim(mod,  ev(ID=4, amt=1866.83, cmt="SOSTSC", time=0, ii=24*14, addl=25), end=8784, delta=2*24,Req=c("P1NPsim","CTXsim","lsBMDsimSCLER","thBMDsimSCLER"))
  
  
  dat<-as.data.frame(out) %>% mutate(label= "PBO Q2W")
  dat2<-as.data.frame(out2) %>% mutate(label= "180 mg Q4W")
  dat3<-as.data.frame(out3) %>% mutate(label= "180 mg Q2W")
  dat4<-as.data.frame(out4) %>% mutate(label= "270 mg Q2W")
  
  allsim1<-rbind(dat,dat2,dat3,dat4)
  allsim1$evid<-0
  
  # dosing records # 
  
  id1<-data.frame(time=seq(0,8064,24*28),ID=1)
  id2<-data.frame(time=seq(0,8064,24*28),ID=2)
  id3<-data.frame(time=seq(0,8064,24*14),ID=3)
  id4<-data.frame(time=seq(0,8064,24*14),ID=4)
  
  dosing<-rbind(id1,id2,id3,id4)
  dosing$evid<-1
  
  sost54sim<-merge(allsim1,dosing, all=T)
  
  
  sost54dat<-read.table(file="./R/sost_54.csv", header=TRUE, sep=",", as.is=T, na.strings = '.')
  sost54dat$evid<-0
  
  alldata<-merge(sost54sim,sost54dat,by=c("ID","time","evid","label"),all=T)
  
  
  alldata$label<-with(alldata,reapply(label,ID,locf))
  alldata$label <- factor(alldata$label, levels = c("PBO Q2W","180 mg Q4W", "180 mg Q2W", "270 mg Q2W"))
  sost54dat$label <- factor(sost54dat$label, levels = c("180 mg Q4W", "180 mg Q2W", "270 mg Q2W","PBO Q2W"))
  
  alldata$thBMDub[alldata$thBMDub==100.00]<-NA
  
  ## P1NP
  
  cols <- c("Simulated"="blue","Data"="firebrick")
  P1NP_plot<- ggplot(data=subset(alldata,evid==0), aes(x=time/24/28)) + 
    facet_wrap(~label,scales="free_x") + 
    scale_x_continuous(name= "time (months)", breaks=c(0,3,6,12)) + ylab("% of baseline") + 
    theme(strip.text.x = element_text(size=11), axis.text=element_text(colour="grey20",size=12),
          legend.position="top", legend.text=element_text(size=15), axis.title=(element_text(size=15))) +
    geom_line(aes(y=P1NPsim, group=1,colour="Simulated"), size = 1) +
    geom_point(aes(y=P1NP,colour="Data"), size=3,na.rm=TRUE) +
    geom_rug(data=subset(alldata,evid==1),mapping = aes(x=time/24/28), color='blue') +
    labs(title = "P1NP +/- 95% CI") + theme(plot.title=element_text(size=15, vjust=1.5)) +
    geom_errorbar(data=sost54dat,aes(x=time/24/28, ymax= P1NPub,ymin = P1NPlb),width=0.2) +
    scale_colour_manual(name="",values=cols) 
  
  ## CTX
  CTX_plot<- ggplot(data=subset(alldata,evid==0), aes(x=time/24/28)) +
    facet_wrap(~label,scales="free_x") + 
    scale_x_continuous(name= "time (months)", breaks=c(0,3,6,12)) + ylab("% of baseline") + 
    theme(strip.text.x = element_text(size=11), axis.text=element_text(colour="grey20",size=12),
          legend.position="top", legend.text=element_text(size=15), axis.title=(element_text(size=15))) +
    geom_line(aes(y=CTXsim,group=1,colour="Simulated"), lwd=1) +
    geom_point(aes(y=CTX, colour="Data"), size=3,na.rm=TRUE) +
    geom_rug(data=subset(alldata,evid==1),mapping = aes(x=time/24/28), color='blue') +
    labs(title = "CTx +/- 95% CI") + theme(plot.title=element_text(size=15, vjust=1.5)) +
    geom_errorbar(data=sost54dat,aes(x=time/24/28,ymax= CTXub,ymin = CTXlb),width=0.2) +
    scale_colour_manual(name="",values=cols) 
  
  ## lsBMD
  
  lsBMD_plot <- ggplot(data=subset(alldata, evid==0), aes(x=time/24/28)) +
    facet_wrap(~label,scales="free_x") + 
    scale_x_continuous(name= "time (months)", breaks=c(0,3,6,12)) + ylab("% of baseline") + 
    theme(strip.text.x = element_text(size=11), axis.text=element_text(colour="grey20",size=12),
          legend.position="top", legend.text=element_text(size=15), axis.title=(element_text(size=15))) +
    geom_line(aes(y=lsBMDsimSCLER, colour="Simulated"), lwd = 1) +
    geom_point(aes(y=lsBMD, colour="Data"), size=3, na.rm=TRUE) +
    geom_rug(data=subset(alldata,evid==1),mapping = aes(x=time/24/28), color='blue') +
    labs(title = "Lumbar Spine BMD +/- 95% CI") + theme(plot.title=element_text(size=15, vjust=1.5)) +
    geom_errorbar(data=sost54dat,aes(x=time/24/28,ymax= lsBMDub,ymin = lsBMDlb),width=0.2) +
    scale_colour_manual(name="",values=cols) 
  
  ## thBMD
  
  thBMD_plot <-ggplot(data=subset(alldata, evid==0), aes(x=time/24/28)) +
    facet_wrap(~label,scales="free_x") +   facet_wrap(~label,scales="free_x") + 
    scale_x_continuous(name= "time (months)", breaks=c(0,3,6,12)) + ylab("% of baseline") + 
    theme(strip.text.x = element_text(size=11), axis.text=element_text(colour="grey20",size=12),
          legend.position="top", legend.text=element_text(size=15), axis.title=(element_text(size=15))) +
    geom_line(aes(y=thBMDsimSCLER, colour="Simulated"), size = 1) +
    geom_point(aes(y=thBMD, colour="Data"), size=3, na.rm=TRUE) +
    geom_rug(data=subset(alldata,evid==1),mapping = aes(x=time/24/28), color='blue') +
    labs(title = "Total Hip BMD +/- 95% CI") + theme(plot.title=element_text(size=15, vjust=1.5)) +
    geom_errorbar(data=sost54dat,aes(x=time/24/28,ymax= thBMDub,ymin = thBMDlb),width=0.2) +
    scale_colour_manual(name="",values=cols) 
  
    if (endpoint == "P1NP") return ( P1NP_plot )
    if (endpoint == "CTX") return ( CTX_plot )
    if (endpoint == "lsBMD") return ( lsBMD_plot )
    if (endpoint == "thBMD") return ( thBMD_plot )
}

##' This is a function to simulate combination arms (ACOP 2015 poster) ## 
##' @export
##' @return plot of combination arms 

sim_combo_arms<- function() {
  
  mod <- cabone()
  
  deno<-read.table(file="./R/deno_FINAL_check(5).csv", header=TRUE, sep=",", as.is=TRUE, na.strings='.')
  deno7<-subset(deno, ID==7)
  teri<-read.table(file="./R/teri_data8.csv", header=TRUE, sep=",", as.is=T, na.strings = '.')
  combo_dat<-read.table(file="./R/Leder_study_combination_arm.csv", header=TRUE, sep=",", as.is=T, na.strings = '.')
  teri34<-subset(teri, ID==34)
  
  
  events1 <- ev(amt = 4856.962,cmt="TERISC",ii=24, addl= 671) 
  events2 <- ev(amt = 6E7,cmt="DENSC",ii=4032, addl= 3)
  events3 <-  events1 + events2
  
  
  ## this code will simulate a distribution but it takes a long time. ##
  
  # idata <- within(data.frame(ID=1:100), {
  #   koutBMDlsCOMBO <- exp(log(0.000137) + rnorm(100, 0,sqrt(0.0277))) # ~15% CV
  #   gamOClsDEN_TERI <- exp(log(0.102) + rnorm(100,0,sqrt(0.0277))) # ~15% CV 
  #   DEN_TERI_COMBO<-1
  # })
  # 
  
  # idata_split<-split(idata,idata$ID)
  # out <- mclapply(seq_along(idata_split), mc.cores=8, function(i) {
  #   events <- events3 
  #   out <- try(mod %>% ev(events) %>% idata_set(idata_split[[i]]) %>% mrgsim(end=16128, Req="lsBMDsimCOMBO", carry.out=names(idata_split[[i]])))
  #   if(inherits(out, "try-error")) return(list(irep=i))
  #   #out <- label(out, irep=i)
  #   #out[out$time < 1008 | out$time > 0,] #
  #   as.tbl(out) %>% mutate(irep=i)
  # }) %>% bind_rows
  # 
  ## Simulate combination arm ##
  events<-events3
  out2 <- mod %>% ev(events) %>% param(DEN_TERI_COMBO=1) %>% mrgsim(end=16128, Req="lsBMDsimCOMBO", delta=24*7) 
  sim_mean_combo<-as.data.frame(out2)
  
  ### calculate 95% CI of simulated data ### 
  # mat <- tapply(out$lsBMDsimCOMBO, INDEX=out$time, stats::quantile, probs = c(0.05,0.5,0.95), simplify=T)
  # 
  # new_mat<-do.call("rbind",mat)
  # new_mat<-as.data.frame(new_mat)
  # time<-out$time[out$ID==1]
  # new_mat$time<-time[3:length(time)]
  # names(new_mat)<-c("LB","MED","UB","time")
  # new_mat_combo<-new_mat
  # sims_combo<-out
  
  ## simulate teriparatide only arm ## 
  out1 <- mod %>% ev(events1) %>% mrgsim(end=16128, Req="lsBMDsimTERI", delta=24*7)
  sim_mean_teri<-as.data.frame(out1)
  ## simulate denosumab only arm ## 
  out3 <- mod %>% ev(events2) %>% mrgsim(end=16128,Req="lsBMDsimDEN", delta=24*7)
  den_output2<-as.data.frame(out3)
  
  ## plot ## 
  dodge <- position_dodge(.1)
  ltyp <- c("combination"="#56B4E9","teriparatide"="#000000","denosumab"="#E69F00")
  ggplot() +
    geom_line(aes(x=time/24/28, y=lsBMDsimCOMBO, group=1, colour="combination"), data=subset(sim_mean_combo,!is.na(lsBMDsimCOMBO)), size=0.5) +
    theme(strip.text.x = element_text(size=11), axis.text=element_text(colour="grey20",size=12),
          legend.position="top", legend.text=element_text(size=15), axis.title=(element_text(size=15))) +
    scale_x_continuous(name= "time (months)", breaks=c(0,3,6,12,18,24,36)) + ylab("% of baseline") +
   
    geom_point(aes(x=time/24/28, y=lsBMD, group=1), data=subset(combo_dat,ID %in% unique(combo_dat$ID[!is.na(combo_dat$lsBMD)])), 
               size=2,colour="#56B4E9",position = position_jitter(w = 0), na.rm=TRUE) +
    labs(title = "lumbar spine BMD  +/- 90% CI \n 2 year treatment") + theme(plot.title=element_text(size=15, vjust=1.5)) +
    geom_errorbar(aes(x=time/24/28, ymax= lsBMDub,ymin = lsBMDlb),data=subset(combo_dat,ID %in% unique(combo_dat$ID[!is.na(combo_dat$lsBMD)])),width=0.2,na.rm=TRUE) +
    #geom_ribbon(aes(ymin=LB, ymax=UB), data=new_mat_combo, fill="#56B4E9", alpha=0.4) 
    scale_colour_manual(name="",values=ltyp) +
    
    ## teri
    geom_line(aes(x=time/24/28, y=lsBMDsimTERI, group=2, colour="teriparatide"), data=sim_mean_teri, size=0.5) +
    geom_point(aes(x=time/24/28, y=lsBMD, group=2, colour="teriparatide"), data=subset(teri34,ID %in% unique(teri34$ID[!is.na(teri34$lsBMD)])), size=2,na.rm=TRUE) +
    geom_errorbar(aes(x=time/24/28, ymax= lsBMDub,ymin = lsBMDlb),data=subset(teri34,ID %in% unique(teri34$ID[!is.na(teri34$lsBMD)])),width=0.2,na.rm=TRUE) +
   # geom_ribbon(aes(ymin=LB, ymax=UB), data=new_mat_teri, fill="#000000", alpha=0.4) +
    ## deno
    geom_line(aes(x=time/24/28, y=lsBMDsimDEN, group=3, colour="denosumab"), data=subset(den_output2,!is.na(lsBMDsimDEN)), size=0.5) +
    geom_point(aes(x=time/24/28, y=BMD.lumbar.spine_dv, group=3, colour="denosumab"), data=subset(deno7,ID %in% unique(deno7$ID[!is.na(deno7$BMD.lumbar.spine_dv)])), size=2,na.rm=TRUE) +
    geom_errorbar(aes(x=time/24/28, ymax= lsBMDub,ymin = lsBMDlb),data=subset(deno7,ID %in% unique(deno7$ID[!is.na(deno7$BMD.lumbar.spine_dv)])),width=0.2,na.rm=TRUE)
   # geom_ribbon(aes(ymin=LB, ymax=UB), data=new_mat_deno2, fill="#E69F00", alpha=0.4) 
  }
