$PROB
  
  //* Repository UUID: 18a99bec-3f78-43e9-9f4c-9c5f5529cddc
  //* Revision: 602
  //* Last Changed Rev: 596
  //* Last Changed Date: 2015-09-03 11:12:00 -0400 (Thu, 03 Sep 2015)
  //* Text Last Updated: 2016-10-07 16:23:01 -0400 (Fri, 07 Oct 2016)
  //* Checksum: bec0ec2a05e24199c69c59e306ca449dc736a6ab
  
  
  // 08-12-15 KTB
  //   * Converting to newest mrgsolve format
  // 08-24-15 KTB 
  //   * Fixed baseline Ctriol / Ctriol50 bug (see double INparenCtriol = ...)   
  
$GLOBAL 
#define max(a,b) ((a) > (b) ? (a) : (b))
#define F11 T85
#define PicOCkin Pic0
  
#define SETINIT if(NEWIND <=1) 
  
#define PKCP (PKCENT/PKVC)
#define DENCP (DENCENT/DENVC)
  
#define CaConc0 (P0/V1)
#define PTHconc0 (PTH0/V1)
#define OB (OBfast + OBslow)
  
#define PTHconc (PTH/V1)
#define CaConc (P/V1)
#define C1 (P/V1)
#define C2 (ECCPhos/V1)
#define C8 (B/V1)
#define D  ROB1
#define Osteoclast OC
#define OB0 (OBfast0 + OBslow0)
#define Calcitriol0 (B_0/V1)
  
  // KTB: Aug 12, 2015
  //   * These are aliases for different initial conditions
  //   * They should be eliminated at some point (e.g. find Q0 and replace with Q_0)
  //   
#define Q0 Q_0
#define RNK0 RNK_0
#define L0 L_0
#define ROB10 ROB1_0
#define Qbone0 Qbone_0
#define O0 O_0
#define RX20 RX2_0
#define CREB0 CREB_0
#define M0 M_0
#define TGFBact0 TGFBact_0
#define TGFB0 TGFB_0
#define PREPTH0 PREPTH_0
#define PTH0 PTH_0
#define B0 B_0 // M Riggs added 21Aug2015
#define ECCPhos0 ECCPhos_0
#define OBfast0 OBfast_0
#define OC0 OC_0
#define P0 P_0
#define T0 T_0
#define BMDfn0 BMDfn_0
#define BMDls0 BMDls_0
#define BMDlsDEN0 BMDlsDEN_0
#define EST0 EST_0
// #define GFR0 GFR_0
#define GFRdt0 GFRdt_0
#define GFRmin0 GFRmin_0  
#define OBslow0 OBslow_0
#define IntraPO0 IntraPO_0
  
$MAIN 
double V5878 = VI ; // THETA1*exp(ETA(1));
double CL5878 = CLI ; //THETA2*exp(ETA(2));
double KA5878 = KAI ; // THETA3*exp(ETA(3));
double K5878 = CL5878/V5878;

D_GUT5878 = D1I ; // THETA4*exp(ETA(4)) 

// _ALAG(1) = THETAx; //Abs lag

FGF23_0 = BFGF23;
FGFR_0 = 1.0 ;

RPHOS_0 = 1.0 ;
//PHOSPHATE_0 = BPO4; 
//PTH_0 = BPTH ;
//CA_0 = BCA ;
//VITD_0 = BVD ; 




  TGFB_0 = Pic0*1000;
  TGFBact_0 = Pic0;
  OBfast_0 = OBtot0*FracOBfast;
  OBslow_0 = OBtot0*(1-FracOBfast);
  M_0 = k3*RNK_0*L_0/k4;
  N_0 = k1*O_0*L_0/k2;
  AOH_0 = B_0/10;

  IntraPO_0 = IntraPOss ; 
  
  //double week = 0; 
  // if(TIME>131490.0) week = (TIME-131490.0)/(7*24) ;  // 131490.0 = hours per 15 years
  // double SevTabs = 0 ;
  // if(TIME>131490.0)  SevTabs = 2.0 + SevEffMaxTab*pow(week,SevTabsGam)/(pow(week,SevTabsGam) + pow(SevTabsET50,SevTabsGam));
  
  
  
$PARAM 
  SevTabs = 0.0 
  T49 = 51.8
  IntraPOss = 1500 // 3226.0

  SevEff = 1.0

  SevEffMaxTab = 6.5 
  SevTabsGam = 0.8 

  T59ckdGam = 0.15
  PTHckdGam = 1.33

  //perC2 = 0.7

// double week = SOLVERTIME*7*24; 
// double SevTabs = 2.0 + SevEffMaxTab*pow(week,SevTabsGam)/(pow(week,SevTabsGam) + pow(SevTabsET50,SevTabsGam)); 
  
  CLI = 25.2  
  VI = 95.5
  KAI = 4.21
  D1I = 0.64
  
  CMTFLAG=1.0
  
  BFGF23 = 45.0
  BPO4 = 3.4
  BPTH  = 31.0
  BCA  = 8.8
  BVD =  53.3

  // Fit 7
  // ctriolSTIMpoMax FGFRIC50renal    kFGF23 FGFRrenalMAX FGF23IC50vitD
  // 18.15961        15.0738       0.1866604    0.5332671      18.56396

  kFGF23 = 0.1866604         // from fit 7
  kFGF= 0.08         

  FGF23IC50 = 6.0
  ctriolSTIMgam = 0.6
  // ctriolSTIMmax = 3.020653 // 
  // ctriolSTIMec50=1.489593  // 
  phosSTIMfgf23Gam=1.0 
  htrMTT=3.0
  ftrMTT=2.0          
  T71=0.01

  FGFRIC50renal = 15.0738 // from fit 7
  FGFRrenalMAX = 0.5332671    // from fit 7

  FGFRIC50bone = 0.3
  FGFRboneMAX = 0.6

  koutRPHOS = 0.02
  ctriolSTIMpoMax = 18.15961 // from fit 7

  FGF23STIMMAX =  1.0
  FGF23STIMMAXvitD =  4.0
  
  FGF23IC50vitD=18.56396    // from fit 7
  
  CaDay = 88.0



// EXISTING PARAMETERS
$FIXED
  GFRtau=8.0 //## years over which GFR declines
  GFRdelta=0 // 95 //## ml/min 

  molRatPCa = 0.464 // 0.5988 // hydroxyap Ca/P ratio reported to be 1.67

  OBtot0 = 0.00501324
  k1 = 0.00000624
  k2 = 0.112013
  k3 = 0.00000624
  k4 = 0.112013
  V1= 14

  FracJ14 = 0.107763
  J14OCmax = 0.543488
  J14OCgam = 1.6971
  FracJ15 = 0.114376
  kinRNKgam = 0.151825
  koutRNK = 0.00323667
  MOCratioGam = 0.603754
  Da = 0.7/24
  OBtgfGAM = 0.0111319
  koutTGF0 = 0.0000298449
  koutTGFGam = 0.919131
  OCtgfGAM = 0.593891
  EmaxPicROB = 3.9745
  PicROBgam = 1.80968
  FracPicROB = 0.883824
  PicOBgam = 0.122313
  FracPicOB = 0.000244818
  EmaxPicOB = 0.251636
  E0Meff = 0.388267
  EmaxMeffOC = 3.15667
  kinOCgam = 8.53065
  EmaxPicOC = 1.9746
  FracPicOC = 0.878215
  PicOCgam = 1.0168
  E0RANKL = 3.80338
  EmaxL = 0.469779
  T16 = 1.06147
  T64 = 0.05
  T65 = 6.3
  T67 = 1.54865
  AlphOHgam =  0.111241
  k14a =  0.0000244437
  HApMRT = 3.60609
  koutL = 0.00293273
  OsteoEffectGam = 0.173833
  TESTPOWER = 1
  opgPTH50 = 3.85
  IO = 0
  RX2Kout0 = 0.693
  E0rx2Kout = 0.125
  EmaxPTHRX2x = 5
  E0crebKin = 0.5
  EmaxPTHcreb = 3.39745
  crebKout = 0.00279513
  bcl2Kout = 0.693
  ScaEffGam = 0.9
  PhosEff0 = 1.52493
  PhosEff50 = 1.3021
  PhosEffGam = 8.25229
  PO4inhPTHgam = 0
  T69 = 0.10
  Reabs50 = 1.57322
  T7 = 2
  T9 = 90
  T70 = 0.01
  T33 = 0.003
  T34 = 0.037
  T35 = 90
  CaPOgam = 1
  T46 = 1.142
  T52 =  0.365
  OralPhos = 10.5/24 // MW PO4 = 94.9714 g/mol; 10.5 mmol / day ~ 1000 mg / day
  F12 = 0.7
  // T49 = 51.8
  // T55 = 0.019268
  PicOBgamkb = 2.92375
  MultPicOBkb = 3.11842
  FracPic0kb = 0.764028
  E0RUNX2kbEffFACT = 1.01
  RUNkbGAM = 3.81644
  T43 = 1.03856
  T45 = 3.85
  T84 = 0.25
  RUNkbMaxFact = 0.638114
  RUNX20 = 10
  Frackb = 0.313186
  T81 = 0.75
  T87 = 0.0495
  T28 = 0.9
  OralCa = 24.055/24
  T310 = 0.105929
  T77 = 0.909359
  T80 = 4
  CtriolPTgam = 12.5033
  CtriolMax = 4.1029
  CtriolMin = 0.9
  PTout = 0.0001604
  kout = 100/14
  T58 = 6249.09
  T59 = 11.7387
  T61 = 96.25
  IPTHint = 0
  IPTHinf = 0
  Pic0 = 0.228142
  LsurvOCCgam = 3.0923
  EmaxLpth = 1.30721
  kO = 15.8885
  kb = 0.000605516
  LsurvOCgam =3.09023
  FracOBfast = 0.797629
  TERIKA = 10.4
  TERIVC = 94.4
  TERIVD  = 7
  TERICL = 62.2
  TERIF = 1
  PKKA = 0
  PKVC = 10
  PKQ1 = 0
  PKQ2 = 0
  PKVP1 = 1
  PKVP2 = 1
  PKCL = 0
  PKVMAX=0
  PKKM=1
  PPTHKS = 0.3
  PPTHKIN=3
  CAEC50=1.8183
  CAGAM=11.7387
  baseBMDfnBAS = 0.7788
  baseBMDfnETH = 0.08629
  baseBMDfnBMI = 0.01137
  koutBMDls = 0.000397
  koutBMDlsDEN = 0.000145
  koutBMDfnBAS = 0.000005651
  koutBMDfnETH = -0.3781
  koutBMDfnBMI = 0.01141
  //      gamOB = 0.0739
  gamOB = 0.0793
  gamOCls = 0.14
  gamOClsDEN = 0.0679
  gamOCfnBAS = 0.3101
  gamOCfnETH = 0.2319
  gamOCfnBMI = -0.01957
  ETHN=0
  BMI=28
  kdenosl = 1.98411e-06
  E2scalePicB1 = 0.0000116832
  
  
  ESTON = 0
  koutEST=0.05776227
  menoDUR= 14551.56 //as.hour(as.year(1.66))
  ageGAM = -2.3
  age50 = 0.64
  ageENTER = 359406 // as.hour(as.year(41))
  ageDONE = 447066 // as.hour(as.year(51))
  tgfbGAM = 0.0374 //updated 6/20/12 KTB
  tgfbactGAM = 0.045273
  robGAM = 0.16
  obGAM = 0.000012
  maxTmESTkid = 0.923737
  //         ##          gamE2frac = 2.1 ## 2.05 # 2.1
  //         ##          e50E2frac = 0.11 ## 0.135 #0.15
  //         ##          E2FRAC0 = 0.3
  //         ##          E2FRACmax = 0.63 ## 0.65
  //         ##      baseHazard = -5.42
  //         ##     HazBMDcov = -7.95
  //GFRtau=10 //## years over which GFR declines
  //GFRdelta=0 //## ml/min 


$CMTN 
GUT5878
CENTRAL5878
FGF23
FGFR
FGFRbone
RPHOS


$INIT 
  GUT5878 = 0.0
  CENTRAL5878 = 0.0
  FGF23 = 0.0
  FGFR = 1.0
  FGFRbone = 1.0
  RPHOS = 1.0

  HTR =  1.0
  HTR1 = 1.0
  HTR2 = 1.0
  HTR3 = 1.0
  HTR4 = 1.0
  HTR5 = 1.0
  HTR6 = 1.0
  HTR7 = 1.0
  HTR8 = 1.0
  

  PTH = 53.90
  S = 0.5
  PTmax = 1.00
  B = 1260.0
  SC = 0.0
  P = 32.90
  ECCPhos = 16.8
  T = 1.58471
  R = 0.50
  HAp = 1.00
  PhosGut = 0.839
  IntraPO = 3226.0
  OC = 0.00115398
  ROB1 = 0.00104122
  L = 0.4
  RNK= 10.0
  O = 4.0
  Q = 100.0
  Qbone = 24900.0
  RX2 = 10.0
  CREB = 10.0
  BCL2 = 100.0
  TERISC = 0
  TERICENT=0
  PKGUT=0, PKCENT=0,PKPER1 = 0,PKPER2 = 0
  
  UCA=0
  UPO = 0
  //VALUE1=0
  //VALUE2=0
  //VALUE3=0
  TGFB=0
  TGFBact=0
  OBfast=0
  OBslow=0
  M=0
  N=0
  AOH=126
  EST  = 1
  GFR = 100/16.667
  GFRdt = 85.0/16.667
  GFRmin = 15.0/16.667
// VALUE3 = 0

$ODE 
  // ASP5878 specific ODEs
  
  
  dxdt_GUT5878 = -KA5878*GUT5878;
  dxdt_CENTRAL5878 = KA5878*GUT5878 - K5878*CENTRAL5878;


  double JFGF23 = BFGF23*kFGF; 
  double FGF23INH = 1.0 ;
  double FGFRINH = 1.0 ;
  double FGFRINHbone = 1.0 ;
  double FGF23STIM = 1.0 ;
  double FGF23STIMvitD = 1.0 ;
  
  double PhosSTIM = 1.0;
  
  double CP5878 = (1000.0*CENTRAL5878/V5878) ; 

  if(CP5878>0.01) FGF23INH = 1.0-(CP5878/(CP5878+FGF23IC50)); else FGF23INH = 1.0 ;
  if(CP5878>0.01) FGFRINH = 1.0-((FGFRrenalMAX*CP5878)/(CP5878+FGFRIC50renal)); else FGFRINH = 1.0 ;
  if(CP5878>0.01) FGFRINHbone = 1.0-((FGFRboneMAX*CP5878)/(CP5878+FGFRIC50bone)); else FGFRINHbone = 1.0 ;
  if(CP5878>0.01) FGF23STIM = 1.0+(FGF23STIMMAX*CP5878/(CP5878+FGF23IC50)); else FGF23STIM = 1.0 ; // effect on PTH
  if(CP5878>0.01) FGF23STIMvitD = 1.0+(FGF23STIMMAXvitD*CP5878/(CP5878+FGF23IC50vitD)); else FGF23STIMvitD = 1.0 ;
  // if(CP>0.01) FGF23INH = 0.7; else FGF23INH = 1.0 ;
  if(ECCPhos>16.800) PhosSTIM = pow((ECCPhos/ECCPhos0),phosSTIMfgf23Gam) ;

  double ctriolSTIM = pow((C8/Calcitriol0),ctriolSTIMgam) ;
  
  // double VDnull = (C8/Calcitriol0) - 1.0 ;
  // if(C8/Calcitriol0 < 1.0) VDnull = 0.0;
  // double ctriolSTIM = 1.0 + VDnull*ctriolSTIMmax / (VDnull + ctriolSTIMec50) ; 
  
  dxdt_FGF23 = JFGF23*FGFRINHbone*PhosSTIM*ctriolSTIM - kFGF*FGF23*FGFRbone; 

  
  // effect of FGF23 on renal PO4 excretion
    
  dxdt_FGFR = kFGF23*FGFRINH - kFGF23*FGFR; 
  dxdt_FGFRbone = kFGF*FGFRINHbone - kFGF*FGFRbone; 
  
    
  // transit compartment for Vitamin D (c'triol) effect... to represent transcriptional effect as
  // suggested in: Wöhrle, S., Bonny, O., Beluch, N., Gaulis, S., Stamm, C., Scheibler, M., … Graus-Porta, D. (2011). 
  // FGF receptors control vitamin D and phosphate homeostasis by mediating renal FGF-23 signaling and regulating FGF-23 
  // expression in bone. Journal of Bone and Mineral Research, 26(10), 2486–2497. http://doi.org/10.1002/jbmr.478
  
  double htrKtr = (1+8)/htrMTT;
  dxdt_HTR1 =  htrKtr*FGF23STIMvitD - HTR1*htrKtr;
  dxdt_HTR2 =  htrKtr*(HTR1-HTR2);
  dxdt_HTR3 =  htrKtr*(HTR2-HTR3);
  dxdt_HTR4 =  htrKtr*(HTR3-HTR4);
  dxdt_HTR5 =  htrKtr*(HTR4-HTR5);
  dxdt_HTR6 =  htrKtr*(HTR5-HTR6);
  dxdt_HTR7 =  htrKtr*(HTR6-HTR7);
  dxdt_HTR8 =  htrKtr*(HTR7-HTR8);
  dxdt_HTR =   htrKtr*HTR8 - htrKtr*HTR;
  
  
  
  //***************************************************************
  //   Calcium / bone model algebraic relationships and
  //   differential equations
  //****************************************************************/
  
  
  //  /*******************************************/
  //  /*  START EDITING HERE                     */
  //  /*******************************************/
  //  /* CHANGES - KTB 2/22/2012 */
  //  /* Pass in initial calcium conc (P0)... already psssing in PTH0 */
  //  /* 2.35 = 32.9/14 = initial calcium concentration */
  //  /* 3.85 = 53.9 / 14 = initial PTH concentration */
  //  /* Substitute CaConc0 for 2.35 */
  //  /* Substitute PTHconc0 for 3.85 */
  // c("Q", "RNK","L", "ROB1",
  //                 "Qbone","O","RX2","CREB",
  //                 "M", "TGFBact","TGFB","PREPTH",
  //                 "PTH","OBfast","OBslow","OC","P","T",
  //                 "EST","GFR"
  //                 )
  
  
  double PhosEffect = 0;
  double J48 = 0;
  double J27 = 0;
  double RUNX2 = 0;
  double kinEST = 0;
  
  //* parameters derived from SS initial conditions */
  double T13 = (CaDay/24)/Q0;
  
  double T15 = CaDay/(CaConc0*V1*24);
  
  double J14OC50= exp(log((J14OCmax*pow(OC0,J14OCgam)/T13) - pow(OC0,J14OCgam))/J14OCgam);
  
  double OCeqn = (J14OCmax*pow(Osteoclast,J14OCgam))/(pow(Osteoclast,J14OCgam) + pow(J14OC50,J14OCgam));
  
  double kinRNK = (koutRNK*RNK0 + k3*RNK0*L0 - k4*M0) / pow(TGFBact0,kinRNKgam) ;
  
  double MOCratio = M/Osteoclast;
  
  double MOCratio0 = M0/OC0;
  
  double MOCratioEff = pow((MOCratio/MOCratio0), MOCratioGam);
  
  double J14OCdepend = OCeqn*Q0*FracJ14*MOCratioEff;
  
  double J14 = T13*Q0*(1-FracJ14) + J14OCdepend;
  
  
  //* 0.464, reported as the molar ratio of P / Ca in hydroxyapatite. */
  // double J41 = 0.464*J14;
  double J41 = molRatPCa*J14;
  
  double bigDb = kb*OB0*Pic0/ROB10;
  
  double kinTGF = koutTGF0*TGFB0;
  
  double koutTGF = koutTGF0*(pow((TGFB/TGFB0),koutTGFGam));
  
  double koutTGFact = koutTGF0*1000;
  
  double koutTGFeqn = koutTGF*TGFB*(pow((Osteoclast/OC0), OCtgfGAM));
  
  double E0PicROB = FracPicROB*Pic0;
  
  double EC50PicROBparen= (EmaxPicROB*pow(TGFBact0,PicROBgam) / (Pic0 - E0PicROB)) - pow(TGFBact0,PicROBgam);
  
  double EC50PicROB = exp(log(EC50PicROBparen)/PicROBgam);
  
  double Dr = kb*OB0/Pic0;
  
  double PicROB = E0PicROB + EmaxPicROB*pow(TGFBact,PicROBgam)/(pow(TGFBact,PicROBgam) + pow(EC50PicROB,PicROBgam));
  
  double ROBin = Dr*PicROB;
  
  double E0PicOB = FracPicOB*Pic0;
  
  double EC50PicOBparen = (EmaxPicOB*pow(TGFBact0,PicOBgam)/(Pic0 - E0PicOB)) - pow(TGFBact0,PicOBgam);
  
  double EC50PicOB = exp(log(EC50PicOBparen)/PicOBgam);
  
  double PicOB = E0PicOB + EmaxPicOB*pow(TGFBact,PicOBgam) / (pow(TGFBact,PicOBgam) + pow(EC50PicOB,PicOBgam));
  
  double KPT =1*(bigDb/PicOB);
  
  double EC50MeffOC = exp(log(pow(M0, kinOCgam)*EmaxMeffOC/(1-E0Meff) - pow(M0, kinOCgam))/kinOCgam);
  
  double MeffOC = E0Meff + (EmaxMeffOC * pow(M, kinOCgam)/(pow(M, kinOCgam) + pow(EC50MeffOC,kinOCgam)));
  
  double kinOC2 = Da*PicOCkin*MeffOC*OC0;
  
  double E0PicOC = FracPicOC*Pic0;
  
  double EC50PicOCparen = (EmaxPicOC*pow(TGFBact0, PicOCgam)/(Pic0 - E0PicOC)) - pow(TGFBact0, PicOCgam);
  
  double EC50PicOC = exp(log(EC50PicOCparen)/PicOCgam);
  
  double PicOC = E0PicOC + ((EmaxPicOC*pow(TGFBact, PicOCgam))/(pow(TGFBact, PicOCgam) + pow(EC50PicOC, PicOCgam)));
  
  double PiL0 = (k3/k4)*L0;
  
  double PiL = M/10;
  
  double EC50survInPar = (E0RANKL - EmaxL)*(pow(PiL0, LsurvOCgam)/(E0RANKL - 1)) - pow(PiL0, LsurvOCgam);
  
  double EC50surv = exp(log(EC50survInPar)/LsurvOCgam);
  
  double LsurvOC = E0RANKL - (E0RANKL - EmaxL)*(pow(PiL, LsurvOCgam)/(pow(PiL, LsurvOCgam) + pow(EC50surv, LsurvOCgam)));
  
  double KLSoc = Da*PicOC*LsurvOC;
  
  double T66 = (pow(T67, AlphOHgam) + pow(PTHconc0, AlphOHgam) )/pow(PTHconc0, AlphOHgam) ;
  
  double k15a = k14a*Qbone0/Q0 ;
  
  double J14a = k14a*Qbone;
  
  double J15a = k15a*Q ;
  
  //* Hydroxy-apatite */
  double kLShap = 1/HApMRT;
  
  double kHApIn = kLShap/OB0;
  
  //* Calcium flux from plasma into bone */
  double J15 = (T15*P*(1-FracJ15) + T15*P*FracJ15*HAp);
  
  //* 0.464, reported as the molar ratio of P / Ca in hydroxyapatite. */
  // double J42 = 0.464*J15;
  double J42 = molRatPCa*J15;
  
  double Osteoblast = OBfast + OBslow;
  
  double kinLbase = koutL*L0;
  
  double OsteoEffect = pow((Osteoblast/OB0), OsteoEffectGam) ;
  
  double PTH50 = EmaxLpth*PTHconc0 - PTHconc0 ;
  
  double LpthEff = EmaxLpth*(PTHconc) / ((PTH50*pow(OsteoEffect,TESTPOWER)) + (PTHconc)) ;
  
  double kinL = kinLbase*(OsteoEffect)*LpthEff;
  
  double pObase = kO*O0;
  
  double pO = pObase*(D/ROB10)*((PTHconc+(opgPTH50*(D/ROB10)))/(2*PTHconc))+ IO;
  
  double RX2Kin = RX2Kout0*RX20;
  
  double EC50PTHRX2x = ((EmaxPTHRX2x*PTHconc0)/(RX2Kout0 - E0rx2Kout)) - PTHconc0;
  
  double RX2Kout = E0rx2Kout + EmaxPTHRX2x*PTHconc/(PTHconc+EC50PTHRX2x);
  
  
  //*******************************************************/
  //* START CREB-RELATED EQUATIONS                        */
  //*******************************************************/
  double EC50PTHcreb = ((EmaxPTHcreb*PTHconc0)/(1-E0crebKin)) -  PTHconc0;
  
  double crebKin0= crebKout*CREB0;
  
  double crebKin = crebKin0* (E0crebKin + EmaxPTHcreb*PTHconc/(PTHconc+EC50PTHcreb));
  
  double bcl2Kin = RX2*CREB*bcl2Kout;
  //*******************************************************/
  
  //*******************************************************/
  //* START PHOS-RELATED EQUATIONS                        */
  //*******************************************************/
  
  //* QUESTION:
  //   Is 1.2 just the initial extracellular phosphate concentration?
  //   16.8 mmol /14 L = 1.2 mM*/
  
  //* C2 is extracellular phosphate concentration */
  double PO4inhPTH = pow((C2/1.2),PO4inhPTHgam);
  
  double PhosEffTop = (PhosEff0 - 1)*( pow(1.2, PhosEffGam) + pow(PhosEff50, PhosEffGam) );
  
  double PhosEffBot =PhosEff0 * pow(1.2, PhosEffGam);
  
  double PhosEffMax =  PhosEffTop / PhosEffBot;
  
  double PhosEff = PhosEff0 - (PhosEffMax*PhosEff0 * pow(C2, PhosEffGam) /(pow(C2, PhosEffGam)  + pow(PhosEff50, PhosEffGam)));
  
  if (C2 > 1.2) PhosEffect = PhosEff ; else PhosEffect = 1;
  
  double T68 = T66*pow(PTHconc, AlphOHgam)/(pow(T67, AlphOHgam)*PO4inhPTH+pow(PTHconc, AlphOHgam)) ;
  
  // double SE = T65*T68*PhosEffect*HTR;
  
  // The PhosEffect was built in here to represent, in a sense, the feedback from FGF23
  // since FGF23 was not previosuly in the model. Now the HTR effect makes the PhosEffect
  // redundant -- and misleading -- since the FGF23 and PO4 increases are actually reflective 
  // of FGFR inhibition and so AOH (and so calcitriol) production should increase rather
  // decline (as would be predicted solely by FGF23 or PO4 increase) as previously modeled.
  double SE = T65*T68*HTR;
  
  
  //* Equations relating to calcitriol-dependent calcium absorption */
  double T36 = T33 + (T34-T33)*(pow(C8,CaPOgam)/(pow(T35,CaPOgam)+ pow(C8,CaPOgam)));
  double T37 = T34 - (T34-T33)*(pow(C8,CaPOgam)/(pow(T35,CaPOgam)+ pow(C8,CaPOgam)));
  
  //* ======================================
  //   RENAL CALCIUM HANDLING
  //   ======================================
  //   Calcium filtration rate in kidney ;
  //   We assume that 50% of filtered calcium is reabsorbed in a PTH-independent manner;
  //   ... and 50% is reabsorbed in a PTH-dependent manner
  //   Fraction unbound assumed to be 0.6
  //
  
  // double GFR = GFRdt + GFRmin ; // use this GFR if eliciting CKD
  
  double CaFilt = 0.6*0.5*GFR*CaConc;
  
  //* Maximum calcium reabsorption in the kidney - PTH sensitive*/
  double mtmEST = (1-maxTmESTkid)/(1-0.1);  //*(1-maxEST)/(1-minEST)*/
  double tmEST = 1 - mtmEST + mtmEST*EST;
  
  double ReabsMax = tmEST * (0.3*GFR*CaConc0 - 0.149997)*(Reabs50 + CaConc0) / CaConc0;
  
  //* Effect of PTH on calcium reabsorption */
  double T17 = PTHconc0*T16 - PTHconc0;
  
  double ReabsPTHeff = (T16*PTHconc)/(PTHconc + T17);
  
  //* PTH-sensitive calcium reabsorption in kidney */
  //* Reabs50 = 1.573 = H(4-u)-delta */
  double CaReabsActive =  (ReabsMax*C1/(Reabs50 + C1))*ReabsPTHeff;
  
  double T20 = CaFilt - CaReabsActive;
  
  double T10 = T7*C8/(C8+T9);
  
  //* Temporary calcium excretion rate */
  double J27a = (2-T10)*T20; 
  //double J27a = T20 ; // (2-T10)*T20; ; 20Jan2016 M Riggs \url{http://onlinelibrary.wiley.com/store/10.1111/vec.12036/asset/vec12036.pdf?v=1&t=ijncre8x&s=ccf360908765b7ba28ec6f341b4ac4e17ee7d974}
  //                                        Table 1: C'triol effects urine Ca only when Ca levels too high or low
  //                                       ... as parameterized in (2-T10)*T20, C'triol effects Ca regardless of Ca conc = over-estimating C'triol influence on J27
  //                                       ... for now, leave it's effect out
  
  
  //* J27 will be the flux of calcium out of the plasma via the kidney */
  if (J27a<0)  J27 = 0 ; else  J27 = J27a;
  
  double ScaEff = pow( (CaConc0/CaConc), ScaEffGam);
  
  double T72 = 90.0 * ScaEff;
  
  double T73 = T71 * (C8 - T72);
  
  double T74 = (exp(T73) - exp(-T73)) / (exp(T73) + exp(-T73));
  
  double T75 = T70 * (0.85 * (1.0 + T74) + 0.15) ;
  
  double T76 = T70 * (0.85 * (1.0 - T74) + 0.15);
  
  //* phosphate renal excretion */
  double T47 = T46*0.88*GFR;
  
  //double J48a = 0.88*GFR*C2 - T47;
  double J48a = (((0.88*GFR*C2*FGFR) - T47)); // INHrenal) / FGF23INHrenal ;
  //double J48a = (0.88*GFR*C2*FGFR - T47)*pow(((GFR*16.667)/100),0.45); // use this for CKD
  
  if (J48a < 0) J48 = 0 ; else J48 = J48a;
  
  /* phosphate oral absorption */
  // double J53 = T52*PhosGut;
  // Base model did not include influence of Vitamin D on Phosphate absorption
  // evidence that Ctriol can nearly double phosphate absorption... 
  //       https://www.ncbi.nlm.nih.gov/pmc/articles/PMC388157/pdf/pnas00057-0034.pdf
  // T83 (=R/0.5) describes this for calcium, start by using similar expression for PO4... RPHOS
  double J53 = T52*PhosGut*RPHOS*SevEffon ; // pow(T83,phosFctriol);
  
  
  // double J54 = T49*C2;
  
  double T55 = T49*(1.2 / IntraPO0) ; // T49*C2_0 / intraPO_0
  //double J56 = T55*IntraPO;
  double J56 = T55*IntraPO; //*pow((IntraPOss/IntraPO),4);
  double J54 = T49*C2*pow((IntraPOss/IntraPO),4) ; // slows how much PO4 can go into IntraPO
  
  
  //* Parameters describing TGF-beta effects on Osteoblast and clast differentiation and apoptosis */
  double E0PicOBkb = MultPicOBkb*Pic0;
  
  double EmaxPicOBkb = FracPic0kb*Pic0;
  
  double EC50PicOBparenKb = ((E0PicOBkb - EmaxPicOBkb)*pow(TGFBact0,PicOBgamkb)) / (E0PicOBkb - Pic0)  - pow(TGFBact0,PicOBgamkb);
  
  double EC50PicOBkb = exp(log(EC50PicOBparenKb)/PicOBgamkb);
  
  double PicOBkb = E0PicOBkb - (E0PicOBkb  - EmaxPicOBkb)*pow(TGFBact,PicOBgamkb) / (pow(TGFBact,PicOBgamkb) + pow(EC50PicOBkb,PicOBgamkb));
  
  //* MMR 01-Mar-2012 adding estrogen effect that propogates through to OB apoptosis*/
  double PicOBkbEff = (PicOBkb/Pic0)*(1/(pow(EST,E2scalePicB1))) ;
  
  //* Parameters describing osteoblast apoptosis as affected by PTH (continuous vs intermitent) */
  double E0RUNX2kbEff= E0RUNX2kbEffFACT*kb;
  
  if (BCL2 > 105)  RUNX2 = BCL2 - 90 ; else  RUNX2 = 10;
  
  double RUNkbMax = E0RUNX2kbEff*RUNkbMaxFact;
  
  double INparen = (RUNkbMax * pow(RUNX20,RUNkbGAM)) / (E0RUNX2kbEff - kb) - pow(RUNX20,RUNkbGAM);
  
  double RUNkb50 = exp(log(INparen)/RUNkbGAM);
  
  double RUNX2kbPrimeEff = RUNkbMax*pow(RUNX2,RUNkbGAM) / (pow(RUNX2,RUNkbGAM) + pow(RUNkb50,RUNkbGAM));
  
  double kbprime = E0RUNX2kbEff*PicOBkbEff - RUNX2kbPrimeEff;
  
  double kbslow = kbprime*Frackb;
  
  double kbfast = (kb*OB0 + kbslow*OBfast0 - kbslow*OB0) / OBfast0 ;
  
  double Frackb2 = kbfast/kbprime;
  
  //***********************************************************/
  //* Equations relating to calcium movement to/from the gut */
  //*********************************************************/
  double T29 = (T28*T0 - T310*T0)/T310;
  
  double T31 = T28*T/(T+T29);
  
  //* R is calcitriol-dependent gut Ca2+ absorption */
  double T83 = R/0.5;

    
  //* J40 = calcium flux from gut to plasma */
  double J40 = T31*T*T83/(T + T81) + T87*T;
  
  //* T85 relates to extent of absorption of orally-administered dose */
  
  double T85Rpart = pow(R, T80)/(pow(R,T80) + pow(T81,T80));
  double T85 = T77*T85Rpart;
  
  //*****************************/
  //* Calcitriol equations     */
  //***************************/
  
  double INparenCtriol =((CtriolMax - CtriolMin) * pow(Calcitriol0, CtriolPTgam)) / (CtriolMax - 1)- pow(Calcitriol0,CtriolPTgam);

  double Ctriol50 = exp(log(INparenCtriol) / CtriolPTgam) ;
  
  double CtriolPTeff = CtriolMax - (CtriolMax - CtriolMin) * pow(C8, CtriolPTgam) / (pow(C8, CtriolPTgam) + pow(Ctriol50, CtriolPTgam));
  
  double PTin = PTout * CtriolPTeff;
  
  double FCTD = (S / 0.5) * PTmax;
  
  //*  Not used with precursor pool model */
  double INparenCa =(T58 - T61) * pow(CaConc0, T59) / (T58 - 385) - pow(CaConc0, T59);
  double T60 = exp(log(INparenCa) / T59) ;
  double T63 =  T58 - (T58 - T61) * pow((CaConc), T59) / (pow((CaConc), T59) + pow(T60, T59));
  
  //* Comment out precursor pool */
  //* double EMAXCA = pow(CaConc,CAGAM)/(pow(CaConc,CAGAM)+pow(CAEC50,CAGAM)); */
  
  //* Zero-order precursor production rate */
  //*  double PPTHR0 = PREPTH0*PPTHKS + PTH0*kout; */
  
  
  /// BEGIN, use if mimic CKD...
  //   double T59ckd = T59*pow(FCTD,T59ckdGam) ; // chosen to raise T59 (greater sigmoidicity) when FCTD raised thru full CKD
  // double PTHckd = pow(FCTD,PTHckdGam) ; //  chosen to raise PTH at CaConc of 2.35 when FCTD raised thru CKD
  
  // // Line below includes shape change with FCTD change... (as PT gland size changes, so does max 
  // // and shape (gamma, T59))
  // double INparenCa =(FCTD*T58 - T61) * pow(CaConc0, T59ckd) / (FCTD*T58 - 385*PTHckd) - pow(CaConc0, T59ckd);
  // double T60 = exp(log(INparenCa) / T59ckd) ;
  
  // double INparenCaNoCKD =(FCTD*T58 - T61) * pow(CaConc0, T59ckd) / (FCTD*T58 - 385*PTHckd) - pow(CaConc0, T59ckd);
  // double T60noCKD = exp(log(INparenCaNoCKD) / T59) ;
  
  // double rhoNativeCa = pow((CaConc), T59ckd) / (pow((CaConc), T59ckd) + pow(T60, T59ckd)); // *(1-rhoAMG)
  // double rhoCaNoCKD = pow((CaConc), T59) / (pow((CaConc), T59) + pow(T60noCKD, T59)); // *(1-rhoAMG)
  
  // double T63 = FCTD*T58 - (FCTD*T58 - T61)*rhoNativeCa ; // T58 - (T58 - T61)*rhoNativeCa;
  
  //double EPTH = T63 ; // *FCTD;// PPTHKIN*PREPTH*(1-EMAXCA)*FCTD;
  // ... END CKD section
  
  
  
  //* PTH production rate from precursor */
  // double EPTH = T63*FCTD;
   double EPTH = T63*FCTD*FGF23STIM ;
  
  //* Infused and subcutaneously administered PTH */
  double IPTH= 0.693*SC + IPTHinf;
  
  //* Total PTH input rate */
  double SPTH = EPTH + IPTH;
  
  //* Teriparatide pk */
  double TERIPKIN = TERISC*TERICL/TERIVC;
  
  //* Plasma PTH (pmol)
  //   SPTH = PTH input rate
  //   kout = PTH first order elimination rate
  //   TERIPKIN = first order rate from tpar subq dosing into plasma
  //
  dxdt_PTH = SPTH - kout*PTH + TERIPKIN;
  
  //* PTH precursor pool */
  //* dxdt_PREPTH =  PPTHR0 - PPTHKS*PREPTH - EPTH; */
  
  dxdt_S = (1 - S) * T76 - (S* T75);
  
  //* PT gland max capacity */
  dxdt_PTmax = PTin - PTout * PTmax;
  
  //dxdt_B = AOH - T69 * B;
  dxdt_B = AOH - T69 * B;
  
  //* Subcutaneous PTH administration */
  //* NOT USED... use TERISC instead */
  dxdt_SC = IPTHint - 0.693*SC;
  
  dxdt_AOH = SE - T64*AOH ;
  
  //* J14 = nu(12-4) calcium flux from bone into plasma
  //   J15 = nu(4-12) calcium flux from plasma into bone
  //   J27 = nu(4-u)  calcium flux from plasma to urine
  //   J40 = nu(1-4)  calcium flux from gut to plasma
  //
  dxdt_P = J14 - J15- J27 + J40; // - JperiphCa; // use JperiCa in CKD
  
  
  
  //double SevEffon = 1.0 ;
  double SevEffon = 1.0 - 0.08*SevTabs;
  
  //if (GFR < 0.9251283)  SevEffon = SevEff ; else  SevEffon = 1.0;
  
  //* Extracelluar phosphate (mmol) */
  //dxdt_ECCPhos = J41  - J42 - J48 + J53*SevEffon - J54 + J56 ; // - molRatPCa*JperiphCa;
  dxdt_ECCPhos = J41  - J42 - J48 + J53 ; // - J54 + J56;
  
  //* Oral calcium */
  //* CMT: T  UNITS: mmol */
  //* J40 --> flux from gut to plasma */
  //* F11 == T85 by definition */
  dxdt_T = OralCa*T85 - J40;
  
  //* Calcitriol-dependent Ca2+ absorption */
  //* CMT: R */
  dxdt_R = T36*(1- R) - T37*R;

  //* Calcitriol-dependent PO4 absorption */
  //* CMT: RPHOS */
  // from: ctriolSTIMpo = C8*ctriolSTIMpoMax / (C8 + ctriolSTIMpoEC50)  
  // at initial conditions:
  // 1 = Calcitriol0*max / (Calcitriol0 + ec50)
  // Calcitriol0 + ec50 = Calcitriol0*max
  // ec50 = Calcitriol0*(max - 1)
  double ctriolSTIMpoEC50 = Calcitriol0*(ctriolSTIMpoMax - 1) ;
  double ctriolSTIMpo = C8*ctriolSTIMpoMax / (C8 + ctriolSTIMpoEC50) ; 

  dxdt_RPHOS = koutRPHOS*ctriolSTIMpo - koutRPHOS*RPHOS;
  
  //* Hydroxyapatite */
  dxdt_HAp = kHApIn*Osteoblast - kLShap*HAp;
  
  //*******************/
  //* Estrogen piece*/
  //*******************/
  
  double AGE = ageENTER + T0;
  
  
  
  double ageONSET = ageDONE-menoDUR;
  
  if(AGE < ageONSET) kinEST = koutEST * pow((AGE/ageENTER),ageGAM);
  
  if(AGE >= ageONSET) kinEST = koutEST * pow((AGE/ageENTER),ageGAM) * (1 - age50 * (pow((AGE-ageONSET),2)/(pow((menoDUR/2),2) + pow((AGE-ageONSET),2))));
  
  dxdt_EST = (kinEST - koutEST * EST)*ESTON;
  
  //*dxdt_E2FRAC= 0;*/
  
  //* Osteoblasts were considered to exist as two populations: fast and slow.
  //   Fast and slow refer to removal rates (kbfast, kbslow); input assumed to
  //   be the same for each.  Total osteoblasts = OBfast + OBslow */
  //*Estrogen effect added: E2dosePicB1 --> PicOBkbEff --> kbprime --> kbslow, kbfast and Frackb2*/
  dxdt_OBfast = (bigDb/PicOB)*D*FracOBfast*Frackb2  - kbfast*OBfast; /* */
  
  dxdt_OBslow = (bigDb/PicOB)*D*(1-FracOBfast)*Frackb - kbslow*OBslow;
  
  //* d/dt(extracellular phosphate) = J41 -  J42 - J48 + J53 - J54 + J56
  //   d/dt(intracellular phosphate) = J54  -  J56
  
  //      The exchange fluxes of PO4 between ECF and bone (J41 and J42)
  //      set same as respective Ca fluxes but multiplied by the stoichiometric
  //      factor of 0.464, reported as the molar ratio of P / Ca in hydroxyapatite.
  // 
  //      d/dt(dietary phosphate) = OralPhos*F12 -J53
  //      Influx set ~ 1 g (10.5 mmol) of PO4 daily.
  //      Bioavailability assumed = 0.7
  
  // dxdt_PhosGut = OralPhos*F12*pow(T83,phosFctriol) - J53;
  dxdt_PhosGut = OralPhos*F12*RPHOS*SevEffon - J53;
  
  dxdt_IntraPO = J54 - J56;
  
  //* OC: Active Osteoclasts */
  dxdt_OC = kinOC2 - KLSoc*OC;
  
  //* D = ROB1; Responding Osteoblasts */
  dxdt_ROB1 = ROBin * pow(1/EST,robGAM) - KPT*ROB1;
  
  //* Latent TGF-beta pool, production dependent on osteoblast function */
  dxdt_TGFB = kinTGF*(pow((Osteoblast/OB0),OBtgfGAM)) * pow(1/EST,tgfbGAM) - koutTGFeqn * pow(EST,tgfbactGAM);
  
  //* active TGF-beta pool, production dependent on osteoclast function
  //   koutTGFeqn = koutTGF*TGFB*(pow((Osteoclast/OC0), OCtgfGAM))
  //
  dxdt_TGFBact = koutTGFeqn * pow(EST,tgfbactGAM) - koutTGFact*TGFBact;
  

  //****************************************************/
  //*Treatment effects*/
  //****************************************************/
  
  //* L: RANK-L */
  dxdt_L = kinL- koutL*L - k1*O*L + k2*N - k3*RNK*L + k4*M ; 
  
  //* RNK: RANK */
  dxdt_RNK = kinRNK*pow(TGFBact,kinRNKgam) - koutRNK*RNK - k3*RNK*L  + k4*M;
  
  //* M: RANK - RANK-L complex */
  dxdt_M = k3*RNK*L - k4*M;
  
  //* N: OPG - RANK-L complex */
  dxdt_N = k1*O*L - k2*N;
  
  //* O: OPG */
  dxdt_O = pO - k1*O*L + k2*N - kO*O;
  
  //* d/dt(Q) = J15-J14+ J14a-J15a  ;Q=exchangeable bone Ca
  //   d/dt(Qbone) = -J14a + J15a    ;Qbone=non(immediately)-exchangeable bone Ca
  //   ~ 99% of the total Ca stored in bone; approximately 100 mmol of 25000 - 30000 mmol
  //   of total skeletal Ca considered immediately exchangeable with plasma Ca
  //
  dxdt_Q = J15 - J14 + J14a - J15a;
  
  dxdt_Qbone = J15a - J14a;
  
  //*
  //   d/dt(BCL2) = bcl2Kin - bcl2Kout*BCL2
  //   d/dt(RX2) = RX2Kin - RX2Kout*RX2
  //   d/dt(CREB) = crebKin - crebKout*CREB
  
  //      Initial conditions set empirically:
  //      both RX2 and CREB started at 10, BCL2 was the product (10*10=100).
  // 
  //      BCL2 assumed half-life of 1 hour: rate const set to 0.693 (bcl2Kout)
  // 
  //      BCL2 affected osteoblast survival:
  //      decreases elimination rate constant for OB (kbprime)
  
  dxdt_RX2 = RX2Kin - RX2Kout*RX2 ;
  
  dxdt_CREB = crebKin - crebKout*CREB;
  
  dxdt_BCL2 = bcl2Kout*CREB*RX2 - bcl2Kout*BCL2;
  
  //* Teriparatide PK info */
  //* TERIPKIN = TERISC * TERICL/TERIVC */
  dxdt_TERISC = -TERIPKIN;
  dxdt_TERICENT  = TERIPKIN - TERICENT*TERIKA;
  
  // GENERAL PK COMPARTMENT 
  double GFR0 = GFRdt0 + GFRmin0 ; 
  double GFRend = GFR0 - GFRdelta/16.667;
  double GFRtau_ = GFRtau*8766;
  
  double kGFR = -log(GFRend/GFRdt0)/GFRtau_;
  
  dxdt_GFRdt = -kGFR*GFRdt;
  dxdt_GFRmin = 0 ; //  M Riggs added 01-Aug-2015
  
  
  //* NONLINEAR PIECE */
  double PKCLNL = PKVMAX/(PKKM+PKCP);
  
  //* PKGUT */
  dxdt_PKGUT = -PKKA*PKGUT;
  
  //* PKCENT */
  dxdt_PKCENT = PKKA*PKGUT + PKQ1*PKPER1/PKVP1 + PKQ2*PKPER2/PKVP2 - (PKQ1+PKQ2+PKCL+PKCLNL)*PKCENT/PKVC;
  
  //* PKPER1 */
  dxdt_PKPER1 = PKQ1*PKCENT/PKVC - PKQ1*PKPER1/PKVP1;
  
  //* PKPER2 */
  dxdt_PKPER2 = PKQ2*PKCENT/PKVC - PKQ2*PKPER2/PKVP2;
  
  
  //* Collects calcium in urine - cumulative rate */
  //* This is a differential equation */
  dxdt_UCA = J27;
  dxdt_UPO = J48;
  
  //* VALUE1 */
  //  dxdt_VALUE1 = 0;
  
  //* VALUE2 */
  //  dxdt_VALUE2 = 0;
  
  //* VALUE3 */
   // dxdt_VALUE3 = 0;

   $TABLE 
     
   double VitDpcfb = 100*(B/V1)/Calcitriol0 ; // calcitriol % change from baseline
   double PHOSpcfb = 100*ECCPhos / ECCPhos0; // phosphate % change from baseline

   double PTHpm = (PTH/V1);          // PTH conc in pM 
   double PTHpg = (PTH/V1*9.4) ;     // PTH conc in pg/mL
   double PTHpcfb = 100*PTH/PTH0 ;   // PTH % change from baseline 
   double calcitriol = C8 ;          // calcitriol conc ()
   double phosphate = C2 ;           // phosphate conc ()
   double OBchange = 100*OB/OB0;         // Bone specific alkaline phosphate (percent of baseline)
   double OCchange  = 100*OC/OC0;         // serum CTx (percent of baseline)
   double CAchange = (100*P/P0) ;    // Calcium (total, serum), %change from baseline 
   double IPRED = CP5878 ; 
   double convPhos = 3.095975232198142; // # converts PO4 mM to phosphorous mg/dL
   // # for mass conversion: 1 mmol phosphorus = 31 mg phosphorus
   double convCtriol = 34.61538461538461 / 90.0 ; // # http://www.endmemo.com/medical/unitconvert/Vitamin_D.php
   
   if(CMTFLAG==3) IPRED = (convPhos*ECCPhos/V1); // PHOSPHATE (mg/dL)
   if(CMTFLAG==4) IPRED = FGF23; // FGF23
   if(CMTFLAG==5) IPRED = (PTH/V1*9.4); // PTH (pg/mL)
   if(CMTFLAG==6) IPRED = (4*P/V1); // calcium (mg/dL)
   if(CMTFLAG==7) IPRED = (convCtriol*B/V1); // calcitriol (pg/mL)
   
   double IPREDpcfb = CP5878;
   if(CMTFLAG==3) IPREDpcfb = (100*ECCPhos / ECCPhos0 - 100); // PHOSPHATE
   if(CMTFLAG==4) IPREDpcfb = (100*FGF23 / BFGF23 - 100); // FGF23
   if(CMTFLAG==5) IPREDpcfb = (100*PTH / PTH0 - 100); // PTH
   if(CMTFLAG==6) IPREDpcfb = (100*P/P0 - 100); // calcium
   if(CMTFLAG==7) IPREDpcfb = (100*(B/V1)/Calcitriol0 - 100); // calcitriol  
   
//  table(OB) = OBfast + OBslow;
//  table(OCfrac) = OC/OC_0;
//  table(OBfrac) = (OBfast+OBslow)/(OBfast_0 + OBslow_0);


$CAPTURE 
  FGF23STIM FGF23INH CtriolPTeff EPTH CP5878 J41 J42 J48 J53 J54 J56 SE T68 PhosEffect
  EPTH T63 FCTD ScaEff T72 ctriolSTIM T83 J40 FGFRINH T36 T37 ctriolSTIMpo IPRED IPREDpcfb
  J14 J15 J27 J40 J14OCdepend T13 FracJ14 PTHconc CaConc
  
  

