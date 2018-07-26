[PROB]

// URL: https://github.com/metrumresearchgroup/OpenBoneMin
  
  // 13-Oct-2017 MMR
  
[GLOBAL]
#define max(a,b) ((a) > (b) ? (a) : (b))
#define F11 T85
#define PicOCkin Pic0
  
#define SETINIT if(NEWIND <=1) 
  
#define PKCP (PKCENT/PKVC)
#define DENCP (DENCENT/DENVC)
  
#define CaConc0 (P_0/V1)
#define PTHconc0 (PTH_0/V1)
#define OB (OBfast + OBslow)
  
#define PTHconc (PTH/V1)
#define CaConc (P/V1)
#define C1 (P/V1)
#define C2 (ECCPhos/V1)
#define C8 (B/V1)
#define D  ROB1
#define Osteoclast OC
#define OB0 (OBfast_0 + OBslow_0)
#define Calcitriol0 (B_0/V1)
  // denosumab concentration (mol)
#define DENMOL (DENCENT/DENVC/150000)*1000*1000


[MAIN]

TGFB_0 = Pic0*1000.0;
TGFBact_0 = Pic0;
OBfast_0 = OBtot0*FracOBfast;
OBslow_0 = OBtot0*(1-FracOBfast);
M_0 = k3*RNK_0*L_0/k4;
N_0 = k1*O_0*L_0/k2;
AOH_0 = B_0/10.0;
F_DENSC = DENF*1E6;
F_TERISC = TERIF;

  
[PARAM]
  
  OBtot0 = 0.00501324
  k1 = 0.00000624
  k2 = 0.112013
  k3 = 0.00000624
  k4 = 0.112013
  V1= 14.0
  CaDay = 88.0
  FracJ14 = 0.107763
  J14OCmax = 0.543488
  J14OCgam = 1.6971
  FracJ15 = 0.114376
  kinRNKgam = 0.151825
  koutRNK = 0.00323667
  MOCratioGam = 0.603754
  Da = 0.7/24.0
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
  opgPTH50 = 3.85
  RX2Kout0 = 0.693
  E0rx2Kout = 0.125
  EmaxPTHRX2x = 5.0
  E0crebKin = 0.5
  EmaxPTHcreb = 3.39745
  crebKout = 0.00279513
  bcl2Kout = 0.693
  ScaEffGam = 0.9
  PhosEff0 = 1.52493
  PhosEff50 = 1.3021
  PhosEffGam = 8.25229
  PO4inhPTHgam = 0.0
  T69 = 0.10
  Reabs50 = 1.57322
  T7 = 2.0
  T9 = 90.0
  T70 = 0.01
  T71 = 0.03
  T33 = 0.003
  T34 = 0.037
  T35 = 90.0
  CaPOgam = 1.0
  T46 = 1.142
  T52 =  0.365
  OralPhos = 10.5/24
  F12 = 0.7
  T49 = 51.8
  T55 = 0.019268
  PicOBgamkb = 2.92375
  MultPicOBkb = 3.11842
  FracPic0kb = 0.764028
  E0RUNX2kbEffFACT = 1.01
  RUNkbGAM = 3.81644
  T43 = 1.03856
  T45 = 3.85
  T84 = 0.25
  RUNkbMaxFact = 0.638114
  RUNX20 = 10.0
  Frackb = 0.313186
  T81 = 0.75
  T87 = 0.0495
  T28 = 0.9
  OralCa = 24.055/24
  T310 = 0.105929
  T77 = 0.909359
  T80 = 4.0
  CtriolPTgam = 12.5033
  CtriolMax = 4.1029
  CtriolMin = 0.9
  PTout = 0.0001604
  kout = 100/14
  T58 = 6249.09
  T59 = 11.7387
  T61 = 96.25
  IPTHint = 0.0
  IPTHinf = 0.0
  Pic0 = 0.228142
  LsurvOCCgam = 3.0923
  EmaxLpth = 1.30721
  kO = 15.8885
  kb = 0.000605516
  LsurvOCgam =3.09023
  FracOBfast = 0.797629
  TERIKA = 10.4
  TERIVC = 94.4
  TERIVD  = 7.0
  TERICL = 62.2
  TERIF = 1.0
  PKKA = 0.0
  PKVC = 10.0
  PKQ1 = 0.0
  PKQ2 = 0.0
  PKVP1 = 1.0
  PKVP2 = 1.0
  PKCL = 0.0
  PKVMAX=0.0
  PKKM=1.0

  koutBMDls = 0.000397
  koutBMDlsDEN = 0.000145
  koutBMDfnBAS = 0.000005651
  gamOB = 0.0793
  gamOCls = 0.14
  gamOClsDEN = 0.0679
  gamOCfnBAS = 0.3101

  ETHN=0.0
  BMI=28.0
  kdenosl = 2.0e-06
  E2scalePicB1 = 0.0000116832
  

  // Denosumab PK parameters from...
  //  ##  M. Peterson B. Stouch D. Chen S. Baughman D. Holloway P. Bekker and S. Martin.
  //  ##  A population pk/pd model describes the rapid profound and sustained suppression of urinary
  //  ##  n-telopeptide following administration of amg 162 a fully human monoclonal antibody against
  //  ##  rankl to healthy postmenopausal women.
  //  ## The AAPS Journal 24(6 Abstract W4340) 2004.
  // 
  //  ## DENK(1112) = 0.0141
  //  ## DENK(1211) = 0.00798
  DENVMAX = 3110.0
  DENKM = 188.0
  DENVC = 2340.0
  DENVP = 1324.0  // = Q/K(1211)
  DENCL = 2.75
  DENQ = 18.67 // = K(1211) * VC
  DENKA = 0.00592
  DENF = 0.729
  
  // This code block feeds in to control estrogen decline during menopause transition
  ESTON = 0.0
  koutEST=0.05776227
  menoDUR = 8736*1.66//as.hour(as.year(1.66))
  ageGAM = -2.3
  age50 = 0.64
  ageENTER = 8736*41//as.hour(as.year(41))
  ageDONE = 8736*51//as.hour(as.year(51))
  tgfbGAM = 0.0374 //
  tgfbactGAM = 0.045273
  robGAM = 0.16
  obGAM = 0.000012
  maxTmESTkid = 0.923737

  // To use to get GFR to decline with time
  GFRtau=10.0 //## years over which GFR declines
  GFRdelta=0.0 //## ml/min 

$CMTN DENSC, TERISC

$INIT 
  
  DENSC=0.0
  PTH = 53.90
  S = 0.5
  PTmax = 1.00
  B = 1260.0
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
  TERISC = 0.0
  TERICENT=0.0
  PKGUT=0.0, PKCENT=0.0,PKPER1 = 0.0,PKPER2 = 0.0
  
  DENCENT=0.0,DENPER = 0.0
  UCA=0.0
  TGFB=0.0
  TGFBact=0.0
  OBfast=0.0
  OBslow=0.0
  M=0.0
  N=0.0
  AOH=126.0
  EST  = 1.0
  BMDls = 1.0
  BMDlsDEN = 1.0
  BMDfn = 1.0
  GFR = 100.0/16.667

[ODE]
  //***************************************************************
  //   Calcium / bone model algebraic relationships and
  //   differential equations
  //****************************************************************/
  
  double PhosEffect = 0.0;
  double J48 = 0.0;
  double J27 = 0.0;
  double RUNX2 = 0.0;
  double kinEST = 0.0;
  
  //* parameters derived from SS initial conditions */
  double T13 = (CaDay/24)/Q_0;
  
  double T15 = CaDay/(CaConc0*V1*24);
  
  double J14OC50= exp(log((J14OCmax*pow(OC_0,J14OCgam)/T13) - pow(OC_0,J14OCgam))/J14OCgam);
  
  double OCeqn = (J14OCmax*pow(Osteoclast,J14OCgam))/(pow(Osteoclast,J14OCgam) + pow(J14OC50,J14OCgam));
  
  double kinRNK = (koutRNK*RNK_0 + k3*RNK_0*L_0 - k4*M_0) / pow(TGFBact_0,kinRNKgam) ;
  
  double MOCratio = M/Osteoclast;
  
  double MOCratio0 = M_0/OC_0;
  
  double MOCratioEff = pow((MOCratio/MOCratio0), MOCratioGam);
  
  double J14OCdepend = OCeqn*Q_0*FracJ14*MOCratioEff;
  
  double J14 = T13*Q_0*(1-FracJ14) + J14OCdepend;
  
  
  //* 0.464, reported as the molar ratio of P / Ca in hydroxyapatite. */
  double J41 = 0.464*J14;
  
  double bigDb = kb*OB0*Pic0/ROB1_0;
  
  double kinTGF = koutTGF0*TGFB_0;
  
  double koutTGF = koutTGF0*(pow((TGFB/TGFB_0),koutTGFGam));
  
  double koutTGFact = koutTGF0*1000;
  
  double koutTGFeqn = koutTGF*TGFB*(pow((Osteoclast/OC_0), OCtgfGAM));
  
  double E0PicROB = FracPicROB*Pic0;
  
  double EC50PicROBparen= (EmaxPicROB*pow(TGFBact_0,PicROBgam) / (Pic0 - E0PicROB)) - pow(TGFBact_0,PicROBgam);
  
  double EC50PicROB = exp(log(EC50PicROBparen)/PicROBgam);
  
  double Dr = kb*OB0/Pic0;
  
  double PicROB = E0PicROB + EmaxPicROB*pow(TGFBact,PicROBgam)/(pow(TGFBact,PicROBgam) + pow(EC50PicROB,PicROBgam));
  
  double ROBin = Dr*PicROB;
  
  double E0PicOB = FracPicOB*Pic0;
  
  double EC50PicOBparen = (EmaxPicOB*pow(TGFBact_0,PicOBgam)/(Pic0 - E0PicOB)) - pow(TGFBact_0,PicOBgam);
  
  double EC50PicOB = exp(log(EC50PicOBparen)/PicOBgam);
  
  double PicOB = E0PicOB + EmaxPicOB*pow(TGFBact,PicOBgam) / (pow(TGFBact,PicOBgam) + pow(EC50PicOB,PicOBgam));
  
  double KPT =1*(bigDb/PicOB);
  
  double EC50MeffOC = exp(log(pow(M_0, kinOCgam)*EmaxMeffOC/(1-E0Meff) - pow(M_0, kinOCgam))/kinOCgam);
  
  double MeffOC = E0Meff + (EmaxMeffOC * pow(M, kinOCgam)/(pow(M, kinOCgam) + pow(EC50MeffOC,kinOCgam)));
  
  double kinOC2 = Da*PicOCkin*MeffOC*OC_0;
  
  double E0PicOC = FracPicOC*Pic0;
  
  double EC50PicOCparen = (EmaxPicOC*pow(TGFBact_0, PicOCgam)/(Pic0 - E0PicOC)) - pow(TGFBact_0, PicOCgam);
  
  double EC50PicOC = exp(log(EC50PicOCparen)/PicOCgam);
  
  double PicOC = E0PicOC + ((EmaxPicOC*pow(TGFBact, PicOCgam))/(pow(TGFBact, PicOCgam) + pow(EC50PicOC, PicOCgam)));
  
  double PiL0 = (k3/k4)*L_0;
  
  double PiL = M/10;
  
  double EC50survInPar = (E0RANKL - EmaxL)*(pow(PiL0, LsurvOCgam)/(E0RANKL - 1)) - pow(PiL0, LsurvOCgam);
  
  double EC50surv = exp(log(EC50survInPar)/LsurvOCgam);
  
  double LsurvOC = E0RANKL - (E0RANKL - EmaxL)*(pow(PiL, LsurvOCgam)/(pow(PiL, LsurvOCgam) + pow(EC50surv, LsurvOCgam)));
  
  double KLSoc = Da*PicOC*LsurvOC;
  
  double T66 = (pow(T67, AlphOHgam) + pow(PTHconc0, AlphOHgam) )/pow(PTHconc0, AlphOHgam) ;
  
  double k15a = k14a*Qbone_0/Q_0 ;
  
  double J14a = k14a*Qbone;
  
  double J15a = k15a*Q ;
  
  //* Hydroxy-apatite */
  double kLShap = 1/HApMRT;
  
  double kHApIn = kLShap/OB0;
  
  //* Calcium flux from plasma into bone */
  double J15 = (T15*P*(1-FracJ15) + T15*P*FracJ15*HAp);
  
  //* 0.464, reported as the molar ratio of P / Ca in hydroxyapatite. */
  double J42 = 0.464*J15;
  
  double Osteoblast = OBfast + OBslow;
  
  double kinLbase = koutL*L_0;
  
  double OsteoEffect = pow((Osteoblast/OB0), OsteoEffectGam) ;
  
  double PTH50 = EmaxLpth*PTHconc0 - PTHconc0 ;
  
  double LpthEff = EmaxLpth*(PTHconc) / ((PTH50*OsteoEffect) + (PTHconc)) ;
  
  double kinL = kinLbase*(OsteoEffect)*LpthEff;
  
  double pObase = kO*O_0;
  
  double pO = pObase*(D/ROB1_0)*((PTHconc+(opgPTH50*(D/ROB1_0)))/(2*PTHconc)) ;
  
  double RX2Kin = RX2Kout0*RX2_0;
  
  double EC50PTHRX2x = ((EmaxPTHRX2x*PTHconc0)/(RX2Kout0 - E0rx2Kout)) - PTHconc0;
  
  double RX2Kout = E0rx2Kout + EmaxPTHRX2x*PTHconc/(PTHconc+EC50PTHRX2x);
  
  
  //*******************************************************/
  //* START CREB-RELATED EQUATIONS                        */
  //*******************************************************/
  double EC50PTHcreb = ((EmaxPTHcreb*PTHconc0)/(1-E0crebKin)) -  PTHconc0;
  
  double crebKin0= crebKout*CREB_0;
  
  double crebKin = crebKin0* (E0crebKin + EmaxPTHcreb*PTHconc/(PTHconc+EC50PTHcreb));
  
  double bcl2Kin = RX2*CREB*bcl2Kout;
  //*******************************************************/
  
  //*******************************************************/
  //* START PHOS-RELATED EQUATIONS                        */
  //*******************************************************/
  
  //* C2 is extracellular phosphate concentration */
  double PO4inhPTH = pow((C2/1.2),PO4inhPTHgam);
  
  double PhosEffTop = (PhosEff0 - 1)*( pow(1.2, PhosEffGam) + pow(PhosEff50, PhosEffGam) );
  
  double PhosEffBot =PhosEff0 * pow(1.2, PhosEffGam);
  
  double PhosEffMax =  PhosEffTop / PhosEffBot;
  
  double PhosEff = PhosEff0 - (PhosEffMax*PhosEff0 * pow(C2, PhosEffGam) /(pow(C2, PhosEffGam)  + pow(PhosEff50, PhosEffGam)));
  
  if (C2 > 1.2) PhosEffect = PhosEff ; else PhosEffect = 1;
  
  double T68 = T66*pow(PTHconc, AlphOHgam)/(pow(T67, AlphOHgam)*PO4inhPTH+pow(PTHconc, AlphOHgam)) ;
  
  double SE = T65*T68*PhosEffect;
  
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
  double CaFilt = 0.6*0.5*GFR*CaConc;
  
  //* Maximum calcium reabsorption in the kidney - PTH sensitive*/
  double mtmEST = (1-maxTmESTkid)/(1-0.1);  
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
  
  //* J27 will be the flux of calcium out of the plasma via the kidney */
  if (J27a<0)  J27 = 0 ; else  J27 = J27a;
  
  double ScaEff = pow( (CaConc0/CaConc), ScaEffGam);
  
  double T72 = 90 * ScaEff;
  
  double T73 = T71 * (C8 - T72);
  
  double T74 = (exp(T73) - exp(-T73)) / (exp(T73) + exp(-T73));
  
  double T75 = T70 * (0.85 * (1 + T74) + 0.15) ;
  
  double T76 = T70 * (0.85 * (1 - T74) + 0.15);
  
  //* phosphate renal excretion */
  double T47 = T46*0.88*GFR;
  
  double J48a = 0.88*GFR*C2 - T47;
  
  if (J48a < 0) J48 = 0 ; else J48 = J48a;
  
  /* phosphate oral absorption */
  double J53 = T52*PhosGut;
  
  double J54 = T49*C2;
  
  double J56 = T55*IntraPO;
  
  //* Parameters describing TGF-beta effects on Osteoblast and clast differentiation and apoptosis */
  double E0PicOBkb = MultPicOBkb*Pic0;
  
  double EmaxPicOBkb = FracPic0kb*Pic0;
  
  double EC50PicOBparenKb = ((E0PicOBkb - EmaxPicOBkb)*pow(TGFBact_0,PicOBgamkb)) / (E0PicOBkb - Pic0)  - pow(TGFBact_0,PicOBgamkb);
  
  double EC50PicOBkb = exp(log(EC50PicOBparenKb)/PicOBgamkb);
  
  double PicOBkb = E0PicOBkb - (E0PicOBkb  - EmaxPicOBkb)*pow(TGFBact,PicOBgamkb) / (pow(TGFBact,PicOBgamkb) + pow(EC50PicOBkb,PicOBgamkb));
  
  //* Estrogen effect that propogates through to OB apoptosis*/
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
  
  double kbfast = (kb*OB0 + kbslow*OBfast_0 - kbslow*OB0) / OBfast_0 ;
  
  double Frackb2 = kbfast/kbprime;
  
  //***********************************************************/
  //* Equations relating to calcium movement to/from the gut */
  //*********************************************************/
  double T29 = (T28*T_0 - T310*T_0)/T310;
  
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
  
  // S is the PTH gland pool cmt
  double FCTD = (S / 0.5) * PTmax;
  
  double INparenCa =(T58 - T61) * pow(CaConc0, T59) / (T58 - 385) - pow(CaConc0, T59);
  double T60 = exp(log(INparenCa) / T59) ;
  double T63 =  T58 - (T58 - T61) * pow((CaConc), T59) / (pow((CaConc), T59) + pow(T60, T59));
  
  //* Total PTH input (production) rate */
  double SPTH = T63*FCTD;
  
  //* Teriparatide pk */
  double TERIPKIN = TERISC*TERICL/TERIVC;
  
  //* Plasma PTH (pmol)
  //   SPTH = PTH input rate
  //   kout = PTH first order elimination rate
  //   TERIPKIN = first order rate from tpar subq dosing into plasma
  //
  dxdt_PTH = SPTH - kout*PTH + TERIPKIN;
  
  // PTH gland pool
  dxdt_S = (1 - S) * T76 - (S* T75);
  
  //* PT gland max capacity */
  dxdt_PTmax = PTin - PTout * PTmax;
  
  dxdt_B = AOH - T69 * B;
  
  // 1-alpha-hydroxylase (AOH) cmt 
  dxdt_AOH = SE - T64*AOH ;
  
  //* J14 = nu(12-4) calcium flux from bone into plasma
  //   J15 = nu(4-12) calcium flux from plasma into bone
  //   J27 = nu(4-u)  calcium flux from plasma to urine
  //   J40 = nu(1-4)  calcium flux from gut to plasma
  //
  dxdt_P = J14 - J15- J27 + J40;
  
  //* Extracelluar phosphate (mmol) */
  //* d/dt(extracellular phosphate) = J41 -  J42 - J48 + J53 - J54 + J56
  //   d/dt(intracellular phosphate) = J54  -  J56
  
  //      The exchange fluxes of PO4 between ECF and bone (J41 and J42)
  //      set same as respective Ca fluxes but multiplied by the stoichiometric
  //      factor of 0.464, reported as the molar ratio of P / Ca in hydroxyapatite.
  // 
  //      d/dt(dietary phosphate) = OralPhos*F12 -J53
  
  dxdt_PhosGut = OralPhos *F12 - J53;
  dxdt_IntraPO = J54 - J56;  
  dxdt_ECCPhos = J41  - J42 - J48 + J53 - J54 + J56;
  
  //* Oral calcium */
  //* CMT: T  UNITS: mmol */
  //* J40 --> flux from gut to plasma */
  //* F11 == T85 by definition */
  dxdt_T = OralCa*T85 - J40;
  
  //* Calcitriol-dependent Ca2+ absorption */
  //* CMT: R */
  dxdt_R = T36*(1- R) - T37*R;
  
  //* Hydroxyapatite */
  dxdt_HAp = kHApIn*Osteoblast - kLShap*HAp;
  
  //*******************/
  //* Estrogen piece*/
  //*******************/
  
  double AGE = ageENTER + SOLVERTIME;
  double ageONSET = ageDONE-menoDUR;
  
  if(AGE < ageONSET) kinEST = koutEST * pow((AGE/ageENTER),ageGAM);
  
  if(AGE >= ageONSET) kinEST = koutEST * pow((AGE/ageENTER),ageGAM) * (1 - age50 * (pow((AGE-ageONSET),2)/(pow((menoDUR/2),2) + pow((AGE-ageONSET),2))));
  
  dxdt_EST = (kinEST - koutEST * EST)*ESTON;
  
  //* Osteoblasts were considered to exist as two populations: fast and slow.
  //   Fast and slow refer to removal rates (kbfast, kbslow); input assumed to
  //   be the same for each.  Total osteoblasts = OBfast + OBslow */

  // Estrogen effect added: E2dosePicB1 --> PicOBkbEff --> kbprime --> kbslow, kbfast and Frackb2*/
  dxdt_OBfast = (bigDb/PicOB)*D*FracOBfast*Frackb2  - kbfast*OBfast; 
  
  dxdt_OBslow = (bigDb/PicOB)*D*(1-FracOBfast)*Frackb - kbslow*OBslow;
  
  //* OC: Active Osteoclasts */
  dxdt_OC = kinOC2 - KLSoc*OC;
  
  //* D = ROB1; Responding Osteoblasts */
  dxdt_ROB1 = ROBin * pow(1/EST,robGAM) - KPT*ROB1;
  
  //* Latent TGF-beta pool, production dependent on osteoblast function */
  dxdt_TGFB = kinTGF*(pow((Osteoblast/OB0),OBtgfGAM)) * pow(1/EST,tgfbGAM) - koutTGFeqn * pow(EST,tgfbactGAM);
  
  //* active TGF-beta pool, production dependent on osteoclast function
  //   koutTGFeqn = koutTGF*TGFB*(pow((Osteoclast/OC_0), OCtgfGAM))
  //
  dxdt_TGFBact = koutTGFeqn * pow(EST,tgfbactGAM) - koutTGFact*TGFBact;
  
  //**********************************/
  //* BMD BONE MINERAL DENSITY */
  //**********************************/
  //* Solve for kinBMD based on OC=OC_0, OB=OB0, and BMD=BMD0 */
  
  //*Lumbar spine*/
  float kinBMDls =  koutBMDls*BMDls_0;
  
  dxdt_BMDls = kinBMDls * pow(OB/OB0,gamOB) - koutBMDls * pow(OC/OC_0,gamOCls) * BMDls;
  
  //*Lumbar spine with DENOSUMAB*/
  float kinBMDlsDEN =  koutBMDlsDEN*BMDlsDEN_0;
  
  dxdt_BMDlsDEN = kinBMDlsDEN * pow(OB/OB0,gamOB) - koutBMDlsDEN * pow(OC/OC_0,gamOClsDEN) * BMDlsDEN;
  
  //*Femoral neck*/
  double gamOCfn = gamOCfnBAS;     
  double koutBMDfn = koutBMDfnBAS; 
  
  float kinBMDfn =  koutBMDfn*BMDfn_0;
  
  dxdt_BMDfn = kinBMDfn * pow(OB/OB0,gamOB) - koutBMDfn * pow(OC/OC_0,gamOCfn) * BMDfn;
  
  //****************************************************/
  //*Treatment effects*/
  //****************************************************/
  
  //* L: RANK-L */
  dxdt_L = kinL- koutL*L - k1*O*L + k2*N - k3*RNK*L + k4*M -  kdenosl*DENMOL*L;
  
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
  
  //* GENERAL PK COMPARTMENT */
  
  //* NONLINEAR PIECE */
  double PKCLNL = PKVMAX/(PKKM+PKCP);
  
  //* PKGUT */
  dxdt_PKGUT = -PKKA*PKGUT;
  
  //* PKCENT */
  dxdt_PKCENT = PKKA*PKGUT + PKQ1*PKPER1/PKVP1 + PKQ2*PKPER2/PKVP2 - 
    (PKQ1+PKQ2+PKCL+PKCLNL)*PKCENT/PKVC;
  
  //* PKPER1 */
  dxdt_PKPER1 = PKQ1*PKCENT/PKVC - PKQ1*PKPER1/PKVP1;
  
  //* PKPER2 */
  dxdt_PKPER2 = PKQ2*PKCENT/PKVC - PKQ2*PKPER2/PKVP2;
  
  
  //* DENOSUMAB PK */
  double  DENCLNL =  (DENVMAX/(DENKM+DENCP));
  dxdt_DENSC =  -DENKA*DENSC;
  dxdt_DENCENT = DENKA*DENSC + DENQ*DENPER/DENVP - 
    (DENQ+DENCL+DENCLNL)*DENCENT/DENVC ;
  dxdt_DENPER = DENQ*DENCENT/DENVC - DENQ*DENPER/DENVP;
  //* DENOSUMAB PK */

  
  double GFRend = GFR_0 - GFRdelta/16.667;
  double GFRtau_ = GFRtau*8766;
  
  double kGFR = -log(GFRend/GFR_0)/GFRtau_;
  
  dxdt_GFR = -kGFR*GFR;
  
    
  //* Collects calcium in urine - cumulative rate */
  //* This is a differential equation */
  dxdt_UCA = J27;


[TABLE]
capture PTHpm = (PTH/V1);   //* PTH conc in pM */
capture PTHc = (PTH/V1*9.4) ;     //* PTH conc in pg/mL */
capture CAchange = (100*P/P_0) ;   /*/ Calcium (total, serum), %change from baseline */
capture OBchange = 100*OB/OB0;  //* bone-specific alkaline phosphatase, %change from baseline */
capture OCchange = 100*OC/OC_0;  //* serum CTx, %change from baseline */
capture CaC =  P/V1 ;   //* Calcium (total, serum) in mM */
capture BMDlsDENchange = (BMDlsDEN-1)*100;
capture OBtot = OBfast + OBslow;
capture OCfrac = OC/OC_0;
capture OBfrac = (OBfast+OBslow)/(OBfast_0 + OBslow_0);

$CAPTURE DENMOL T43 DENCP
  
