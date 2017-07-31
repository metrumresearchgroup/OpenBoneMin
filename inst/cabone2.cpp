[PROB]

// 07-28-17 RJB
// Converting 'thesis' model to newest mrgsolve format

[GLOBAL]
#define max(a,b) ((a) > (b) ? (a) : (b))
#define F11 T85
#define PicOCkin Pic0

#define SETINIT if(NEWIND <=1)


// Definitions: Initialize Compartments
  
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
#define OBfast0 OBfast_0
#define OC0 OC_0
#define P0 P_0
#define T0 T_0
#define BMDfn0 BMDfn_0
#define BMDls0 BMDls_0
#define BMDlsDEN0 BMDlsDEN_0
#define EST0 EST_0
#define GFR0 GFR_0
#define OBslow0 OBslow_0
#define OCY0 OCY_0
#define TOL0 TOL_0

#define PKCP (PKCENT/PKVC)
#define DENCP (DENCENT/DENVC)

#define CaConc0 (P0/V1)
#define PTHconc0 (PTH0/V1)
#define OB (OBfast*trans + OBslow)

#define PTHconc (PTH/V1)
#define CaConc (P/V1)
#define C1 (P/V1)
#define C2 (ECCPhos/V1)
#define C8 (B/V1)
#define D  ROB1
#define Osteoclast OC
#define OB0 (OBfast0 + OBslow0)
#define Osteoblast (OBfast*trans + OBslow)
#define Calcitriol0 (B_0/V1)

#define BMDlsSCLER0 BMDlsSCLER_0
#define BMDfnSCLER0 BMDfnSCLER_0
#define BMDthSCLER0 BMDthSCLER_0

#define BMDlsTERI0 BMDlsTERI_0
#define BMDfnTERI0 BMDfnTERI_0
#define BMDthTERI0 BMDthTERI_0

#define BMDlsDEN0 BMDlsDEN_0
#define BMDfnDEN0 BMDfnDEN_0
#define BMDthDEN0 BMDthDEN_0


  // denosumab concentration (mol) 
#define DENMOL (DENCENT/DENVC/150000)*1000*1000

 // sclerostin (nmol/L) and sclerostin ab concentration (nmol) 
#define SOSTCP (SOSTCENT/SOSTVC)
 
 [MAIN]
 TGFB_0 = Pic0*1000;
 TGFBact_0 = Pic0;
 OBfast_0 = OBtot0*FracOBfast;
 OBslow_0 = OBtot0*(1-FracOBfast);
 M_0 = k3*RNK_0*L_0/k4;
 N_0 = k1*O_0*L_0/k2;
 AOH_0 = B_0/10;
 F_DENSC = DENF*1E6;
 F_TERISC = TERIF;
 F_SOSTSC = SCLERF;

 

// PARAMETER VALUES 

[PARAM]
  OBtot0 = 0.00501324
  k1 = 0.00000624
  k2 = 0.112013
  k3 = 0.00000624
  k4 = 0.112013
  V1= 14
  CaDay = 88
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
  FracPicROB =  0.883824  
  PicOBgam = 0.122313
  FracPicOB = 0.000244818
  EmaxPicOB = 0.251636
  E0Meff = 0.388267
  EmaxMeffOC = 3.15667
  kinOCgam = 8.53065, EmaxPicOC = 1.9746, FracPicOC = 0.878215, PicOCgam = 1.0168
  E0RANKL = 3.80338, EmaxL = 0.469779
  T16 = 1.06147
  T64 = 0.05
  T65 = 6.3
  T67 = 1.54865
  AlphOHgam =  0.111241
  k14a =  0.0000244437
  HApMRT = 3.60609
  koutL = 0.00293273
  
  TotOsteoEffectGam = 0.173833
  

  TESTPOWER = 1
  opgPTH50 = 3.85
  IO = 0
  RX2Kout0 = 0.693, E0rx2Kout = 0.125
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
  T71 = 0.03
  T33 = 0.003
  T34 = 0.037
  T35 = 90
  CaPOgam = 1
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
  TERIKA = 10.4, TERIVC = 94.4, TERIVD  = 7, TERICL = 62.2, TERIF = 0.95 //Teriparatide PK
  PKKA = 0, PKVC = 10, PKQ1 = 0, PKQ2 = 0    //Generic PK params
  PKVP1 = 1, PKVP2 = 1, PKCL = 0, PKVMAX=0, PKKM=1 //Generic PK params cont
  PPTHKS = 0.3
  PPTHKIN=3
  CAEC50=1.8183
  CAGAM=11.7387
  baseBMDfnBAS = 0.7788
  baseBMDfnETH = 0.08629
  baseBMDfnBMI = 0.01137
  koutBMDls = 0.000397
  koutBMDlsDEN = 7.374996e-05 //0.000145 
  koutBMDfnDEN = 0.0001186424 //0.000145 
  koutBMDthDEN = 0.0001080459 //0.000145 
  koutBMDfnBAS = 0.000005651
  koutBMDfnETH = -0.3781
  koutBMDfnBMI = 0.01141
  gamOB = 0.0793
  gamOCls = 0.14
  gamOClsDEN = 7.912201e-02 //0.0679
  gamOCfnDEN = 0.0514712439 //0.0679
  gamOCthDEN = 0.0552337198 //0.0679
  gamOCfnBAS = 0.3101
  gamOCfnETH = 0.2319
  gamOCfnBMI = -0.01957
  ETHN=0
  BMI=28
  kdenosl = 1.98411e-06 //2.612931e-07
  E2scalePicB1 = 0.0000116832
  FracOBE = 20
  
  // Denosumab Params from Peterson, et al.The AAPS Journal, 24(6 Abstract W4340), 2004 //
  DENVMAX = 3110, DENKM = 188, DENVC = 2340
  DENVP = 1324  // = Q/K(12,11) 
  DENCL = 2.75
  DENQ = 18.67 // = K(12,11) * VC 
  DENKA = 0.00592, DENF = 0.729
  
////DENVMAX = 3660,  // Estimates pulled from Berkeley Model  
////DENKM = 164,  
//#DENVC = 2380,  
//#DENVP = 1870,  
//#DENQ = 64,  
//#DENCL = 2.75,	
//#DENKA = 0.00551,	
//#DENF = 0.75	
  
  
  
  ESTON = 0
  koutEST=0.05776227
  menoDUR=8736*1.66//as.hour(as.year(1.66))
  ageGAM = -2.3
  age50 = 0.64
  ageENTER = 8736*41//as.hour(as.year(41))
  ageDONE = 8736*51//as.hour(as.year(51))
  tgfbGAM = 0.0374 // updated 6/20/12 KTB
  tgfbactGAM = 0.045273
  robGAM = 0.16
  obGAM = 0.000012
  maxTmESTkid = 0.923737
  GFRtau=10 // years over which GFR declines
  GFRdelta=0

  // Sclerostin Piece //
   SOSTVMAX = 5.87/24, SOSTKM = 0.453, SOSTVC = 2.9, SOSTVP = 3.29 
  SOSTCL = 0.254/24, SOSTQ = 0.467/24, SOSTKA = 0.187/24, SCLERF=0.904
  
  // Sclerostin indirect response model // 
  SOSTKIN=3.725/24, SOSTKOUT= 25/24, SOSTK0 = 0.197/24
  
  // P1NP indirect response model : last changed 05/28/14 // 
  
  FracOCY = 0.5 // 20% of OB become OCY from Holmen, 2005
  EMAXSCLER = 4.670795836 // #updated 04/09
  gammaDr = 0.044584465 //0.0703
  koutPRE = 4
  
  // gammascler = 0.5 #
  gammaOCY = 0.276280938  // updated 04/09
  gammaOPG = 1.597073748 // updated 04/09
  SCLEROBgam =  0.162250232 // updated 04/09
  kout_T = 0.006073441 //
  
  koutBMDSCLER = 0.000145
  gamOBSCLERls =  0.75766752 // 4/22 change
  gamOCSCLER = 0.06530286 // change 04/22
  
  gamOBSCLERth = 0.22459622 //4/22 change
  gamOBSCLERfn = 0.13145334 // 4/22 change
  
  gamOBSCLER = 0.099
  gamOCsSCLER = 0.109
  
  kout_BMDdel = 0.00245716 
  kout_BMDdelTERI = 0.001
  
  TYPE = 2

  koutBMDlsTERI =  0.00055370 
  koutBMDthTERI =  0.0001394248 
  koutBMDfnTERI =  0.000066284   
  
  gamOCfnTERI =   0.21199 
  gamOClsTERI = 0.016916 
  gamOCthTERI = 0.13118   
  
  koutBMDlsCOMBO = 0.0001368262 
  koutBMDthCOMBO = 0.0001048986 
  koutBMDfnCOMBO = 0.0001275736 
  
  gamOClsDEN_TERI = 0.1015274694 
  gamOCthDEN_TERI = 0.0704668698 //# updated 6/8/15
  gamOCfnDEN_TERI = 0.0671132952 //updated 6/8/15
  
  
  gamOBTERIfn =  0.495529  // updated as of 03/30
  
  gamOBTERIls =  0.271226 //  updated as of 04/2
  
  gamOBTERIth =  0.29803 //updated 03/30
  
  thTERIfactor = 0.35775

  SMAX = 8.690800116 //updated 04/09
  
  kout_TOL = 0.001901105 //
  kant = -2
  
  baseHazard = 0.058265425   // from fracture510.summary.csv
  HazBMDCov = 1.3259
  //  HazBMDcfbCov =  4.83   // ffrom fracture313.summary.csv
  //  HazBMDbaseCov =  0.39   // from fracture313.summary.csv
  HazpostMenoAgeCov = 0.02690102025  //from fracture510.summary.csv
  HazradFracCov = -0.20419  // from fracture510.summary.csv
  HazBMICov = -0.021603069564   // from fracture510.summary.csv
  agepostmeno= 20  
  radFracInc= 0
  
  BmdbaseRef = 0.783
  postMenoAgeRef = 20
  ageLastMenPeriodRef = 51.7
  bmiRef = 27.1              
  lsBMDbase = 0.8
  Scler_BMD_type=0,Deno_BMD_type=0,Teri_BMD_type=0,DEN_TERI_COMBO=0,Combo_BMD_type=0, SCLER_DEN_SEQ=0
  betaDrug_DENO=-1.73,betaDrug_TERI=-0.898,betaDrug_SCLER=-0.9,betaDrug_COMBO=-1.213 //combo ave deno&bisphos


$CMTN 
SOSTSC, DENSC, TERISC

// INITIAL CONDITIONS 

$INIT

    PTH = 53.90     // (pmol)
    S = 0.5         // PTH gland pool    
    PTmax = 1.00    // PT gland max capacity
    B = 1260.0      // Circulating calicitriol (pmol)
    SC = 0.0        // Subcu PTH compartment (pmol)
    P = 32.90       // Extracellular calcium (mmol)
    ECCPhos = 16.8  // Extracellular phosphate (mmol)
    T = 1.58471     // Oral calcium (mmol)
    R = 0.50        // Calcitriol dependent ca absorption
    HAp = 1.00      // Hydroxyapetite conc
    PhosGut = 0.839  // Oral phosphate (mmol)
    IntraPO = 3226.0 // Intracellular phosphate (mmol)
    OC = 0.00115398  // Osteoclast population 
    ROB1 = 0.00104122 // Responding osteoblasts
    L = 0.4         // RANKL concentration
    RNK= 10.0       // RANK concentration  
    O = 4.0         // OPG
    Q = 100.0       // Immediate-exchangeable bone calcium (mmol)
    Qbone = 24900.0 // Non-immediate exchangeable bone calcium (mmol)
    RX2 = 10.0      // RunX2
    CREB = 10.0     // Creb
    BCL2 = 100.0    // Bcl-2
    TERISC = 0, TERICENT=0   // Teriparatide PK
    PKGUT=0, PKCENT=0, PKPER1 = 0, PKPER2 = 0  // Generic PK
    DENCENT=0, DENPER = 0,  DENSC=0 //Denosumab PK 
    UCA=0           // Urine Calcium (pmol)
    TGFB=0          // Latent TGF beta
    TGFBact=0       // Active TGF beta
    OBfast=0        // Fast differentiating osteoblasts
    OBslow=0        // Slow differentiating Osteoblasts
    M=0             // RANK-RANKL complex
    N=0             // OPG-RANKL complex
    AOH=126         // 1-alpha hydroxylase (pmol)        
    EST  = 1        // estrogen
    BMDls = 1       // lumbar spine BMD for combination therapy
  //  BMDlsDEN = 1    // lumbar spine BMD after Denos Rx
    BMDfn = 1       // femoral neck BMD
    GFR = 100/16.667  // GFR
    SOSTSC = 0, SOSTCENT=0, SOSTPER=0 //Sclerostin Ab PK
    SCLER = 0.149         // Sclerostin compartment  (nmol/L)
    //P1NP = 100  #      # P1NP (% of baseline)#
    OCY = 0.0709 //0.0117          // Osteocytes
    trans = 1 // translation compartment
    
    BMDlsSCLER = 1       // Sclerostin lumbar spine BMD
    BMDfnSCLER = 1       // Sclerostin femoral neck spine BMD
    BMDthSCLER = 1       // Sclerostin total hip BMD
    
    BMDlsDEN = 1     // Denosumab lumbar spine BMD
    BMDfnDEN = 1    // Denosumab femoral neck spine BMD
    BMDthDEN = 1    // Denosumab total hip BMD
    
    delBMDls = 1 // delay compartment for lumbar spine BMD
    delBMDth = 1 // delay compartment for total hip BMD
    delBMDfn = 1 // delay compartment for femoral neck BMD
    
    BMDlsTERI = 1       // TERI lumbar spine BMD //
    BMDfnTERI = 1       // TERI femoral neck spine BMD //
    BMDthTERI = 1       // TERI lumbar spine BMD
    delBMDlsTERI = 1 // delay compartment for lumbar spine BMD with TERI
    delBMDfnTERI = 1 // delay compartment for femoral neck BMD with TERI
    delBMDthTERI = 1 // delay compartment for total hip BMD with TERI
    
    BMDlsCOMBO = 1 // DEN/TERI combo lumbar spine BMD //
    BMDthCOMBO = 1 // DEN/TERI combo total hip BMD //
    BMDfnCOMBO = 1 // DEN/TERI combo femoral neck BMD //
    
    TOL = 1
    cumHazard = 0 

 
//%


// TOL_0 = kin_TOL/kout_TOL; 
 
//END_main



///
//   Calcium / bone model algebraic relationships and
//   differential equations

[ODE]
  //
  //  START EDITING HERE                     
  //
  // CHANGES - KTB 2/22/2012 
  // Pass in initial calcium conc (P0)... already psssing in PTH0 
  // 2.35 = 32.9/14 = initial calcium concentration 
  // 3.85 = 53.9 / 14 = initial PTH concentration 
  // Substitute CaConc0 for 2.35 
  // Substitute PTHconc0 for 3.85 
// c("Q", "RNK","L", "ROB1",
//                 "Qbone","O","RX2","CREB",
//                 "M", "TGFBact","TGFB","PREPTH",
//                 "PTH","OBfast","OBslow","OC","P","T",
//                 "BMDfn","BMDls","BMDlsDEN","EST","GFR"
//                 )

  double PhosEffect = 0;
  double J48 = 0;
  double J27 = 0;
  double RUNX2 = 0;
  double kinEST = 0;

  // parameters derived from SS initial conditions 
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
  
  
  // 0.464, reported as the molar ratio of P / Ca in hydroxyapatite. 
  double J41 = 0.464*J14;
  
  double bigDb = (kb*OB0*Pic0/ROB10);
  
  double kinTGF = koutTGF0*TGFB0;
  
  double koutTGF = koutTGF0*(pow((TGFB/TGFB0),koutTGFGam));
  
  double koutTGFact = koutTGF0*1000;
  
  double koutTGFeqn = koutTGF*TGFB*(pow((Osteoclast/OC0), OCtgfGAM));
  
  double E0PicROB = (FracPicROB*Pic0);
  
  double EC50PicROBparen = (EmaxPicROB*pow(TGFBact0,PicROBgam) / (Pic0 - E0PicROB)) - pow(TGFBact0,PicROBgam);
  
  double EC50PicROB = (exp(log(EC50PicROBparen)/PicROBgam));
  
  double Dr = (kb*OB0/Pic0) ;
  
  double PicROB = (E0PicROB + EmaxPicROB*pow(TGFBact,PicROBgam)/(pow(TGFBact,PicROBgam) + pow(EC50PicROB,PicROBgam)));
  
  double ROBin = (Dr*PicROB);
  
  double SC50 = (SMAX - 1);
  
  double E0PicOB = FracPicOB*Pic0;
  
  double EC50PicOBparen = (EmaxPicOB*pow(TGFBact0,PicOBgam)/(Pic0 - E0PicOB)) - pow(TGFBact0,PicOBgam);
  
  double EC50PicOB = exp(log(EC50PicOBparen)/PicOBgam);
  
  double PicOB = E0PicOB + EmaxPicOB*pow(TGFBact,PicOBgam) / (pow(TGFBact,PicOBgam) + pow(EC50PicOB,PicOBgam));
  
  double KPT = (bigDb/PicOB);
  
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
  
  // Hydroxy-apatite 
  double kLShap = 1/HApMRT;
  
  double kHApIn = kLShap/OB0;
  
  // Calcium flux from plasma into bone 
  double J15 = (T15*P*(1-FracJ15) + T15*P*FracJ15*HAp);
  
  // 0.464, reported as the molar ratio of P / Ca in hydroxyapatite. 
  double J42 = 0.464*J15;
  
  double kinLbase = koutL*L0;
  
  double fracOBEffect = FracOBE/TOL ;         //   /TOL 20 #5 #4.5 #2 #5.8 #14 #25 
  
  double OsteoEffect = pow((Osteoblast/OB0), (TotOsteoEffectGam/fracOBEffect)) ; 
  
  double OsteoCYEffect = pow((OCY/OCY0), (TotOsteoEffectGam*(1-(1/fracOBEffect)))) ; 
  
  double PTH50 = EmaxLpth*PTHconc0 - PTHconc0 ;
  
  double LpthEff = EmaxLpth*(PTHconc) / ((PTH50*pow((OsteoCYEffect*OsteoEffect),TESTPOWER)) + (PTHconc)) ;
  
  double kinL = kinLbase*(OsteoCYEffect)*(OsteoEffect)*LpthEff;
  
  double pObase = kO*O0;
  
  double pO = pObase*(D/ROB10)*((PTHconc+(opgPTH50*(D/ROB10)))/(2*PTHconc))+ IO;
  
  double RX2Kin = RX2Kout0*RX20;
  
  double EC50PTHRX2x = ((EmaxPTHRX2x*PTHconc0)/(RX2Kout0 - E0rx2Kout)) - PTHconc0;
  
  double RX2Kout = E0rx2Kout + EmaxPTHRX2x*PTHconc/(PTHconc+EC50PTHRX2x);
  
  
  //
  // START CREB-RELATED EQUATIONS                        
  //
  double EC50PTHcreb = ((EmaxPTHcreb*PTHconc0)/(1-E0crebKin)) -  PTHconc0;
  
  double crebKin0= crebKout*CREB0;
  
  double crebKin = crebKin0* (E0crebKin + EmaxPTHcreb*PTHconc/(PTHconc+EC50PTHcreb));
  
  double bcl2Kin = RX2*CREB*bcl2Kout;
  //
  
  //
  // START PHOS-RELATED EQUATIONS                        
  //
  
  // QUESTION:
  //     Is 1.2 just the initial extracellular phosphate concentration?
  //     16.8 mmol /14 L = 1.2 mM
  
  // C2 is extracellular phosphate concentration 
  double PO4inhPTH = pow((C2/1.2),PO4inhPTHgam);
  
  double PhosEffTop = (PhosEff0 - 1)*( pow(1.2, PhosEffGam) + pow(PhosEff50, PhosEffGam) );
  
  double PhosEffBot =PhosEff0 * pow(1.2, PhosEffGam);
  
  double PhosEffMax =  PhosEffTop / PhosEffBot;
  
  double PhosEff = PhosEff0 - (PhosEffMax*PhosEff0 * pow(C2, PhosEffGam) /(pow(C2, PhosEffGam)  + pow(PhosEff50, PhosEffGam)));
  
  double PhosEffect;
  if (C2 > 1.2) PhosEffect = PhosEff ; else PhosEffect = 1;
  
  double T68 = T66*pow(PTHconc, AlphOHgam)/(pow(T67, AlphOHgam)*PO4inhPTH+pow(PTHconc, AlphOHgam)) ;
  
  double SE = T65*T68*PhosEffect;
  
  // Equations relating to calcitriol-dependent calcium absorption 
  double T36 = T33 + (T34-T33)*(pow(C8,CaPOgam)/(pow(T35,CaPOgam)+ pow(C8,CaPOgam)));
  double T37 = T34 - (T34-T33)*(pow(C8,CaPOgam)/(pow(T35,CaPOgam)+ pow(C8,CaPOgam)));
  
  // ======================================
  //     RENAL CALCIUM HANDLING
  //     ======================================
  //     Calcium filtration rate in kidney ;
  //     We assume that 50% of filtered calcium is reabsorbed in a PTH-independent manner;
  //     ... and 50% is reabsorbed in a PTH-dependent manner
  //     Fraction unbound assumed to be 0.6
  
  double CaFilt = 0.6*0.5*GFR*CaConc;
  
  // Maximum calcium reabsorption in the kidney - PTH sensitive
  double mtmEST = (1-maxTmESTkid)/(1-0.1);  //(1-maxEST)/(1-minEST)
  double tmEST = 1 - mtmEST + mtmEST*EST;
  
  double ReabsMax = tmEST * (0.3*GFR*CaConc0 - 0.149997)*(Reabs50 + CaConc0) / CaConc0;
  
  // Effect of PTH on calcium reabsorption 
  double T17 = PTHconc0*T16 - PTHconc0;
  
  double ReabsPTHeff = (T16*PTHconc)/(PTHconc + T17);
  
  // PTH-sensitive calcium reabsorption in kidney 
  // Reabs50 = 1.573 = H(4-u)-delta 
  double CaReabsActive =  (ReabsMax*C1/(Reabs50 + C1))*ReabsPTHeff;
  
  double T20 = CaFilt - CaReabsActive;
  
  double T10 = T7*C8/(C8+T9);
  
  // Temporary calcium excretion rate 
  double J27a = (2-T10)*T20;
  
  // J27 will be the flux of calcium out of the plasma via the kidney 
  if (J27a<0)  J27 = 0 ; else  J27 = J27a;
  
  double ScaEff = pow( (CaConc0/CaConc), ScaEffGam);
  
  double T72 = 90 * ScaEff;
  
  double T73 = T71 * (C8 - T72);
  
  double T74 = (exp(T73) - exp(-T73)) / (exp(T73) + exp(-T73));
  
  double T75 = T70 * (0.85 * (1 + T74) + 0.15) ;
  
  double T76 = T70 * (0.85 * (1 - T74) + 0.15);
  
  // phosphate renal excretion 
  double T47 = T46*0.88*GFR;
  
  double J48a = 0.88*GFR*C2 - T47;
  
  if (J48a < 0) J48 = 0 ; else J48 = J48a;
  
  // phosphate oral absorption 
  double J53 = T52*PhosGut;
  
  double J54 = T49*C2;
  
  double J56 = T55*IntraPO;
  
  // Parameters describing TGF-beta effects on Osteoblast and clast differentiation and apoptosis 
  double E0PicOBkb = MultPicOBkb*Pic0; 
  
  double EmaxPicOBkb = FracPic0kb*Pic0;
  
  double EC50PicOBparenKb = ((E0PicOBkb - EmaxPicOBkb)*pow(TGFBact0,PicOBgamkb)) / (E0PicOBkb - Pic0)  - pow(TGFBact0,PicOBgamkb);
  
  double EC50PicOBkb = exp(log(EC50PicOBparenKb)/PicOBgamkb);
  
  double PicOBkb = E0PicOBkb - (E0PicOBkb  - EmaxPicOBkb)*pow(TGFBact,PicOBgamkb) / (pow(TGFBact,PicOBgamkb) + pow(EC50PicOBkb,PicOBgamkb));
  
  // MMR 01-Mar-2012 adding estrogen effect that propogates through to OB apoptosis 
  
  double E0RUNX2kbEff =(E0RUNX2kbEffFACT*kb);
  
  double PicOBkbEff = (PicOBkb/Pic0)*(1/(pow(EST,E2scalePicB1))) ;
  
  double kbprime =  (E0RUNX2kbEff*PicOBkbEff - RUNX2kbPrimeEff);
  double kbslow = (kbprime*Frackb);
  double kbfast = ((kb*OB0 + kbslow*OBfast0 - kbslow*OB0) / OBfast0 );
  
  
  
  // Parameters describing osteoblast apoptosis as affected by PTH (continuous vs intermitent) 
  // 4/14/15 change this to only function when TERI or exogenous PTH is dosed  
  
  if (BCL2 > 105) RUNX2 = BCL2 - 90 ; else RUNX2 = 10;
  
  double RUNkbMax = E0RUNX2kbEff*RUNkbMaxFact;
  
  double INparen = (RUNkbMax * pow(RUNX20,RUNkbGAM)) / (E0RUNX2kbEff - kb) - pow(RUNX20,RUNkbGAM);
  
  double RUNkb50 = exp(log(INparen)/RUNkbGAM);
  
  double RUNX2kbPrimeEff = RUNkbMax*pow(RUNX2,RUNkbGAM) / (pow(RUNX2,RUNkbGAM) + pow(RUNkb50,RUNkbGAM));
  
  
  // Sclerostin Effect - prop to OCY apoptosis rate 
  double SCLEREFF = ((SCLER/SCLER_0));
  double PTHEFFECT = (pow(PTH/PTH_0,20)); 
  
  double SCLEREFF_TOL= (SCLER/SCLER_0);
  double SCLER_TOL = (SMAX*SCLEREFF_TOL/(SC50 + SCLEREFF_TOL));
  
  double kin_TOL = 1 * kout_TOL;
  
  dxdt_TOL = kin_TOL*(SCLER_TOL) - kout_TOL*TOL;
  
  double koutOCY = (OB0*FracOCY)*(pow((SCLEREFF),gammaOCY));
  
  double Frackb2 = kbfast/kbprime;
  
  //
  // Equations relating to calcium movement to/from the gut 
  //
  double T29 = (T28*T0 - T310*T0)/T310;
  
  double T31 = T28*T/(T+T29);
  
  // R is calcitriol-dependent gut Ca2+ absorption 
  double T83 = R/0.5;
  
  // J40 = calcium flux from gut to plasma 
  double J40 = T31*T*T83/(T + T81) + T87*T;
  
  // T85 relates to extent of absorption of orally-administered dose 
  
  double T85Rpart = pow(R, T80)/(pow(R,T80) + pow(T81,T80));
  double T85 = T77*T85Rpart;
  
  //
  // Calcitriol equations     
  //
  double INparenCtriol =((CtriolMax - CtriolMin) * pow(Calcitriol0, CtriolPTgam)) / (CtriolMax - 1)- pow(Calcitriol0,CtriolPTgam);
  
  double Ctriol50 = exp(log(INparenCtriol) / CtriolPTgam) ;
  
  double CtriolPTeff = CtriolMax - (CtriolMax - CtriolMin) * pow(C8, CtriolPTgam) / (pow(C8, CtriolPTgam) + pow(Ctriol50, CtriolPTgam));
  
  double PTin = PTout * CtriolPTeff;
  
  double FCTD = (S / 0.5) * PTmax;
  
  //  Not used with precursor pool model 
  double INparenCa =(T58 - T61) * pow(CaConc0, T59) / (T58 - 385) - pow(CaConc0, T59);
  double T60 = exp(log(INparenCa) / T59) ;
  double T63 =  T58 - (T58 - T61) * pow((CaConc), T59) / (pow((CaConc), T59) + pow(T60, T59));
  
  // Comment out precursor pool 
  // double EMAXCA = pow(CaConc,CAGAM)/(pow(CaConc,CAGAM)+pow(CAEC50,CAGAM)); 
  
  // Zero-order precursor production rate 
  //  double PPTHR0 = PREPTH0*PPTHKS + PTH0*kout; 
  
  // PTH production rate from precursor 
  double EPTH = T63*FCTD;// PPTHKIN*PREPTH*(1-EMAXCA)*FCTD; 
  
  // Infused and subcutaneously administered PTH 
  double IPTH= 0.693*SC + IPTHinf;
  
  // Total PTH input rate 
  double SPTH = EPTH + IPTH;
  
  // Teriparatide pk 
  double TERIPKIN = TERISC*TERICL/TERIVC;
  
  // Plasma PTH (pmol)
  //     SPTH = PTH input rate
  //     kout = PTH first order elimination rate
  //     TERIPKIN = first order rate from tpar subq dosing into plasma
  
  dxdt_PTH = SPTH - kout*PTH + TERIPKIN;
  
  //  if(DEN_TERI_COMBO==1) dxdt_RX2 = 0; 
  
  
  // PTH precursor pool 
  // dxdt_PREPTH =  PPTHR0 - PPTHKS*PREPTH - EPTH; 
  
  dxdt_S = (1 - S) * T76 - (S* T75);
  
  // PT gland max capacity 
  dxdt_PTmax = PTin - PTout * PTmax;
  
  dxdt_B = AOH - T69 * B;
  
  // Subcutaneous PTH administration 
  // NOT USED... use TERISC instead 
  dxdt_SC = IPTHint - 0.693*SC;
  
  dxdt_AOH = SE - T64*AOH ;
  
  // J14 = nu(12-4) calcium flux from bone into plasma
  //     J15 = nu(4-12) calcium flux from plasma into bone
  //     J27 = nu(4-u)  calcium flux from plasma to urine
  //     J40 = nu(1-4)  calcium flux from gut to plasma
  
  dxdt_P = J14 - J15- J27 + J40;
  
  // Extracelluar phosphate (mmol) 
  dxdt_ECCPhos = J41  - J42 - J48 + J53 - J54 + J56;
  
  // Oral calcium 
  // CMT: T  UNITS: mmol 
  // J40 --> flux from gut to plasma 
  // F11 == T85 by definition 
  dxdt_T = OralCa*T85 - J40;
  
  // Calcitriol-dependent Ca2+ absorption 
  // CMT: R 
  dxdt_R = T36*(1- R) - T37*R;
  
  // Hydroxyapatite 
  dxdt_HAp = kHApIn*Osteoblast - kLShap*HAp;
  
  //
  // Estrogen piece
  //
  
  double AGE = ageENTER + T_0;
  
  double ageONSET = ageDONE-menoDUR;
  
  if(AGE < ageONSET) kinEST = koutEST * pow((AGE/ageENTER),ageGAM);
  
  if(AGE >= ageONSET) kinEST = koutEST * pow((AGE/ageENTER),ageGAM) * (1 - age50 * (pow((AGE-ageONSET),2)/(pow((menoDUR/2),2) + pow((AGE-ageONSET),2))));
  
  dxdt_EST = (kinEST - koutEST * EST)*ESTON;
  
  //dxdt_E2FRAC= 0;
  
  // Osteoblasts were considered to exist as two populations: fast and slow.
  //     Fast and slow refer to removal rates (kbfast, kbslow); input assumed to
  //     be the same for each.  Total osteoblasts = OBfast + OBslow 
  //Estrogen effect added: E2dosePicB1 --> PicOBkbEff --> kbprime --> kbslow, kbfast and Frackb2
  dxdt_OBfast = (bigDb/PicOB)*D*FracOBfast*Frackb2  - kbfast*OBfast; // 
  
  dxdt_OBslow = (bigDb/PicOB)*D*(1-FracOBfast)*Frackb - kbslow*OBslow;
  
  // translation compartment to delay SCLER effect on OB 
  
  double kin_T = kout_T;
  
  double EC50SCLER =  (exp(log(EMAXSCLER-1)/SCLEROBgam));
  
  dxdt_trans = kin_T*(EMAXSCLER*pow(SCLEREFF,SCLEROBgam)/(pow(EC50SCLER,SCLEROBgam)+pow(SCLEREFF,SCLEROBgam))) - kout_T*trans;
  
  //P1NPdum = ((2050 * pow((BSAPdum), 1.8) / ((pow(467, 1.8) + pow((BSAPdum),1.8))))-20); 
  
  // d/dt(extracellular phosphate) = J41 -  J42 - J48 + J53 - J54 + J56
  //     d/dt(intracellular phosphate) = J54  -  J56
  //
  //     The exchange fluxes of PO4 between ECF and bone (J41 and J42)
  //     set same as respective Ca fluxes but multiplied by the stoichiometric
  //     factor of 0.464, reported as the molar ratio of P / Ca in hydroxyapatite.
  //
  //     d/dt(dietary phosphate) = OralPhos*F12 -J53
  //     Influx set ~ 1 g (10.5 mmol) of PO4 daily.
  //     Bioavailability assumed = 0.7
  
  dxdt_PhosGut = OralPhos *F12 - J53;
  
  dxdt_IntraPO = J54 - J56;
  
  // OC: Active Osteoclasts 
  dxdt_OC = kinOC2 - KLSoc*OC;
  
  // D = ROB1; Responding Osteoblasts 
  dxdt_ROB1 = ROBin * pow(1/EST,robGAM) - pow((SCLEREFF),gammaDr)*KPT*ROB1 ;
  
  // Latent TGF-beta pool, production dependent on osteoblast function 
  dxdt_TGFB = kinTGF*(pow((Osteoblast/OB0),OBtgfGAM)) * pow(1/EST,tgfbGAM) - koutTGFeqn * pow(EST,tgfbactGAM);
  
  // active TGF-beta pool, production dependent on osteoclast function
  //  koutTGFeqn = koutTGF*TGFB*(pow((Osteoclast/OC0), OCtgfGAM))
  
  dxdt_TGFBact = koutTGFeqn * pow(EST,tgfbactGAM) - koutTGFact*TGFBact;
  
  // OCY: Active Osteocytes 
  // in amt is some fraction of osteoblasts leaving the system 
  
  dxdt_OCY = OB*FracOCY*OCY0 - koutOCY*OCY;
  
  
  //
  // BMD BONE MINERAL DENSITY 
  //
  // Solve for kinBMD based on OC=OC0, OB=OB0, and BMD=BMD0 
  
  //Lumbar spine
  double kinBMDls =  koutBMDls*BMDls0;
  
  dxdt_BMDls = kinBMDls * pow(OB/OB0,gamOB) - koutBMDls * pow(OC/OC0,gamOCls) * BMDls;
  
  //Lumbar spine with DENOSUMAB
  double kinBMDlsDEN =  koutBMDlsDEN*BMDlsDEN0;
  
  dxdt_BMDlsDEN = kinBMDlsDEN * pow(OB/OB0,gamOB) - koutBMDlsDEN * pow(OC/OC0,gamOClsDEN) * BMDlsDEN;
  
  //Femoral neck
  double baseBMDfn = baseBMDfnBAS; // (1+baseBMDfnETH*ETHN)*(1+baseBMDfnBMI*(BMI-27));
  double gamOCfn = gamOCfnBAS;     // (1+gamOCfnETH*ETHN)*(1+gamOCfnBMI*(BMI-27)); 
  double koutBMDfn = koutBMDfnBAS; // (1+koutBMDfnETH*ETHN)*(1+koutBMDfnBMI*(BMI-27));
  
  double kinBMDfn =  koutBMDfn*BMDfn0;
  
  dxdt_BMDfn = kinBMDfn * pow(OB/OB0,gamOB) - koutBMDfn * pow(OC/OC0,gamOCfn) * BMDfn;
  
  //Femoral neck with DENOSUMAB
  
  double kinBMDfnDEN =  koutBMDfnDEN*BMDfnDEN0;
  
  dxdt_BMDfnDEN = kinBMDfnDEN * pow(OB/OB0,gamOB) - koutBMDfnDEN * pow(OC/OC0,gamOCfnDEN) * BMDfnDEN;
  
  
  //Total hip with DENOSUMAB
  
  double kinBMDthDEN =  koutBMDthDEN*BMDthDEN0;
  
  dxdt_BMDthDEN = kinBMDthDEN * pow(OB/OB0,gamOB) - koutBMDthDEN * pow(OC/OC0,gamOCthDEN) * BMDthDEN;
  
  
  // BMD for anabolic therapies 
  
  // SCLEROSTIN    
  
  
  //Lumbar spine
  
  double kinBMDlsSCLER =  koutBMDSCLER*BMDlsSCLER0;
  
  double kin_BMDdel = kout_BMDdel;
  
  dxdt_delBMDls = kin_BMDdel * pow(OB/OB0,gamOBSCLERls) - kout_BMDdel*delBMDls;
  
  dxdt_BMDlsSCLER = kinBMDlsSCLER * delBMDls - koutBMDSCLER * pow(OC/OC0,gamOCSCLER) * BMDlsSCLER;
  
  //Femoral neck
  
  double kinBMDfnSCLER =  koutBMDSCLER*BMDfnSCLER0;
  
  dxdt_delBMDfn = kin_BMDdel* pow(OB/OB0,gamOBSCLERfn) - kout_BMDdel*delBMDfn;
  
  dxdt_BMDfnSCLER = kinBMDfnSCLER * delBMDfn  - pow(OC/OC0,gamOCSCLER) * koutBMDSCLER * BMDfnSCLER;
  
  //Total hip
  
  double kinBMDthSCLER =  koutBMDSCLER*BMDthSCLER0;
  
  
  dxdt_delBMDth = kin_BMDdel * pow(OB/OB0,gamOBSCLERth) - kout_BMDdel*delBMDth;
  
  dxdt_BMDthSCLER = kinBMDthSCLER * delBMDth - koutBMDSCLER * pow(OC/OC0,gamOCSCLER) * BMDthSCLER;
  
  // 
  //     TERIPARATIDE      
  
  //Lumbar spine
  
  
  
  double kinBMDlsTERI =  koutBMDlsTERI*BMDlsTERI0;
  
  double kin_BMDdelTERI = kout_BMDdelTERI;
  
  dxdt_delBMDlsTERI = kin_BMDdelTERI * pow(OB/OB0,gamOBTERIls) - kout_BMDdelTERI*delBMDlsTERI;
  
  dxdt_BMDlsTERI = kinBMDlsTERI * delBMDlsTERI - koutBMDlsTERI * pow(OC/OC0,gamOClsTERI) * BMDlsTERI;
  
  //Femoral neck
  
  double kinBMDfnTERI =  koutBMDfnTERI*BMDfnTERI0;
  
  dxdt_delBMDfnTERI = kin_BMDdelTERI* pow(OB/OB0,gamOBTERIfn) - kout_BMDdelTERI*delBMDfnTERI;
  
  dxdt_BMDfnTERI = kinBMDfnTERI * delBMDfnTERI  - pow(OC/OC0,gamOCfnTERI) * koutBMDfnTERI * BMDfnTERI;
  
  //Total hip
  
  double kinBMDthTERI =  koutBMDthTERI*BMDthTERI0;
  
  dxdt_delBMDthTERI = kin_BMDdelTERI* pow(OB/OB0,gamOBTERIth) - kout_BMDdelTERI*delBMDthTERI;
  
  dxdt_BMDthTERI = kinBMDthTERI * delBMDthTERI  - pow(OC/OC0,gamOCthTERI) * koutBMDthTERI * BMDthTERI;
  
  
  // BMD with combination therapy /
  
  // DEN/TERI  
  double kinBMDlsCOMBO = koutBMDlsCOMBO;
  
  dxdt_BMDlsCOMBO = 0;
  if(DEN_TERI_COMBO==1) dxdt_BMDlsCOMBO = kinBMDlsCOMBO * pow(OB/OB0,gamOB) - koutBMDlsCOMBO* pow(OC/OC0,gamOClsDEN_TERI)*BMDlsCOMBO;
  
  double kinBMDthCOMBO = koutBMDthCOMBO;
  dxdt_BMDthCOMBO = 0;
  if(DEN_TERI_COMBO==1) dxdt_BMDthCOMBO = kinBMDthCOMBO * pow(OB/OB0,gamOB) - koutBMDthCOMBO* pow(OC/OC0,gamOCthDEN_TERI)*BMDthCOMBO;
  
  double kinBMDfnCOMBO = koutBMDfnCOMBO;
  dxdt_BMDfnCOMBO = 0;
  if(DEN_TERI_COMBO==1) dxdt_BMDfnCOMBO = kinBMDfnCOMBO * pow(OB/OB0,gamOB) - koutBMDfnCOMBO* pow(OC/OC0,gamOCfnDEN_TERI)*BMDfnCOMBO;
  
  
  // FRACTURE PROBABILITY 
  //
  // BMDhat = 0.8 g/cm^2 
  // postMenoAgehat = 20 y 
  // BMIhat = 27.1 kg/m^2 
  //req dataset items = agepostmeno (at baseline), radFracInc, BMI, lsBMDbase 
  // code in switches: 
  // DEN_TERI_COMBO = 0 : no drug
  // DEN_TERI_COMBO = 1 : combinination TERI/DEN
  // DEN_TERI_COMBO = 2 : TERI only
  // DEN_TERI_COMBO = 3 : DEN only
  // DEN_TERI_COMBO = 4 : SCLER only 
  
  
  double PBO_BMD = 1;
  double betaDrug = 0;
  if (DEN_TERI_COMBO==0) { 
    PBO_BMD = BMDls; 
  }   
  
  double COMBO_BMD = 1;
  if (DEN_TERI_COMBO==1) { 
    COMBO_BMD = BMDlsCOMBO; 
    betaDrug = betaDrug_COMBO;
  }  
  double TERI_BMD = 1;
  double TERI_BMDth = 1;
  double TERI_BMDfn = 1;
  if (DEN_TERI_COMBO==2) { 
    TERI_BMD = BMDlsTERI; 
    betaDrug = betaDrug_TERI;
    TERI_BMDth = BMDthTERI;
    TERI_BMDfn = BMDfnTERI;
  }     
  
  double DEN_BMD = 1; 
  double DEN_BMDth = 1; 
  double DEN_BMDfn = 1; 
  if (DEN_TERI_COMBO==3) { 
    DEN_BMD = BMDlsDEN; 
    betaDrug = betaDrug_DENO;
    DEN_BMDth=BMDthDEN;
    DEN_BMDfn=BMDfnDEN;
  } 
  
  double SCLER_BMD = 1; 
  double SCLER_BMDth = 1; 
  double SCLER_BMDfn = 1; 
  
  if (DEN_TERI_COMBO==4) { 
    SCLER_BMD = BMDlsSCLER; 
    betaDrug = betaDrug_SCLER;
    SCLER_BMDth=BMDthSCLER;
    SCLER_BMDfn=BMDfnSCLER;
  } 
  
  
  
  if(SCLER_DEN_SEQ==1){
    SCLER_BMD = BMDlsSCLER; 
    SCLER_BMDth=BMDthSCLER;
    SCLER_BMDfn=BMDfnSCLER;
    DEN_BMD = BMDlsDEN;
    DEN_BMDth=BMDthDEN;
    DEN_BMDfn=BMDfnDEN;
  }
  
  
  double lsBMDtot = (PBO_BMD + COMBO_BMD + TERI_BMD + DEN_BMD + SCLER_BMD)-4; //(PBO_BMD + COMBO_BMD + TERI_BMD + DEN_BMD + SCLER_BMD);
  
  double thBMDtot = (1 + TERI_BMDth + DEN_BMDth + SCLER_BMDth)-3; // added 09/25
  
  double fnBMDtot = (1 +  TERI_BMDfn + DEN_BMDfn + SCLER_BMDfn)-3; // added 09/25
  
  
  // Hazard Model Stuff // 
  
  double drug = betaDrug;
  
  double lsBMD = lsBMDtot*lsBMDbase; //nominal BMD
  
  double lsBMDCFB = (lsBMDtot*lsBMDbase) - lsBMDbase; //CFB BMD
  
  double postMenoAge = agepostmeno + T_0/365.25/24; 
  
  
  double Hazard = baseHazard*exp(HazBMDCov*log(lsBMD/lsBMDbase) + 
    HazpostMenoAgeCov*(postMenoAge-postMenoAgeRef )+ 
    HazradFracCov*radFracInc+HazBMICov*(BMI-bmiRef) + betaDrug);
  //double Hazard = exp(baseHazard *(1+HazBMDcov*(BMDfn-0.8)));
  // double Hazard = exp(baseHazard) * pow(BMDfn/0.8,HazBMDcov); 
  
  dxdt_cumHazard = Hazard / 365.25/24; 
  
  double Survival = exp(-cumHazard);
  
  //Treatment effects
  //
  
  // L: RANK-L 
  dxdt_L = kinL- koutL*L - k1*O*L + k2*N - k3*RNK*L + k4*M -  kdenosl*DENMOL*L;
  
  // RNK: RANK 
  dxdt_RNK = kinRNK*pow(TGFBact,kinRNKgam) - koutRNK*RNK - k3*RNK*L  + k4*M;
  
  // M: RANK - RANK-L complex 
  dxdt_M = k3*RNK*L - k4*M;
  
  // N:  - RANK-L complex 
  dxdt_N = k1*O*L - k2*N;
  
  // O: OPG 
  dxdt_O = pO*pow((SCLEREFF),gammaOPG) - k1*O*L + k2*N - kO*O;
  
  // *pow((SCLEREFF+1),gammaOPG) 
  
  // d/dt(Q) = J15-J14+ J14a-J15a  ;Q=exchangeable bone Ca
  //       d/dt(Qbone) = -J14a + J15a    ;Qbone=non(immediately)-exchangeable bone Ca
  //       ~ 99% of the total Ca stored in bone; approximately 100 mmol of 25000 - 30000 mmol
  //       of total skeletal Ca considered immediately exchangeable with plasma Ca
  
  dxdt_Q = J15 - J14 + J14a - J15a;
  
  dxdt_Qbone = J15a - J14a;
  
  //
  //     d/dt(BCL2) = bcl2Kin - bcl2Kout*BCL2
  //     d/dt(RX2) = RX2Kin - RX2Kout*RX2
  //     d/dt(CREB) = crebKin - crebKout*CREB
  //
  //     Initial conditions set empirically:
  //     both RX2 and CREB started at 10, BCL2 was the product (10*10=100).
  //
  //     BCL2 assumed half-life of 1 hour: rate const set to 0.693 (bcl2Kout)
  //
  //     BCL2 affected osteoblast survival:
  //     decreases elimination rate constant for OB (kbprime)
  
  dxdt_RX2 = RX2Kin - RX2Kout*RX2 ;
  
  dxdt_CREB = crebKin - crebKout*CREB;
  
  dxdt_BCL2 = bcl2Kout*CREB*RX2 - bcl2Kout*BCL2;
  
  // Teriparatide PK info 
  // TERIPKIN = TERISC * TERICL/TERIVC 
  dxdt_TERISC = -TERIPKIN;
  dxdt_TERICENT  = TERIPKIN - TERICENT*TERIKA;
  
  // GENERAL PK COMPARTMENT 
  
  double GFRend = GFR0 - GFRdelta/16.667;
  double GFRtau_ = GFRtau*8766;
  
  double kGFR = -log(GFRend/GFR0)/GFRtau_;
  
  dxdt_GFR = -kGFR*GFR;
  
  // NONLINEAR PIECE 
  double PKCLNL = PKVMAX/(PKKM+PKCP);
  
  // PKGUT 
  dxdt_PKGUT = -PKKA*PKGUT;
  
  // PKCENT 
  dxdt_PKCENT = PKKA*PKGUT + PKQ1*PKPER1/PKVP1 + PKQ2*PKPER2/PKVP2 - (PKQ1+PKQ2+PKCL+PKCLNL)*PKCENT/PKVC;
  
  // PKPER1 
  dxdt_PKPER1 = PKQ1*PKCENT/PKVC - PKQ1*PKPER1/PKVP1;
  
  // PKPER2 
  dxdt_PKPER2 = PKQ2*PKCENT/PKVC - PKQ2*PKPER2/PKVP2;
  
  
  // DENOSUMAB PK 
  double  DENCLNL =  (DENVMAX/(DENKM+DENCP));
  dxdt_DENSC =  -DENKA*DENSC;
  dxdt_DENCENT = DENKA*DENSC + DENQ*DENPER/DENVP - (DENQ+DENCL+DENCLNL)*DENCENT/DENVC ;
  dxdt_DENPER = DENQ*DENCENT/DENVC - DENQ*DENPER/DENVP;
  
  
  // SCLEROSTIN AB PK AND SCLEROSTIN 
  
  double SOSTCLNL = SOSTVMAX/(SOSTKM+SOSTCP);
  dxdt_SOSTSC = -SOSTKA*SOSTSC;
  dxdt_SOSTCENT = SOSTKA*SOSTSC - SOSTCLNL*SOSTCENT-(SOSTCL*(SOSTCP))-(SOSTQ*(SOSTCP)) + (SOSTQ*(SOSTPER/SOSTVP));
  dxdt_SOSTPER = (SOSTQ*(SOSTCENT/SOSTVC))-(SOSTQ*SOSTPER/SOSTVP);
  dxdt_SCLER = SOSTKIN-SOSTKOUT*SCLER-(SOSTK0-SOSTKOUT)*SCLER*(SOSTCENT/SOSTVC)/(SOSTKM+(SOSTCENT/SOSTVC));
  
  // BSAP: % of baseline 
  
  double BSAP = (OB/OB0*100); 
  
  
  // Collects calcium in urine - cumulative rate 
  // This is a differential equation 
  dxdt_UCA = J27;
  


[TABLE]
capture BSAP_OB = (OB/OB0*100);
capture SCLERmeas = (100*OCY/OCY_0);
capture inrate = (((kbfast*OBfast)+(kbslow*OBslow)) * FracOCY);
capture P1NPsim = ((2050 * pow((BSAP), 1.8) / ((pow(467, 1.8) + pow((BSAP),1.8))))-20.4181); 
capture ROB1out = (pow((SCLEREFF),gammaDr)*KPT*ROB1);
capture picROB = ((E0PicROB + EmaxPicROB*pow(TGFBact,PicROBgam)/(pow(TGFBact,PicROBgam) + pow(EC50PicROB,PicROBgam))));
capture CTXsim = (OC/OC0*100);
 
capture lsBMDsimSCLER = (BMDlsSCLER*100);
capture lsBMDsimTERI = (BMDlsTERI*100);
capture lsBMDsimDEN = (BMDlsDEN*100);
capture thBMDsimSCLER = (BMDthSCLER*100);
capture thBMDsimTERI = (BMDthTERI*100);
capture thBMDsimDEN = (BMDthDEN*100);
capture fnBMDsimSCLER = (BMDfnSCLER*100);
capture fnBMDsimTERI = (BMDfnTERI*100);
capture fnBMDsimDEN = (BMDfnDEN*100);
capture lsBMDsimCOMBO = (BMDlsCOMBO*100);
capture thBMDsimCOMBO = (BMDthCOMBO*100);
capture fnBMDsimCOMBO = (BMDfnCOMBO*100);


double DV = ((2050 * pow((BSAP), 1.8) / ((pow(467, 1.8) + pow((BSAP),1.8))))-20.4181);
if(TYPE==1) DV = (OC/OC0*100);

double SclerBMDpred;
if(Scler_BMD_type==0)  SclerBMDpred = (BMDlsSCLER*100); 
if(Scler_BMD_type==1)  SclerBMDpred = (BMDthSCLER*100);
if(Scler_BMD_type==2)  SclerBMDpred = (BMDfnSCLER*100);

double DenoBMDpred;
if(Deno_BMD_type==0)  DenoBMDpred = (BMDlsDEN*100); 
if(Deno_BMD_type==1)  DenoBMDpred = (BMDthDEN*100);
if(Deno_BMD_type==2)  DenoBMDpred = (BMDfnDEN*100);

double TeriBMDpred;
if(Teri_BMD_type==0)  TeriBMDpred = (BMDlsTERI*100); 
if(Teri_BMD_type==1)  TeriBMDpred = (BMDthTERI*100);
if(Teri_BMD_type==2)  TeriBMDpred = (BMDfnTERI*100);

double ComboBMDpred;
if(Combo_BMD_type==0)  ComboBMDpred = (BMDlsCOMBO*100); 
if(Combo_BMD_type==1)  ComboBMDpred = (BMDthCOMBO*100);
if(Combo_BMD_type==2)  ComboBMDpred = (BMDfnCOMBO*100);


capture PTHpm = (PTH/V1);   //* PTH conc in pM */
capture PTHc = (PTH/V1*9.4) ;     //* PTH conc in pg/mL */
capture CAchange = (100*P/P0) ;   /*/ Calcium (total, serum), %change from baseline */
capture OBchange = 100*OB/OB0;  //* bone-specific alkaline phosphatase, %change from baseline */
capture OCchange = 100*OC/OC0;  //* serum CTx, %change from baseline */
capture CaC =  P/V1 ;   //* Calcium (total, serum) in mM */
capture BMDlsDENchange = (BMDlsDEN-1)*100;
capture OBtot = OBfast + OBslow;
capture OCfrac = OC/OC_0;
capture OBfrac = (OBfast+OBslow)/(OBfast_0 + OBslow_0);

$CAPTURE DENMOL T43 DENCP SOSTCP DenoBMDpred TeriBMDpred ComboBMDpred SclerBMDpred
