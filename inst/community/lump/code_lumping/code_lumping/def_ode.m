
% (1) Define second (derived) parameter values
% (2) Define nonlinear ODEs

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%  Lumped version of
%    Peterson MC, Riggs MM (2010) Bone 46:49-63
%                        +
%    Peterson MC, Riggs MM (2012) CPT Pharmacometrics Syst Pharmacol 1:e14
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

%% (1) Define second (derived) parameter values

Osteoblast = y(29-2) + y(30-2);

%----------------------------------------------------------------------
% kin21 (zero-order)
%----------------------------------------------------------------------
kinRNK = (koutRNK*IC0(21-2) + k3*IC0(21-2)*IC0(22-2) - k4*IC0(24-2)) / IC0(20-2)^kinRNKgam ;
kin21 = kinRNK*y(20-2)^kinRNKgam;

%----------------------------------------------------------------------
% kout21 (first-order)
%----------------------------------------------------------------------
% koutRNK

%----------------------------------------------------------------------
% kin04
%----------------------------------------------------------------------
k0412 = CaDay/24;
T13 = k0412/IC0(12);             % typical kout12
J14OC50= exp(log((J14OCmax*IC0(18-2)^J14OCgam/T13) - IC0(18-2)^J14OCgam)/J14OCgam);            % d181204 in the paper
Osteoclast =  y(18-2);
% H181204 in dydt(4)
OCeqn = (J14OCmax*Osteoclast^J14OCgam)/(Osteoclast^J14OCgam + J14OC50^J14OCgam);

% last product in v1204
MOCratio = y(24-2)/Osteoclast;
MOCratio0 = IC0(24-2)/IC0(18-2);
MOCratioEff = (MOCratio/MOCratio0)^MOCratioGam       ;

J14OCdepend = IC0(12) * FracJ14 * OCeqn * MOCratioEff;   

% v1204 in the paper
J14 = k0412*(1-FracJ14) + J14OCdepend;
% kin04 = J14;

%----------------------------------------------------------------------
% kout04
%----------------------------------------------------------------------
T15 = k0412/IC0(4);              % typical kout04
% v0412 in the paper, transit to bone Ca depends on HAp amount (y14)
J15 = T15*(1-FracJ15) + T15*FracJ15*y(14);     % product with y(4) in ODE
kout0412 = J15;

% for v040u, 2-H0604
C8 = y(6)/Vva;
C8_0 = IC0(6)/Vva;
% H0604 in dydt(4)
T10 = 2 - T7*C8/(C8+T9);

% for v040u, 0.3*GFR*A(4)/Vva
C1 = y(4)/Vva ;
C1_0 = IC0(4)/Vva;
CaFilt = 0.6*0.5*GFR/Vva;                          % divide by Vva so A can be used with coefficient

% for v040u, H040u
ReabsMax = (0.3*GFR*C1_0 - 0.149997)*(Reabs50 + C1_0) / C1_0;
% kout in H040u in dydt(4)
k040u = ReabsMax/(Reabs50 + C1)/Vva;            % divide by Vva so A can be used with coefficient

% for v040u, H07040u
C4 = y(7)/Vva;
C7_0 = IC0(7)/Vva;
T17 = C7_0*T16 - C7_0;           % d07040u in the paper
% H07040u in dydt(4)
ReabsPTHeff = T16*C4/(C4+T17) ;

CaReabsActive = k040u*ReabsPTHeff ;
% 2nd parenthesis
T20 = CaFilt - CaReabsActive;             % coefficient!!!!!!!!!!!!!!!!!

% v040u (without y4) in dydt(4)
J27a = T10*T20;                       % final coefficient!!!!!!!!!!!!!!!!!
if J27a*y(4)<0
    J27 = 0;
else
    J27 = J27a;
end
kout040u = J27;

% kout04 = kout0412 + kout040u;

%----------------------------------------------------------------------
% kin02
%----------------------------------------------------------------------
% H+0602 in dydt(2)
T36 = T33 + (T34-T33)*(C8^CaPOgam/(T35^CaPOgam+ C8^CaPOgam));
% H-0602 in dydt(2)
T37 = T34 - (T34-T33)*(C8^CaPOgam/(T35^CaPOgam+ C8^CaPOgam));

% dydt(2) = T36*(1-y2) - T37*y2 = "T36 - (T36+T37)*y2"           !!!!!!!!!!!!!!!!!!!!!!
kin02 = T36;

%----------------------------------------------------------------------
% kout02
%----------------------------------------------------------------------
kout02 = T36 + T37;

%----------------------------------------------------------------------
% kin05
%----------------------------------------------------------------------
J41 = 0.464*J14  ;               % v1405 in the paper
kin1405 = J41;

% a part of v050u without y5 is included to kin05 (zero-order)
if y(5)/Vva < T46
    phi050u = y(5)/Vva;
else
    phi050u = T46;
end
T47 = 0.88*GFR*phi050u;
kin050u = T47;

kin05 = kin1405 + kin050u;

%----------------------------------------------------------------------
% kout05
%----------------------------------------------------------------------
J42 = 0.464*kout0412;                % kout in v0514 in the paper
kout0514 = J42;

kout050u = 0.88*GFR/Vva;

kb0305 = T52;
% J53 = T52*y(3);

kb0508 = T49/Vva;
% J54 = T49*C2;

kb0805 = T55;
% J56 = T55*y(8);

%----------------------------------------------------------------------
% kin19
%----------------------------------------------------------------------
kinTGF = koutTGF0*IC0(19-2);

kin19 = kinTGF*(Osteoblast/(IC0(29-2)+IC0(30-2)))^OBtgfGAM;

%----------------------------------------------------------------------
% kout19
%----------------------------------------------------------------------
koutTGF = koutTGF0;
ratio19 = ((y(19-2)/IC0(19-2))^koutTGFGam);
koutTGFeqn = koutTGF * ratio19 * ((Osteoclast/IC0(18-2))^OCtgfGAM);

%----------------------------------------------------------------------
% kin20
%----------------------------------------------------------------------
% koutTGFeqn*y(19-2)

%----------------------------------------------------------------------
% kout20
%----------------------------------------------------------------------
koutTGFact = 1000*koutTGF0;

%----------------------------------------------------------------------
% kin16
%----------------------------------------------------------------------
kin16 = kb*(IC0(29-2)+IC0(30-2));              % typical kin16

Dr = kin16/Pic0  	  ;          % for the 1st term in dydt(16-1)

E0PicROB = FracPicROB*Pic0;                % r2016 in H+2016
EC50PicROBparen= (EmaxPicROB*IC0(20-2)^PicROBgam / (Pic0 - E0PicROB)) - IC0(20-2)^PicROBgam;
EC50PicROB = exp(log(EC50PicROBparen)/PicROBgam);
% H+2016 in dydt(16-1)
PicROB = E0PicROB + EmaxPicROB*y(20-2)^PicROBgam/(y(20-2)^PicROBgam + EC50PicROB^PicROBgam);

ROBin2 = Dr*PicROB;
ROBin = ROBin2;

%----------------------------------------------------------------------
% kout16 (also used in dydt(29-2) (17a)
%----------------------------------------------------------------------
bigDb = kb*(IC0(29-2)+IC0(30-2))*Pic0/IC0(16-1);     % for the 1st term in dydt(29-2) (17a)

E0PicOB = FracPicOB*Pic0;
EC50PicOBparen = (EmaxPicOB*IC0(20-2)^PicOBgam/(Pic0 - E0PicOB)) - IC0(20-2)^PicOBgam;
EC50PicOB = exp(log(EC50PicOBparen)/PicOBgam);
% H+2017 in dydt(16-1)
PicOB = E0PicOB + EmaxPicOB*y(20-2)^PicOBgam / (y(20-2)^PicOBgam + EC50PicOB^PicOBgam);

KPT =1*(bigDb/PicOB)  ;

%----------------------------------------------------------------------
% kin18
%----------------------------------------------------------------------
kin18 = Da*IC0(18-2);             % typical kin18

PicOCkin = Pic0;

% EC50MeffOC = exp(log(IC0(24-2)^kinOCgam*EmaxMeffOC/(1-E0Meff) - IC0(24-2)^kinOCgam)/kinOCgam);
% H+2418S in dydt(18-2)

% gamma>8, range of y(24-2) is less than EC50 == 0.0002, so effect is set as zero
MeffOC = E0Meff; % + (EmaxMeffOC * y(24-2)^kinOCgam/(y(24-2)^kinOCgam + EC50MeffOC^kinOCgam));

kinOC2 = kin18*PicOCkin*MeffOC;

%----------------------------------------------------------------------
% kout18
%----------------------------------------------------------------------
E0PicOC = FracPicOC*Pic0;
EC50PicOCparen = (EmaxPicOC*IC0(20-2)^PicOCgam/(Pic0 - E0PicOC)) - IC0(20-2)^PicOCgam;
EC50PicOC = exp(log(EC50PicOCparen)/PicOCgam);         
% H+2018D in dydt(18-2)
PicOC = E0PicOC + ((EmaxPicOC*y(20-2)^PicOCgam)/(y(20-2)^PicOCgam + EC50PicOC^PicOCgam));

% use y(24-2) rather than y(22-2)
PiL0 = (k3/k4)*IC0(22-2);                % IC0(24-2)
PiL = y(24-2)/10;
EC50survInPar = (E0RANKL - EmaxL)*(PiL0^LsurvOCgam/(E0RANKL - 1)) - PiL0^LsurvOCgam;
EC50surv = exp(log(EC50survInPar)/LsurvOCgam);
% H-2218D in dydt(18-2)
LsurvOC = E0RANKL - (E0RANKL - EmaxL)*(PiL^LsurvOCgam/(PiL^LsurvOCgam + EC50surv^LsurvOCgam));
   
KLSoc = Da*PicOC*LsurvOC;              % #*3

%----------------------------------------------------------------------
% Ca transit between IC (12) and non-IC (13)
%----------------------------------------------------------------------
k15a = k14a*IC0(13)/IC0(12);                    % k1213 in dydt(12), ok

%----------------------------------------------------------------------
% kout14
%----------------------------------------------------------------------
kLShap = 1/HApMRT;

%----------------------------------------------------------------------
% kin14
%----------------------------------------------------------------------
kHApIn = kLShap/(IC0(29-2)+IC0(30-2));

kHApIn2 = kHApIn*Osteoblast;

%----------------------------------------------------------------------
% kin22
%----------------------------------------------------------------------
kinLbase = koutL*IC0(22-2);                % typical kin22
   
OsteoEffect = (Osteoblast/(IC0(29-2)+IC0(30-2)))^OsteoEffectGam ;

PTH50 = EmaxLpth*C7_0 - C7_0 ;
PTHconc = C4;
LpthEff = EmaxLpth*(PTHconc) / ((PTH50*OsteoEffect^TESTPOWER) + (PTHconc))  ;

kinL = kinLbase*OsteoEffect*LpthEff;

%----------------------------------------------------------------------
% kout22 (binding to RANK (21-2) and OPG (23-2) should also be considered
%----------------------------------------------------------------------
% koutL

%----------------------------------------------------------------------
% kin23
%----------------------------------------------------------------------
pObase = kO*IC0(23-2);                     % typical kin23

D = y(16-1);
pO = pObase*(D/IC0(16-1))*((PTHconc+(opgPTH50*(D/IC0(16-1))))/(2*PTHconc))+ IO;

%----------------------------------------------------------------------
% kout23 (binding to RANKL (22-2) should also be considered
%----------------------------------------------------------------------
% kO

%----------------------------------------------------------------------
% kin26
%----------------------------------------------------------------------
RX2Kin = RX2Kout0*IC0(26-2);                    % k26S in the paper

%----------------------------------------------------------------------
% kout26
%----------------------------------------------------------------------
EC50PTHRX2x = ((EmaxPTHRX2x*C7_0)/(RX2Kout0 - E0rx2Kout)) - C7_0;
% H+0726D in dydt(26-2)
RX2Kout = E0rx2Kout + EmaxPTHRX2x*PTHconc/(PTHconc+EC50PTHRX2x);

%----------------------------------------------------------------------
% kin27
%----------------------------------------------------------------------
crebKin0= crebKout*IC0(27-2);        % typical kin27

EC50PTHcreb = ((EmaxPTHcreb*C7_0)/(1-E0crebKin)) -  C7_0;
% H+0727S
H0727S = E0crebKin + EmaxPTHcreb*PTHconc/(PTHconc+EC50PTHcreb);

crebKin = crebKin0 * H0727S;

%----------------------------------------------------------------------
% kout27
%----------------------------------------------------------------------
% crebKout

%----------------------------------------------------------------------
% kin28
%----------------------------------------------------------------------
bcl2Kin = y(26-2)*y(27-2)*bcl2Kout;

%----------------------------------------------------------------------
% kout28
%----------------------------------------------------------------------
% bcl2Kout

%----------------------------------------------------------------------
% kin09
%----------------------------------------------------------------------
C2 = y(5)/Vva ;
C2_0 = IC0(5)/Vva;
PO4inhPTH = (C2/C2_0)^PO4inhPTHgam  ;
T66 = (T67^AlphOHgam + C7_0^ AlphOHgam )/C7_0^ AlphOHgam ;
% H0709 in dydt(9)
T68 = T66*C4^AlphOHgam/(T67^AlphOHgam*PO4inhPTH+C4^AlphOHgam) ;

PhosEffTop = (PhosEff0 - 1)*( C2_0^PhosEffGam + PhosEff50^PhosEffGam ) ;
PhosEffBot = PhosEff0 * C2_0^PhosEffGam ;
PhosEffMax =  PhosEffTop / PhosEffBot;
% H-0509 in dydt(9)
PhosEff = PhosEff0 - (PhosEffMax*PhosEff0 * C2^PhosEffGam /(C2^ PhosEffGam  + PhosEff50^PhosEffGam));
if y(5) >= IC0(5)
        PhosEffect = PhosEff;
else
        PhosEffect = 1;
end

SE = T65*T68*PhosEffect;

%----------------------------------------------------------------------
% kout09
%----------------------------------------------------------------------
% T64

%----------------------------------------------------------------------
% kin10
%----------------------------------------------------------------------
CaConc0 = IC0(4)/Vva;
CaConc = y(4)/Vva;
ScaEff =  (CaConc0/CaConc)^ ScaEffGam ;

T72 = C8_0 * ScaEff  ;                    % C8_0 = dT0604 in the paper = 90

T73 = T71 * (C8 - T72) ;

T74 = (exp(T73) - exp(-T73)) / (exp(T73) + exp(-T73));

% T+0604 in dydt(10)
Tp0604 = 1 + T74;
T75 = T70 * (0.85 * Tp0604 + 0.15)  ;
% T-0604 in dydt(10)
Tm0604 = 1 - T74;
T76 = T70 * (0.85 * Tm0604 + 0.15);

% dydt(10) = (1-y10)*T76 - y10*T75 = "T76 - (T76+T75)*y10"           !!!!!!!!!!!!!!!!!!!!!!
kin10 = T76;

%----------------------------------------------------------------------
% kout10
%----------------------------------------------------------------------
kout10 = T76+T75;

%----------------------------------------------------------------------
% kin30 (17b)
%----------------------------------------------------------------------
% No zero-order input

%----------------------------------------------------------------------
% kout30 (17b)
%----------------------------------------------------------------------
E0PicOBkb = MultPicOBkb*Pic0;                  % a2017D
EmaxPicOBkb = FracPic0kb*Pic0;                 % r2017D
EC50PicOBparenKb = ((E0PicOBkb - EmaxPicOBkb)*IC0(20-2)^PicOBgamkb) / (E0PicOBkb - Pic0)  - IC0(20-2)^PicOBgamkb;
EC50PicOBkb = exp(log(EC50PicOBparenKb)/PicOBgamkb);
% H-2017D in dydt(29-2) (17a)
PicOBkb = E0PicOBkb - (E0PicOBkb  - EmaxPicOBkb)*y(20-2)^PicOBgamkb / (y(20-2)^PicOBgamkb + EC50PicOBkb^PicOBgamkb);

PicOBkbEff = PicOBkb/Pic0;

E0RUNX2kbEff= E0RUNX2kbEffFACT*kb;
RUNkbMax = E0RUNX2kbEff*RUNkbMaxFact;          % a2817D
INparen = (RUNkbMax * RUNX20^RUNkbGAM) / (E0RUNX2kbEff - kb) - RUNX20^RUNkbGAM;
RUNkb50 = exp(log(INparen)/RUNkbGAM);
if y(28-2) > 105
    RUNX2 = y(28-2) - 90;
else
    RUNX2 = 10;
end
% H2817D in dydt(29-2) (17a)
RUNX2kbPrimeEff = RUNkbMax*RUNX2^RUNkbGAM / (RUNX2^RUNkbGAM + RUNkb50^RUNkbGAM);

kbprime = E0RUNX2kbEff*PicOBkbEff - RUNX2kbPrimeEff;

kbslow = kbprime*Frackb;

kb30 = KPT*(1-FracOBfast)*Frackb;

%----------------------------------------------------------------------
% kin29 (17a)
%----------------------------------------------------------------------
% No zero-order input

%----------------------------------------------------------------------
% kout29 (17a)
%----------------------------------------------------------------------
kbfast = (kb*(IC0(29-2)+IC0(30-2)) + kbslow*IC0(29-2)*FracOBfast - kbslow*(IC0(29-2)+IC0(30-2))) / IC0(29-2);

Frackb2 = kbfast/kbprime;

kb29 = KPT*FracOBfast*Frackb2;

%----------------------------------------------------------------------
% kin01
%----------------------------------------------------------------------
T85Rpart = y(2)^T80/(y(2)^T80 + T81^T80) ;
% H+0201 in dydt(1)
T85 = T77*T85Rpart;
F11 = T85;

kin01 = OralCa*F11;

%----------------------------------------------------------------------
% kout01
%----------------------------------------------------------------------
T29 = (T28*T0 - 0.17533*T0)/0.17533;           % d0104
% H0104 in dydt(1)
T31 = T28*y(1)/(y(1)+T29);

T83 = y(2)/IC0(2);

kout0104 = T31*T83/(y(1)+T81) + T87;
% J40 = T31*y(1)*T83/(y(1) + T81) + T87*y(1);

%----------------------------------------------------------------------
% kin11
%----------------------------------------------------------------------
% use C8_0 rather than C8
INparenCtriol =((CtriolMax - CtriolMin) * C8_0^CtriolPTgam) / (CtriolMax - 1) - C8_0^CtriolPTgam ;
Ctriol50 = exp(log(INparenCtriol) / CtriolPTgam) ;
% H-0611
CtriolPTeff = CtriolMax - (CtriolMax - CtriolMin) * C8^CtriolPTgam / (C8^CtriolPTgam + Ctriol50^CtriolPTgam) ;

PTin= PTout * CtriolPTeff;

%----------------------------------------------------------------------
% kout11
%----------------------------------------------------------------------
% PTout

%----------------------------------------------------------------------
% kin07
%----------------------------------------------------------------------
INparenCa = (T58 - T61) * CaConc0^T59 / (T58 - C7_0*100) - CaConc0^T59  ;
T60 = exp(log(INparenCa) / T59)  ;
% H-04107 in dydt(7)
T63 =  T58 - (T58 - T61) * (CaConc)^T59 / ((CaConc)^T59 + T60^T59);

FCTD = (y(10) / IC0(10)) * y(11);

EPTH = T63 * FCTD;

% IPTH = 0.693*SC +  IPTHinf;
IPTH = 0;

SPTH = EPTH + IPTH ;

%----------------------------------------------------------------------
% kout07
%----------------------------------------------------------------------
kout = T57/14;                       % k7D in the paper

%----------------------------------------------------------------------
% kin03
%----------------------------------------------------------------------
kin03 = OralPhos*F12;

%----------------------------------------------------------------------
% kout03
%----------------------------------------------------------------------
% kb0305

%----------------------------------------------------------------------
% denosumab PK (y32 = total, C = free)
%----------------------------------------------------------------------
C = 0.5*((y(32-2)/V1 - y(34-2) - Kss) + sqrt((y(32-2)/V1 - y(34-2) - Kss)^2 + 4*Kss*y(32-2)/V1));
CLtot = CL + kint*V1*y(34-2)/(Kss+C);
drug = (kint-koutL)*C/(Kss+C);
mic34 = (kint-kdeg)*C/(Kss+C);

%----------------------------------------------------------------------
% BMD
%----------------------------------------------------------------------
BSAP = y(29-2) + y(30-2);
SCTX = y(18-2);

Kin_BMD  = 0.000146*IC_BMD0*(BSAP/(IC0(29-2)+IC0(30-2)))^0.0739;
Kout_BMD = 0.000146*(SCTX/IC0(18-2))^0.0779;    

%% (2) Define nonlinear ODEs

nm = length(IC0) + length(IC_BMD0);
dydt=zeros(nm,1);
dydt(1)  = [kin01 - kout0104*y(1)];
dydt(2)  = [kin02 - kout02*y(2)];
dydt(3)  = [kin03 - kb0305*y(3)];
dydt(4)  = [J14 - (kout0412+kout040u)*y(4) + kout0104*y(1)];
dydt(5)  = [kin05 - kout0514*y(4) - (kout050u+kb0508)*y(5) + kb0305*y(3) + kb0805*y(8)];
dydt(6)  = [y(9) - T69*y(6)];
dydt(7)  = [SPTH - kout*y(7)];
dydt(8)  = [kb0508*y(5) - kb0805*y(8)];
dydt(9)  = [SE - T64*y(9)];
dydt(10) = [kin10 - kout10*y(10)];
dydt(11) = [PTin  - PTout *y(11)];
dydt(12) = [kout0412*y(4) - J14 + k14a*y(13) - k15a*y(12)];
dydt(13) = [k15a*y(12) - k14a*y(13)];
dydt(14) = [kHApIn2  - kLShap*y(14)];                        % [v0514-v1405+f1514-f1415];
dydt(16-1) = [ROBin - KPT*y(16-1)];
dydt(18-2) = [kinOC2 - KLSoc*y(18-2)];                                                            % osteoclast
dydt(19-2) = [kin19 - koutTGFeqn*y(19-2)];                                                        % latent TGF
dydt(20-2) = [koutTGFeqn*y(19-2) - koutTGFact*y(20-2)];                                           % active TGF
dydt(21-2) = [kin21 - koutRNK*y(21-2) - k3*y(21-2)*y(22-2)  + k4*y(24-2)];                        % RANK
dydt(22-2) = [kinL - (koutL+k1*y(23-2)+k3*y(21-2)+drug)*y(22-2) + k2*y(25-2) + k4*y(24-2)];       % RANKL
dydt(23-2) = [pO - (kO+k1*y(22-2))*y(23-2) + k2*y(25-2)];                                         % OPG             
dydt(24-2) = [k3*y(21-2)*y(22-2) - k4*y(24-2)];                                                   % RANK-RANKL complex
dydt(25-2) = [k1*y(23-2)*y(22-2) - k2*y(25-2)];                                                   % OPG-RANKL complex
dydt(26-2) = [RX2Kin - RX2Kout*y(26-2)];
dydt(27-2) = [crebKin - crebKout*y(27-2)];
dydt(28-2) = [bcl2Kin - bcl2Kout*y(28-2)];
dydt(29-2) = [kb29*y(16-1) - kbfast*y(29-2)];
dydt(30-2) = [kb30*y(16-1) - kbslow*y(30-2)];
dydt(31-2) = [-ka*y(31-2)];
dydt(32-2) = [ka*y(31-2)-CLtot*C-Q*(C-y(33-2)/V2)];                % total drug = free drug (C) + complex (RC)
dydt(33-2) = [Q*(C-y(33-2)/V2)];
dydt(34-2) = [ksyn - (kdeg+mic34)*y(34-2)];    % ng/mL, total RANKL = free RANKL (R) + complex (RC)
dydt(35-2) = [Kin_BMD - Kout_BMD*y(35-2)];
