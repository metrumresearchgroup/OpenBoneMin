
% Function for getting inductively-calculated K0 and K (rate constants vector/matrix)
% for bone system

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%  Lumped version of
%    Peterson MC, Riggs MM (2010) Bone 46:49-63
%                        +
%    Peterson MC, Riggs MM (2012) CPT Pharmacometrics Syst Pharmacol 1:e14
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

%% 

function [K0,K] = Kmat_induc_PD(y0,nm,IC0,C)

% dydt(1)  = [kin01 - kout0104*y(1)];   % y(5)_0 = iss for all time points
% dydt(2)  = [kin02 - kout02*y(2)];
% dydt(3)  = [kin03 - kb0305*y(3)];
% dydt(4)  = [J14 - (kout0412+kout040u)*y(4) + kout0104*y(1)];
% dydt(5)  = [kin05 - kout0514*y(4) - (kout050u+kout0508)*y(5) + kb0305*y(3) + kb0805*y(8)];
% dydt(6)  = [y(9) - T69*y(6)];
% dydt(7)  = [SPTH - kout*y(7)];
% dydt(8)  = [kout0508*y(5) - kb0805*y(8)];
% dydt(9)  = [SE - T64*y(9)];
% dydt(10) = [kin10 - kout10*y(10)];
% dydt(11) = [PTin - PTout*y(11)];
% dydt(12) = [kout0412*y(4) - J14 + k14a*y(13) - k15a*y(12)];
% dydt(13) = [k15a*y(12) - k14a*y(13)];
% dydt(14) = [kHApIn2  - kLShap*y(14)];                        % [v0514-v1405+f1514-f1415];
% dydt(15) = [0];                                              % [f1415-f1514];
% dydt(16-1) = [ROBin - KPT*y(16)];
% dydt(17) = [0];
% dydt(18-2) = [kinOC2 - KLSoc*y(18-2)];                                                                   % osteoclast
% dydt(19-2) = [kin19 - koutTGFeqn*y(19-2)];                                                               % latent TGF
% dydt(20-2) = [koutTGFeqn*y(19-2) - koutTGFact*y(20-2)];                                                    % active TGF
% dydt(21-2) = [kin21 - koutRNK*y(21-2) - k3*y(21-2)*y(22-2)  + k4*y(24-2)];                                     % RANK
% dydt(22-2) = [kinL - (koutL+k1*y(23-2)+k3*y(21-2)+drug)*y(22-2) + k2*y(25-2) + k4*y(24-2)];       % RANKL
% dydt(23-2) = [pO - kO*y(23-2) - k1*y(23-2)*y(22-2) + k2*y(25-2)];                                              % OPG             
% dydt(24-2) = [k3*y(21-2)*y(22-2) - k4*y(24-2)];                                                              % RANK-RANKL complex
% dydt(25-2) = [k1*y(23-2)*y(22-2) - k2*y(25-2)];                                                              % OPG-RANKL complex
% dydt(26-2) = [RX2Kin - RX2Kout*y(26-2)];
% dydt(27-2) = [crebKin - crebKout*y(27-2)];
% dydt(28-2) = [bcl2Kin - bcl2Kout*y(28-2)];
% dydt(29-2) = [kb29*y(16) - kbfast*y(29-2)];
% dydt(30-2) = [kb30*y(16) - kbslow*y(30-2)];

%% hyperbolic functions used for K-matrix

Model_parameter_values0

K0 = zeros(nm,1);
K  = zeros(nm,nm);

%----------------------------------------------------------------------
% kin01
%----------------------------------------------------------------------
% no zero order

T85Rpart = y0(2)^T80/(y0(2)^T80 + T81^T80) ;
% H+0201 in dydt(1)
T85 = T77*T85Rpart;
F11 = T85;
% kb0201 = OralCa*F11;
K0(1) = OralCa*F11;

%----------------------------------------------------------------------
% kout01
%----------------------------------------------------------------------
T29 = (T28*T0 - 0.17533*T0)/0.17533;           % d0104
% H0104 in dydt(1)
T31 = T28*y0(1)/(y0(1)+T29);
T83 = y0(2)/IC0(2);
kb0104 = T31*T83/(y0(1)+T81) + T87;

%----------------------------------------------------------------------
% kin02
%----------------------------------------------------------------------
C8 = y0(6)/Vva;
% H+0602 in dydt(2)
% T36 = T33;               %  + (T34-T33)*(C8^CaPOgam/(T35^CaPOgam+ C8^CaPOgam));
T36_2 = T33 + (T34-T33)*(C8^CaPOgam/(T35^CaPOgam+ C8^CaPOgam));
% H-0602 in dydt(2)
T37 = T34 - (T34-T33)*(C8^CaPOgam/(T35^CaPOgam+ C8^CaPOgam));

% dydt(2) = T36*(1-y2) - T37*y2 = "T36 - (T36+T37)*y2"           !!!!!!!!!!!!!!!!!!!!!!
K0(2) = T36_2;

%----------------------------------------------------------------------
% kout02
%----------------------------------------------------------------------
kout02 = T36_2 + T37;

% kb0602 = (T34-T33)/Vva*(C8^(CaPOgam-1)/(T35^CaPOgam+ C8^CaPOgam));   

%----------------------------------------------------------------------
% kin03
%----------------------------------------------------------------------
K0(3) = OralPhos*F12;

%----------------------------------------------------------------------
% kout03
%----------------------------------------------------------------------
% kb0305

%----------------------------------------------------------------------
% kin04
%----------------------------------------------------------------------
k0412 = CaDay/24;
T13 = k0412/IC0(12);             % typical kout12
J14OC50= exp(log((J14OCmax*IC0(18-2)^J14OCgam/T13) - IC0(18-2)^J14OCgam)/J14OCgam);            % d181204 in the paper
Osteoclast =  y0(18-2);
% H181204 in dydt(4)
OCeqn = (J14OCmax*Osteoclast^J14OCgam)/(Osteoclast^J14OCgam + J14OC50^J14OCgam);

% last product in v1204
MOCratio = y0(24-2)/Osteoclast;
MOCratio0 = IC0(24-2)/IC0(18-2);
MOCratioEff = (MOCratio/MOCratio0)^MOCratioGam       ;

J14OCdepend = IC0(12) * FracJ14 * OCeqn * MOCratioEff;   

% v1204 in the paper
J14 = k0412*(1-FracJ14) + J14OCdepend;
K0(4) = J14;

%----------------------------------------------------------------------
% kout04
%----------------------------------------------------------------------
T15 = k0412/IC0(4);              % typical kout04
% v0412 in the paper, transit to bone Ca depends on HAp amount (y14)
J15 = T15*(1-FracJ15) + T15*FracJ15*y0(14);     % product with y(4) in ODE
kout0412 = J15;

% for v040u, 2-H0604
T10 = 2 - T7*C8/(C8+T9);

% for v040u, 0.3*GFR*A(4)/Vva
CaFilt = 0.6*0.5*GFR/Vva;                       % divide by Vva so A can be used with coefficient

C1 = y0(4)/Vva ;
C1_0 = IC0(4)/Vva;
% for v040u, H040u
ReabsMax = (0.3*GFR*C1_0 - 0.149997)*(Reabs50 + C1_0) / C1_0;
% kout in H040u in dydt(4)
k040u = ReabsMax/(Reabs50 + C1)/Vva;            % divide by Vva so A can be used with coefficient

% for v040u, H07040u
C4 = y0(7)/Vva;
C7_0 = IC0(7)/Vva;
T17 = C7_0*T16 - C7_0;           % d07040u in the paper
% H07040u in dydt(4)
ReabsPTHeff = T16*C4/(C4+T17) ;

CaReabsActive = k040u*ReabsPTHeff ;
% 2nd parenthesis
T20 = CaFilt - CaReabsActive;             % coefficient!!!!!!!!!!!!!!!!!

% v040u (without y4) in dydt(4)
J27a = T10*T20;                       % final coefficient!!!!!!!!!!!!!!!!!
if J27a*y0(4)<0
    J27 = 0;
else
    J27 = J27a;
end
kout040u = J27;

%----------------------------------------------------------------------
% kin05
%----------------------------------------------------------------------
J41 = 0.464*J14  ;               % v1405 in the paper
kin1405 = J41;

% a part of v050u without y5 is included to kin05 (zero-order)
if y0(5)/Vva < T46
    phi050u = y0(5)/Vva;
else
    phi050u = T46;
end
T47 = 0.88*GFR*phi050u;
kin050u = T47;

K0(5) = kin1405 + kin050u;       % - v050u;   %  -kout0514*y0(4);

%----------------------------------------------------------------------
% kout05
%----------------------------------------------------------------------
J42 = 0.464*kout0412 ;                % kout in v0514 in the paper
kout0514 = J42;

kout050u = 0.88*GFR/Vva;

kb0305 = T52;
% J53 = T52*y(3);

kb0508 = T49/Vva;

kb0805 = T55;
% J56 = T55*y(8);

%----------------------------------------------------------------------
% kin06
%----------------------------------------------------------------------
% no zero order

% K0(6) = y0(9);

%----------------------------------------------------------------------
% kin07
%----------------------------------------------------------------------
% no zero order

CaConc0 = IC0(4)/Vva;
CaConc = y0(4)/Vva;

INparenCa = (T58 - T61) * CaConc0^T59 / (T58 - C7_0*100) - CaConc0^T59  ;
T60 = exp(log(INparenCa) / T59)  ;
% H-04107 in dydt(7)
T63 =  T58 - (T58 - T61) * (CaConc)^T59 / ((CaConc)^T59 + T60^T59);

FCTD = (y0(10) / IC0(10));           %  * y0(11);

EPTH = T63 * FCTD;
% IPTH = 0.693*SC +  IPTHinf;
IPTH = 0;
SPTH = EPTH + IPTH ;

kb1107 = SPTH;
% K0(7) = SPTH;

%----------------------------------------------------------------------
% kout07
%----------------------------------------------------------------------
kout = T57/14;                       % k7D in the paper

%----------------------------------------------------------------------
% kin09
%----------------------------------------------------------------------
% no zero order

C2 = y0(5)/Vva ;
C2_0 = IC0(5)/Vva;
PO4inhPTH = (C2/C2_0)^PO4inhPTHgam  ;
T66 = (T67^AlphOHgam + C7_0^ AlphOHgam )/C7_0^ AlphOHgam ;
% H0709 in dydt(9)
T68 = T66*C4^(AlphOHgam)/(T67^(AlphOHgam)*PO4inhPTH+C4^AlphOHgam) ;        % T66/Vva

PhosEffTop = (PhosEff0 - 1)*( C2_0^PhosEffGam + PhosEff50^PhosEffGam ) ;
PhosEffBot = PhosEff0 * C2_0^PhosEffGam ;
PhosEffMax =  PhosEffTop / PhosEffBot;
% H-0509 in dydt(9)
PhosEff = PhosEff0 - (PhosEffMax*PhosEff0 * C2^PhosEffGam /(C2^ PhosEffGam  + PhosEff50^PhosEffGam));
if y0(5) >= IC0(5)
        PhosEffect = PhosEff;
else
        PhosEffect = 1;
end

SE = T65*T68*PhosEffect;

% kb0709 = SE;
K0(9) = SE;

%----------------------------------------------------------------------
% kout09
%----------------------------------------------------------------------
% T64

%----------------------------------------------------------------------
% kin10
%----------------------------------------------------------------------
ScaEff =  (CaConc0/CaConc)^ ScaEffGam ;

C8_0 = IC0(6)/Vva;
T72 = C8_0 * ScaEff  ;                    % C8_0 = dT0604 in the paper = 90

T73 = T71 * (C8 - T72) ;

T74 = (exp(T73) - exp(-T73)) / (exp(T73) + exp(-T73));

% T+0604 in dydt(10)
Tp0604 = 1 + T74;
T75 = T70 * (0.85 * Tp0604 + 0.15);
% T-0604 in dydt(10)
Tm0604 = 1 - T74;
T76 = T70 * (0.85 * Tm0604 + 0.15);

% dydt(10) = (1-y10)*T76 - y10*T75 = "T76 - (T76+T75)*y10"           !!!!!!!!!!!!!!!!!!!!!!
K0(10) = T76;

%----------------------------------------------------------------------
% kout10
%----------------------------------------------------------------------
kout10 = T76+T75;
     
%----------------------------------------------------------------------
% kin11
%----------------------------------------------------------------------
% use C8_0 rather than C8!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INparenCtriol =((CtriolMax - CtriolMin) * C8_0^CtriolPTgam) / (CtriolMax - 1) - C8_0^CtriolPTgam ;
Ctriol50 = exp(log(INparenCtriol) / CtriolPTgam) ;               % ok
% H-0611
CtriolPTeff = CtriolMax - (CtriolMax - CtriolMin) * C8^CtriolPTgam / (C8^CtriolPTgam + Ctriol50^CtriolPTgam) ;

K0(11) = PTout * CtriolPTeff;
% K0(11) = PTout * CtriolMax;
% kb0611 = -PTout * ((CtriolMax - CtriolMin)/Vva) * C8^(CtriolPTgam-1) / (C8^CtriolPTgam + Ctriol50^CtriolPTgam) ;

%----------------------------------------------------------------------
% kout11
%----------------------------------------------------------------------
% PTout

%----------------------------------------------------------------------
% kin12
%----------------------------------------------------------------------
K0(12) = -J14 + k14a*y0(13);       % + kout0412*y0(4) + k14a*y0(13);

%----------------------------------------------------------------------
% Ca transit between IC (12) and non-IC (13)
%----------------------------------------------------------------------
k15a = k14a*IC0(13)/IC0(12);                    % k1213 in dydt(12), ok

%----------------------------------------------------------------------
% kin14
%----------------------------------------------------------------------
% no zero-order 

%----------------------------------------------------------------------
% kout14
%----------------------------------------------------------------------
kLShap = 1/HApMRT;

kHApIn = kLShap/(IC0(29-2)+IC0(30-2));

Osteoblast = y0(29-2) + y0(30-2);
% K0(14) = kHApIn*Osteoblast;

%----------------------------------------------------------------------
% kin16
%----------------------------------------------------------------------
kin16 = kb*(IC0(29-2)+IC0(30-2));              % typical kin16

Dr = kin16/Pic0  	  ;          % for the 1st term in dydt(16-1)

E0PicROB = FracPicROB*Pic0;                % r2016 in H+2016
EC50PicROBparen= (EmaxPicROB*IC0(20-2)^PicROBgam / (Pic0 - E0PicROB)) - IC0(20-2)^PicROBgam;
EC50PicROB = exp(log(EC50PicROBparen)/PicROBgam);       
% H+2016 in dydt(16-1)
% PicROB = E0PicROB + EmaxPicROB*y0(20-2)^PicROBgam/(y0(20-2)^PicROBgam + EC50PicROB^PicROBgam);

% ROBin2 = Dr*PicROB;
% K0(16-1) = ROBin2;
K0(16-1) = Dr*E0PicROB;      % + Dr*EmaxPicROB*y0(20-2)^(PicROBgam-1)/(y0(20-2)^PicROBgam + EC50PicROB^PicROBgam);

% kb2016 = Dr*EmaxPicROB*y0(20-2)^(PicROBgam-1)/(y0(20-2)^PicROBgam + EC50PicROB^PicROBgam);      % PicROBgam-1
kb2016 = Dr*EmaxPicROB*1/(y0(20-2)^PicROBgam + EC50PicROB^PicROBgam);

%----------------------------------------------------------------------
% kout16 (also used in dydt(29-2) (17a)
%----------------------------------------------------------------------
bigDb = kb*(IC0(29-2)+IC0(30-2))*Pic0/IC0(16-1);     % for the 1st term in dydt(29-2) (17a)

E0PicOB = FracPicOB*Pic0;
EC50PicOBparen = (EmaxPicOB*IC0(20-2)^PicOBgam/(Pic0 - E0PicOB)) - IC0(20-2)^PicOBgam;
EC50PicOB = exp(log(EC50PicOBparen)/PicOBgam);
% H+2017 in dydt(16-1)
PicOB = E0PicOB + EmaxPicOB*y0(20-2)^PicOBgam / (y0(20-2)^PicOBgam + EC50PicOB^PicOBgam);

KPT =1*(bigDb/PicOB)  ;

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
PicOBkb = E0PicOBkb - (E0PicOBkb  - EmaxPicOBkb)*y0(20-2)^PicOBgamkb / (y0(20-2)^PicOBgamkb + EC50PicOBkb^PicOBgamkb);

PicOBkbEff = PicOBkb/Pic0;

E0RUNX2kbEff= E0RUNX2kbEffFACT*kb;
RUNkbMax = E0RUNX2kbEff*RUNkbMaxFact;          % a2817D
INparen = (RUNkbMax * RUNX20^RUNkbGAM) / (E0RUNX2kbEff - kb) - RUNX20^RUNkbGAM;
RUNkb50 = exp(log(INparen)/RUNkbGAM);
if y0(28-2) > 105
    RUNX2 = y0(28-2) - 90;
else
    RUNX2 = 10;
end
% H2817D in dydt(29-2) (17a)
RUNX2kbPrimeEff = RUNkbMax*RUNX2^RUNkbGAM / (RUNX2^RUNkbGAM + RUNkb50^RUNkbGAM);

kbprime = E0RUNX2kbEff*PicOBkbEff - RUNX2kbPrimeEff;

kbslow = kbprime*Frackb;

kb30 = KPT*(1-FracOBfast)*Frackb;
% K0(30-2) = KPT*(1-FracOBfast)*Frackb*y0(16-1);

%----------------------------------------------------------------------
% kin29 (17a)
%----------------------------------------------------------------------
% No zero-order input

%----------------------------------------------------------------------
% kout29 (17a)
%----------------------------------------------------------------------
kbfast = (kb*(IC0(29-2)+IC0(30-2)) + kbslow*IC0(29-2)*FracOBfast - kbslow*(IC0(29-2)+IC0(30-2))) / IC0(29-2) ;          % add "*FracOBfast"!!!!!!!!!!!!!

Frackb2 = kbfast/kbprime;

kb29 = KPT*FracOBfast*Frackb2;
% K0(29-2) = KPT*FracOBfast*Frackb2*y0(16-1);

%----------------------------------------------------------------------
% kin18
%----------------------------------------------------------------------
kin18 = Da*IC0(18-2);             % typical kin18

PicOCkin = Pic0;

% EC50MeffOC = exp(log(IC0(24-2)^kinOCgam*EmaxMeffOC/(1-E0Meff) - IC0(24-2)^kinOCgam)/kinOCgam);
% H+2418S in dydt(18-2)
% MeffOC = E0Meff + (EmaxMeffOC * y0(24-2)^kinOCgam/(y0(24-2)^kinOCgam + EC50MeffOC^kinOCgam));

K0(18-2) = kin18*PicOCkin*E0Meff;                % MeffOC;

% gamma>8, range of y(24-2) is less than EC50 == 0.0002, so effect is set as zero
% kb2418 =   kin18*PicOCkin*(EmaxMeffOC * y0(24-2)^(kinOCgam-1)/(y0(24-2)^kinOCgam + EC50MeffOC^kinOCgam));

%----------------------------------------------------------------------
% kout18
%----------------------------------------------------------------------
E0PicOC = FracPicOC*Pic0;
EC50PicOCparen = (EmaxPicOC*IC0(20-2)^PicOCgam/(Pic0 - E0PicOC)) - IC0(20-2)^PicOCgam;
EC50PicOC = exp(log(EC50PicOCparen)/PicOCgam);              % ok
% H+2018D in dydt(18-2)
PicOC = E0PicOC + ((EmaxPicOC*y0(20-2)^PicOCgam)/(y0(20-2)^PicOCgam + EC50PicOC^PicOCgam));

% use y(24-2) rather than y(22-2)
PiL0 = (k3/k4)*IC0(22-2);                % IC0(24-2)
PiL = y0(24-2)/10;
EC50survInPar = (E0RANKL - EmaxL)*(PiL0^LsurvOCgam/(E0RANKL - 1)) - PiL0^LsurvOCgam;
EC50surv = exp(log(EC50survInPar)/LsurvOCgam);
% H-2218D in dydt(18-2)
% LsurvOC = E0RANKL - (E0RANKL - EmaxL)*(PiL^LsurvOCgam/(PiL^LsurvOCgam + EC50surv^LsurvOCgam));
   
KLSoc = Da*PicOC*E0RANKL;              % #*3

kb2418 = Da*PicOC*(E0RANKL - EmaxL)*(1/10)/(PiL^LsurvOCgam + EC50surv^LsurvOCgam)*y0(18-2)*PiL^(LsurvOCgam-1);

%----------------------------------------------------------------------
% kin19
%----------------------------------------------------------------------
kinTGF = koutTGF0*IC0(19-2);

K0(19-2) = kinTGF*(Osteoblast/(IC0(29-2)+IC0(30-2)))^OBtgfGAM;

%----------------------------------------------------------------------
% kout19
%----------------------------------------------------------------------
koutTGF = koutTGF0;
ratio19 = ((y0(19-2)/IC0(19-2))^koutTGFGam);
koutTGFeqn = koutTGF * ratio19 * ((Osteoclast/IC0(18-2))^OCtgfGAM);

%----------------------------------------------------------------------
% kin20
%----------------------------------------------------------------------
% koutTGFeqn*y(19-2)

koutTGFeqn2 = koutTGF * ratio19 / (IC0(18-2)^OCtgfGAM) * y0(19-2) * y0(18-2)^(OCtgfGAM-1);

%----------------------------------------------------------------------
% kout20
%----------------------------------------------------------------------
koutTGFact = 1000*koutTGF0;

%----------------------------------------------------------------------
% kin21 (zero-order)
%----------------------------------------------------------------------
kinRNK = (koutRNK*IC0(21-2) + k3*IC0(21-2)*IC0(22-2) - k4*IC0(24-2)) / IC0(20-2)^kinRNKgam ;
K0(21-2) = kinRNK*y0(20-2)^kinRNKgam;           %  + k4*y0(24-2);

%----------------------------------------------------------------------
% kout21 (first-order)
%----------------------------------------------------------------------
% koutRNK

% kb2021 = kinRNK*y0(20-2)^(kinRNKgam-1);

%----------------------------------------------------------------------
% kin22
%----------------------------------------------------------------------
kinLbase = koutL*IC0(22-2);                % typical kin22
   
OsteoEffect = (Osteoblast/(IC0(29-2)+IC0(30-2)))^OsteoEffectGam ;
OsteoEffect0 = Osteoblast^(OsteoEffectGam-1)/((IC0(29-2)+IC0(30-2))^OsteoEffectGam) ;

PTH50 = EmaxLpth*C7_0 - C7_0 ;
PTHconc = C4;
LpthEff = EmaxLpth*(PTHconc) / ((PTH50*OsteoEffect^TESTPOWER) + (PTHconc))  ;
% LpthEff = EmaxLpth*(1) / ((PTH50*OsteoEffect^TESTPOWER) + (PTHconc)) /Vva ;

% K0(22-2) = kinLbase*OsteoEffect*LpthEff;
% kb0722 = kinLbase*OsteoEffect*LpthEff;
kb1822 = kinLbase*OsteoEffect0*LpthEff;

%----------------------------------------------------------------------
% kout22 (binding to RANK (21-2) and OPG (23-2) should also be considered
%----------------------------------------------------------------------
% koutL

%----------------------------------------------------------------------
% kin23
%----------------------------------------------------------------------
pObase = kO*IC0(23-2);                     % typical kin23

D = y0(16-1);
K0(23-2) = pObase*(D/IC0(16-1))*((PTHconc+(opgPTH50*(D/IC0(16-1))))/(2*PTHconc))+ IO;   %  + k2*y0(25-2);

%----------------------------------------------------------------------
% kout23 (binding to RANKL (22-2) should also be considered
%----------------------------------------------------------------------
% kO

% kb1623 = pObase*(1/IC0(16-1))*((PTHconc+(opgPTH50*(D/IC0(16-1))))/(2*PTHconc))+ IO;

%----------------------------------------------------------------------
% kin24
%----------------------------------------------------------------------
% K0(24-2) = k3*y0(21-2)*y0(22-2);

%----------------------------------------------------------------------
% kin25
%----------------------------------------------------------------------
% K0(25-2) = k1*y0(23-2)*y0(22-2);

%----------------------------------------------------------------------
% kin26
%----------------------------------------------------------------------
K0(26-2) = RX2Kout0*IC0(26-2);                    % k26S in the paper

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

K0(27-2) = crebKin0 * H0727S;

% kb0727 =   crebKin0 * EmaxPTHcreb*(1/Vva)/(PTHconc+EC50PTHcreb);

%----------------------------------------------------------------------
% kout27
%----------------------------------------------------------------------
% crebKout

%----------------------------------------------------------------------
% kin28
%----------------------------------------------------------------------
% K0(28-2) = y0(26-2)*y0(27-2)*bcl2Kout;

%----------------------------------------------------------------------
% kout28
%----------------------------------------------------------------------
% bcl2Kout

%----------------------------------------------------------------------
% denosumab PK (y32 = total, C = free)
%----------------------------------------------------------------------
drug = (kint-koutL)*C/(Kss+C);

%% Define the parameter matrix ('K' matrix) of the ODEs in the original model

K(1,1) = -kb0104;
% K(1,2) = kb0201;

K(2,2) = -kout02;
% K(2,6) = kb0602;

K(3,3) = -kb0305;

K(4,1) = kb0104;
K(4,4) = -(kout0412+kout040u);

K(5,3) = kb0305;
K(5,4) = -kout0514;
K(5,5) = -(kout050u+kb0508);            % -kb0508;       % 
K(5,8) = kb0805;

K(6,6) = -T69;
K(6,9) = 1;

K(7,7) = -kout;
K(7,11) = kb1107;

K(8,5) = kb0508;
K(8,8) = -kb0805;

% K(9,7) = kb0709;
K(9,9) = -T64;

K(10,10) = -kout10;

% K(11,6) = kb0611;
K(11,11) = -PTout;

K(12,4)  = kout0412;
K(12,12) = -k15a;
% K(12,13) =  k14a;                    %%%%%%%%%%%%%%%%

K(13,12) =  k15a;
K(13,13) = -k14a;

K(14,14) = -kLShap;
K(14,29-2) = kHApIn;
K(14,30-2) = kHApIn;

K(16-1,16-1) = -KPT;
K(16-1,20-2) = kb2016;

K(18-2,18-2) = -KLSoc;
% K(18-2,20-2) = -kb2018;
K(18-2,24-2) = kb2418;

K(19-2,19-2) = -koutTGFeqn;

K(20-2,18-2) =  koutTGFeqn2;
K(20-2,20-2) = -koutTGFact;

% K(21-2,20-2) = kb2021;
K(21-2,21-2) = -(koutRNK + k3*y0(22-2));
% K(21-2,22-2) = -k3*y0(21-2);
K(21-2,24-2) = k4;                          % from complex

% K(22-2,7)    = kb0722;
% K(22-2,21-2) = -k3*y0(22-2);
K(22-2,29-2) = kb1822;
K(22-2,30-2) = kb1822;
K(22-2,22-2) = -(koutL + k3*y0(21-2) + k1*y0(23-2) + drug);      
% K(22-2,23-3) = -k1*y0(22-2);
K(22-2,24-2) = k4;
K(22-2,25-2) = k2;

% K(23-2,16-1) = kb1623;
% K(23-2,22-2) = -k1*y0(23-2);
K(23-2,23-2) = -(kO + k1*y0(22-2));
K(23-2,25-2) = k2;                          % from complex

% K(24-2,21-2) = k3*y0(22-2);
K(24-2,22-2) = k3*y0(21-2);
K(24-2,24-2) = -k4;

K(25-2,22-2) = k1*y0(23-2);
% K(25-2,23-2) = k1*y0(22-2);
K(25-2,25-2) = -k2;

K(26-2,26-2) = -RX2Kout;

% K(27-2,7)    = kb0727;
K(27-2,27-2) = -crebKout;

K(28-2,27-2) = y0(26-2)*bcl2Kout;
K(28-2,28-2) = -bcl2Kout;

K(29-2,16-1) = kb29;
K(29-2,29-2) = -kbfast;

K(30-2,16-1) = kb30;
K(30-2,30-2) = -kbslow;
