
% Define parameter values

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%  Lumped version of
%    Peterson MC, Riggs MM (2010) Bone 46:49-63
%                        +
%    Peterson MC, Riggs MM (2012) CPT Pharmacometrics Syst Pharmacol 1:e14
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

%% 

k1 = 0.00000624;                               % k2124 in dydt(21)
k2 = 0.112013;                                 % k2421 in dydt(21)
k3 = 0.00000624;                               % k2124 in dydt(21)
k4 = 0.112013;                                 % k2421 in dydt(21)
kinRNKgam = 0.151825 ;                         % g2021 in dydt(21)
koutRNK = 0.00323667;                          % k21D in dydt(21)
Vva= 14;                                       % Vascular volume
CaDay = 88;                                    % k0412*24 in the paper, Ca exchanged between bone and plasma (mmol/day)
FracJ14 = 0.107763;                            % phi1204 in dydt(4)
J14OCmax = 0.543488;                           % a181204 in dydt(4)
J14OCgam = 1.6971;                             % g181204 in dydt(4)
FracJ15 = 0.114376 ;                           % phi0412 in dydt(4)
MOCratioGam = 0.603754;                        % g241204 in dydt(4)
Da = 0.7/24; 	                               % k18D in dydt(18)
OBtgfGAM = 0.0111319;                          % g1719S in dydt(19) 
koutTGF0 = 0.0000298449;                       % k1920 in dydt(19)
koutTGFGam = 0.919131;                         % g1920 in dydt(19)

OCtgfGAM = 0.593891;                           % g181920 in dydt(19)

EmaxPicROB = 3.9745;                           % a2016

PicROBgam = 1;                                 % g2016 in dydt(16)

FracPicROB = 0.883824;                         % used to calculate r2016
PicOBgam = 0.122313;                           % g2017 in dydt(16) and (17)
FracPicOB = 0.000244818;                       % used to calculate r2017
EmaxPicOB = 0.251636;                          % a2017 in dydt(16) and (17)
E0Meff = 0.388267;                             % r2418S in dydt(18)
EmaxMeffOC = 3.15667;                          % a2418S

kinOCgam = 8.53065;                            % g2418S in dydt(18)

EmaxPicOC = 1.9746;                            % a2018D
FracPicOC = 0.878215;                          % used to calculate r2018D

PicOCgam = 1.0168;                             % g2018D in dydt(18)

E0RANKL = 3.80338;                             % a2218D in dydt(18)
EmaxL = 0.469779;                              % r2218D in dydt(18)
GFR = 100/16.667;                              % 6.0
T16 = 1.06147;                                 % a07040u in the paper
T64 = 0.05;                                    % k090D in dydt(9)
T65 = 6.3;                                     % k090S in dydt(9)
T67 = 1.54865;                                 % d0709 in dydt(9)
AlphOHgam =  0.111241;                         % g0709 in dydt(9)
PhosEff0 = 1.52493;                            % a0509 in dydt(9)
PhosEff50 = 1.3021;                            % d0509 in dydt(9)
PhosEffGam = 8.25229;                          % g0509 in dydt(9)
PO4inhPTHgam = 0;                       
k14a =  0.0000244437;                          % k1312 in dydt(12)
HApMRT = 3.60609;
koutL = 0.00293273;                            % k22D in the paper
IO = 0;                              
RX2Kout0 = 0.693;                              % used to calculate k26S
E0rx2Kout = 0.125;                             % r0726D in dydt(26)
EmaxPTHRX2x = 5;                               % a0726D-r0726D in dydt(26)
E0crebKin = 0.5 ;                              % r0727S in dydt(27)
EmaxPTHcreb = 3.39745;                         % a0727S-r0727S in dydt(27)              
crebKout = 0.00279513;                         % k27D in dydt(27)
bcl2Kout = 0.693;                              % k28D in dydt(28)
ScaEffGam = 0.9;                               % g0410 in dydt(10)
T70 = 0.01;                                    % a10 in dydt(10)
T71 = 0.03;                                    % bT0604 in dydt(10)
T69 = 0.10;                                    % k6D in dydt(6)
Reabs50 = 1.57322;                             % d040u in dydt(4)
T7 = 2;                                        % a0604 in dydt(4)
T9 = 90;                                       % d0604 in dydt(4) = IC0(6)/Vva
T33 = 0.003;                                   % r0602 in dydt(2)
T34 = 0.037;                                   % a0602 in dydt(2)
T35 = 90;                                      % d0602 in dydt(2) = IC0(6)/Vva
CaPOgam = 1;                                   % g0602 in dydt(2)
T46 = 1.142;	                               % phi050u in dydt(5)
T52 =  0.365;                                  % k0305 in dydt(5)
T49 = 51.8;                                    % k0508 in dydt(5)
T55 = 0.019268;                                % k0805 in dydt(5)
OralPhos = 10.5/24;                            % D3 in dydt(3)
F12 = 0.7;                                     % F3 in dydt(3)
PicOBgamkb = 2.92375;                          % g2017D in dydt(30) (17b)
MultPicOBkb = 3.11842;                         % used to calculate a2017D
FracPic0kb = 0.764028;                         % used to calculate r2017D
E0RUNX2kbEffFACT = 1.01;                       % phi28k17D in the paper
RUNkbGAM = 3.67798;                            % g2817D
RUNkbMaxFact = 0.638114;                       % used to calculate a2817D
RUNX20 = 10;                                   % IC0(26)
Frackb = 0.313186;                             % phik17D in the paper
T81 = 0.75;                                    % d0201 in dydt(1)
T77 = 0.909359;                                % a0201 in dydt(1)
T80 = 4;                                       % g0201 in dydt(1)
T87 = 0.0495;                                  % k0104 in dydt(1)
T0 = 1.58471;                                  % IC0(1)
T28 = 0.9;                                     % a0104 in dydt(1)
OralCa = 24.055/24;                            % D1 in dydt(1)
CtriolPTgam = 12.5033;                         % g0611 in dydt(11)
CtriolMax = 4.1029;                            % a0611 in dydt(11)
CtriolMin = 0.9;                               % r0611 in dydt(11)
PTout = 0.0001604;                             % k0101 in dydt(11)
OsteoEffectGam = 0.173833;                     % g1722 in dydt(22)
EmaxLpth = 1.30721;                            % a0722 in dydt(22)
TESTPOWER = 1;                       
T57 = 100;                                     % used to calculate k7D
T58 = 6249.09;                                 % a04107 in dydt(7)
T59 = 11.7387;                                 % g04107 in dydt(7)
T61 = 96.25;                                   % r04107 in dydt(7)
IPTHint = 0;                                   % ok
IPTHinf = 0;                                   % ok
Pic0 = 0.228142;                               % ok
opgPTH50 = 3.85;                               % d0723 in dydt(23)
kO = 15.8885;                                  % k23D in dydt(23)
kb = 0.000605516;                              % k17D in dydt(29) (17a)
FracOBfast = 0.797629;                         % phi17a in dydt(29) (17a)

LsurvOCgam = 3.09023;                          % g2218D in dydt(18)

% for denosumab PK (Sutjandra et al, CPK, 2011)
CL = 3.06/1000;                        % L/h/66kg
V1 = 2490/1000;                        % L/66kg
Q = 37.9/1000;                         % L/h/66kg
V2 = 1360/1000;                        % L/66kg
ka = 0.212/24*(66/71.5)^-0.577;        % /h;
F1 = 0.638;
Kss = 138*10;                          % ng/mL
kint = 0.00795;                        % /h

R0 = 614;
kdeg = 0.00148;
ksyn = R0*kdeg;
