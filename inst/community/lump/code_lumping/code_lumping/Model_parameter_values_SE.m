
% Parameter values with uncertainty in the original model

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%  Lumped version of
%    Peterson MC, Riggs MM (2010) Bone 46:49-63
%                        +
%    Peterson MC, Riggs MM (2012) CPT Pharmacometrics Syst Pharmacol 1:e14
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

%% 

% nsim = 250;

rng(1,'twister');
s = rng;
rng(s);

% for denosumab PK (Sutjandra et al, CPK, 2011)
CLall = exp(normrnd(log(3.06/1000),0.00987,[1 nsim]));       % L/h/66kg
V1all = exp(normrnd(log(2490/1000),0.01750,[1 nsim]));       % L/66kg
Qall  =  exp(normrnd(log(37.9/1000),0.00010,[1 nsim]));      % L/h/66kg
V2all = exp(normrnd(log(1360/1000),0.00010,[1 nsim]));       % L/66kg
kaall = exp(normrnd(log(0.212/24*(66/71.5)^-0.577),0.0701,[1 nsim]));          % /h;
% F1all = exp(normrnd(log(0.638),0.00010,[1 nsim]));
Kssall= exp(normrnd(log(138*10),0.00010,[1 nsim]));          % ng/mL
kintall = exp(normrnd(log(0.00795),0.000101,[1 nsim]));      % /h

R0all = exp(normrnd(log(614),0.0182,[1 nsim]));
kdegall = exp(normrnd(log(0.00148),0.00010,[1 nsim]));
% ksynall = R0*kdeg;

OralCa_all = exp(normrnd(log(24.055/24),0.5,[1 nsim]));
FracJ14_all = exp(normrnd(log(0.107763),4,[1 nsim]));
EmaxLpth_all = exp(normrnd(log(1.30721),3,[1 nsim]));
EmaxPicROB_all = exp(normrnd(log(3.9745),0.05,[1 nsim]));
EmaxPicOC_all = exp(normrnd(log(1.9746),0.10,[1 nsim]));

Param_SE = [CLall
            V1all
            Qall
            V2all
            kaall
            Kssall
            kintall
            R0all
            kdegall
            OralCa_all
            FracJ14_all
            EmaxLpth_all
            EmaxPicROB_all
            EmaxPicOC_all];
