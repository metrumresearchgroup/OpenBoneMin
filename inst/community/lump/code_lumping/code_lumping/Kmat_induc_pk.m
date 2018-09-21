
% Function for getting inductively-calculated K0 and K (rate constants vector/matrix)
% for denosumab PK

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%  Lumped version of
%    Peterson MC, Riggs MM (2010) Bone 46:49-63
%                        +
%    Peterson MC, Riggs MM (2012) CPT Pharmacometrics Syst Pharmacol 1:e14
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

%% 

function [K0,K] = Kmat_induc_pk(y0,nm)

% dydt(31-2) = [-ka*y(31-2)];
% dydt(32-2) = [ka*y(31-2)-CLtot*C-Q*(C-y(33-2)/V2)];                % total drug = free drug (C) + complex (RC)
% dydt(33-2) = [Q*(C-y(33-2)/V2)];
% dydt(34-2) = [ksyn - (kdeg+mic34)*y(34-2)];    % ng/mL, total RANKL = free RANKL (R) + complex (RC)

%% hyperbolic functions used for K-matrix

Model_parameter_values0

Crest = 0.5*((0        - y0(4) - Kss) + sqrt((y0(2)/V1 - y0(4) - Kss)^2 + 4*Kss*y0(2)/V1));
C     = 0.5*((y0(2)/V1 - y0(4) - Kss) + sqrt((y0(2)/V1 - y0(4) - Kss)^2 + 4*Kss*y0(2)/V1));

CLtot = CL + kint*V1*y0(4)/(Kss+C);                % y(34-2) should be ng/mL
mic34 = (kint-kdeg)*C/(Kss+C);

K0 = zeros(nm,1);
K0(2) = -CLtot*Crest - Q*Crest;
K0(3) = Q*Crest;
K0(4) = ksyn;

% Define the parameter matrix ('K' matrix) of the ODEs in the original model
K = [-ka                 0                    0          0
      ka     -(CLtot*0.5/V1 + Q*0.5/V1)      Q/V2        0
      0               Q*0.5/V1              -Q/V2        0
      0                  0                    0     -(kdeg + mic34)];
