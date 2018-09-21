
% Function to solve original nonlinear system with uncertainty

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%  Lumped version of
%    Peterson MC, Riggs MM (2010) Bone 46:49-63
%                        +
%    Peterson MC, Riggs MM (2012) CPT Pharmacometrics Syst Pharmacol 1:e14
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

%%

function dydt = pkpdfun_SE(t,y,IC0,IC_BMD0,j,Param_SE)

nm = length(IC0) + length(IC_BMD0);

% Definition for parameters without precision (SE)
Model_parameter_values0

% Definition for parameters with precision (SE)
CL = Param_SE(1,j);
V1 = Param_SE(2,j);
Q  = Param_SE(3,j);
V2 = Param_SE(4,j);
ka = Param_SE(5,j);
Kss= Param_SE(6,j);
kint = Param_SE(7,j);
R0 = Param_SE(8,j);
kdeg = Param_SE(9,j);

ksyn = R0*kdeg;

OralCa = Param_SE(10,j);
EmaxLpth = Param_SE(12,j);
if EmaxLpth<1
    EmaxLpth=1;
end
EmaxPicROB = Param_SE(13,j);
EmaxPicOC = Param_SE(14,j);

def_ode
