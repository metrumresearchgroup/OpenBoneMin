
% Function to solve original nonlinear system

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%  Lumped version of
%    Peterson MC, Riggs MM (2010) Bone 46:49-63
%                        +
%    Peterson MC, Riggs MM (2012) CPT Pharmacometrics Syst Pharmacol 1:e14
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

%% 

function dydt = pkpdfun(t,y,IC0,IC_BMD0)

% nm = length(IC0) + length(IC_BMD0);

Model_parameter_values0

def_ode
