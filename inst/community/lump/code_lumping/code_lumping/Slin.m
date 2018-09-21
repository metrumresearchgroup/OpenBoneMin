
% Setting for linearization

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%  Lumped version of
%    Peterson MC, Riggs MM (2010) Bone 46:49-63
%                        +
%    Peterson MC, Riggs MM (2012) CPT Pharmacometrics Syst Pharmacol 1:e14
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

%% 

% Time range of interest
global incr_t;
year = 4;
maxt = 52*7*year;          % days
incr_t = 0.5;              % days

% Parameter values
Model_parameter_values0 

% Initial condition
IC_setting
IC0 = IC;
IC_BMD = 100;
IC_BMD0 = IC_BMD;

% Dosing information
numd = 2*year;
tDosed = [];
for i=0:numd-1
    tDosed = [tDosed;52*7*(i/2)];    %#ok<AGROW> % dosing stop at year-0.5
end
Dose0 = 60;
Dose = F1*Dose0*1000*ones(length(tDosed),1); % ug
nd = length(Dose);

% For inductive linearization
global max_n incr_ME npk;
max_n=80;       % Maximum number of iterations for inductive approximation
incr_ME = incr_t*2;     % step size in matrix exponential solution
cri_thre = 0.001;

%% End of the code
