
% Setting for lumping

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%  Lumped version of
%    Peterson MC, Riggs MM (2010) Bone 46:49-63
%                        +
%    Peterson MC, Riggs MM (2012) CPT Pharmacometrics Syst Pharmacol 1:e14
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

%% 

% Number of replications for VPC
nsim = 250;

% Indicate the output state (e.g. output = [2] meant the 2nd state for observation)
% See the file IC_setting.m for state #
output = [16;27;28];

% Constraint for lumping
constraint = [16     % state #16 remains unlumped
              15     % state #15 remains unlumped
              27     % state #27 remains unlumped
              28];   % state #28 remains unlumped
% e.g. if you additionally want to force states #1 to #3 to group together,
% constraint = [16 0  0     % state #16 remains unlumped
%               15 0  0     % state #15 remains unlumped
%               27 0  0     % state #27 remains unlumped
%               28 0  0     % state #28 remains unlumped
%               1  2  3];   % states #1 to 3 are forced to group together

% Increment for m (number of states)
incr = 1;

% Initial grouping for Simulated Annealing
% length = number of states in the original model
% Each number represents the state # in the lumped model
% Vx = [5 5 5 5 5 5 5 5 5 5 5 5 5 5 4 3 5 6 5 5 5 5 5 5 5 5 1 2]; % 6 states
Vx = [8 8 8 8 8 8 8 8 8 8 5 8 8 8 4 3 8 6 8 5 8 7 5 8 8 8 1 2]; % 8 states

mf = max(Vx);     % minimum number of states tested as a lumped model
ml = 8;          % maximum number of states tested as a lumped model

% Simulated Annealing (SA) Scenario
maxit = 5000*4/10;
num_dis = maxit/4;
temp_ini = 100000; % 1000, 10000, 100000 %temperature change here
decayRate = 0.999;

%% End of the code
