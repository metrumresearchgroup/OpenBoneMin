
% Function for solving ODEs using matrix exponential with Pade approximation

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%  Lumped version of
%    Peterson MC, Riggs MM (2010) Bone 46:49-63
%                        +
%    Peterson MC, Riggs MM (2012) CPT Pharmacometrics Syst Pharmacol 1:e14
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

%% 

function A = ME_solution(K0,K,IC,TT)  

% for liner ODEs: dy/dt=K*y --> transformed into matrix exponential: dy/dt=exp(t*K)*y(0)

A_last = zeros(length(K0),length(TT));      % state*time
for j=1:length(TT)
    ext = expmdemo1(K*TT(j));
    A1 = ext*IC;
    A2 = (ext-eye(size(ext))) / K * K0;       
    A_last(:,j) = A1 + A2;
end
A = A_last';
