
% Initial condition

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%  Lumped version of
%    Peterson MC, Riggs MM (2010) Bone 46:49-63
%                        +
%    Peterson MC, Riggs MM (2012) CPT Pharmacometrics Syst Pharmacol 1:e14
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

%% 

IC = zeros(nmpd+4,1);

IC(1) = 1.2375;    % 1, Gut Ca, 1.2375 in R, 1.585 in the paper
IC(2) = 0.5;       % 2, Gut Ca abs
IC(3) = 0.839;     % 3, Gut PO4
IC(4) = 32.9;      % 4, Vas Ca
IC(5) = 16.8;      % 5, Vas PO4
IC(6) = 1260;      % 6, Vas Calcitriol
IC(7) = 53.9;      % 7, Vas PTH
IC(8) = 3226;      % 8, Int PO4
IC(9) = 126;       % 9, 1-alpha-OH
IC(10)= 0.5;       % 10, PT pool
IC(11)= 1;         % 11, PT max
IC(12)= 100;       % 12, Ca IC
IC(13)= 24900;     % 13, Ca non-IC
IC(14)= 1;         % 14, HAp (not PO4 IC, assume same with Ca)
IC(16-1)= 0.00104122; % 16, responding OB
IC(18-2)= 0.001154;   % 18, OC
IC(19-2)= 228.1;      % 19, latent TGF beta
IC(20-2)= 0.2281;     % 20, active TGF beta
IC(21-2)= 10;         % 21, RANK
IC(22-2)= 0.4;        % 22, RANKL
IC(23-2)= 4;          % 23, OPG
IC(24-2)= k1*IC(21-2)*IC(22-2)/k2;      % 24, RANK-RANKL
IC(25-2)= k1*IC(23-2)*IC(22-2)/k2;      % 25, OPG-RANKL
IC(26-2)= 10;         % 26, Runx2
IC(27-2)= 10;         % 27, CREB
IC(28-2)= 100;        % 28, Bcl-2       % 100 in the paper, 10 in the original code
IC(29-2)= 0.00501324*FracOBfast;        % 29, OB fast (17a)  % 0.666757; 
IC(30-2)= 0.00501324*(1-FracOBfast);    % 30, OB slow (17b)         

% for denosumab PK
IC(31-2) = 0;         % abs
IC(32-2) = 0;         % central (total denosumab)
IC(33-2) = 0;         % peripheral
IC(34-2) = R0;        % ng/mL, free RANKL (R) + complex (RC)
