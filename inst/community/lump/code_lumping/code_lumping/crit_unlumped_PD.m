
% Solve original nonlinear system

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%  Lumped version of
%    Peterson MC, Riggs MM (2010) Bone 46:49-63
%                        +
%    Peterson MC, Riggs MM (2012) CPT Pharmacometrics Syst Pharmacol 1:e14
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

%% 

ti = [tDosed
      maxt];          % add time after the last dose

IC_all = [IC;IC_BMD];
IC_all0 = IC_all;

nm = length(IC_all);

y = [];
c = [];
for dos = 1:nd
%      dos
     IC_all(31-2) = IC_all(31-2) + Dose(dos);
     TT0=ti(dos):incr_t:ti(dos+1); % d
     TT = (TT0-TT0(1))*24;         % h 
     
     options=odeset('RelTol',1e-5,'NonNegative',[1:nm]);
     [t,y1] = ode15s(@(t,y1) pkpdfun(t, y1, IC0, IC_BMD0), TT, IC_all, options);
     
     y = [y;y1];   %#ok<AGROW>
     c = [c;TT0']; %#ok<AGROW>   % d
     for m=1:nm
         IC_all(m)=y1(end,m);
     end
     if(dos<nd)
      y(end,:) = [];
      c(end) = [];
     end
     if dos==nd
         IC_all = IC_all0;           % IC returned to the original one
     end
end            

C = 0.5*((y(:,32-2)/V1 - y(:,34-2) - Kss) + sqrt((y(:,32-2)/V1 - y(:,34-2) - Kss).^2 + 4*Kss*y(:,32-2)/V1));

OB = y(:,29-2) + y(:,30-2);
OB0 = IC0(29-2) + IC0(30-2);
CFB_OB = OB/OB0 * 100;
OC = y(:,18-2);
CFB_OC = OC/IC0(18-2) * 100;

BMD = y(:,35-2);
CFB_BMD = (BMD-IC_all0(35-2))/IC_all0(35-2)*100;
