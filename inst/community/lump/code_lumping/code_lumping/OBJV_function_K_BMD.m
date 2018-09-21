
% Function for calculating SS for each lumped model

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%  Lumped version of
%    Peterson MC, Riggs MM (2010) Bone 46:49-63
%                        +
%    Peterson MC, Riggs MM (2012) CPT Pharmacometrics Syst Pharmacol 1:e14
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

%% 

function [ ARD,CFB_BMD_hat,pp,K0_L,K_L,y_hat,r_BMD ] = OBJV_function_K_BMD(M,K0_all,K_all,IC,TTpw,output,CFB_BMD,year)

global incr_t incr_ME;

y_hat = [];
ccc = [];
BMD = [];

IC_hat = M*IC;
nm = length(IC_hat);

K0_L = zeros(length(IC_hat),length(TTpw));         
K_L = zeros(length(IC_hat),length(IC_hat),length(TTpw));        

IC_BMD = 100;
IC_BMD0 = IC_BMD;

% choose the output state for the calculation of OBJV
% first need to find which row the output state is in now
output_row = zeros(length(output),1);
for i=1:length(output)
 output_row(i) = find(M(:,output(i))==1);
end

TT = TTpw(1):incr_ME:TTpw(end);
     for j=2:length(TT)
         
        tt = TT(j-1):incr_t:TT(j);    % interval = incr_ME
        tt2 = tt - TT(j-1);           % t0 converted to 0
        tt2 = tt2*24;                 % h
        
        b = tt';

        %------------------------
        % Solution for each state in the lumped model
        %------------------------
        if mod(TT(j-1),52*7/2) == 0 && TT(j-1) < 52*7*year
            a0 = TT(j-1);
        else
            a0 = (TT(j-1)+TT(j))/2;
        end
        a = find(abs(TTpw-a0)<0.0001);
        
        % K_hat at a particular time is calculated based on K at a particular time
        K0_hat = M*K0_all(:,a);
        K_hat = M*K_all(:,:,a)*pinv(M);
        
        if mod(TT(j-1),52*7/2) == 0 && TT(j-1) < 52*7*year
            for ppj=a:a+1
             K0_L(:,ppj) = K0_hat;
             K_L(:,:,ppj) = K_hat;
            end
        elseif j==length(TT)
            for ppj=a-1:a+1
             K0_L(:,ppj) = K0_hat;
             K_L(:,:,ppj) = K_hat;
            end
        else
            for ppj=a-1:a
             K0_L(:,ppj) = K0_hat;
             K_L(:,:,ppj) = K_hat;
            end
        end
        
        y_hat0 = ME_solution(K0_hat,K_hat,IC_hat,tt2);
        
        %------------------------
        % Solution for BMD from the lumped model
        %------------------------
        [r_BMD,~] = size(y_hat0);
        if mod(TT(j-1),52*7/2) == 0 && TT(j-1) < 52*7*year
            aa = 1;
        else
            aa = (1+r_BMD)/2;      % + tspan(1,1);
        end

        % fast and slow OB unlumped
        BSAP = y_hat0(aa,output_row(2)) + y_hat0(aa,output_row(3));
        SCTX = y_hat0(aa,output_row(1));
        
        K0_BMD = 0.000146*IC_BMD0*(BSAP/(IC(output(2))+IC(output(3)))).^0.0739;
        K_BMD = -0.000146*(SCTX/IC(output(1))).^0.0779;
        
        BMD0 = ME_solution(K0_BMD,K_BMD,IC_BMD,tt2);
        
        y_hat = [y_hat;y_hat0];   %#ok<AGROW>
        ccc = [ccc; b];           %#ok<AGROW>
        BMD = [BMD;BMD0];         %#ok<AGROW>
        if(j<length(TT))
         for m=1:nm
          IC_hat(m,1) = y_hat(end,m);
         end
         IC_BMD = BMD0(end);
         y_hat(end,:) = [];
         ccc(end,:) = [];
         BMD(end,:) = [];
        end
     end         % end of step size ME

pp = -9999;

CFB_BMD_hat = real((BMD-IC_BMD0)/IC_BMD0*100);
dif2  = (CFB_BMD(:,1) - CFB_BMD_hat(:,1)).^2;
ARD = sqrt(sum(dif2));
if isnan(ARD)
 ARD = 99999;
end
  
end
