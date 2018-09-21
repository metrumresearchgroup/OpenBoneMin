
% This is the run file

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%  Original model is based on:
%    Peterson MC, Riggs MM (2010) Bone 46:49-63
%                        +
%    Peterson MC, Riggs MM (2012) CPT Pharmacometrics Syst Pharmacol 1:e14
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

% Step 1. linearize the original nonlinear model
%  deliverable: linearized original model (y)

%                 dy/dt = K0_all + K_all * y

% Step 2. lump the linearized original model
%  deliverable : linearized lumped model (y_hat) & lumping matrix (L)

%                       y_hat = L * y   --->   y = L_inv * y_hat
%                                             (L_inv: pseudo inverse of L)
%                 dy_hat/dt = L * K0_all + (L * K_all * L_inv) * y

diary results.txt
clc
clear variables
close all

% Define each state
comexp_num0 = [];
for i=1:28
    if i<10
        comexp_num0 = [comexp_num0; '' num2str(i) '      '];      %#ok<AGROW>
    else
        comexp_num0 = [comexp_num0; '' num2str(i) '     '];       %#ok<AGROW>
    end
end
comexp_num = transpose(cellstr(comexp_num0));
nmpd = length(comexp_num);
comexp_char0 = ['Gut Ca ';'Gut abs';'Gut PO4';'Va Ca  ';'Va PO4 ';'Va Calt';'Va PTH ';
                'Int PO4';'1-a-OH ';'PT pool';'PT max ';'Ca IC  ';'Ca nIC ';'HAp    ';
                          'rOB    ';          'OC     ';'LA TGFb';'Ac TGFb';'RANK   ';
                'RANKL  ';'OPG    ';'RANK-L ';'OPG-L  ';'Runx2  ';'CREB   ';'Bcl-2  ';
                'OB fast';'OB slow'];
comexp_char = transpose(cellstr(comexp_char0));
comexp_original = [comexp_num;comexp_char];

Slin
Slum

%% 
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% 1. linearize the original nonlinear model
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

%% Prediction from original nonlinear model

disp(strcat('<Runtime to get prediction from original nonlinear model>'))
tic

crit_unlumped_PD

toc

fulltime = c;
CFB_full = zeros(length(fulltime),3);
CFB_full(:,1) = CFB_OB;
CFB_full(:,2) = CFB_OC;
CFB_full(:,3) = CFB_BMD;

end_ori = find(c==maxt);           % find(c==52*7*year);
cori = c(1:end_ori,1);               % time
Cori = C(1:end_ori,1);               % Conc

OBori = OB(1:end_ori,1);
CFB_OBori = CFB_OB(1:end_ori,1);
OCori = OC(1:end_ori,1);
CFB_OCori = CFB_OC(1:end_ori,1);

BMDori = BMD(1:end_ori,1);
CFB_BMDori = CFB_BMD(1:end_ori,1);

%% Inductive linearization for denosumab PK

% initial condition
ICpk = IC0(nmpd+1:end);
ICpk0 = ICpk;
nm = length(ICpk);

disp(strcat('<Runtime to get linearized PK model>'))
tic

ttt = [0:incr_t:maxt]';   % d
C_induc = zeros(length(ttt),max_n);

yy = zeros(nm,1);         % for y[n-1]

ARD = zeros(max_n,1);
ARDc = zeros(max_n,1);
n = 0;
cri = 9999;
while cri>=cri_thre && n<max_n
    n = n + 1;

    ICpk = IC0(nmpd+1:end);
    
    A_1_0 = [];             % Prediction from 4 PK-related states
    C_1_0 = [];             % Derived from 4 PK-related states
    cc = [];
    
    for dos = 1:nd
%      dos
     ICpk(1) = ICpk(1) + Dose(dos);
     TT=ti(dos):incr_ME:ti(dos+1);    % d
     for j=2:length(TT)
        tt = TT(j-1):incr_t:TT(j);    % interval = incr_ME
        tt2 = tt - TT(j-1);           % t0 converted to 0
        tt2 = tt2*24;                 % h
        
        b = tt';       % + tspan(1,1);
        
        if dos==1 && j==2
            a0 = TT(j-1);
            a = find(abs(c-a0)<0.0001);
        else
            a00 = (TT(j-1)+TT(j))/2;      % + tspan(1,1);
            a = find(abs(c-a00)<0.0001);
        end
             
        % set y[n-1]
        for i=1:nm
         if(n==1)
          yy(i) = ICpk0(i);
         else
          yy(i) = y_1(a,i);
         end
        end
        
        [K0,K] = Kmat_induc_pk(yy,nm);  
        A_induc0 = ME_solution(K0,K,ICpk,tt2);
    
        C_induc0 = 0.5*((A_induc0(:,2)/V1 - A_induc0(:,4) - Kss) + sqrt((A_induc0(:,2)/V1 - A_induc0(:,4) - Kss).^2 + 4*Kss*A_induc0(:,2)/V1));
        
        A_1_0 = [A_1_0;A_induc0];   %#ok<AGROW>
        C_1_0 = [C_1_0;C_induc0];   %#ok<AGROW>
        cc = [cc; b];               %#ok<AGROW>
        if(j<length(TT))
         for m=1:nm
          ICpk(m,1) = A_induc0(end,m);
         end
         A_1_0(end,:) = [];
         C_1_0(end,:) = [];
         cc(end,:) = [];
        end
     end         % end of step size ME
     if(dos<nd)
      A_1_0(end,:) = [];
      C_1_0(end,:) = [];
      cc(end,:) = [];
     end
    end          % end of dosing
    
    y_1 = A_1_0;
    C_induc(:,n) = C_1_0;
    
    % difference between reference and linearisation
    reldif = abs(Cori - C_induc(:,n)) ./ C_induc(:,n);
    ARD(n) = max(reldif);
    
    if n>1
     % difference between successive iterations
     reldif1 = abs(C_induc(:,n) - C_induc(:,n-1)) ./ C_induc(:,n-1);
     cri = max(reldif1);
     ARDc(n) = cri;
    end
    
end  % end of inductive approximation

npk = n;
        
toc        

figure(1)
semilogy(cc/365,Cori,'k-','LineWidth',1); 
hold on
semilogy(cc/365,C_induc(:,round(n/20*1)),'r--','LineWidth',2);
semilogy(cc/365,C_induc(:,round(n/20*3.5)),'g--','LineWidth',2);
semilogy(cc/365,C_induc(:,round(n/20*20)),'b--','LineWidth',2);
n1 = {['n=' num2str(round(n/20*1)) '']};
text(0.9,2.5,n1,'fontsize',11);
n2 = {['n=' num2str(round(n/20*3.5)) '']};
text(2.025,25,n2,'fontsize',11);
n3 = {['n=' num2str(round(n/20*20)) '']};
text(3,11000,n3,'fontsize',11);
xlim([0 year])
ylim([0.1 100000])
set(gca,'fontsize',16)
xlabel('Time (year)','fontsize',20)
ylabel('Free denosumab (ng/mL)','fontsize',20)
hold off
print(gcf,'-dtiff','-r600','PK_inductive.tif')

figure(2)
semilogy(1:max_n,ARDc(:,1),'k^-');
xlabel('Number of inductive approximations','fontsize',18)
ylabel('E_{S} for free denosumab','fontsize',18)
set(gca,'fontsize',14)
set(gca,'xtick',[0:5:max_n]);
ylim([10^-4 10^2])
print(gcf,'-dtiff','-r600','ES_PK.tif')

%% Inductive linearization for bone model

disp(strcat('<Runtime to get linearized bone model>'))
tic

aPointsb = zeros(length(cc),max_n);
aPointsc = zeros(length(cc),max_n);

ICPD = IC0(1:nmpd);
ICPD0 = ICPD;
nm = length(ICPD);

yy = zeros(nm,1);         % for y[n-1]

K0_all = zeros(nm,length(cc));
K_all = zeros(nm,nm,length(cc));

ARD_b = zeros(max_n,1);
ARDc_b = zeros(max_n,1);
ARD_c = zeros(max_n,1);
ARDc_c = zeros(max_n,1);
n = 0;
cri = 9999;
while cri>=cri_thre && n<max_n
    n = n + 1;
    
    K0_all = zeros(nm,length(cc));
    K_all = zeros(nm,nm,length(cc));
    
    if rem(n,max_n/5)==0
     n
    end
    
    ICPD = IC0(1:nmpd);
    
    y_1_0 = [];
    ccc = [];
    
     TT = cc(1):incr_ME:cc(end);
     for j=2:length(TT)
        tt = TT(j-1):incr_t:TT(j);    % interval = incr_ME
        tt2 = tt - TT(j-1);           % t0 converted to 0
        tt2 = tt2*24;                 % h
        
        b = tt';
        
        if  mod(TT(j-1),52*7/2) == 0 && TT(j-1) < 52*7*year
            a0 = TT(j-1);
        else
            a0 = (TT(j-1)+TT(j))/2;
        end
        a = find(abs(c-a0)<0.0001);
        
        Cindj = C_induc(a,npk);       % denosumab conc.
     
        % set y[0]
        for i=1:nm
         if(n==1)
          yy(i) = ICPD0(i);
         else
          yy(i) = y_1(a,i);
         end
        end
        
        [K0,K] = Kmat_induc_PD(yy,nm,ICPD0,Cindj);  
        
        if mod(TT(j-1),52*7/2) == 0 && TT(j-1) < 52*7*year
            for ppj=a:a+1
             K0_all(:,ppj) = K0;
             K_all(:,:,ppj) = K;
            end
        elseif j==length(TT)
            for ppj=a-1:a+1
             K0_all(:,ppj) = K0;
             K_all(:,:,ppj) = K;
            end
        else
            for ppj=a-1:a
             K0_all(:,ppj) = K0;
             K_all(:,:,ppj) = K;
            end
        end
        
        y_induc0 = ME_solution(K0,K,ICPD,tt2);
        
        y_1_0 = [y_1_0;y_induc0];   %#ok<AGROW>
        ccc = [ccc; b];             %#ok<AGROW>
        if(j<length(TT))
         for m=1:nm
          ICPD(m,1) = y_induc0(end,m);
         end
         y_1_0(end,:) = [];
         ccc(end,:) = [];
        end
     end         % end of step size ME
    
    y_1 = y_1_0;
    aPointsb(:,n) = y_1(:,29-2) + y_1(:,30-2);
    aPointsc(:,n) = y_1(:,18-2);
    
    % difference between reference and linearisation
    reldifb = abs(OBori - aPointsb(:,n)) ./ aPointsb(:,n);
    ARD_b(n) = max(reldifb);
    reldifc = abs(OCori - aPointsc(:,n)) ./ aPointsc(:,n);
    ARD_c(n) = max(reldifc);
    
    if n>1        
     % difference between successive iterations
     reldif1b = abs(aPointsb(:,n) - aPointsb(:,n-1)) ./ aPointsb(:,n-1);
     crib = max(reldif1b);
     ARDc_b(n) = crib;
     reldif1c = abs(aPointsc(:,n) - aPointsc(:,n-1)) ./ aPointsc(:,n-1);
     cric = max(reldif1c);
     ARDc_c(n) = cric;
     
     cri = max(crib,cric);
     
    end
    
end  % end of inductive approximation

toc        
 
pairscb = -9999*ones(1,max_n);
   for i=3:2:n
      pairscb(i) = (ARDc_b(i)+ARDc_b(i-1))/2;
   end
pairscc = -9999*ones(1,max_n);
   for i=3:2:n
      pairscc(i) = (ARDc_c(i)+ARDc_c(i-1))/2;
   end
pairscb3 = pairscb(pairscb>=0);
pairscc3 = pairscc(pairscc>=0);

CFB_OB_induc = aPointsb./OB0 * 100;
CFB_OC_induc = aPointsc./IC0(18-2) * 100;

figure(21)
plot(cc/365,CFB_OBori,'k-','LineWidth',1); 
hold on
plot(cc/365,CFB_OB_induc(:,round(n/5*1)),'r--','LineWidth',2);
plot(cc/365,CFB_OB_induc(:,round(n/5*2)),'g--','LineWidth',2);
plot(cc/365,CFB_OB_induc(:,round(n/5*5)),'b--','LineWidth',2);
n1 = {['n=' num2str(round(n/5*1)) '']};
text(0.8,30,n1,'fontsize',11);
n2 = {['n=' num2str(round(n/5*2)) '']};
text(3.5,60,n2,'fontsize',11);
n3 = {['n=' num2str(round(n/5*5)) '']};
text(3.5,85,n3,'fontsize',11);
xlim([0 year])
ylim([-50 150])
set(gca,'fontsize',16)
xlabel('Time (year)','fontsize',20)
ylabel('% of baseline for OB','fontsize',20)
hold off
print(gcf,'-dtiff','-r600','OB_inductive.tif')

figure(22)
semilogy(1:max_n,ARDc_b(:,1),'k^-');
hold on
semilogy(3:2:n,pairscb3,'k--','LineWidth',2);
hold off
xlabel('Number of inductive approximations','fontsize',18)
ylabel('E_{S} for OB','fontsize',18)
set(gca,'fontsize',14)
ylim([10^-4 10^1])
xlim([3 20])
legend('Original','Average');
legend('boxoff')
print(gcf,'-dtiff','-r600','ES_OB.tif')

figure(23)
plot(cc/365,CFB_OCori,'k-','LineWidth',1); 
hold on
plot(cc/365,CFB_OC_induc(:,round(n/5*1)),'r--','LineWidth',2);
plot(cc/365,CFB_OC_induc(:,round(n/5*2)),'g--','LineWidth',2);
plot(cc/365,CFB_OC_induc(:,round(n/5*5)),'b--','LineWidth',2);
n1 = {['n=' num2str(round(n/5*1)) '']};
text(1.7,6,n1,'fontsize',11);
n2 = {['n=' num2str(round(n/5*2)) '']};
text(3.55,25,n2,'fontsize',11);
n3 = {['n=' num2str(round(n/5*5)) '']};
text(3.3,40,n3,'fontsize',11);
xlim([0 year])
ylim([-50 150])
set(gca,'fontsize',16)
xlabel('Time (year)','fontsize',20)
ylabel('% of baseline for OC','fontsize',20)
hold off
print(gcf,'-dtiff','-r600','OC_inductive.tif')

figure(24)
semilogy(1:max_n,ARDc_c(:,1),'k^-');
hold on
semilogy(3:2:n,pairscc3,'k--','LineWidth',2);
hold off
xlabel('Number of inductive approximations','fontsize',18)
ylabel('E_{S} for OC','fontsize',18)
set(gca,'fontsize',14)
set(gca,'xtick',[0:5:max_n]);
ylim([10^-4 10^1])
xlim([3 20])
legend('Original','Average');
legend('boxoff')
print(gcf,'-dtiff','-r600','ES_OC.tif')

CFB_last = zeros(length(cc),3);
CFB_last(:,1) = CFB_OB_induc(:,n);
CFB_last(:,2) = CFB_OC_induc(:,n);

IC_setting

%-------------------------------------------------------------------;
%-------------------------------------------------------------------;
%-------------------------------------------------------------------;
% K0_all & K_all
%      dy/dt = K0_all(time-varying) + K_all(time-varying)*y
%       ---> dy^/dt = M*K0_all + M*K_all*pinv(M)*y^
%-------------------------------------------------------------------;
%-------------------------------------------------------------------;
%-------------------------------------------------------------------;

%% BMD prediction from linearized model

BMD = [];
ccc = [];

     TT = cc(1):incr_ME:cc(end);
     for j=2:length(TT)
        tt = TT(j-1):incr_t:TT(j);    % interval = incr_ME
        tt2 = tt - TT(j-1);           % t0 converted to 0
        tt2 = tt2*24;                 % h
        
        b = tt';       
        
        if TT(j-1)==0 || TT(j-1)==52*7/2
            a0 = TT(j-1);
        else
            a0 = (TT(j-1)+TT(j))/2;    
        end
        a = find(abs(c-a0)<0.0001);
        
        BSAP = aPointsb(a,n);
        SCTX = aPointsc(a,n);
     
        K0 = 0.000146*IC_BMD0*(BSAP/(ICPD0(29-2)+ICPD0(30-2)))^0.0739;
        K = -0.000146*(SCTX/ICPD0(18-2))^0.0779;
        
        BMD0 = ME_solution(K0,K,IC_BMD,tt2);
        
        BMD = [BMD;BMD0];   %#ok<AGROW>
        ccc = [ccc; b];     %#ok<AGROW>
        if(j<length(TT))
         IC_BMD = BMD0(end);
         BMD(end,:) = [];
         ccc(end,:) = [];
        end
     end         % end of step size ME
     
     CFB_BMD = (BMD-IC_BMD0)/IC_BMD0*100;
     
%-------------------------------------------------------------------;
%-------------------------------------------------------------------;
% CFB_BMD = prediction from original model
%-------------------------------------------------------------------;
%-------------------------------------------------------------------;
     
figure(31)
plot(cc/30,CFB_BMD,'b--','LineWidth',1); 
hold on
plot(cc/30,CFB_BMDori,'k-','LineWidth',1);
set(gca,'fontsize',16)
xlabel('Month','fontsize',20)
ylabel('LSBMD,% of baseline','fontsize',20)
legend('Inductive','ode15s','Location','southeast')
legend('boxoff')
print(gcf,'-dtiff','-r600','BMD.tif')

IC_BMD = IC_BMD0;

%% 
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% 2. lump the linearized model
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

%% Prediction from original model with uncertainty (for VPC)

%Parameter values with uncertainty in original paper
Model_parameter_values_SE
 
disp(strcat(num2str(nsim),' replications with SE'))
tic

%Simulate nsim times
simulationpk = zeros(length(ttt),nsim);
simulationb = zeros(length(ttt),nsim);
simulationc = zeros(length(ttt),nsim);
simulationBMD = zeros(length(ttt),nsim);
for j=1:nsim
    
    if rem(j,nsim/5)==0
      j
    end
            
crit_unlumped_PD_SE %Prediction of original unlumped system

simulationpk(:,j)= C;
simulationb(:,j) = CFB_OB;
simulationc(:,j) = CFB_OC;
simulationBMD(:,j) = CFB_BMD;
 
end

toc

% percentile
p2_5 = zeros(length(ttt),4);
p97_5 = zeros(length(ttt),4);
for i=1:length(ttt)
    
    p2_5(i,1) =  quantile(simulationpk(i,:),0.025);
    p97_5(i,1) = quantile(simulationpk(i,:),0.975);
    
    p2_5(i,2) =  quantile(simulationb(i,:),0.025);
    p97_5(i,2) = quantile(simulationb(i,:),0.975);
    
    p2_5(i,3) =  quantile(simulationc(i,:),0.025);
    p97_5(i,3) = quantile(simulationc(i,:),0.975);
    
    p2_5(i,4) =  quantile(simulationBMD(i,:),0.025);
    p97_5(i,4) = quantile(simulationBMD(i,:),0.975);
    
end

% K returned to the typical value
Model_parameter_values0

IC = ICPD0;

%% SA to find the best lumped model within m

% Initial lumping matrix produced from Vx
M = vector2matrix(Vx,mf);
[r,c0]=size(M);

Vbest = zeros(length(IC),c0);
Vlast = zeros(length(IC),c0);
ARDbest = zeros(length(IC),1);
ARDlast = zeros(length(IC),1);
itisbest = zeros(length(IC),1);
CFB_hat_best = zeros(length(cc),1,length(IC));
CFB_hat_last = zeros(length(cc),1,length(IC));
ARDHist = zeros(length(IC),length(maxit));
ARDbestHist = zeros(length(IC),length(maxit));
ppbest = zeros(length(IC),1);
pplast = zeros(length(IC),1);
vpcall = -99*ones(length(IC),1);

V = matrix2vector(M);              % for all states in the original model, get the number of state in the lumped model

while r <= ml
    
  disp(strcat('m=',num2str(r)))
  tic
    
  nrows = r;
      
  if r==mf      
   [ARD,CFB_hat,pp,K0_L,K_L,y_hat,r_BMD] = OBJV_function_K_BMD(M,K0_all,K_all,IC,cc,output,CFB_BMD,year);
   Vbest(r,:) = V;
   Mbest = M;
   ARDbest(r,1) = ARD;
   for jh=1:1
    CFB_hat_best(:,jh,r) = CFB_hat(:,jh); % keep simulation result for best model
   end
   ppbest(r,1) = pp;
  else
   Vbest(r,:) = Vbest(r-incr,:);
   ARDbest(r,1) = ARDbest(r-incr,1);
   for jh=1:1
    CFB_hat_best(:,jh,r) = CFB_hat_best(:,jh,r-incr); % keep simulation result for best model
   end
   ppbest(r,1) = ppbest(r-incr,1);
  end
  
  % SA initialisation
  itis = 1;
  temp = temp_ini; 
  itisHist = zeros(1,length(maxit));
  tempHist = zeros(1,length(maxit));
  while itis <= maxit
    for staten = 1:length(IC)
      if nrows==1
          VNew = V;
      else
          VNew = nearestNeighbour(V, nrows, constraint, staten);
      end
      MNew = vector2matrix(VNew,nrows);
      [ARDnew,CFB_hat_new,pp_new,K0_L,K_L,y_hat,r_BMD] = OBJV_function_K_BMD(MNew,K0_all,K_all,IC,cc,output,CFB_BMD,year);
            
      if ARDnew >= ARD 
        p = exp(-(ARDnew-ARD)/temp); % Boltzman factor
        if p > rand
          V = VNew;
          ARD = ARDnew;
          CFB_hat = CFB_hat_new;
          pp = pp_new;
        end
      end
      if ARDnew < ARD
          V = VNew;
          ARD = ARDnew;
          CFB_hat = CFB_hat_new;
          pp = pp_new;
      end
      if ARDnew < ARDbest(r,1)
          V = VNew;
          M = MNew;
          ARD = ARDnew;
          CFB_hat = CFB_hat_new;
          pp = pp_new;
          Vbest(r,:) = V;
          Mbest = M;
          ARDbest(r,1) = ARD;
          for jh=1:1
           CFB_hat_best(:,jh,r) = CFB_hat(:,jh); % keep simulation result for best model
          end
          ppbest(r,1) = pp;
          itisbest(r) = itis;
      end
    end
    
    itisHist(itis) = itis;
    ARDHist(r,itis) = ARD;
    ARDbestHist(r,itis) = ARDbest(r,1);
    tempHist(itis) = temp;
    
    temp=temp*decayRate;
    itis = itis+1;
    
    if rem(itis,num_dis)==0
      r
      itis
    end
  end
  
  %--------------------------------------------------;
  % for last result in SA
  %--------------------------------------------------;
  Vlast(r,:) = V;
  Mlast = M;
  [rlast,clast] = size(ARDHist);
  ARDlast(r,1) = ARDHist(r,clast);
  for jh=1:1
   CFB_hat_last(:,jh,r) = CFB_hat(:,jh); % keep simulation result for best model
  end
  pplast(r,1) = pp;

  [ vpcall(r,:) ] = OBJV_function_vpc(CFB_hat_best(:,:,r),p2_5(:,4),p97_5(:,4));
  
  % Create a new M with additonal row(s) of zeros
  if r<c0             % if before reaching original model
      V2 = V;
      for iii=1:incr
          V2(iii) = r + iii;
      end
      V = V2;
%       M2 = [Mbest;zeros(incr,c0)]; 
%       M = M2;
  end
%   [r,~]=size(M);
  r = max(V);
  
toc

end

for rep3 = mf:incr:ml
    
disp(strcat('m=',num2str(rep3)))

Vbestc = num2cell(Vbest(rep3,:));
comexp_original2 = [comexp_original;Vbestc];
disp(comexp_original2)

figure(200+rep3)
subplot(1,2,1)
semilogy(itisHist,ARDHist(rep3,:),'b-','LineWidth',2);
xlim([0 maxit])
ylim([10^-11 10^4])
xlabel('Number of iterations','fontsize',12)
ylabel('Sum of squares (last)','fontsize',12)

subplot(1,2,2)
semilogy(itisHist,ARDbestHist(rep3,:),'b-','LineWidth',2);
xlim([0 maxit])
ylim([10^-1 10^4])
xlabel('Number of iterations','fontsize',12)
ylabel('Sum of squares','fontsize',12)
eval(['print -dtiff -r600 Diagnostics_m' num2str(rep3) '.tif'])
    
%---------------------------------------------------------------
% VPC
%---------------------------------------------------------------
figure(300)
plot(ttt/30,real(CFB_hat_best(:,1,rep3)),'r-','LineWidth',1);
hold on
plot(ttt/30,CFB_BMDori,'k--','LineWidth',1); % original
plot(ttt/30,p2_5(:,4),'b--','LineWidth',1);
plot(ttt/30,p97_5(:,4),'b--','LineWidth',1);
text(10,13,[num2str(round(vpcall(rep3,1),2)) '% outside 95% CI'],'fontsize',9);
% xlim([0 12])
ylim([0 20])
set(gca,'fontsize',9)
xlabel('Time (month)','fontsize',16)
ylabel('% of baseline BMD','fontsize',16)
hold off
legend(['Reduced model (m=' num2str(rep3) ')'],'Original model','95%CI from original model','Location','northwest')
legend('boxoff')
% set(gcf,'PaperPositionMode','Auto')
eval(['print -dtiff -r600 VPC' num2str(rep3) '.tif'])

% Save SS for calculation of composite criterion
str = ['ARDbest', num2str(rep3), ' = ARDbest(rep3);'];
eval(str);
str2 = ['save ARDbest', num2str(rep3), '.mat ARDbest', num2str(rep3)];
eval(str2);

end

IC = ICPD0;

L = Mbest;
L_inv = pinv(L);

%%

diary off

%% End of the code
