
% Calculate criterion values for lumped model selection

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%  Lumped version of
%    Peterson MC, Riggs MM (2010) Bone 46:49-63
%                        +
%    Peterson MC, Riggs MM (2012) CPT Pharmacometrics Syst Pharmacol 1:e14
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

%% Composite Criterion

mfc = 6;
mlc = 28;

% Import SS for each m to calculate criterion values
ARDbest2 = zeros(length(IC),1);
for i=mf:incr:ml
    str = ['load ARDbest', num2str(i), ';'];
    eval(str);
    str2 = ['ARDbest2(i) = ARDbest', num2str(i), ';'];
    eval(str2);
end    
ARDbest2 = real(ARDbest2);

if(ml==length(IC)-1)
    lasto = 1;
else
    lasto = 0;
end;
m = [mf:incr:ml+lasto]';            % include original model

penal = [1:length(IC)]';

ia = 1;

obj3 = zeros(length(mf:incr:ml+lasto),4);
for p0=50:50                    % 30:20:90
    
    alpha = p0/100;
    
  CRa = (ARDbest-0)/(ARDbest(mfc)-0);
  CRb = (penal-mfc)/(length(IC)-mfc);
  obj = alpha*CRa + (1-alpha)*CRb;

CRa2 = CRa(mf:incr:ml+lasto);
CRb2 = CRb(mf:incr:ml+lasto);
obj2 = obj(mf:incr:ml+lasto);

obj3(:,ia) = obj2;

figure(2000)
plot(m,CRa2,'-bo','LineWidth',1);
hold on
plot(m,CRb2,'-b^','LineWidth',1);
plot(m,obj2,'--g*','LineWidth',2);
xlim([mf ml])
ylim([0 1])
set(gca,'fontsize',12)
set(gca,'xtick',mf:1:length(IC)+1)
xlabel('Number of states','fontsize',12)
ylabel('Criteria','fontsize',12)
hold off
legend('Criterion for deviation','Criterion for penalty','OBJ')
legend('boxoff')
eval(['print -dtiff -r600 Criteria_' num2str(p0) 'percent.tif'])

ia = ia + 1;

end

ia = 1;

%% End of the code
