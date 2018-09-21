
% Calculate percentage outside CI

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%  Lumped version of
%    Peterson MC, Riggs MM (2010) Bone 46:49-63
%                        +
%    Peterson MC, Riggs MM (2012) CPT Pharmacometrics Syst Pharmacol 1:e14
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

%%

function [ vpc ] = OBJV_function_vpc(C_hat_best,p2_5,p97_5)

[~,c_BMD] = size(C_hat_best);

count=zeros(c_BMD,1);
for j=1:length(C_hat_best)
    for p=1:c_BMD
     if C_hat_best(j,p)<p2_5(j,p) || p97_5(j,p)<C_hat_best(j,p) || isnan(C_hat_best(j,p))
        count(p) = count(p)+1;
     end
    end
end

vpc = count/length(C_hat_best)*100;                    % percent of deviation

end
