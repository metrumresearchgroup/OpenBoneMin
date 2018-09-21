
% Function for converting M to V 

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%  Lumped version of
%    Peterson MC, Riggs MM (2010) Bone 46:49-63
%                        +
%    Peterson MC, Riggs MM (2012) CPT Pharmacometrics Syst Pharmacol 1:e14
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

%% 

function [ V ] = matrix2vector( M )

[m,n] = size(M);
V=zeros(1,n);
for i=1:n
  for j=1:m
    if (M(j,i)==1)
      V(i) = j;
      break;
    end
  end
end

end
