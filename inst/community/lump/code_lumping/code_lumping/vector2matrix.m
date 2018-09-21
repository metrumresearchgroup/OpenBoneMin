
% Function to convert V to M (lumping matrix)

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%  Lumped version of
%    Peterson MC, Riggs MM (2010) Bone 46:49-63
%                        +
%    Peterson MC, Riggs MM (2012) CPT Pharmacometrics Syst Pharmacol 1:e14
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

%%

function [ M ] = vector2matrix( V, rows)

M = zeros(rows, length(V));

 for j=1:length(V)
   M(V(j),j) = 1;
 end

end
