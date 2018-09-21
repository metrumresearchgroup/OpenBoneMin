
% Function for providing different M in SA

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%  Lumped version of
%    Peterson MC, Riggs MM (2010) Bone 46:49-63
%                        +
%    Peterson MC, Riggs MM (2012) CPT Pharmacometrics Syst Pharmacol 1:e14
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

%% 

function [ Vnew ] = nearestNeighbour( V , rows , constraint , staten)

j = V(staten);

Vnew = V;

b = nonzeros(constraint);
c = 0;
for i=1:length(b)
    if staten == b(i)
        c = c + 1;
    end
end

if c ~= 1                  % for free states

  V2 = V(V==j);
  
  if length(V2)>1               % only when this state consists of more than one original state, not to make all zeros row
      
    [l,~] = size(constraint);
    population0 = l+1:rows;
    jnew2 = population0(population0 ~= j);
    jnew = jnew2(randsample(1:length(jnew2),1));
%     Vnew = V;
    Vnew(staten) = jnew;
      
  end 

end

while length(unique(Vnew)) ~= rows
    num = 1:length(V);
    staten2 = num(num ~= staten);
    staten3 = randsample(staten2,1);
    Vnew(staten3) = randsample(1:rows,1);
end

end
