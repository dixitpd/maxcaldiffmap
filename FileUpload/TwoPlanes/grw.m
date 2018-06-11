function [ kij prb ] = grw( Delta )

nS = size(Delta,1);
nr = sum(Delta');prb = nr/sum(nr);
nr = repmat(nr,nS,1);
kij = Delta./nr';


end

