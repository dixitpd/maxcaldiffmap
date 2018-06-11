function [ kij prb ] = grw( Delta,nS )

nr = sum(Delta');nr = repmat(nr,nS,1);
kij = Delta./nr';
[v e] = eigs(sparse(kij'),5);
e = diag(e);ii = abs(e-1);ii = find(ii==min(ii));
prb = v(:,ii);prb = prb/sum(prb);

end

