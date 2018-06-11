function [ kij prb eta ] = merw( Delta )

[v e] = eigs(sparse(Delta),5);
e = diag(e);ii = find(e==max(e));
e
nu  = v(:,ii);nu = nu/norm(nu);eta = e(ii);
kij = diag(1./nu)*Delta*(diag(nu))/eta;
%
[v e] = eigs(sparse(kij'),5);
e = diag(e);
ii = abs(e-1);ii = find(ii==min(ii));
prb = v(:,ii);prb = prb/sum(prb);

end

