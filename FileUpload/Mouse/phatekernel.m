function [ delt ] = phatekernel( dists,neib,beta )

nS = size(dists,1);
for i=1:nS
    a = dists(i,:);m = sort(a);
    sg(i) = m(neib+1);
end
%
f1 = repmat(sg,nS,1);f2 = f1';
%
delt = exp(-(dists./f1).^(beta)) + exp(-(dists./f2).^(beta));


end

