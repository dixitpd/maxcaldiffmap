function [ km ] = get_discrete_mat( wmat,prbs )

nS = size(wmat,1);
prbs = prbs/sum(prbs);
lamb = ones(nS,1);
bets = ones(nS,1);lbx = lamb + 1e5;
while norm ( (lamb-lbx)./lamb ) > 1e-5
    lbx = lamb;
    aa   = wmat*lamb;
    bets = prbs./aa;
    aa   = wmat*bets;
    lamb = prbs./aa;
end
%
c = bets./prbs;
km = diag(c)*wmat*diag(lamb);



end