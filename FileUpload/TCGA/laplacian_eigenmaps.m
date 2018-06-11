clc
clear
%
load tcga_raw
'done'
%
mnx = mean(tcga);stx = std(tcga);
ii = find(stx==0);tcga(:,ii) = [];
tcga = zscore(tcga);
nS = size(tcga,1);
f = tcga.*(tcga);f(isnan(f)) = 0;f = sum(f')';f = f/mean(f);
pf = exp(-6*f);
% %
pdt = squareform(pdist(tcga));
aa  = reshape(pdt,nS*nS,1);aa(aa==0) = [];gM = prctile(aa,10);
W = exp(-pdt.*pdt/(2*gM*gM));
kD = grw(W,nS);
L = sparse(eye(nS))-kD;
[vD eD] = eigs(L,5,'SR');eD = diag(eD);
[a b] = sort(eD);vD = vD(:,b([2 3]));
kP = get_discrete_mat(W,pf);
L = sparse(eye(nS))-kP;
[vP eP] = eigs(L,5,'SR');eP = diag(eP);
[a b] = sort(eP);vP = vP(:,b([2 3]));
% % 
i0 = find(labels==0);
i1 = find(labels==1);
i2 = find(labels==2);
i3 = find(labels==3);
i4 = find(labels==4);
%
subplot(1,3,1)
hold on
plot(vP(i0,1),vP(i0,2),'ko')
plot(vP(i1,1),vP(i1,2),'ro')
plot(vP(i2,1),vP(i2,2),'bo')
plot(vP(i3,1),vP(i3,2),'go')
plot(vP(i4,1),vP(i4,2),'co')
subplot(1,3,2)
hold on
plot(vD(i0,1),vD(i0,2),'ko')
plot(vD(i1,1),vD(i1,2),'ro')
plot(vD(i2,1),vD(i2,2),'bo')
plot(vD(i3,1),vD(i3,2),'go')
plot(vD(i4,1),vD(i4,2),'co')
%
indx{1} = i0;indx{2} = i1;indx{3} = i2;indx{4} = i3;indx{5} = i4;
for i=1:5
    for j=1:5
        d1(i,j)    = distdist(tcga(indx{i},:),tcga(indx{j},:));
        d2(i,j)    = distdist(vP(indx{i},:),vP(indx{j},:));
        d3(i,j)    = distdist(vD(indx{i},:),vD(indx{j},:));
    end
end
for i=1:5
    for j=1:5
        dplain(i,j) = d1(i,j)/sqrt(d1(i,i)*d1(j,j));
        dvp(i,j)    = d2(i,j)/sqrt(d2(i,i)*d2(j,j));
        dvd(i,j)    = d3(i,j)/sqrt(d3(i,i)*d3(j,j));
    end
end
%
dplain = dplain-diag(diag(dplain));
a1 = reshape(dplain,25,1);a1(a1==0) = [];
dvp = dvp-diag(diag(dvp));
a2 = reshape(dvp,25,1);a2(a2==0) = [];
dvd = dvd-diag(diag(dvd));
a3 = reshape(dvd,25,1);a3(a3==0) = [];
[mean(a1) mean(a2) mean(a3) 100*(mean(a2)-mean(a3))/mean(a3)]
sD = silhouette(vD,labels,'Euclidean');
sP = silhouette(vP,labels,'Euclidean');
subplot(1,3,3)
plot(sP,sD,'ko')
hold on
plot([-1 1],[-1 1])