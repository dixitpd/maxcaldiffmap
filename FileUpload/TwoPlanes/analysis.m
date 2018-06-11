clc
clear
%
% Create synthetic data
%
nSamp = 500; % number of sample points in each plane
nDim = 10; % total number of dimensions
nClust = 2;inds{1} = 1:nSamp;inds{2} = nSamp+1:2*nSamp;
labels(inds{1}) = 0;labels(inds{2}) = 1;
for i=1:nDim-1
    d1(:,i) = randn(nSamp,1);
    d2(:,i) = randn(nSamp,1);
end
% coordinates of the two planes
d1 = [zeros(nSamp,1) d1];  % 
d2 = [2*ones(nSamp,1) d2]; %
data = [d1;d2];
% pairwise distance in the data and selection of epsilon
pdt = pdist(data);epsilon = prctile(pdt,10);
pdt = squareform(pdt);
% constructing the kernel
W = exp(-pdt.*pdt/(2*epsilon*epsilon));
nS = size(pdt,1);
%
subplot(2,3,1)
hold on
plot3(d1(:,1),d1(:,2),d1(:,3),'ko')
plot3(d2(:,1),d2(:,2),d2(:,3),'ro')
%
subplot(2,3,2)
hold on
plot3(d1(:,2),d1(:,3),d1(:,4),'ko')
plot3(d2(:,2),d2(:,3),d2(:,4),'ro')
subplot(2,3,3)
hold on
a = pdt(1:nSamp,1:nSamp);a = reshape(a,nSamp*nSamp,1);
a(a==0) = [];
b = pdt(nSamp+1:2*nSamp,nSamp+1:2*nSamp);b = reshape(b,nSamp*nSamp,1);
b(b==0) = [];x1 = [a;b];
x2 = reshape(pdt(1:nSamp,nSamp+1:2*nSamp),nSamp*nSamp,1);
[b a] = hist(x1,0:0.1:10);b = b/sum(b);
plot(a,b,'b')
[b a] = hist(x2,0:0.1:10);b = b/sum(b);
plot(a,b,'b--')
%
% Standard Markov chain 
%
[kD pD] = grw(W);
%
% Pathe entropy maximized Markov chain 
%
[kP pP] = merw(W);
%
% Construct diffusion maps
%
[vD eD] = get_n_vecs(kD,3,1);
[vP eP] = get_n_vecs(kP,3,1);
%
subplot(2,3,4)
hold on
plot(vP(1:nSamp,1),vP(1:nSamp,2),'ko')
plot(vP(nSamp+1:2*nSamp,1),vP(nSamp+1:2*nSamp,2),'ro')
%
subplot(2,3,5)
hold on
plot(vD(1:nSamp,1),vD(1:nSamp,2),'ko')
plot(vD(nSamp+1:2*nSamp,1),vD(nSamp+1:2*nSamp,2),'ro')
%
nClust = 2;inds{1} = 1:nSamp;inds{2} = nSamp+1:2*nSamp;
clear d1 d2
for i=1:nClust
    for j=1:nClust
        d1(i,j)    = distdist(data(inds{i},:),data(inds{j},:));
        d2(i,j)    = distdist(vP(inds{i},:),vP(inds{j},:));
        d3(i,j)    = distdist(vD(inds{i},:),vD(inds{j},:));
    end
end
for i=1:nClust
    for j=nClust
        dplain(i,j) = d1(i,j)/sqrt(d1(i,i)*d1(j,j));
        dvp(i,j)    = d2(i,j)/sqrt(d2(i,i)*d2(j,j));
        dvd(i,j)    = d3(i,j)/sqrt(d3(i,i)*d3(j,j));
    end
end
%
dplain = dplain-diag(diag(dplain));
a1 = reshape(dplain,nClust*nClust,1);a1(a1==0) = [];
dvp = dvp-diag(diag(dvp));
a2 = reshape(dvp,nClust*nClust,1);a2(a2==0) = [];
dvd = dvd-diag(diag(dvd));
a3 = reshape(dvd,nClust*nClust,1);a3(a3==0) = [];
[mean(a1) mean(a2) mean(a3) 100*(mean(a2)-mean(a3))/mean(a3)]
sD = silhouette(vD,labels,'Euclidean');
sP = silhouette(vP,labels,'Euclidean');
subplot(2,3,6)
plot(sP,sD,'o')
hold on
plot([-1 1],[-1 1])

%