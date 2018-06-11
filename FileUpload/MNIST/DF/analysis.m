clc
clear
%
load ../../data
nsk     = 50;
lbls    = data(1:nsk:end,1);
data    = data(1:nsk:end,2:end);
ds = [1 2 3 4 5];
inds = [];
for i=1:length(ds)
    inds = [inds;find(lbls==ds(i))];
    ix{i} = find(lbls==ds(i));
end
data = data(inds,:);
lbls = lbls(inds);
inds = [];
for i=1:length(ds)
    inds = [inds;find(lbls==ds(i))];
    ix{i} = find(lbls==ds(i));
end
nS = size(data,1);
mappedX = tsne(full(data), lbls, 2, 30, []);
load mappedXDF
%Xdata = [];
Xdata = [Xdata;mappedX];
size(Xdata)
%
nItem = length(ds);
for i=1:nItem
    indx{i} = find(lbls==ds(i));
end
for i=1:nItem
    for j=1:nItem
        d3(i,j)    = distdist(mappedX(indx{i},:),mappedX(indx{j},:));
    end
end
for i=1:nItem
    for j=1:nItem
        dvd(i,j)    = d3(i,j)/sqrt(d3(i,i)*d3(j,j));
    end
end
%
a3 = reshape(dvd,nItem*nItem,1);a3(a3==0) = [];
%aa = [];ss = [];
aa = [aa;mean(a3)]
s = silhouette(mappedX,lbls);
ss = [ss;s];
save mappedXDF Xdata aa ss
