clc
clear
%
load data
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
load DF/mappedXDF
subplot(1,3,1)
sD = reshape(ss,599,5);sD = mean(sD');
st = 599*4+1;en = 599*5;
aD = aa;
f = Xdata(st:en,:);
subplot(1,3,1)
hold on
plot(f(ix{1},1),f(ix{1},2),'ko')
plot(f(ix{2},1),f(ix{2},2),'ro')
plot(f(ix{3},1),f(ix{3},2),'bo')
plot(f(ix{4},1),f(ix{4},2),'go')
plot(f(ix{5},1),f(ix{5},2),'co')
%
load PD/mappedXPD
aa = aa(1:5);
aP = aa;ss = ss(1:599*5);
sP = reshape(ss,599,5);sP = mean(sP');
st = 599*2+1;en = 599*3;
f = Xdata(st:en,:);
subplot(1,3,2)
hold on
plot(f(ix{1},1),f(ix{1},2),'ko')
plot(f(ix{2},1),f(ix{2},2),'ro')
plot(f(ix{3},1),f(ix{3},2),'bo')
plot(f(ix{4},1),f(ix{4},2),'go')
plot(f(ix{5},1),f(ix{5},2),'co')
%
subplot(1,3,3)
plot(sP,sD,'ko')
hold on
plot([-1 1],[-1 1],'r--')
