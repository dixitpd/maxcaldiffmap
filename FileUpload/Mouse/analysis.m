clc
clear

%
% Perform t-test to show better separation
% 

%
% 1:124 ->CLP 125:248 ->GMP 249:372 ->HSC 373:496 ->LMPP 497:597 ->PreMegE
%
load cleaned_up_data
data = Dct;nS = size(data,1);data(data==-2) = NaN;dt = zeros(nS,nS);
%
for i=1:nS-1
    for j=i+1:nS
        dt = naneucdist(data(i,:),data(j,:));
        pdt(i,j) = dt;pdt(j,i) = dt;
    end
end
W = phatekernel(pdt,5,8);W = W - diag(diag(W));
D = sum(W);D = diag(1./D);
[kD pD] = grw(D*W*D,nS);
[vD eD] = get_n_vecs(kD,3,1);
% % %
enx = -data.*log(data);enx(isnan(enx)) = 0;
enx = sum(enx')';enx = 1./(1+exp(6*enx));
enx = enx/sum(enx);
[kP] = get_discrete_mat(W,enx);
[vP eP] = get_n_vecs(kP,3,1);
% %
i1 = find(lbls==1);i2 = find(lbls==2);
i3 = find(lbls==3);i4 = find(lbls==4);
i5 = find(lbls==5);
subplot(1,2,1)
mx = max(max(vP));mn = min(min(vP));
vP = (vP-mn)/(mx-mn);
hold on
plot3(vP(i1,1),vP(i1,2),vP(i1,3),'bo')
plot3(vP(i2,1),vP(i2,2),vP(i2,3),'ro')
plot3(vP(i3,1),vP(i3,2),vP(i3,3),'go')
plot3(vP(i4,1),vP(i4,2),vP(i4,3),'ko')
plot3(vP(i5,1),vP(i5,2),vP(i5,3),'co')
mx = max(max(vP));mn = min(min(vP));
xlim([mn mx])
ylim([mn mx])
zlim([mn mx])
view(170,45)
% %
subplot(1,2,2)
hold on
mx = max(max(vD));mn = min(min(vD));
vD = (vD-mn)/(mx-mn);
plot3(vD(i1,1),vD(i1,2),vD(i1,3),'bo')
plot3(vD(i2,1),vD(i2,2),vD(i2,3),'ro')
plot3(vD(i3,1),vD(i3,2),vD(i3,3),'go')
plot3(vD(i4,1),vD(i4,2),vD(i4,3),'ko')
plot3(vD(i5,1),vD(i5,2),vD(i5,3),'co')
mx = max(max(vD));mn = min(min(vD));
xlim([mn mx])
ylim([mn mx])
zlim([mn mx])
view(170,60)
for i=1:5
    for j=1:5
        m1 = vP(find(lbls==i),:);m2 = vP(find(lbls==j),:);
        dP(i,j) = distdist(m1,m2);
        m1 = vD(find(lbls==i),:);m2 = vD(find(lbls==j),:);
        dD(i,j) = distdist(m1,m2);
    end
end
for i=1:5
    for j=1:5
        d1(i,j) = dP(i,j)/sqrt(dP(i,i)*dP(j,j));
        d2(i,j) = dD(i,j)/sqrt(dD(i,i)*dD(j,j));
    end
end
dPd = sqrt(diag(1./diag(dP)));dP = dPd*dP*dPd;dP = dP - diag(diag(dP));
dDd = sqrt(diag(1./diag(dD)));dD = dDd*dD*dDd;dD = dD - diag(diag(dD));
aa = reshape(dP,25,1);aa(aa==0) = [];
bb = reshape(dD,25,1);bb(bb==0) = [];
[h p] = ttest(aa-bb)
mean((aa-bb)./bb)