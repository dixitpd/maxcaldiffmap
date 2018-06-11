function d = distdist( c1,c2 )
d = 0;
d1 = 0;
d2 = 0;
%
l1 = length(c1);
l2 = length(c2);
d1 = mean(sum((c1.*c1)'));
d2 = mean(sum((c2.*c2)'));
d3 = 2*(sum(mean(c1).*mean(c2)));
d = d1+d2-d3;
d = sqrt(d);
% takes square root
end
