function [ vecs e] = get_n_vecs( kmat,n,t )

[v e] = eigs(sparse(kmat),n+1);
e = diag(e);[a b] = sort(e,'descend');
e = a(2:end);e = e/e(1);
for i=1:n
    tt = v(:,b(i+1));tt = tt / norm(tt);
    vecs(:,i) = tt;(e(i)^t);
end


end

