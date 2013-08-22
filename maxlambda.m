
function maxl=maxlambda(adj)

[V,D]=eig(adj);
maxl=max(diag(D));