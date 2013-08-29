% The eigenvalues of the Laplacian of the graph
% INPUTs: adjacency matrix
% OUTPUTs: laplacian eigenvalues, sorted

function s=graph_spectrum_adjmatr_unsort(adj)

[v,D]=eig(adj);
s=diag(D); % sort in decreasing order