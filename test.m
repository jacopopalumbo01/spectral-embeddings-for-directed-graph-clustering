W = [0 1 0 0
     0 0 1 0
     0 0 0 1
     1 0 0 0];

W = [0 1 0 0 0 0
     0 0 1 1 0 0
     1 0 0 0 0 0
     0 0 0 0 1 0
     0 0 1 0 0 1
     0 0 0 1 0 0];


G = digraph(W);
plot(G);

[clusters, eigvec, eigval] = BCS(W, 2)
