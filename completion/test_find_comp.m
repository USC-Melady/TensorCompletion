n1 = 100;
n2 = 200;
p = 0.2;
seed = 2;
d = 3;r=3;

G = rand(n1+n2, n1+n2);
G(G > p) =1;
G(G <=p ) = 0;

submat_idx = find_comp([n1,n2], sparse(G), seed, r, d);

