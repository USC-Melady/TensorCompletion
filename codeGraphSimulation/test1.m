n = 100;
r = 2;
d = 3;
c = log(n);
p = c* r /n;


rho = 0.05;
[ G, subIndexes ] = generateGraph( n,rho,p );

nRuns = 1;
[ pre1,rec1,pre2,rec2 ] = getResult( G, d, r ,n, nRuns);