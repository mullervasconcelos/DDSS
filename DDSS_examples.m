% Generate random covariance matrices

global n data

%set dimension
n = 2;

%2 dimensional examples in the paper

S1 = [1 0; 0 1];
S2 = [1 0.4; 0.4 1];

M1 = [0;0];
M2 = [4;2];

mu = [M1';M2'];

sigma = cat(3,S1,S2);

gm = gmdistribution(mu,sigma,[0.75 0.25])

%T = [];

%for N = [10 100 1000 10000 100000]

%data = random(gm,N);

%tic

%main_unicast

%test_subgradient

%options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');

%[x,J]=fminunc('cost_unicast',zeros(1,n),options);

%[x,J]=fminunc('cost_broadcast',zeros(1,2*n*(n-1)),options);
 
%T = [T (toc)]

%end

