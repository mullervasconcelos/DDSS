% Generate random covariance matrices

clear all
clc

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

% 
% 
% H1 = 0.25*randn(n);
% 
% S1 = H1*H1';
% 
% eig(S1);
% 
% M1 = -1 + rand(n,1);
% 
% %x1 = mvnrnd(M1,S1)
% 
% H2 = 0.3*randn(n);
% 
% S2 = H2*H2';
% 
% eig(S2);
% 
% M2 = rand(n,1);
% 
% H3 = 0.5*randn(n);
% 
% S3 = H3*H3';
% 
% eig(S3);
% 
% M3 = 1 + 0.5*rand(n,1);
% 
% mu = [M1';M2';M3'];
% 
% sigma = cat(3,S1,S2,S3);
% 
% gm = gmdistribution(mu,sigma,[0.3 0.4 0.3]);

T = [];

for N = [10 100 1000 10000 100000];

data = random(gm,N);

tic
%main_unicast
%test_subgradient
options = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton');
% %options = optimoptions(@fminunc,'Algorithm','quasi-newton');
%[x,J]=fminunc('cost_unicast',zeros(1,n),options);
 [x,J]=fminunc('cost_broadcast',zeros(1,2*n*(n-1)),options);
T = [T (toc)]


end




%var(data)

%hist(data(:,3),100)


%x = mvnrnd(M,S)