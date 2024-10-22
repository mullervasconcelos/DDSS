clear all

clc

global n data theta

%Bivariate Gaussian mixture example in the paper

S1 = [1 0; 0 1];
S2 = [1 0.4; 0.4 1];

M1 = [0;0];
M2 = [4;2];

mu = [M1';M2'];

sigma = cat(3,S1,S2);

gm = gmdistribution(mu,sigma,[0.75 0.25]);

N = 1000;

data = random(gm,N);

[N,n] = size(data);

%%%%%%%%%% CCP Unicast

tic

mean_data = mean(data);

x0 = mean_data;

theta = x0;

J0 = cost_unicast(x0);

g = zeros(N,n);

for k = 1:N
    
    g(k,:) = subgradient_unicast(data(k,:));
    
end

thetanew = 0.5*mean(g)+mean_data;

J1 = cost_unicast(thetanew);

delta = J0 - J1;

if delta < 0 
    
    disp('ERROR!')
    
end

while delta>=10^-4
    
theta = thetanew;

g = zeros(N,n);

for k = 1:N
    
    g(k,:) = subgradient_unicast(data(k,:));
    
end

thetanew = 0.5*mean(g)+mean_data;

J1 = cost_unicast(theta);

J2 = cost_unicast(thetanew);

delta =  J1 - J2;

if delta < 0 
    
    disp('ERROR!')
    
    break
    
end

end

toc

thetastar = thetanew

Jstar = cost_unicast(thetastar)

M = 100

U = [];

T=[];

for m =1:M

   m
   
data = random(gm,N);

theta = mean(data);

% PROBLEM = createOptimProblem('fmincon','objective','cost_unicast','x0',mean(data));
% GS = GlobalSearch('Display','off');
% [theta_m,Jstar_m]=run(GS,PROBLEM);

tic
options = optimoptions(@patternsearch,'Display','off','FunctionTolerance', 1e-4);
[thetastar_m,Jstar_m]=patternsearch(@cost_unicast,thetastar,[],[],[],[],[],[],[],options);
T = [T; (toc)];

U = [U; cost_unicast(thetastar) - Jstar_m];
 
end

alpha = 0.05;
 
Ubar = mean(U);
 
epsilon = tinv(1-alpha,M-1)*sqrt(var(U))/sqrt(M)
 
gap = Ubar + epsilon

sum(U<0)

