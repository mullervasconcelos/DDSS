
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

N = 100000;

data = random(gm,N);

[N,n] = size(data);

%%%%%%%%%% CCP Broadcast


d = 2*n*(n-1);

mean_data = mean(data);

theta0 = zeros(d,1);

% 
% tic
% options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton')
% [thetastar,Jstar]=fminunc('cost_broadcast',theta0,options)
% toc


mean_data_squared = mean(data.^2);

A = sparse(d,d); 

%tic
for i = 1:n
    
    Ai = zeros(2);
   
    Ai(1,1) = mean_data_squared(i);
    Ai(1,2) = mean_data(i);
    Ai(2,1) = mean_data(i);
    Ai(2,2) = 1;
    
    for j=1:n-1
    
    A(2*(n-1)*(i-1)+1 + 2*(j-1)  :2*(n-1)*(i-1)+2+ 2*(j-1),2*(n-1)*(i-1)+1+ 2*(j-1):2*(n-1)*(i-1)+2+ 2*(j-1)) = Ai;
    
    end
    
end
%toc

A = 2*A;

b = zeros(2,n*(n-1));

%tic

for j = 1:n
    
    for i=1:j-1
        
        b(:, (n-1)*(j-1) +i) =  [mean(data(:,i).*data(:,j)); mean_data(i)];
        
    end
    
    for i=j+1:n
        
       
        
        b(:, (n-1)*(j-1) + i-1) =  [mean(data(:,i).*data(:,j)); mean_data(i)];
        
    end
           
        
end
    
b = reshape(b,d,1);

%toc


b = 2*b;


%tic

tic

theta = theta0;

J0 = cost_broadcast(theta0);

g = zeros(N,d);

for k = 1:N
    
    g(k,:) = subgradient_broadcast(data(k,:)');
    
end

g = mean(g)';

%tic
thetanew = A\(g+b);
%toc

J1 = cost_broadcast(thetanew);

delta = J0 - J1;

if delta < 0 
    
    disp('ERROR!')
    
end

while delta>=10^-4

theta = thetanew;

g = zeros(N,d);

%tic
for k = 1:N
    
    g(k,:) = subgradient_broadcast(data(k,:)');
    
end
%toc


g = mean(g)';

%tic
thetanew = A\(g+b);
%toc

J1 = cost_broadcast(theta);

J2 = cost_broadcast(thetanew);

delta =  J1 - J2;

if delta < 0 
    
    disp('ERROR!')
    
end

end

thetastar = thetanew

Jstar = J2


toc
%%%%%%%%%%%%%%%%%%%%%%%%

% VALIDATION

M = 100;

U = [];

tic

T=[];

for m =1:M
    
    m

data = random(gm,N);

%xhat = mean(data);

%PROBLEM = createOptimProblem('fmincon','objective','cost_broadcast','x0',theta0);
%GS = GlobalSearch('Display','off');
%[theta_m,Jstar_m]=run(GS,PROBLEM);

tic
options = optimoptions(@patternsearch,'Display','off','FunctionTolerance', 1e-4);
[thetastar_m,Jstar_m]=patternsearch(@cost_broadcast,thetastar,[],[],[],[],[],[],[],options);
T = [T; (toc)];


U = [U; cost_broadcast(thetastar) - Jstar_m];
 
end

alpha = 0.05;
 
Ubar = mean(U);
 
epsilon = tinv(1-alpha,M-1)*sqrt(var(U))/sqrt(M)
 
gap= Ubar + epsilon

toc
 
sum(U<0)

