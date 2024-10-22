clear all

clc

global n theta data

DDSS_examples

N = 10000;
 
data = random(gm,N);

mean_data = mean(data);
 
theta0 = zeros(1,n);

theta = theta0;

J0 = cost_unicast(theta0);

g = zeros(N,n);

for k = 1:N
    
    g(k,:) = subgradient_unicast(data(k,:));
    
end

thetanew = 0.5*mean(g)+mean_data;

J1 = cost_unicast(thetanew);

delta = abs(J0 - J1);

while delta>=10^-4
    
theta = thetanew;

g = zeros(N,n);

for k = 1:N
    
    g(k,:) = subgradient_unicast(data(k,:));
    
end

thetanew = 0.5*mean(g)+mean_data;

J1 = cost_unicast(theta);

J2 = cost_unicast(thetanew);

delta =  abs(J1 - J2);

end

thetastar = thetanew;

Jstar = cost_unicast(thetastar)
 
Jblind = sum(var(data)) - max(var(data))
 

