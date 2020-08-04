global n xhat data


%DDSS_examples

%x

%xhat

%subgradient_unicast(x,xhat)

% delta = 1e-6;
% 
% N = 10000;
% 
% data = random(gm,N);
% 

mean_data = mean(data);
% 
x0 = zeros(1,n);
% 


%x0 = zeros(1,n);

%N = 1000;

%data = random(gm,N);

%tic
%xhat = x0;
%options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton')
%[thetastar,Jstar]=fminunc('cost_unicast',xhat,options)
%toc





xhat = x0;

J0 = cost_unicast(x0);

g = zeros(N,n);

for k = 1:N
    
    g(k,:) = subgradient_unicast(data(k,:));
    
end

xhatnew = 0.5*mean(g)+mean_data;

J1 = cost_unicast(xhatnew);

delta = abs(J0 - J1);

while delta>=10^-4
    
xhat = xhatnew;

g = zeros(N,n);

for k = 1:N
    
    g(k,:) = subgradient_unicast(data(k,:));
    
end

xhatnew = 0.5*mean(g)+mean_data;

J1 = cost_unicast(xhat);

J2 = cost_unicast(xhatnew);

delta =  abs(J1 - J2);

end



% 
 xhatstar = xhatnew;
% 
 Jstar = cost_unicast(xhatstar);
% 
% Jblind = sum(var(data)) - max(var(data))
% 
% Jblind - Jstar

% 
% norm(xhatstar1-xhatstar2,2)
% 
% Jstar1
% Jstar2


% M = 100;
% 
% G = [];
% 
% for m = 1:M
%    
% 
% data = random(gm,N);
% 
% mean_data = mean(data);
% 
% x0 = mean_data;
% 
% g = zeros(N,n);
% 
% for k = 1:N
%     
%     g(k,:) = subgradient_unicast(data(k,:));
%     
% end
% 
% xhatnew = 0.5*mean(g)+mean_data;
% 
% while norm(xhat-xhatnew,2)>=delta
% 
% xhat = xhatnew;
% 
% g = zeros(N,n);
% 
% for k = 1:N
%     
%     g(k,:) = subgradient_unicast(data(k,:));
%     
% end
% 
% xhatnew = 0.5*mean(g)+mean_data;
% 
% end
% 
% %toc
% 
% xhatstar_test = xhatnew;
% 
% Jstar_test = cost_unicast(xhatstar_test);
% 
% G = [G; cost_unicast(xhatstar) - Jstar_test];
% 
% end
% 
% hist(G,30)
% 
% alpha = 0.01;
% 
% Gbar = mean(G);
% 
% epsilon = tinv(1-alpha,M-1)*sqrt(var(G))/sqrt(M)
% 
% Gbar + epsilon
% 
% 
% %vbar = mean(V)
% 
% %sqrt(varbar)
% 
% 
