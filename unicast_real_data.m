%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UNICAST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

global n theta data

load('train_batch.mat')

data = train_batch;

[N,n] = size(data);

tic

mean_data = mean(data);

x0 = mean_data;

theta = x0;

J0 = cost_unicast(theta)

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
    
end

end

thetastar = thetanew
 
Jstar = cost_unicast(thetanew)
 
toc
 
% COMPARISON WITH BLIND SCHEDULING 

Jblind = sum(var(data)) - max(var(data))

((Jblind - Jstar)/Jblind)*100

% VALIDATION ANALYSIS

load('test_batches.mat')

M = 20;

U = [];

T = [];


for m = 1:M
   
    m
    
data = test_batches((m-1)*N+1:m*N,:);

tic
options = optimoptions(@patternsearch,'Display','off','FunctionTolerance', 1e-4);
[thetastar_m,Jstar_m]=patternsearch(@cost_unicast,x0,[],[],[],[],[],[],[],options);
T = [T; (toc)];

U = [U; cost_unicast(thetastar) - Jstar_m]

end

alpha = 0.05;

Ubar = mean(U);

epsilon = tinv(1-alpha,M-1)*sqrt(var(U))/sqrt(M)

Ubar + epsilon

sum(U<0)

T

%x0 = mean_data;

%tic
%xhat = x0;
%options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');
%[xhatstar,Jstar]=ga(@cost_unicast,n);
%xhatstar
%Jstar
%toc





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BROADCAST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%test_batches((3)*s+1:(3+1)*s,1)

%figure

%hist(train_batch(:,1),100)


% %xhat = mean_data;
% 
% xhat = mean_data;

% J0 = cost_unicast(xhat);
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
% J1 = cost_unicast(xhatnew);
% 
% delta = abs(J0 - J1)/(1+abs(J0));
% 
% while delta>=10^-6
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
% J1 = cost_unicast(xhat);
% 
% J2 = cost_unicast(xhatnew);
% 
% delta =  abs(J1 - J2)/(1+abs(J1));
% 
% end
%  
%  xhatstar = xhatnew
% 
%  Jstar = cost_unicast(xhatstar)
%  
%  
%  d = 2*n*(n-1);
% 
% mean_data = mean(data);
% 
% % 
% % tic
% % options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton')
% % [thetastar,Jstar]=fminunc('cost_broadcast',theta0,options)
% % toc
% 
% 
% mean_data_squared = mean(data.^2);
% 
% A = sparse(d,d); 
% 
% for i = 1:n
%     
%     Ai = zeros(2);
%     
%     
%    
%     Ai(1,1) = mean_data_squared(i);
%     Ai(1,2) = mean_data(i);
%     Ai(2,1) = mean_data(i);
%     Ai(2,2) = 1;
%     
%     for j=1:n-1
%     
%     A(2*(n-1)*(i-1)+1 + 2*(j-1)  :2*(n-1)*(i-1)+2+ 2*(j-1),2*(n-1)*(i-1)+1+ 2*(j-1):2*(n-1)*(i-1)+2+ 2*(j-1)) = Ai;
%     
%     end
%     
% end
% 
% A = 2*A;
% 
% b = zeros(2,n*(n-1));
% 
% for j = 1:n
%     
%     for i=1:j-1
%         
%         b(:, (n-1)*(j-1) +i) =  [mean(data(:,i).*data(:,j)); mean_data(i)];
%         
%     end
%     
%     for i=j+1:n
%         
%         b(:, (n-1)*(j-1) + i-1) =  [mean(data(:,i).*data(:,j)); mean_data(i)];
%         
%     end
%            
%         
% end
%     
% b = reshape(b,d,1);
% 
% theta = zeros(1,d);
% 
% %toc
% 
% b = 2*b;
% 
% %tic
% 
% J0 = cost_broadcast(theta)
% 
% g = zeros(N,d);
% 
% for k = 1:N
%     
%     g(k,:) = subgradient_broadcast(data(k,:)');
%     
% end
% 
% g = mean(g)';
% 
% thetanew = A\(g+b);
% 
% J1 = cost_broadcast(thetanew);
% 
% delta = abs(J0 - J1)/(1+abs(J0));
% 
% while delta>=10^-6
% 
% theta = thetanew;
% 
% g = zeros(N,d);
% 
% for k = 1:N
%     
%     g(k,:) = subgradient_broadcast(data(k,:)');
%     
% end
% 
% g = mean(g)';
% 
% thetanew = A\(g+b);
% 
% J1 = cost_broadcast(theta);
% 
% J2 = cost_broadcast(thetanew);
% 
% delta =  abs(J1 - J2)/(1+abs(J1));
% 
% end
% 
% Jstar_broadcast = J2
%  
% %tic
% %xhat = x0;
% %options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');
% %[thetastar,Jstar]=fminunc('cost_broadcast',thetanew,options)
% %toc
% 
% 
% 
%  
%  