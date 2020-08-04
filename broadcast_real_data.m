%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BROADCAST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

global n theta data

load('train_batch.mat')

data = train_batch;

[N,n] = size(data);

d = 2*n*(n-1);

mean_data = mean(data);

mean_data_squared = mean(data.^2);

A = sparse(d,d); 

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

A = 2*A;

b = zeros(2,n*(n-1));

for j = 1:n
    
    for i=1:j-1
        
        b(:, (n-1)*(j-1) +i) = [mean(data(:,i).*data(:,j)); mean_data(i)];
        
    end
    
    for i=j+1:n
       
        b(:, (n-1)*(j-1) + i-1) = [mean(data(:,i).*data(:,j)); mean_data(i)];
        
    end
           
        
end
    
b = reshape(b,d,1);

b = 2*b;

%%%%%%%



theta0 = zeros(2,n*(n-1));

cost_broadcast(zeros(d,1));

for j = 1:n
    
    for i=1:j-1
        
        theta0(:, (n-1)*(j-1) +i) =  [0; mean_data(i)];
        
    end
    
    for i=j+1:n
       
        theta0(:, (n-1)*(j-1) + i-1) =  [0; mean_data(i)];
        
    end
           
        
end

%mean(data)

%theta0
    
theta0 = reshape(theta0,d,1);

%%%%%%

tic

theta = theta0;

J0 = cost_broadcast(theta);

g = zeros(N,d);

for k = 1:N
    
    g(k,:) = subgradient_broadcast(data(k,:)');
    
end

g = mean(g)';

thetanew = A\(g+b);

J1 = cost_broadcast(thetanew);

delta = J0 - J1;

if delta < 0 
    
    disp('ERROR!')
    
end

while delta>=10^-4

theta = thetanew;

g = zeros(N,d);

for k = 1:N
    
    g(k,:) = subgradient_broadcast(data(k,:)');
    
end

g = mean(g)';

thetanew = A\(g+b);

J1 = cost_broadcast(theta);

J2 = cost_broadcast(thetanew);

delta =  J1 - J2;

if delta < 0 
    
    disp('ERROR!')
    
end

end



thetastar = thetanew;

toc

pause

Jstar = J2

%%%%

Jblind = sum(var(data)) - max(var(data))

((Jblind - Jstar)/Jblind)*100


 % VALIDATION ANALYSIS

load('test_batches.mat')

M = 20;

U = [];

T = [];

for m = 1:20
   
    m
    
data = test_batches((m-1)*N+1:m*N,:);

mean_data = mean(data);

mean_data_squared = mean(data.^2);

A = sparse(d,d); 

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

A = 2*A;

b = zeros(2,n*(n-1));

for j = 1:n
    
    for i=1:j-1
        
        b(:, (n-1)*(j-1) +i) =  [mean(data(:,i).*data(:,j)); mean_data(i)];
        
    end
    
    for i=j+1:n
       
        b(:, (n-1)*(j-1) + i-1) =  [mean(data(:,i).*data(:,j)); mean_data(i)];
        
    end
           
        
end
    
b = reshape(b,d,1);

b = 2*b;


% options = optimoptions(@simulannealbnd,'Display','iter','FunctionTolerance',1e-4,'TemperatureFcn',@temperatureboltz,'InitialTemperature',1);
% [x,fval] = simulannealbnd(@cost_broadcast,thetanew,[],[],options)

%opts = optimoptions(@fmincon,'Algorithm','sqp','Display','iter','FunctionTolerance',1e-2,'MaxFunctionEvaluations',Inf);
%problem = createOptimProblem('fmincon','objective',@cost_broadcast,'x0',thetanew,'lb',[],'ub',[],'options',opts);


% opts = optimoptions(@fminunc,'Algorithm','quasi-newton','Display','iter','StepTolerance',1e-4);
% problem = createOptimProblem('fminunc','objective',@cost_broadcast,'x0',thetanew,'options',opts);
% %ms = MultiStart;
% [x,f] = fminunc(problem)

%%%%%%

tic

theta = thetastar;

J0 = cost_broadcast(theta);

g = zeros(N,d);

for k = 1:N
    
    g(k,:) = subgradient_broadcast(data(k,:)');
    
end

g = mean(g)';

thetanew = A\(g+b);

J1 = cost_broadcast(thetanew)

delta = J0 - J1;

if delta < 0 
    
    disp('ERROR!')
    
end

while delta>=10^-4

theta = thetanew;

g = zeros(N,d);

for k = 1:N
    
    g(k,:) = subgradient_broadcast(data(k,:)');
    
end

g = mean(g)';

thetanew = A\(g+b);

J1 = cost_broadcast(theta);

J2 = cost_broadcast(thetanew);

delta =  J1 - J2;

if delta < 0 
    
    disp('ERROR!')
    
end

end

thetastar_m = thetanew;

Jstar_m = J2;

%%%%

T = [T; (toc)]

U = [U; cost_broadcast(thetastar) - Jstar_m]

end

alpha = 0.05;

Ubar = mean(U);

epsilon = tinv(1-alpha,M-1)*sqrt(var(U))/sqrt(M)

Bound = Ubar + epsilon
 
sum(U<0)













% 
% %pause
% 
% G = [];
% 
% for m = 1:20
%     
%     m
%    
% data = test_batches((m-1)*N+1:m*N,:);
% 
% %%%%%%
% 
% theta0 = zeros(d,1);
% 
% %%%%%%%
% 
% mean_data_squared = mean(data.^2);
% mean_data = mean(data);
% 
% A = sparse(d,d); 
% 
% %tic
% for i = 1:n
%     
%     Ai = zeros(2);
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
% %toc
% 
% A = 2*A;
% 
% b = zeros(2,n*(n-1));
% 
% %tic
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
% %toc
% 
% b = 2*b;
% 
% %%%%%%
% 
% theta = theta0;
% 
% J0 = cost_broadcast(theta0);
% 
% g = zeros(N,d);
% 
% 
% for k = 1:N
%     
%     g(k,:) = subgradient_broadcast(data(k,:)');
%     
% end
% 
% g = mean(g)';
% 
% %tic
% thetanew = A\(g+b);
% %toc
% 
% J1 = cost_broadcast(thetanew);
% 
% delta = abs(J0 - J1);
% 
% while delta>=10^-4
% 
% theta = thetanew;
% 
% g = zeros(N,d);
% 
% %tic
% for k = 1:N
%     
%     g(k,:) = subgradient_broadcast(data(k,:)');
%     
% end
% %toc
% 
% g = mean(g)';
% 
% %tic
% thetanew = A\(g+b);
% %toc
% 
% J1 = cost_broadcast(theta);
% 
% J2 = cost_broadcast(thetanew);
% 
% delta =  abs(J1 - J2);
% 
% end
% 
% Jstar_test = J2
% 
% %%%%%%
% 
% % thetastar_test = thetanew;
% 
% G = [G; cost_broadcast(thetastar) - Jstar_test]
% 
% end
% 
% %hist(G,30)
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
% %%%%
% 
% 
% 
% % 
% % 
% % %%%%%%
% % 
% % xhat = x0;
% % 
% % J0 = cost_unicast(x0);
% % 
% % g = zeros(N,n);
% % 
% % for k = 1:N
% %     
% %     g(k,:) = subgradient_unicast(data(k,:));
% %     
% % end
% % 
% % xhatnew = 0.5*mean(g)+mean_data;
% % 
% % J1 = cost_unicast(xhatnew);
% % 
% % delta = abs(J0 - J1)/(1+abs(J0));
% % 
% % while delta>=10^-6
% %     
% % xhat = xhatnew;
% % 
% % g = zeros(N,n);
% % 
% % for k = 1:N
% %     
% %     g(k,:) = subgradient_unicast(data(k,:));
% %     
% % end
% % 
% % xhatnew = 0.5*mean(g)+mean_data;
% % 
% % J1 = cost_unicast(xhat);
% % 
% % J2 = cost_unicast(xhatnew);
% % 
% % delta =  abs(J1 - J2)/(1+abs(J1));
% % 
% % end
% % 
% % xhatstar = xhatnew
% %  
% % Jstar = cost_unicast(xhatnew)
% %  
% % toc
% %  
% % % COMPARISON WITH BLIND SCHEDULING 
% % 
% % Jblind = sum(var(data)) - max(var(data))
% % 
% % ((Jblind - Jstar)/Jblind)*100
% % 
% % % VALIDATION ANALYSIS
% % 
% % load('test_batches.mat')
% % 
% % M = 20;
% % 
% % %pause
% % 
% % G = [];
% % 
% % for m = 1:M
% %    
% % data = test_batches((m-1)*N+1:m*N,:);
% % 
% % mean_data = mean(data);
% % 
% % x0 = zeros(1,n);
% % 
% % %x0 = mean_data;
% % 
% % %%%%
% % xhat = x0;
% % 
% % J0 = cost_unicast(x0);
% % 
% % g = zeros(N,n);
% % 
% % for k = 1:N
% %     
% %     g(k,:) = subgradient_unicast(data(k,:));
% %     
% % end
% % 
% % xhatnew = 0.5*mean(g)+mean_data;
% % 
% % J1 = cost_unicast(xhatnew);
% % 
% % delta = abs(J0 - J1)/(1+abs(J0));
% % 
% % while delta>=10^-6
% %     
% % xhat = xhatnew;
% % 
% % g = zeros(N,n);
% % 
% % for k = 1:N
% %     
% %     g(k,:) = subgradient_unicast(data(k,:));
% %     
% % end
% % 
% % xhatnew = 0.5*mean(g)+mean_data;
% % 
% % J1 = cost_unicast(xhat);
% % 
% % J2 = cost_unicast(xhatnew);
% % 
% % delta =  abs(J1 - J2)/(1+abs(J1));
% % 
% % end
% % %%%%
% % 
% % %xhatstar = ponto;
% % 
% % xhatstar_test = xhatnew;
% % 
% % Jstar_test = cost_unicast(xhatstar_test);
% % 
% % G = [G; cost_unicast(xhatstar) - Jstar_test];
% % 
% % end
% % 
% % %hist(G,30)
% % 
% % alpha = 0.01;
% % 
% % Gbar = mean(G);
% % 
% % epsilon = tinv(1-alpha,M-1)*sqrt(var(G))/sqrt(M)
% % 
% % Gbar + epsilon
% % 
% % 
% % %vbar = mean(V)
% % 
% % %sqrt(varbar)
% % 
% % 
% % 
% %  
% %  %x0 = mean_data;
% % 
% % %tic
% % %xhat = x0;
% % %options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');
% % %[xhatstar,Jstar]=ga(@cost_unicast,n);
% % %xhatstar
% % %Jstar
% % %toc
% % 
% % 
% % 
% % 
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BROADCAST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % 
% % %test_batches((3)*s+1:(3+1)*s,1)
% % 
% % %figure
% % 
% % %hist(train_batch(:,1),100)
% % 
% % 
% % % %xhat = mean_data;
% % % 
% % % xhat = mean_data;
% % 
% % % J0 = cost_unicast(xhat);
% % % 
% % % g = zeros(N,n);
% % % 
% % % for k = 1:N
% % %     
% % %     g(k,:) = subgradient_unicast(data(k,:));
% % %     
% % % end
% % % 
% % % xhatnew = 0.5*mean(g)+mean_data;
% % % 
% % % J1 = cost_unicast(xhatnew);
% % % 
% % % delta = abs(J0 - J1)/(1+abs(J0));
% % % 
% % % while delta>=10^-6
% % %     
% % % xhat = xhatnew;
% % % 
% % % g = zeros(N,n);
% % % 
% % % for k = 1:N
% % %     
% % %     g(k,:) = subgradient_unicast(data(k,:));
% % %     
% % % end
% % % 
% % % xhatnew = 0.5*mean(g)+mean_data;
% % % 
% % % J1 = cost_unicast(xhat);
% % % 
% % % J2 = cost_unicast(xhatnew);
% % % 
% % % delta =  abs(J1 - J2)/(1+abs(J1));
% % % 
% % % end
% % %  
% % %  xhatstar = xhatnew
% % % 
% % %  Jstar = cost_unicast(xhatstar)
% % %  
% % %  
% % %  d = 2*n*(n-1);
% % % 
% % % mean_data = mean(data);
% % % 
% % % % 
% % % % tic
% % % % options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton')
% % % % [thetastar,Jstar]=fminunc('cost_broadcast',theta0,options)
% % % % toc
% % % 
% % % 
% % % mean_data_squared = mean(data.^2);
% % % 
% % % A = sparse(d,d); 
% % % 
% % % for i = 1:n
% % %     
% % %     Ai = zeros(2);
% % %     
% % %     
% % %    
% % %     Ai(1,1) = mean_data_squared(i);
% % %     Ai(1,2) = mean_data(i);
% % %     Ai(2,1) = mean_data(i);
% % %     Ai(2,2) = 1;
% % %     
% % %     for j=1:n-1
% % %     
% % %     A(2*(n-1)*(i-1)+1 + 2*(j-1)  :2*(n-1)*(i-1)+2+ 2*(j-1),2*(n-1)*(i-1)+1+ 2*(j-1):2*(n-1)*(i-1)+2+ 2*(j-1)) = Ai;
% % %     
% % %     end
% % %     
% % % end
% % % 
% % % A = 2*A;
% % % 
% % % b = zeros(2,n*(n-1));
% % % 
% % % for j = 1:n
% % %     
% % %     for i=1:j-1
% % %         
% % %         b(:, (n-1)*(j-1) +i) =  [mean(data(:,i).*data(:,j)); mean_data(i)];
% % %         
% % %     end
% % %     
% % %     for i=j+1:n
% % %         
% % %         b(:, (n-1)*(j-1) + i-1) =  [mean(data(:,i).*data(:,j)); mean_data(i)];
% % %         
% % %     end
% % %            
% % %         
% % % end
% % %     
% % % b = reshape(b,d,1);
% % % 
% % % theta = zeros(1,d);
% % % 
% % % %toc
% % % 
% % % b = 2*b;
% % % 
% % % %tic
% % % 
% % % J0 = cost_broadcast(theta)
% % % 
% % % g = zeros(N,d);
% % % 
% % % for k = 1:N
% % %     
% % %     g(k,:) = subgradient_broadcast(data(k,:)');
% % %     
% % % end
% % % 
% % % g = mean(g)';
% % % 
% % % thetanew = A\(g+b);
% % % 
% % % J1 = cost_broadcast(thetanew);
% % % 
% % % delta = abs(J0 - J1)/(1+abs(J0));
% % % 
% % % while delta>=10^-6
% % % 
% % % theta = thetanew;
% % % 
% % % g = zeros(N,d);
% % % 
% % % for k = 1:N
% % %     
% % %     g(k,:) = subgradient_broadcast(data(k,:)');
% % %     
% % % end
% % % 
% % % g = mean(g)';
% % % 
% % % thetanew = A\(g+b);
% % % 
% % % J1 = cost_broadcast(theta);
% % % 
% % % J2 = cost_broadcast(thetanew);
% % % 
% % % delta =  abs(J1 - J2)/(1+abs(J1));
% % % 
% % % end
% % % 
% % % Jstar_broadcast = J2
% % %  
% % % %tic
% % % %xhat = x0;
% % % %options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');
% % % %[thetastar,Jstar]=fminunc('cost_broadcast',thetanew,options)
% % % %toc
% % % 
% % % 
% % % 
% % %  
% % %  