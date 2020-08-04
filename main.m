clear all
clc

global n data

table = [];

nvar = [2 5 10 20 50]

for n = nvar
    
    n
    
    T = [];

for k=1:10


DDSS_examples

tic
main_unicast
%test_subgradient
%PROBLEM = createOptimProblem('fmincon','objective','cost_unicast','x0',mean(data));
%GS = GlobalSearch('Display','off');
%tic
%[x_hat_m,Jstar_m]=run(GS,PROBLEM);
T = [T; (toc)];


end

table = [table T]


end


errorbar(nvar,mean(table),std(table))


% d = 2*n*(n-1);
% 
% mean_data = mean(data);
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
%        
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
% tic
% %test_subgradient
% options = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton');
% %options = optimoptions(@fminunc,'Algorithm','quasi-newton');
% %[x,J]=fminunc('cost_unicast',mean(data),options);
% [x,J]=fminunc('cost_broadcast',b,options);
% toc

