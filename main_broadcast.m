global n data theta

%n = 4;

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

end

Jstar = J2;



% tic
% options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton')
% [thetastar,Jstar]=fminunc('cost_broadcast',thetanew,options)
% toc



%subgradient_broadcast(data)