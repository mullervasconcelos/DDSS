function[J] = cost_broadcast(theta)

global n data

[w,b] = theta_split(theta);

[N,n] = size(data);

F = zeros(N,1);

%G = zeros(N,1);

for k = 1:N

    X = data(k,:);
    
    f = zeros(n,1);
    
    for j = 1:n
        
        for i = 1:n
            
            f(j) = f(j) + ( X(i) - ( w(i,j)*X(j) + b(i,j) ) )^2;
                  
        end
        
    end
    
    f;
    
    F(k) = min(f);
            
end

J = mean(F);