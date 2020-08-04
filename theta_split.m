function [w,b] = theta_split(theta)

global n

%n = 4;

theta_parsed = reshape(theta,2*(n-1),n);

w = [];
b = [];

for i = 1:2*(n-1)
   
    if mod(i,2) == 1
       
        w = [w; theta_parsed(i,:)];
        
    else
        
        b = [b; theta_parsed(i,:)];
        
    end
        
    
end

wnew = eye(n);

for i = 1:n
   
    for j = 1:n
        
        if i > j 
        
        wnew(i,j) = w(i-1,j);
        
        end
        
        if i < j 
        
        wnew(i,j) = w(i,j);
        
        end
        
    end
    
end

w = wnew;

bnew = zeros(n);

for i = 1:n
   
    for j = 1:n
        
        if i > j 
        
        bnew(i,j) = b(i-1,j);
        
        end
        
        if i < j 
        
        bnew(i,j) = b(i,j);
        
        end
        
    end
    
end

b = bnew;