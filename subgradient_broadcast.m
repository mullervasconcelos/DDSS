function[g] = subgradient_broadcast(x)

global n theta 

d = 2*n*(n-1);

Gstar = -inf;
jstar = 0;

%theta

[w,b] = theta_split(theta);

for j=1:n
    
    G = 0;
    
    for l=1:n
        
        if l~=j 
            
            for i=1:n
                
                if i~=l
                    
                    G = G + (x(i)-w(i,l)*x(l)-b(i,l))^2;
                    
                end
                
            end
            
        end
        
     end
    
    
    
    if G>Gstar
        
        Gstar = G;
        jstar = j;
        
    end
    
end

k = zeros(2,n*(n-1));

%tic

for l = 1:n
    
    if l~=jstar
    
    for i=1:l-1
        
        k(:, (n-1)*(l-1) +i) =  [-2*(x(l)*(x(i)-(w(i,l)*x(l)+b(i,l)))); -2*(x(i)-(w(i,l)*x(l)+b(i,l)))];
        
    end
    
    for i=l+1:n
        
        %(n-1)*(j-1) + i-1
        
        k(:, (n-1)*(l-1) + i-1) =  [-2*(x(l)*(x(i)-(w(i,l)*x(l)+b(i,l)))); -2*(x(i)-(w(i,l)*x(l)+b(i,l)))];
        
    end
    
    end
                
end
    
k = reshape(k,d,1);





% k = [];
% 
% for l = 1:n
%     
%     if l~=jstar
%         
%         for i = 1:n
%             
%             if i~=l
%                 
%                 k = [k; -2*(x(l)*(x(i)-(w(i,l)*x(l)+b(i,l)))); -2*(x(i)-(w(i,l)*x(l)+b(i,l)))];
%                 
%             end
%             
%         end
%         
%     else
%         
%         for i = 1:n
%             
%             if i~=l
%                 
%                 k = [k; 0; 0];
%                 
%             end
%             
%         end
%         
%         
%     end
%     
% end

g = k';

