function[g] = subgradient_unicast(x)

global n xhat

Gstar = -inf;
jstar = 0;



for j=1:n
    
    
    
    G = (x(j)-xhat(j))^2;
    
    if G>Gstar
        
        Gstar = G;
        jstar = j;
        
    end
    
end


g = -2*(x(jstar)-xhat(jstar))*[zeros(jstar-1,1);1;zeros(n-jstar,1)]';

