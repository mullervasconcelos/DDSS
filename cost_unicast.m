function[J] = cost_unicast(xhat)

global data

[N,n] = size(data);

F = zeros(N,1);

G = zeros(N,1);

for k=1:N

F(k) = norm(xhat-data(k,:),2)^2;

G(k) = max((xhat-data(k,:)).^2);

end

J = mean(F-G);