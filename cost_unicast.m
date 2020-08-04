function[J] = cost_unicast(theta)

global data

[N,n] = size(data);

F = zeros(N,1);

G = zeros(N,1);

for k=1:N

F(k) = norm(theta-data(k,:),2)^2;

G(k) = max((theta-data(k,:)).^2);

end

J = mean(F-G);