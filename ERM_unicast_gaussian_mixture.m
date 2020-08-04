global data n

%set dimension
n = 2;

%2 dimensional examples in the paper

S1 = [1 0; 0 1];
S2 = [1 0.4; 0.4 1];

M1 = [0;0];
M2 = [4;2];

mu = [M1';M2'];

sigma = cat(3,S1,S2);

gm = gmdistribution(mu,sigma,[0.75 0.25])

data = random(gm,2);

ezsurf('cost_unicast',[-1,1,-1,1])
