clc

global n data

d = 2*n*(n-1);

theta = rand(d,1);

N = 1000;

data = random(gm,N);

tic
options = optimoptions(@fminunc,'Display','iter')
[thetastar,Jstar]=fminunc(@cost_broadcast,theta,options)
toc

