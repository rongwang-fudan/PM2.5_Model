function [ rou ] = residu( act )

global M Emi X

d2=size(X,1);
alpha=act(1);
beta=act(2);
gama=act(3);
fc=X(:,14)./(X(:,15)+X(:,14)); % observed fraction of particles that have a diameter of 2.5-10 micrometer in PM10 (f10)
Y1=log(X(:,15).*(1-gama*fc)+alpha*X(:,12)+beta*X(:,13));

[b1,bint1,r1,rint1,stats1]=regress(Y1,[ones(d2,1) M Emi]);
rou=-stats1(1,1); % R2 for A

