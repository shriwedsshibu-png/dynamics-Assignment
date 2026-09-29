clc;
clear all;
close all;
zeta=0.4;
omegan=5;
vinital=8;
omegad=omegan.*sqrt(1-zeta^2);
t=linspace(0,5,100);
x=exp(-zeta .*omegan*t) .*vinital .*sin(omegad*t);
plot(t,x);

