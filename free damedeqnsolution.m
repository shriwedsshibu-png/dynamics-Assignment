clc;
clear;
close all;

%% System parameters
m = 1000;       % Mass (kg)
zeta = 0.2;       % Damping ratio

%% Initial conditions
x0 = 0;         % Initial displacement (m)
v0 = 5;         % Initial velocity (m/s)

%% Calculate system properties
wn = 5;                 % Natural frequency (rad/s)
wd = wn*sqrt(1-zeta^2);        % Damped natural frequency (rad/s)

%% Time
t = linspace(0,10,1000);

%% Damped free vibration response
x = (v0/wd) .* exp(-zeta*wn*t) .* sin(wd*t);

%% Plot displacement
plot(t,x);

xlabel('Time (s)');
ylabel('Displacement (m)');
title('Damped Free Vibration: x_0 = 0, v_0 \notequal 0');

grid on;