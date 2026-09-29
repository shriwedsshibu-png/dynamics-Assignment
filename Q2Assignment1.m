clear all
clc
close all

load('elcentro.mat')
ug = p(:) * 0.01;   % The (:) forces this into a vertical column
dt = t(2) - t(1);
time = t(:);        % The (:) forces this into a vertical column
N = length(ug);

mass = [1600 800]; %only mass varies 
k_stiff = [60000 60000];
damp = [0.02 0.02];
%newmark
u_new = zeros(N,2);
v_new = zeros(N,2);
a_new = zeros(N,2);
%Runge-kutta
u_rk = zeros(N,2);
v_rk = zeros(N,2);
a_rk = zeros(N,2);
%dulham
u_dh = zeros(N,2);
v_dh = zeros(N,2);
a_dh = zeros(N,2);

for c=1:2
    m = mass(c);
    k = k_stiff(c);
    z = damp(c);
    
    wn = sqrt(k/m);
    wd = wn*sqrt(1 - z^2);
    cc = 2*m*wn*z;
    
    % newmark refer table 5.4.2
    gam = 0.5; %avg acc
    bet = 0.25; 
    peff = -m*ug;
    
    a_new(1,c) = (peff(1) - cc*v_new(1,c) - k*u_new(1,c))/m;
    khat = k + (gam/(bet*dt))*cc + (1/(bet*dt^2))*m;
    a_const = (1/(bet*dt))*m + (gam/bet)*cc;
    b_const = (1/(2*bet))*m + dt*(gam/(2*bet) - 1)*cc;
    
    for i = 1:N-1
        dp = peff(i+1) - peff(i);
        dp_hat = dp + a_const*v_new(i,c) + b_const*a_new(i,c);
        du = dp_hat/khat;
        
        dv = (gam/(bet*dt))*du - (gam/bet)*v_new(i,c) + dt*(1 - gam/(2*bet))*a_new(i,c);
        da = (1/(bet*dt^2))*du - (1/(bet*dt))*v_new(i,c) - (1/(2*bet))*a_new(i,c);
        
        u_new(i+1,c) = u_new(i,c) + du; %taking the time step and using the last value 
        v_new(i+1,c) = v_new(i,c) + dv;
        a_new(i+1,c) = a_new(i,c) + da;
    end
    
    % rk4 method
    y1 = zeros(N,1); 
    y2 = zeros(N,1);
    
    for i=1:N-1
        ugmid = (ug(i)+ug(i+1))/2;
        
        k1u = y2(i);
        k1v = -(cc/m)*y2(i) - (k/m)*y1(i) - ug(i);
        
        k2u = y2(i) + 0.5*dt*k1v;
        k2v = -(cc/m)*(y2(i) + 0.5*dt*k1v) - (k/m)*(y1(i) + 0.5*dt*k1u) - ugmid;
        
        k3u = y2(i) + 0.5*dt*k2v;
        k3v = -(cc/m)*(y2(i) + 0.5*dt*k2v) - (k/m)*(y1(i) + 0.5*dt*k2u) - ugmid;
        
        k4u = y2(i) + dt*k3v;
        k4v = -(cc/m)*(y2(i) + dt*k3v) - (k/m)*(y1(i) + dt*k3u) - ug(i+1);
        
        y1(i+1) = y1(i) + (dt/6)*(k1u + 2*k2u + 2*k3u + k4u);
        y2(i+1) = y2(i) + (dt/6)*(k1v + 2*k2v + 2*k3v + k4v);
    end
    
    u_rk(:,c) = y1;
    v_rk(:,c) = y2;
    a_rk(:,c) = -2*z*wn.*v_rk(:,c) - (wn^2).*u_rk(:,c) - ug;
    
    % duhamel part 
    A = cumtrapz(time, ug .* exp(z*wn*time) .* cos(wd*time)); % integral using trapezoidal method
    B = cumtrapz(time, ug .* exp(z*wn*time) .* sin(wd*time));
    
    u_dh(:,c) = -(1/wd) .* exp(-z*wn*time) .* (A .* sin(wd*time) - B .* cos(wd*time));
    
    v_dh(:,c) = gradient(u_dh(:,c), dt);
    a_dh(:,c) = -2*z*wn.*v_dh(:,c) - (wn^2).*u_dh(:,c) - ug;
end

% plotting section
figure(1)
subplot(3,1,1)
plot(time, u_new(:,1), 'b')
hold on
grid on
plot(time, u_new(:,2), 'r')
title('NEWMARKS DISPLACEMENT')
legend('case 1', 'case 2')

subplot(3,1,2)
plot(time, v_new(:,1), 'b')
hold on
grid on
plot(time, v_new(:,2), 'r')
title('NEWMARKS VELOCITY')

subplot(3,1,3)
plot(time, a_new(:,1), 'b')
hold on
grid on
plot(time, a_new(:,2), 'r')
title('NEWMARKS ACCELERATION')
xlabel('time (s)')


figure(2)
subplot(3,1,1)
plot(time, u_rk(:,1), 'b')
hold on
grid on
plot(time, u_rk(:,2), 'r')
title('RUNGE KUTTA METHOD DISSPLACEMENT')
legend('case 1', 'case 2')

subplot(3,1,2)
plot(time, v_rk(:,1), 'b')
hold on
grid on
plot(time, v_rk(:,2), 'r')
title('RUNGE KUTTA METHOD VELOCITY')

subplot(3,1,3)
plot(time, a_rk(:,1), 'b')
hold on
grid on
plot(time, a_rk(:,2), 'r')
title('RUNGE KUTTA METHOD ACCELERATION')
xlabel('time (s)')


figure(3)
subplot(3,1,1)
plot(time, u_dh(:,1), 'b')
hold on
grid on
plot(time, u_dh(:,2), 'r')
title('DUHAMELS DISPLACEMENT')
legend('case 1', 'case 2')

subplot(3,1,2)
plot(time, v_dh(:,1), 'b')
hold on
grid on
plot(time, v_dh(:,2), 'r')
title('DUHAMELS VELOCITY')

subplot(3,1,3)
plot(time, a_dh(:,1), 'b')
hold on
grid on
plot(time, a_dh(:,2), 'r')
title('DUMAMELS ACCELERATION')
xlabel('time (s)')
