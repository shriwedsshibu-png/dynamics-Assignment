clear; 
clc; 
close all;

load('elcentro.mat'); %get the data provided by Sir% WE are using Avg Accleration method
ug = p*0.01;  %ground accleration%           
dt = t(2) - t(1);  %TIME STEP%
zeta_array = [0.01, 0.02, 0.03, 0.04]; % It is the Damping ratios i.e 1%, 2%, 3%, 4%
T = 0.01:0.05:4.0;                     % T is Time increment : Natural periods from 0.01 to 4 seconds with 0.05 sec increment upto 4 Sec with the critriea - EQN 5.4.19.b
wn_array = 2 * pi ./ T;                % Row of natural frequency with Time increment
m = 1;                               %unit mass   
Sd = zeros(length(zeta_array), length(T)); %  we are deifining spectral displacement
Sv = zeros(length(zeta_array), length(T)); % velocity spectra
Sa = zeros(length(zeta_array), length(T)); %accleration spectra
gamma = 0.5; % as per given values from textbook
beta = 0.25; % as per given values from textbook

% Generate Response Spectra
for i = 1:length(zeta_array)
    zeta = zeta_array(i);
    for j = 1:length(T)
        wn = wn_array(j);
        k = m * wn^2;              
        c = 2 * m * wn * zeta;      

         % Initialize response as we need it for next time step
        u = zeros(length(ug), 1); %defining displacment
        v = zeros(length(ug), 1); %defining velocity
        a = zeros(length(ug), 1); %defining acc
        
        p_eff = -m * ug; %intial load applied as building direction is oppostise to direction of earth 
        
        a(1) = (p_eff(1) - c * v(1) - k * u(1))/ m; % intial accleration - 0 values as initial is not defined in matlab
        
        k_hat = k + (gamma / (beta * dt)) * c + (1 / (beta * dt^2)) * m; % coeff from Method given in Textbook 
     
    A_const = (1 / (beta * dt)) * m + (gamma / beta) * c; % coeff from Method given in Textbook  
       B_const = (1 / (2 * beta)) * m + dt * (gamma / (2 * beta) - 1) * c; % coeff from Method given in Textbook 
        
       for step = 1:(length(ug) - 1) % values at each step time as given in textbook 
              dp_i = p_eff(step+1) - p_eff(step);
            
           dp_hat_i = dp_i + A_const * v(step) + B_const * a(step);
           du_i = dp_hat_i / k_hat;
            
          dv_i = (gamma / (beta * dt)) * du_i - (gamma / beta) * v(step) ...
                   + dt * (1 - gamma / (2 * beta)) * a(step);
            
          da_i = (1 / (beta * dt^2)) * du_i - (1 / (beta * dt)) * v(step) ...
                   - (1 / (2 * beta)) * a(step);
              u(step+1) = u(step) + du_i;
            v(step+1) = v(step) + dv_i;
            a(step+1) = a(step) + da_i;
        end
        % these values are used in next step of values for next time step 

Sd(i, j) = max(abs(u));          % Displacement
Sv(i, j) = wn * Sd(i, j);        % Velocity
Sa(i, j) = wn^2 * Sd(i, j);      % Acceleration
end
end

% Plot
colors = {'b', 'g', 'k', 'm'}; % colours are defined 
legend_labels = {'\zeta = 1%', '\zeta = 2%', '\zeta = 3%', '\zeta = 4%'};
figure('Name', 'Response Spectra', 'Position', [100, 100, 900, 700]);

% Displacement Response Spectrum
subplot(3,1,1);
hold on; grid on;
for i = 1:length(zeta_array)
    plot(T, Sd(i,:), colors{i}, 'LineWidth', 1.5);
end
title('Displacement Response Spectrum (S_d)');
ylabel('S_d- m ');
legend(legend_labels);

% Velocity Response Spectrum
subplot(3,1,2);
hold on; grid on;
for i = 1:length(zeta_array)
    plot(T, Sv(i,:), colors{i}, 'LineWidth', 1.5);
end
title('Velocity Response Spectrum (S_v)');
ylabel('S_v- m/s');

% Acceleration Response Spectrum
subplot(3,1,3);
hold on; grid on;
for i = 1:length(zeta_array)
    plot(T, Sa(i,:), colors{i}, 'LineWidth', 1.5);
end
title('Acceleration Response Spectrum (S_a)');
xlabel('Natural Period T (s)');
ylabel('S_a - m/s^2');

sgtitle('Newmark Method of Avg Accleration for Response Spectra');
