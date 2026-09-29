clear all
clc
close all

m = 1600;
k = 60000;
zeta = 0.02;

r = 0:0.02:3;

Rd = 1 ./ sqrt((1 - r.^2).^2 + (2 * zeta .* r).^2);

figure(1)
plot(r, Rd, 'b')
hold on
grid on

plot([1 1], [0 30], 'r--')
plot([0 3], [1 1], 'k--')

xlabel('Frequency Ratio')
ylabel('Dynamic Amplification Factor')
title('Harmonic Response for 2% Damping')
legend('DAF Curve', 'Resonance', 'Static Limit')

xlim([0 3])
ylim([0 30])

max_Rd = max(Rd);
disp(max_Rd)
