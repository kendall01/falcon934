clear all;
close all;

%Number of steps
N = 100;

%Run calculations for burn rate
[~, ~, mdot_O2, mdot_f, ringD, centerD, numHoles, t_f] = doubleCircleAreaFun(N);

%Simulate Burn
simulateDouble(ringD, centerD, numHoles, N);

%Calculate Mixture Ratio
phi = mdot_O2 ./ mdot_f;
figure(2)
plot(linspace(0,t_f, N), phi)
title('Mixture Ratio')
xlabel('Burn Time [s]')
ylabel('Mixture Ratio')