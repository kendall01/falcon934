clear all;
close all;

%Number of steps
N = 100;
M = 1;

%Run calculations for burn rate
[CA, SA, mdot_O2, mdot_f, iRingD, oRingD, centerD, numHoles, t_f] = doubleCircleAreaFun(N, M);
%[CA, SA, mdot_O2, mdot_f, iRingD, oRingD, centerW, centerH, ringHoles, t_f] = doubleCircleLineAreaFun(N);

%Simulate Burn
%simulateDouble(iRingD(:,1), oRingD(:,1), centerD(:,1), numHoles, N);

%Calculate Mixture Ratio
phi = mdot_O2 ./ mdot_f;
figure(2)
plot(linspace(0,t_f, N), phi)
title('Mixture Ratio')
xlabel('Burn Time [s]')
ylabel('Mixture Ratio')

figure(3)
plot(1:M-1, centerD(:,1:M-1)')
