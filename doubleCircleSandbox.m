clear all;
close all;

%Number of steps
N = 100;
M = 100;

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


centerD = centerD ./ .0254;
iRingD = iRingD ./ .0254;
oRingD = oRingD ./ .0254;
IR = .0115/.0254;
OR = .019/.0254;
l = linspace(0, 7.375, 99)';
FigHandle = figure(3);
set(FigHandle, 'Position', [100, 100, 1049, 300]);

for i = 1:size(centerD,1)
plot(l, centerD(i,1:M-1)'/2, 'r')
axis manual
axis equal
axis([0 7.375 -.9975 .9975])
xlabel('Length along grain [in]')
ylabel('Distance from Centerline [in]')
title('Lengthwise Cross-section of Fuel Grain Showing Tapered Burn')
hold on
plot(l, -centerD(i,1:M-1)'/2, 'r')
hold on
plot(l, (iRingD(i,1:M-1)'/2)+IR, 'b')
hold on
plot(l, (-iRingD(i,1:M-1)'/2)+IR, 'b')
hold on
plot(l, -((iRingD(i,1:M-1)'/2)+IR), 'b')
hold on
plot(l, -((-iRingD(i,1:M-1)'/2)+IR), 'b')
hold on
plot(l, (oRingD(i,1:M-1)'/2)+OR, 'g')
hold on
plot(l, (-oRingD(i,1:M-1)'/2)+OR, 'g')
hold on
plot(l, -((oRingD(i,1:M-1)'/2)+OR), 'g')
hold on
plot(l, -((-oRingD(i,1:M-1)'/2)+OR), 'g')
hold off
pause(0.03)
end


pause(1)

iRingDrl = [.34,.325];
oRingDrl = [.257, .23];
FigHandle4 = figure(4);
set(FigHandle4, 'Position', [100, 100, 1049, 300]);
plot(l([1,end]), centerD(i,[1,M-1])'/2, 'r')
axis manual
axis equal
axis([0 7.375 -.9975 .9975])
xlabel('Length along grain [in]')
ylabel('Distance from Centerline [in]')
title('Lengthwise Cross-section of Fuel Grain Showing Tapered Burn')
hold on;
plot(l([1,end]), 	   -centerD(N-1,[1,M-1])'/2, 			'r'); hold on;
plot(l([1,end]),  ( 	iRingD( N-1,[1,M-1])'/2)	+IR, 	'b'); hold on;
plot(l([1,end]),  (    -iRingD( N-1,[1,M-1])'/2)	+IR, 	'b'); hold on;
plot(l([1,end]), -(( 	iRingD( N-1,[1,M-1])'/2)	+IR), 	'b'); hold on;
plot(l([1,end]), -((   -iRingD( N-1,[1,M-1])'/2)	+IR), 	'b'); hold on;
plot(l([1,end]), (		oRingD( N-1,[1,M-1])'/2)	+OR, 	'g'); hold on;
plot(l([1,end]), (     -oRingD( N-1,[1,M-1])'/2)	+OR, 	'g'); hold on;
plot(l([1,end]), -((	oRingD( N-1,[1,M-1])'/2)	+OR), 	'g'); hold on;
plot(l([1,end]), -((   -oRingD( N-1,[1,M-1])'/2)	+OR), 	'g'); hold on;
plot(l([1,end]),  ( 	iRingDrl'/2)	+IR, 	'r', 'LineWidth', 1); hold on;
plot(l([1,end]),  (    -iRingDrl'/2)	+IR, 	'r', 'LineWidth', 1); hold on;
plot(l([1,end]), -(( 	iRingDrl'/2)	+IR), 	'r', 'LineWidth', 1); hold on;
plot(l([1,end]), -((   -iRingDrl'/2)	+IR), 	'r', 'LineWidth', 1); hold on;
plot(l([1,end]), (		oRingDrl'/2)	+OR, 	'r', 'LineWidth', 2); hold on;
plot(l([1,end]), (     -oRingDrl'/2)	+OR, 	'r', 'LineWidth', 2); hold on;
plot(l([1,end]), -((	oRingDrl'/2)	+OR), 	'r', 'LineWidth', 2); hold on;
plot(l([1,end]), -((   -oRingDrl'/2)	+OR), 	'r', 'LineWidth', 2)
hold off