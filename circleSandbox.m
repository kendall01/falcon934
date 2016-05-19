clear all;
close all;
N = 50; 

[CA, SA, mdot_O2, mdot_f, ringD, centerD, NUM_RING_HOLES, t_f] = circleAreaFun(N);
tstep = t_f / N;
IN_TO_M = .0254;
R = .018; %m radius of outer ring. .018 works well. how far from the center of the grain does the circle of ring holes lie
OR = 1.995/2 * IN_TO_M; %outer radius of fuel grain

figure(1)
for i = 1:length(centerD)
    %c = (centerD(i)/2);
    %myfun = @(x,y) x^2 + y^2 - c;
    ezpolar(@(x)OR) %Draws exterior of fuel grain
    hold on
    [center(:,1), center(:,2)] = pol2cart(linspace(-2*pi, 2*pi,N)', linspace(centerD(i)/2,centerD(i)/2,N)');
    plot(center(:,1), center(:,2))
    hold on
    [ring(:,1), ring(:,2)] = pol2cart(linspace(-2*pi, 2*pi,N)', linspace(ringD(i)/2,ringD(i)/2,N)');
    for j = 1:NUM_RING_HOLES
        hold on
        x_off = R * cos((2*pi/NUM_RING_HOLES) * j);
        y_off = R * sin((2*pi/NUM_RING_HOLES) * j);
        plot(ring(:,1) + x_off, ring(:,2) + y_off)
    end
%     axis([-OR OR -OR OR])
    axis tight
    hold off
    pause(.0001); 
end

phi = mdot_O2 ./ mdot_f;
figure(2)
plot(linspace(0,t_f, N), phi)
title('Mixture Ratio')