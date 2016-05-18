clear all;
close all;

[CA, SA, mdot_f, ringD, centerD] = circleAreaFun();
NUM_RING_HOLES = 8;
IN_TO_M = .0254;
R = .018; %m radius of outer ring. how far from the center of the grain does the circle of ring holes lie
OR = 1.995/2 * IN_TO_M; %outer radius of fuel grain
N = 200; 

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
        x_off = R * cos((pi/4) * j);
        y_off = R * sin((pi/4) * j);
        plot(ring(:,1) + x_off, ring(:,2) + y_off)
    end
%     axis([-OR OR -OR OR])
    axis tight
    hold off
    pause(.01);
end
