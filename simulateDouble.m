function [] = simulateDouble(ringD, centerD, numHoles, N)

OUT =2;
IN = 1;
IN_TO_M = .0254;
RING_FACTOR = numHoles(OUT)/numHoles(IN); % ratio of radius of outer ring to inner ring. Also, of number of holes since circumference is proportional to radius.
IR = .0115; %m radius of inner ring. .018 works well. how far from the center of the grain does the circle of ring holes lie
R = RING_FACTOR * IR;
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
    %outer ring
    [ring(:,1), ring(:,2)] = pol2cart(linspace(-2*pi, 2*pi,N)', linspace(ringD(i)/2,ringD(i)/2,N)');
    for j = 1:numHoles(OUT)
        hold on
        x_off = R * cos((2*pi/numHoles(OUT)) * j);
        y_off = R * sin((2*pi/numHoles(OUT)) * j);
        plot(ring(:,1) + x_off, ring(:,2) + y_off)
    end
    hold on
    %inner ring
    [inner_ring(:,1), inner_ring(:,2)] = pol2cart(linspace(-2*pi, 2*pi,N)', linspace(ringD(i)/2,ringD(i)/2,N)');
    for j = 1:numHoles(IN)
        hold on
        x_off = IR * cos((2*pi/numHoles(IN)) * j);
        y_off = IR * sin((2*pi/numHoles(IN)) * j);
        plot(inner_ring(:,1) + x_off, inner_ring(:,2) + y_off)
    end
    %     axis([-OR OR -OR OR])
    axis tight
    hold off
    pause(.0001);
end


end